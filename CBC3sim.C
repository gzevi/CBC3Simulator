#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include "TStyle.h"
#include "TMath.h"
#include <TString.h>

using namespace std;


// External inputs
int       NEvents = 100;
int       Nck = 3;
int       Vcth = 60; // DAC units, already pedestal-subtracted
int       HipSuppress = 1; // 0 means don't suppress
int       DLL = 10;
int       verbose = 0; // 1: every event, 2: every clock, 3: every ns
bool      diagnosticHistograms = true; // Make histograms for every event

//////////// Some useful numbers: /////////////////
// 1 Vcth ~ 130 electrons
// 1 MIP ~ 2.5 fC
// Pulse shape peak for 2.5, 5, 7.5, 10 fC is ~ 100, 225, 325, 375 Vcth units 
// MPV of Landau in beam test ~ 144 Vcth ~ 2.5 fC
///////////////////////////////////////////////////


TF1 * Function_Landau() {
  // Currently usina a Landau in Vcth units, but could use fC and then scale the PreampShaper output
  TF1* f = new TF1("f", "landau", 0, 300);
  f->SetParameters(1, 144, 17);
  return f;
}

TF1 * Function_PulseShape(float charge) {
  TF1* f = new TF1("f1", "gaus", 0, Nck*25);
  f->SetParameters(charge, 20, 10);
  return f;
}

float Step0_GetCharge( TF1* f) {
  // Return fC (based on Landau*Gaussian)
  float charge = 2.5; // 1 MIP ~ 2.5 fC
  charge = f->GetRandom();
  return charge; 
}

float Step1_PreampShaper( float charge , int time, TF1* f ) {
  // Return mV at each ns (based on parametrized shape)
  // Pulse shape is longer than a clock cycle
  
  // DLL delay shifts the shape (not sure if forward or backwards)
  int timeCorr = time - DLL;
  float voltage = 0;

  // Test square shape: 0000011111111111111111111100000
  // if ( timeCorr > 5 && timeCorr < 50 )
  //   voltage = charge * 30;

  // External function
  voltage = f->Eval(timeCorr);

  return voltage;
}

bool Step2_Comparator( float voltage , float threshold ) {
  // Return 0 or 1
  // Need to implement hysteresis (will requirer historical knowledge)
  return voltage > threshold;
}

bool Step3_HitDetectSampled( bool comp , int time ) {
  // Return 0 or 1 
  // Also called "40 MHz"
  // [CBC3 Manual] The output from the comparator is sampled using the 40MHz clock from the Delay Locked Loop. Only comparator outputs present on the rising edge of the clock will be captured and the output will only return to zero on the first rising clock edge following the comparators return to zero. The minimum width of output pulse is one clock cycle. Hits following immediately one after another in subsequent clock cycles will be captured even if the channel does not return below the comparator threshold for each hit.
  if (time != 0) return false;
  else return comp;
} 


bool Step3_HitDetectLatched( bool comp , int time, bool compMinusOneNs ) {
  // Also called "Fixed Width"
  // [CBC3 Manual] The output from the Hit comparator is latched for a full 25ns clock period. Non-synchronous comparator output transitions are captured. The output is a fixed 25ns pulse, regardless of the width of the comparator output pulse. Hits following immediately one after another in subsequent clock cycles will be captured provided the channel returns below the comparator threshold for each hit.
  
  // [CBC3 Manual] 2.7 ns "Setup & Hold blind spot" at the beginning of clock cycle
  if (time < 3) return false;

  // Capture transition from 0 to 1
  if ( comp && !compMinusOneNs ) return true;
  else return false;
}


bool Step4_HIPSuppress( bool hit , int &countPreviousHits ) {
  if ( countPreviousHits >= HipSuppress ) {
    return false;
  }
  else return hit;
}



TH1D * RunTest (TFile * f, TString dirName) {

  f->cd();
  TDirectory * testDir = f->mkdir(dirName);
  testDir->cd();


  // Summary Histograms, per run
  TH1D * h_sumHits = new TH1D("SumHits", "SumHits", 5, 0.5, 5.5);
  h_sumHits->GetXaxis()->SetBinLabel(1, "Sampled");
  h_sumHits->GetXaxis()->SetBinLabel(2, "Latched");
  h_sumHits->GetXaxis()->SetBinLabel(3, "OR");
  h_sumHits->GetXaxis()->SetBinLabel(4, "Sampled_HIP");
  h_sumHits->GetXaxis()->SetBinLabel(5, "OR_HIP");
  gStyle->SetOptStat(0);


  TF1 * landau = Function_Landau();

  // Event looper
  for (int ievent = 0; ievent < NEvents; ievent++) {
    if (verbose>=1) cout<<"New Event"<<endl;


    // Diagnostic histograms, per event
    TDirectory * eventDir;
    TH1D* h_pulse;
    TH1D* h_comp;
    TH1D* h_hits;


    if (diagnosticHistograms) {
      eventDir = testDir->mkdir( Form("Event%i", ievent) );
      eventDir->cd();
      h_pulse = new TH1D("PulseShape", "PulseShape", Nck*25, 0, Nck*25);
      h_pulse->GetXaxis()->SetNdivisions(3, kFALSE);
      h_comp  = (TH1D*) h_pulse->Clone();
      h_comp->SetNameTitle("Comparator", "Comparator");
      h_hits = (TH1D*) h_sumHits->Clone();
      h_hits->SetNameTitle("Hits", "Hits");
    }

    float charge = Step0_GetCharge(landau); 

    // Pulse shape function defined right after charge is known
    TF1 * shape = Function_PulseShape(charge);

    // Internal counters
    bool comparatorOutputPreviousNanosecond = false;
    int  HIPcountSampled = 0;
    int  HIPcountOR = 0;

    // Flags for found hits
    bool foundSampledHitThisClock = false;
    bool foundLatchedHitThisClock = false;
    bool foundLatchedHitNextClock = false;

    // Count hits over several clock cycles: 
    int countSampledHits = 0;
    int countLatchedHits = 0;
    int countORHits = 0;
    int countSampledHIPHits = 0;
    int countORHIPHits = 0;


    // Loop over clock cycles
    for (int ick = 0; ick < Nck; ick++) {

      if (verbose>=2) cout<<"Clock cycle "<<ick<<endl;

      foundLatchedHitThisClock = foundLatchedHitNextClock;
      foundSampledHitThisClock = false;
      foundLatchedHitNextClock = false;

      if (verbose>=3) cout << "ns \t Voltage \t Comparator "<<endl;
      // Loop inside a clock cycle
      for (int ins = 0; ins <= 24; ins++) {

        float voltage = Step1_PreampShaper( charge, ins + ick*25, shape );
        bool comp = Step2_Comparator( voltage, Vcth );

        if ( !foundSampledHitThisClock ) foundSampledHitThisClock = Step3_HitDetectSampled( comp, ins );
        if ( !foundLatchedHitNextClock ) foundLatchedHitNextClock = Step3_HitDetectLatched( comp, ins, comparatorOutputPreviousNanosecond );

        if (verbose>=3) cout << ins << " \t " << voltage << "\t\t" << comp << endl;

        if (diagnosticHistograms) {
          h_pulse->SetBinContent(ins+ick*25, voltage);
          h_comp->SetBinContent(ins+ick*25, comp*100); // arbitrary scaling just for easier overlay
        }

        comparatorOutputPreviousNanosecond = comp;
      } // End of clock cycle

      bool foundORHit = ( foundLatchedHitThisClock || foundSampledHitThisClock );
      bool foundSampledHitHIP = Step4_HIPSuppress( foundSampledHitThisClock, HIPcountSampled );
      bool foundORHitHIP      = Step4_HIPSuppress( foundORHit, HIPcountOR );


      // HIP counters
      HIPcountSampled = foundSampledHitThisClock  ? HIPcountSampled+1 : 0;
      HIPcountOR      = foundORHit                ? HIPcountOR+1      : 0;

      // Hit counters
      if (foundSampledHitThisClock)  { countSampledHits++;    h_sumHits->Fill(1); if (diagnosticHistograms) h_hits->Fill(1); }
      if (foundLatchedHitThisClock)  { countLatchedHits++;    h_sumHits->Fill(2); if (diagnosticHistograms) h_hits->Fill(2); } 
      if (foundORHit)                { countORHits++;         h_sumHits->Fill(3); if (diagnosticHistograms) h_hits->Fill(3); }   
      if (foundSampledHitHIP)        { countSampledHIPHits++; h_sumHits->Fill(4); if (diagnosticHistograms) h_hits->Fill(4); }   
      if (foundORHitHIP)             { countORHIPHits++;      h_sumHits->Fill(5); if (diagnosticHistograms) h_hits->Fill(5); } 

      if (verbose>=2) cout<<"Sampled\tLatched\tOR   SampledHIP\tORHIP  [LatchedNextCock]"<<endl;
      if (verbose>=2) cout<<foundSampledHitThisClock<<" \t "<<foundLatchedHitThisClock<<" \t "<<foundORHit<<" \t "<<foundSampledHitHIP<<" \t "<<foundORHitHIP<<" \t "<<foundLatchedHitNextClock << endl;

    } // End of clock

    if (verbose>=1) cout<<"End of event: ";    
    if (verbose>=1) cout<<"Sampled  Latched   OR   SampledHIP    ORHIP "<<endl;
    if (verbose>=1) cout<<"\t\t"<<countSampledHits<<" \t "<<countLatchedHits<<"         "<<countORHits<<"    \t "<<countSampledHIPHits<<"   \t "<<countORHIPHits<<endl;  

    if (diagnosticHistograms) {
      h_pulse->Write();
      h_comp->Write();
      h_hits->Write();
    }

  } // End of event


  testDir->cd();
  h_sumHits->Scale(1./NEvents);
  h_sumHits->GetYaxis()->SetRangeUser(0, h_sumHits->GetMaximum()*1.3);
//  h_sumHits->Draw();
  h_sumHits->Write();
  return h_sumHits;

} // End of RunTest


void CBC3sim () {


  TFile * f = new TFile("TestSim.root", "RECREATE");
  f->cd();

  // Threshold scan
  float vcthStart = 60;
  float vcthRange = 100;
  float vcthStep = 5;
  TH2D * h_vcthScan = new TH2D(Form("ThresholdScan_DLL%i", DLL), "ThresholdScan", (int) vcthRange/vcthStep, vcthStart, vcthStart+vcthRange, 5, 0.5, 5.5);
  for (int i = 0; i <= vcthRange; i += vcthStep) {
    Vcth = vcthStart + i;
    TH1D* h_test = RunTest(f, Form("Vcth_%i", Vcth));
    for (int ibin = 1; ibin <= 5; ibin++ ) {
      h_vcthScan->SetBinContent(h_vcthScan->GetXaxis()->FindBin(Vcth+0.01), ibin, h_test->GetBinContent(ibin));
    }
  }

  // Timing scan
  Vcth = 60; // Reset Vcth to nominal
  float DLLStart = 0;
  float DLLRange = 25;
  float DLLStep = 1;
  TH2D * h_dllScan = new TH2D(Form("DLLScan_Vcth%i", Vcth), "DLLScan", (int) DLLRange/DLLStep, DLLStart, DLLStart+DLLRange, 5, 0.5, 5.5);
  for (int i = 0; i <= DLLRange; i += DLLStep) {
    DLL = DLLStart + i;
    TH1D* h_test = RunTest(f, Form("DLL_%i", DLL));
    for (int ibin = 1; ibin <= 5; ibin++ ) {
      h_dllScan->SetBinContent(h_dllScan->GetXaxis()->FindBin(DLL+0.01), ibin, h_test->GetBinContent(ibin));
    }
  }

  f->cd();
  h_vcthScan->Write();
  h_dllScan->Write();

}


