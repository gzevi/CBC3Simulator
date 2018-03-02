#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TRandom3.h"

using namespace std;

// Currently see no reason to modify these once they are set
const int       NEvents = 1000; // Number of events for each RunTest call
const int       Nck = 4; // Number of clock cycles
const int       HipSuppress = 1; // 0 means don't suppress

// Knobs that are overwritten when performing scans
int       Vcth = 20; // DAC units, already pedestal-subtracted
int       DLL = 10; // DLL delay. If DLL=0, signal pulse starts at t=0. If DLL=10, it starts at 10

// Verbosity etc
int       verbose = 0; // 0: none, 1: every event, 2: every clock, 3: every ns
bool      diagnosticHistograms = false; // Make histograms for every event (reduce NEvents to avoid slowing down)
bool      doVcthScan = true;
bool      doDLLScan = true;
bool      do2DScan = true;

// best guess at shape of flandau 
double defLandauPars[3]={1.0, 144 , 9.0};
double defPulseShapePars[4]={3, 4, 1, 0.0};

//////////// Some useful numbers: /////////////////
// 1 Vcth ~ 130 electrons
// 1 MIP ~ 2.5 fC
// Pulse shape peak for 2.5, 5, 7.5, 10 fC is ~ 100, 225, 325, 375 Vcth units 
// MPV/width of Landau in beam test ~ 144/9 Vcth 
///////////////////////////////////////////////////

TF1 *gaussianNoise() {
  // Use a gaussian around 0 (pedestal subtracted) with a noise of 6 Vcth
  // f(x) = p0*exp(-0.5*((x-p1)/p2)^2)
  TF1 *f = new TF1("f", "gaus", 0, 300);
  f->SetParameters(1, 0, 6);
  f->SetParNames("Normalisation", "Mean (pedestal)", "Sigma (noise)");
  return f;
}
TF1 * Function_Landau(double* pars, TString pFuncName="f", double pXmin=0,double pXmax=300) 
{
  // Currently using a Landau in Vcth units, but could use fC and then scale the PreampShaper output
  TF1* f = new TF1(pFuncName.Data(), "landau", pXmin, pXmax);
  f->SetParameters(pars);//1, 144, 9);
  return f;
}
double StepFunction(double* x , double *par)
{
  double x0 = x[0];
  double tStart = par[0];
  double tStop = par[1];
  double riseTime = par[2];

  return 0.5*(TMath::Erf((x0-tStart)/(sqrt(2.)*riseTime)) - TMath::Erf((x0-tStop)/(sqrt(2.)*riseTime)));
}
TF1 * Function_Comparator(double* pars, TString pFuncName="f", double pXmin=0,double pXmax=300) 
{
  // Currently using a Landau in Vcth units, but could use fC and then scale the PreampShaper output
  TF1* f = new TF1(pFuncName.Data(), StepFunction, pXmin, pXmax, 3);
  f->SetParameters(pars);//1, 144, 9);
  return f;
}

TF1 *Function_Landau() {
  // Currently using a Landau in Vcth units, but could use fC and then scale the PreampShaper output
  TF1 *f = new TF1("f", "landau", 0, 300);
  f->SetParameters(1, 144, 9);
  f->SetParNames("Normalisation", "Most probably value", "Sigma");
  return f;
}

double Double_LanGaus(double *x, double *par) {

  // Fit parameters:
  // par[0] = Width (scale) parameter of Landau density
  // par[1] = Most Probable Value parameter of Landau density  
  // par[2] = Total area (integral -inf to inf, normalization constant)
  // par[3] = Width (sigma) of convoluted Gaussian function

  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  double invSqrd2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift    = -0.22278298;       // Landau maximum location

  // Control constants
  double nCsteps     = 100.0;      // number of convolution steps
  double sigmasGauss =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  double xx;
  double mpc;
  double fLandau;
  double sum = 0.0;
  double step;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  double xlow = x[0] - sigmasGauss * par[3];
  double xupp = x[0] + sigmasGauss * par[3];

  step = (xupp-xlow) / nCsteps;

  // Convolution integral of Landau and Gaussian by sum
  for(double i = 1.0; i <= nCsteps/2; i++) {

    xx = xlow + (i-.5) * step;
    fLandau = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fLandau * TMath::Gaus(x[0], xx, par[3]);

    xx = xupp - (i-.5) * step;
    fLandau = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fLandau * TMath::Gaus(x[0], xx, par[3]);

  }

  return (par[2] * step * sum * invSqrd2pi / par[3]);

}

TF1 * Function_LanGaus() {
  // Currently using a Landau in Vcth units, but could use fC and then scale the PreampShaper output
  TF1* f = new TF1("f", "Double_LanGaus", 0, 300, 4);
  // Initial parameter values to be checked!
  // I set Landau "width" to 9 (is this the same as sigma?)
  // MVP is 144
  // Total area of the convolution is 1 for now, this is a normalization factor
  // Sigma of the gaussian is 6 (CBC3 noise is 6 Vcth units??? To be checked!!!)
  f->SetParameters(9, 144, 1, 6);
  f->SetParNames("Landau width", "Most probably value", "Normalisation", "Gaussian noise");
  return f;
}

double Function_RC(double *x, double *par) {

  double N = par[0];
  double tau = par[1];
  double norm = par[2];
  
  double result =  norm / TMath::Gamma(N+1) * TMath::Power( (x[0]/tau), N) * TMath::Exp(-x[0] / tau);
  if (result < 0) result = 0;
  return result;

}


TF1 * Function_PulseShape(float charge, double* pars, TString pFuncName="f1", double pXmin=-25 , double pXmax=(Nck-1)*25) {
  // TF1* f = new TF1("f1", "gaus", 0, Nck*25);
  // f->SetParameters(charge, 20, 10);

  TF1* f = new TF1(pFuncName.Data(), Function_RC, pXmin, pXmax, 4);
  f->SetParameters(pars);//3, 4, 1,0);

  // Normalize so that peak matches "charge". This way if we draw "140" from the Landau distribution,
  // we will have a PulseShape that has a maximum value of 140
  double max = f->GetMaximum(0, 1000);
  f->SetParameter(2, charge/max);
  //f->SetParameters(3, 4, charge/max,0.0);

  return f;
}

float Step0_GetCharge( TF1* f) {
  // Return fC (based on Landau*Gaussian)
  float charge = 2.5; // 1 MIP ~ 2.5 fC (But this is not used?)
  charge = f->GetRandom(); // (Now charge is in Vcth again)
  return charge; 
}

float Step0_GetNoise( TF1* f) {
  // Return fC (based on Landau*Gaussian)
  //float noise = 2.5; // 1 MIP ~ 2.5 fC
  float noise = f->GetRandom();
  return noise; 
}

float Step1_PreampShaper( int time, TF1* f ) {
  // Return mV at each ns (based on parametrized shape)
  // Pulse shape is longer than a clock cycle
  
  // DLL delay shifts the shape 
  int timeCorr = time - DLL;
  float voltage = 0;

  // External function
  voltage = f->Eval(timeCorr);

  return voltage;
}

float Step1_PreampShaperDoubleHit( int time, TF1* f, TF1* f2 ) {
  // Return mV at each ns (based on parametrized shape)
  // Pulse shape is longer than a clock cycle
  
  // DLL delay shifts the shape 
  int timeCorr = time - DLL;
  float voltage = 0;

  // External function
  voltage = TMath::Max( f->Eval(timeCorr), f2->Eval(timeCorr-25));

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

void labelAxis(TAxis * a) {
  a->SetBinLabel(1, "Sampled");
  a->SetBinLabel(2, "Latched");
  a->SetBinLabel(3, "OR");
  a->SetBinLabel(4, "Sampled_HIP");
  a->SetBinLabel(5, "OR_HIP");
}

vector<TH1D*> RunTest (TDirectory * f, TString dirName) {

  f->cd();
  TDirectory * testDir = f->mkdir(dirName);
  testDir->cd();

  // Summary Histograms, per run
  TH1D * h_sumHits = new TH1D("SumHits", "SumHits", 5, 0.5, 5.5);
  labelAxis(h_sumHits->GetXaxis());
  TH1D * h_sumHitsNextClock = new TH1D("SumHitsNextClock", "SumHitsNextClock", 5, 0.5, 5.5);
  labelAxis(h_sumHitsNextClock->GetXaxis());
  TH1D * h_sumHitsInTwoClocks = new TH1D("SumHitsInTwoClocks", "SumHitsInTwoClocks", 5, 0.5, 5.5);
  labelAxis(h_sumHitsInTwoClocks->GetXaxis());
  gStyle->SetOptStat(0);

  //TF1 * landau = Function_LanGaus(); // This is the convolution!
  TF1 * landau = Function_Landau();
  TF1 * noise = gaussianNoise(); 

  // Event looper
  for (int ievent = 0; ievent < NEvents; ievent++) {
    if (verbose>=1) cout<<"New Event"<<endl;

    // Diagnostic histograms, per event
    TDirectory * eventDir;
    TH1D* h_pulse;
    TH1D* h_comp;
    TH1D* h_hits;
    TH1D* h_hitsNextClock;

    if (diagnosticHistograms) {
      eventDir = testDir->mkdir( Form("Event%i", ievent) );
      eventDir->cd();
      h_pulse = new TH1D("PulseShape", "PulseShape", Nck*25, -25, (Nck-1)*25);
      h_pulse->GetXaxis()->SetNdivisions(4, kFALSE);
      h_comp  = (TH1D*) h_pulse->Clone();
      h_comp->SetNameTitle("Comparator", "Comparator");
      h_hits = (TH1D*) h_sumHits->Clone();
      h_hits->Reset();
      h_hits->SetNameTitle("Hits", "Hits");
      h_hitsNextClock = (TH1D*) h_hits->Clone();
      h_hitsNextClock->SetNameTitle("HitsNextClock", "HitsNextClock");
    }

    // Get the charge randomly from the landau
    float charge = Step0_GetCharge(landau); 
    float charge2 = Step0_GetCharge(landau); 

    // Add a random noise to the charge
    charge += Step0_GetNoise(noise);
    charge2 += Step0_GetNoise(noise);

    // Pulse shape function defined right after charge is known
    TF1 * shape  = Function_PulseShape(charge,defPulseShapePars);
    TF1 * shape2 = Function_PulseShape(charge2,defPulseShapePars);


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
    for (int ick = -1; ick < Nck; ick++) {

      if (verbose>=2) cout<<"Clock cycle "<<ick<<endl;

      foundLatchedHitThisClock = foundLatchedHitNextClock;
      foundSampledHitThisClock = false;
      foundLatchedHitNextClock = false;

      if (verbose>=3) cout << "ns \t Voltage \t Comparator "<<endl;
      // Loop inside a clock cycle
      for (int ins = 0; ins <= 24; ins++) {

        int coarsePlusFineTime = ins + ick*25;

        float voltage = Step1_PreampShaper( coarsePlusFineTime, shape );
        float voltageDoubleHit = Step1_PreampShaperDoubleHit( coarsePlusFineTime, shape, shape2 );
        //voltage = voltageDoubleHit;
        bool comp = Step2_Comparator( voltage, Vcth );

        if ( !foundSampledHitThisClock ) foundSampledHitThisClock = Step3_HitDetectSampled( comp, ins );
        if ( !foundLatchedHitNextClock ) foundLatchedHitNextClock = Step3_HitDetectLatched( comp, ins, comparatorOutputPreviousNanosecond );

        if (verbose>=3) cout << coarsePlusFineTime << " \t " << voltage << "\t\t" << comp << endl;

        if (diagnosticHistograms) {
          int timeBin = h_pulse->GetXaxis()->FindBin(coarsePlusFineTime);
          h_pulse->SetBinContent(timeBin, voltage);
          h_comp->SetBinContent(timeBin, comp*Vcth); 
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
      if (ick == 1) { // "In time" hits (arbitrarily decide that I want hit to be detected in 2nd bunch crossing)
        if (foundSampledHitThisClock)  { countSampledHits++;    h_sumHits->Fill(1); if (diagnosticHistograms) h_hits->Fill(1); }
        if (foundLatchedHitThisClock)  { countLatchedHits++;    h_sumHits->Fill(2); if (diagnosticHistograms) h_hits->Fill(2); } 
        if (foundORHit)                { countORHits++;         h_sumHits->Fill(3); if (diagnosticHistograms) h_hits->Fill(3); }   
        if (foundSampledHitHIP)        { countSampledHIPHits++; h_sumHits->Fill(4); if (diagnosticHistograms) h_hits->Fill(4); }   
        if (foundORHitHIP)             { countORHIPHits++;      h_sumHits->Fill(5); if (diagnosticHistograms) h_hits->Fill(5); } 

        if (verbose>=2) cout<<"Sampled\tLatched\tOR   SampledHIP\tORHIP  [LatchedNextCock]"<<endl;
        if (verbose>=2) cout<<foundSampledHitThisClock<<" \t "<<foundLatchedHitThisClock<<" \t "<<foundORHit<<" \t "<<foundSampledHitHIP<<" \t "<<foundORHitHIP<<" \t "<<foundLatchedHitNextClock << endl;
      }

      if (ick == 2) { // Late and Double hits

        if (foundSampledHitThisClock)  {  
          h_sumHitsNextClock->Fill(1); 
          if (diagnosticHistograms) h_hitsNextClock->Fill(1); 
          if (countSampledHits>0)     h_sumHitsInTwoClocks->Fill(1);
        }
        if (foundLatchedHitThisClock)  {  
          h_sumHitsNextClock->Fill(2); 
          if (diagnosticHistograms) h_hitsNextClock->Fill(2); 
          if (countLatchedHits>0)     h_sumHitsInTwoClocks->Fill(2);
        } 
        if (foundORHit)                {  
          h_sumHitsNextClock->Fill(3); 
          if (diagnosticHistograms) h_hitsNextClock->Fill(3); 
          if (countORHits>0)          h_sumHitsInTwoClocks->Fill(3);
        }   
        if (foundSampledHitHIP)        {  
          h_sumHitsNextClock->Fill(4); 
          if (diagnosticHistograms) h_hitsNextClock->Fill(4); 
          if (countSampledHIPHits>0)  h_sumHitsInTwoClocks->Fill(4);
        }   
        if (foundORHitHIP)             {  
          h_sumHitsNextClock->Fill(5); 
          if (diagnosticHistograms) h_hitsNextClock->Fill(5); 
          if (countORHIPHits>0)       h_sumHitsInTwoClocks->Fill(5);
        } 

      }

    } // End of clock

    if (verbose>=1) cout<<"End of event: ";    
    if (verbose>=1) cout<<"Sampled  Latched   OR   SampledHIP    ORHIP "<<endl;
    if (verbose>=1) cout<<"\t\t"<<countSampledHits<<" \t "<<countLatchedHits<<"         "<<countORHits<<"    \t "<<countSampledHIPHits<<"   \t "<<countORHIPHits<<endl;  

    if (diagnosticHistograms) {
      h_pulse->Write();
      h_comp->Write();
      h_hits->Write();
      h_hitsNextClock->Write();
    }

  } // End of event

  testDir->cd();
  h_sumHits->Scale(1./NEvents);
  h_sumHits->GetYaxis()->SetRangeUser(0, h_sumHits->GetMaximum()*1.3);
  h_sumHitsNextClock->Scale(1./NEvents);
  h_sumHitsNextClock->GetYaxis()->SetRangeUser(0, h_sumHitsNextClock->GetMaximum()*1.3);
  h_sumHitsInTwoClocks->Scale(1./NEvents);
  h_sumHitsInTwoClocks->GetYaxis()->SetRangeUser(0, h_sumHitsInTwoClocks->GetMaximum()*1.3);
 
//  h_sumHits->Draw();
  h_sumHits->Write();
  h_sumHitsNextClock->Write();
  h_sumHitsInTwoClocks->Write();

  vector<TH1D*> v_h_sumHits;
  v_h_sumHits.push_back(h_sumHits);
  v_h_sumHits.push_back(h_sumHitsNextClock);
  v_h_sumHits.push_back(h_sumHitsInTwoClocks);

  return v_h_sumHits;

} // End of RunTest

void CBC3sim () {

  TFile * f = new TFile("TestSim.root", "RECREATE");
  f->cd();

  // Threshold scan
  DLL = 15;
  TDirectory * VcthDir = f->mkdir(Form("ThresholdScan_DLL%i", DLL));
  VcthDir->cd();
  float vcthStart = 20;
  float vcthRange = 180;
  float vcthStep = 2;
  TH2D * h_vcthScan = new TH2D(Form("ThresholdScan_DLL%i", DLL), "ThresholdScan", (int) vcthRange/vcthStep, vcthStart, vcthStart+vcthRange, 5, 0.5, 5.5);
  labelAxis(h_vcthScan->GetYaxis());
  TH2D * h_vcthScanDoubleHits = (TH2D*) h_vcthScan->Clone();
  h_vcthScanDoubleHits->SetNameTitle(Form("ThresholdScanDoubleHits_DLL%i", DLL), "ThresholdScanDoubleHits");
  if (doVcthScan) {
    cout<<"Vcth scan: "<<flush;
    for (int i = 0; i <= vcthRange; i += vcthStep) {
      Vcth = vcthStart + i;
      cout<<Vcth<<" "<<flush;
      vector<TH1D*> v_h_test = RunTest(VcthDir, Form("Vcth_%i", Vcth));
      int vcthBin = h_vcthScan->GetXaxis()->FindBin(Vcth+0.01);
      for (int ibin = 1; ibin <= 5; ibin++ ) {
        h_vcthScan->SetBinContent(vcthBin, ibin, v_h_test.at(0)->GetBinContent(ibin));
        h_vcthScanDoubleHits->SetBinContent(vcthBin, ibin, v_h_test.at(2)->GetBinContent(ibin));
      }
    }
    cout<<endl;
  }

  // Timing scan
  Vcth = 20; // Reset Vcth to nominal
  TDirectory * DLLDir = f->mkdir(Form("DLLScan_Vcth%i", Vcth));
  DLLDir->cd();
  float DLLStart = -25;
  float DLLRange = 100;
  float DLLStep = 1;
  TH2D * h_dllScan = new TH2D(Form("DLLScan_Vcth%i", Vcth), "DLLScan", (int) DLLRange/DLLStep, DLLStart, DLLStart+DLLRange, 5, 0.5, 5.5);
  labelAxis(h_dllScan->GetYaxis());
  TH2D * h_dllScanNextClock = (TH2D*) h_dllScan->Clone();
  TH2D * h_dllScanDoubleHits = (TH2D*) h_dllScan->Clone();
  h_dllScanNextClock->SetNameTitle(Form("DLLScanNextClock_Vcth%i", Vcth), "DLLScanNextClock");
  h_dllScanDoubleHits->SetNameTitle(Form("DLLScanDoubleHits_Vcth%i", Vcth), "DLLScanDoubleHits");
  if (doDLLScan) {
    cout<<"DLL scan: "<<flush;
    for (int i = 0; i <= DLLRange; i += DLLStep) {
      DLL = DLLStart + i;
      cout<<DLL<<" "<<flush;
      vector<TH1D*> v_h_test = RunTest(DLLDir, Form("DLL_%i", DLL));
      int DLLbin = h_dllScan->GetXaxis()->FindBin(DLL+0.01);
      for (int ibin = 1; ibin <= 5; ibin++ ) {
        h_dllScan->SetBinContent(DLLbin, ibin, v_h_test.at(0)->GetBinContent(ibin));
        h_dllScanNextClock->SetBinContent(DLLbin, ibin, v_h_test.at(1)->GetBinContent(ibin));
        h_dllScanDoubleHits->SetBinContent(DLLbin, ibin, v_h_test.at(2)->GetBinContent(ibin));
      }
    }
    cout<<endl;
  }

  // 2D scan
  TDirectory * TwoDDir = f->mkdir("2DScan");
  TwoDDir->cd();

  vcthStart = 20;
  vcthRange = 160;
  vcthStep = 5;
  DLLStart = -25;
  DLLRange = 50;
  DLLStep = 1;

  TH2D * h_2DScanSampled = new TH2D("2DScanSampled", "2DScanSampled", (int) DLLRange/DLLStep, DLLStart, DLLStart+DLLRange, (int) vcthRange/vcthStep, vcthStart, vcthStart+vcthRange);
  TH2D * h_2DScanLatched = (TH2D*) h_2DScanSampled->Clone(); h_2DScanLatched->SetNameTitle("2DScanLatched", "2DScanLatched");
  TH2D * h_2DScanOR = (TH2D*) h_2DScanSampled->Clone(); h_2DScanOR->SetNameTitle("2DScanOR", "2DScanOR");
  TH2D * h_2DScanSampledHIP = (TH2D*) h_2DScanSampled->Clone(); h_2DScanSampledHIP->SetNameTitle("2DScanSampledHIP", "2DScanSampledHIP");
  TH2D * h_2DScanORHIP = (TH2D*) h_2DScanSampled->Clone(); h_2DScanORHIP->SetNameTitle("2DScanORHIP", "2DScanORHIP");
  TH2D * h_2DScanSampledDoubleHits = (TH2D*) h_2DScanSampled->Clone(); h_2DScanSampledDoubleHits->SetNameTitle("2DScanSampledDoubleHits", "2DScanSampledDoubleHits");
  TH2D * h_2DScanLatchedDoubleHits = (TH2D*) h_2DScanSampled->Clone(); h_2DScanLatchedDoubleHits->SetNameTitle("2DScanLatchedDoubleHits", "2DScanLatchedDoubleHits");
  TH2D * h_2DScanORDoubleHits = (TH2D*) h_2DScanSampled->Clone(); h_2DScanORDoubleHits->SetNameTitle("2DScanORDoubleHits", "2DScanORDoubleHits");
  TH2D * h_2DScanSampledHIPDoubleHits = (TH2D*) h_2DScanSampled->Clone(); h_2DScanSampledHIPDoubleHits->SetNameTitle("2DScanSampledHIPDoubleHits", "2DScanSampledHIPDoubleHits");
  TH2D * h_2DScanORHIPDoubleHits = (TH2D*) h_2DScanSampled->Clone(); h_2DScanORHIPDoubleHits->SetNameTitle("2DScanORHIPDoubleHits", "2DScanORHIPDoubleHits");
  if (do2DScan) {
    cout<<"2D scan: "<<flush;
    for (int i = 0; i <= vcthRange; i += vcthStep) {
      for (int j = 0; j <= DLLRange; j += DLLStep) {
        Vcth = vcthStart + i;
        DLL = DLLStart + j;
        cout<<Vcth<<"/"<<DLL<<" "<<flush;
        vector<TH1D*> v_h_test = RunTest(TwoDDir, Form("Vcth%i_DLL%i", Vcth, DLL));
        int vcthBin = h_2DScanSampled->GetYaxis()->FindBin(Vcth+0.01);
        int DLLbin = h_2DScanSampled->GetXaxis()->FindBin(DLL+0.01);
        // cout<<"Ran with DLL "<<DLL<<" and Vcth "<<Vcth<<". Filling 2D scan "<<DLLbin<<" "<<vcthBin<<" with value "<<v_h_test.at(0)->GetBinContent(0)<<endl;
        h_2DScanSampled->SetBinContent(             DLLbin, vcthBin, v_h_test.at(0)->GetBinContent(1));
        h_2DScanLatched->SetBinContent(             DLLbin, vcthBin, v_h_test.at(0)->GetBinContent(2));
        h_2DScanOR->SetBinContent(                  DLLbin, vcthBin, v_h_test.at(0)->GetBinContent(3));
        h_2DScanSampledHIP->SetBinContent(          DLLbin, vcthBin, v_h_test.at(0)->GetBinContent(4));
        h_2DScanORHIP->SetBinContent(               DLLbin, vcthBin, v_h_test.at(0)->GetBinContent(5));
        h_2DScanSampledDoubleHits->SetBinContent(   DLLbin, vcthBin, v_h_test.at(2)->GetBinContent(1));
        h_2DScanLatchedDoubleHits->SetBinContent(   DLLbin, vcthBin, v_h_test.at(2)->GetBinContent(2));
        h_2DScanORDoubleHits->SetBinContent(        DLLbin, vcthBin, v_h_test.at(2)->GetBinContent(3));
        h_2DScanSampledHIPDoubleHits->SetBinContent(DLLbin, vcthBin, v_h_test.at(2)->GetBinContent(4));
        h_2DScanORHIPDoubleHits->SetBinContent(     DLLbin, vcthBin, v_h_test.at(2)->GetBinContent(5));
      }
    }
    cout<<endl;
  }

  f->cd();
  h_vcthScan->Write();
  h_vcthScanDoubleHits->Write();
  h_dllScan->Write();
  h_dllScanNextClock->Write();
  h_dllScanDoubleHits->Write();
  h_2DScanSampled->Write();
  h_2DScanLatched->Write();
  h_2DScanOR->Write();
  h_2DScanSampledHIP->Write();
  h_2DScanORHIP->Write();
  h_2DScanSampledDoubleHits->Write();
  h_2DScanLatchedDoubleHits->Write();
  h_2DScanORDoubleHits->Write();
  h_2DScanSampledHIPDoubleHits->Write();
  h_2DScanORHIPDoubleHits->Write();

}
