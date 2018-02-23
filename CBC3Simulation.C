//#include "simPars.h"
#include "Utils.h"
#include "simSteps.h"

// Parameters for the simulation
#include "simPars.h"

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
    TF1 * shape  = Function_PulseShape(charge, Nck);
    TF1 * shape2 = Function_PulseShape(charge2, Nck);

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

        float voltage = Step1_PreampShaper( coarsePlusFineTime, shape, DLL );
        float voltageDoubleHit = Step1_PreampShaperDoubleHit( coarsePlusFineTime, shape, shape2, DLL );
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
      bool foundSampledHitHIP = Step4_HIPSuppress( foundSampledHitThisClock, HIPcountSampled, HipSuppress );
      bool foundORHitHIP      = Step4_HIPSuppress( foundORHit, HIPcountOR, HipSuppress );

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
      } // End IF clockcycle == 1

      if (ick == 2) { // Late and Double hits

        if (foundSampledHitThisClock)  {  
          h_sumHitsNextClock->Fill(1); 
          if (diagnosticHistograms) h_hitsNextClock->Fill(1); 
          if (countSampledHits>0)     h_sumHitsInTwoClocks->Fill(1); // This means it was in two clockcycles because in clc == 1 it was made bigger than 1 if there was a hit
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

      } // End IF clockcycle == 2

    } // End of clock

    if (verbose>=1) std::cout << "End of event: ";    
    if (verbose>=1) std::cout << "Sampled  Latched   OR   SampledHIP    ORHIP " << std::endl;
    if (verbose>=1) std::cout << "\t\t" << countSampledHits << " \t " << countLatchedHits << "         " << countORHits << "    \t " << countSampledHIPHits << "   \t " << countORHIPHits << std::endl;

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

void CBC3Simulation () {

  TFile * f = new TFile(Form("./simresults/cbc3Simulation_%s.root", cSimFileName), "RECREATE");
  f->cd();

  // Threshold scan
  DLL = 15;
  TDirectory * VcthDir = f->mkdir(Form("ThresholdScan_DLL%i", DLL));
  VcthDir->cd();
  float vcthStart = 20;
  float vcthRange = 450;
  float vcthStep = 2;
  TH2D * h_vcthScan = new TH2D(Form("ThresholdScan_DLL%i", DLL), "ThresholdScan", (int) vcthRange/vcthStep, vcthStart, vcthStart+vcthRange, 5, 0.5, 5.5);
  labelAxis(h_vcthScan->GetYaxis());
  TH2D * h_vcthScanDoubleHits = (TH2D*) h_vcthScan->Clone();
  h_vcthScanDoubleHits->SetNameTitle(Form("ThresholdScanDoubleHits_DLL%i", DLL), "ThresholdScanDoubleHits");
  if (doVcthScan) {
    std::cout << "Vcth scan: " << std::flush;
    for (int i = 0; i <= vcthRange; i += vcthStep) {
      Vcth = vcthStart + i;
      std::cout << Vcth << " " << std::flush;
      vector<TH1D*> v_h_test = RunTest(VcthDir, Form("Vcth_%i", Vcth));
      int vcthBin = h_vcthScan->GetXaxis()->FindBin(Vcth+0.01);
      for (int ibin = 1; ibin <= 5; ibin++ ) {
        h_vcthScan->SetBinContent(vcthBin, ibin, v_h_test.at(0)->GetBinContent(ibin));
        h_vcthScanDoubleHits->SetBinContent(vcthBin, ibin, v_h_test.at(2)->GetBinContent(ibin));
      }
    }
    std::cout << std::endl;
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
    std::cout << "DLL scan: " << std::flush;
    for (int i = 0; i <= DLLRange; i += DLLStep) {
      DLL = DLLStart + i;
      std::cout << DLL << " " << std::flush;
      vector<TH1D*> v_h_test = RunTest(DLLDir, Form("DLL_%i", DLL));
      int DLLbin = h_dllScan->GetXaxis()->FindBin(DLL+0.01);
      for (int ibin = 1; ibin <= 5; ibin++ ) {
        h_dllScan->SetBinContent(DLLbin, ibin, v_h_test.at(0)->GetBinContent(ibin));
        h_dllScanNextClock->SetBinContent(DLLbin, ibin, v_h_test.at(1)->GetBinContent(ibin));
        h_dllScanDoubleHits->SetBinContent(DLLbin, ibin, v_h_test.at(2)->GetBinContent(ibin));
      }
    }
    std::cout << std::endl;
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
    std::cout << "2D scan: " << std::flush;
    for (int i = 0; i <= vcthRange; i += vcthStep) {
      for (int j = 0; j <= DLLRange; j += DLLStep) {
        Vcth = vcthStart + i;
        DLL = DLLStart + j;
        std::cout << Vcth << "/" << DLL << " " << std::flush;
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
    std::cout << std::endl;
  }

  // This is for plotting, a bit ufly though because it is still in this file. To be fixed soon!

//  TCanvas *c0 = new TCanvas();
//  TH1D *p_vcthScanLatched = h_vcthScan->ProjectionX("projectionLatchedVcthScanDDL15", 2, 2, "");
//  p_vcthScanLatched->SetTitle("Threshold scan, latched, DDL15;Threshold [Vcth];Efficiency []");
//  p_vcthScanLatched->Draw("");

//  TCanvas *c1 = new TCanvas();
//  TH1D *p_vcthScanLatchedDiff = (TH1D*) p_vcthScanLatched->Clone(); 
//  p_vcthScanLatchedDiff->SetNameTitle("diffLatchedVcthScanDDL15", "Threshold scan (differential), latched, DDL15;Threshold [Vcth];Efficiency []");
//  for ( int i = 0; i < p_vcthScanLatchedDiff->GetNbinsX()-1; i++ ) {
//    p_vcthScanLatchedDiff->SetBinContent( i, p_vcthScanLatched->GetBinContent(i) - p_vcthScanLatched->GetBinContent(i+1) );
//  }
//  p_vcthScanLatchedDiff->Draw("");
//  
//  TCanvas *c2 = new TCanvas();
//  TH1D *p_vcthScanSampled = h_vcthScan->ProjectionX("projectionSampledVcthScanDDL15", 1, 1, "");
//  p_vcthScanSampled->SetTitle("Threshold scan, sampled, DDL15;Threshold [Vcth];Efficiency []");
//  p_vcthScanSampled->Draw("");

//  TCanvas *c3 = new TCanvas();
//  TH1D *p_vcthScanSampledDiff = (TH1D*) p_vcthScanLatched->Clone();
//  p_vcthScanSampledDiff->SetNameTitle("diffSampledVcthScanDDL15", "Threshold scan (differential), sampled, DDL15;Threshold [Vcth];Efficiency []");
//  for ( int i = 0; i < p_vcthScanSampledDiff->GetNbinsX()-1; i++ ) {
//    p_vcthScanSampledDiff->SetBinContent( i, p_vcthScanSampled->GetBinContent(i) - p_vcthScanSampled->GetBinContent(i+1) );
//  }
//  p_vcthScanSampledDiff->Draw("");

//  TCanvas *c4 = new TCanvas();
//  TH1D *p_vcthScanOR = h_vcthScan->ProjectionX("projectionOrVcthScanDDL15", 3, 3, "");
//  p_vcthScanOR->SetTitle("Threshold scan, OR, DDL15;Threshold [Vcth];Efficiency []");
//  p_vcthScanOR->Draw("");

//  TCanvas *c5 = new TCanvas();
//  TH1D *p_vcthScanORDiff = (TH1D*) p_vcthScanLatched->Clone();
//  p_vcthScanORDiff->SetNameTitle("diffORVcthScanDDL15", "Threshold scan (differential), OR, DDL15;Threshold [Vcth];Efficiency []");
//  for ( int i = 0; i < p_vcthScanORDiff->GetNbinsX()-1; i++ ) {
//    p_vcthScanORDiff->SetBinContent( i, p_vcthScanOR->GetBinContent(i) - p_vcthScanOR->GetBinContent(i+1) );
//  }
//  p_vcthScanORDiff->Draw("");

//  //TF1 *landau = new TF1("f", "landau", 0., 600.);
//  TF1 *langausFitFunction = new TF1("f", Double_LanGaus, 50., 300., 4);  
//  langausFitFunction->SetNpx(500);
//  langausFitFunction->SetParameters(9., 144., 1., 6.);
//  //langausFitFunction->FixParameter(0, 9.);
//  //langausFitFunction->FixParameter(1, 144.);
//  //langausFitFunction->FixParameter(3, 6.);
//  langausFitFunction->SetParLimits(0, 0., 12.);
//  langausFitFunction->SetParLimits(1., 120., 160.);
//  langausFitFunction->SetParLimits(3, 0., 100.);
//  langausFitFunction->SetParNames("width", "MPV", "Norm", "Noise");

//  std::cout << "Latched: " << std::endl;
//  p_vcthScanLatchedDiff->Fit(langausFitFunction, "R");
//  std::cout << "Sampled: " << std::endl;
//  p_vcthScanSampledDiff->Fit(langausFitFunction, "R");
//  std::cout << "Or mode: " << std::endl;
//  p_vcthScanORDiff->Fit(langausFitFunction, "R");

//  c1->SaveAs("./images/latchedDiffFitLandauAddedNoiseWithLangaus.pdf");
//  p_vcthScanLatchedDiff->SaveAs("./rootfiles/latchedDiffFitLandauAddedNoiseWithLangaus.root");
//  c3->SaveAs("./images/sampledDiffFitLandauAddedNoiseWithLangaus.pdf");
//  p_vcthScanSampledDiff->SaveAs("./rootfiles/sampledDiffFitLandauAddedNoiseWithLangaus.root");  
//  c5->SaveAs("./images/orDiffFitLandauAddedNoiseWithLangaus.pdf");
//  p_vcthScanORDiff->SaveAs("./rootfiles/orDiffFitLandauAddedNoiseWithLangaus.root");  

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
