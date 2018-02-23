#include <iostream>

#include "Utils.h"

const char *getFileName (const char * cSimFileName) {
  return Form("cbc3Simulation_%s.root", cSimFileName);
}

void analysis () {

  TFile *cFile = TFile::Open( getFileName("VcthOnly_10000events") );
  TH2* cVcthHist = 0;

  cFile->GetObject("ThresholdScan_DLL15", cVcthHist);
  
  cVcthHist->Draw("colz");
  
  // This is for plotting the 1D histos:

  TCanvas *c0 = new TCanvas();
  TH1D *p_vcthScanLatched = cVcthHist->ProjectionX("projectionLatchedVcthScanDDL15", 2, 2, "");
  p_vcthScanLatched->SetTitle("Threshold scan, latched, DDL15;Threshold [Vcth];Efficiency []");
  p_vcthScanLatched->Draw("");

  TCanvas *c1 = new TCanvas();
  TH1D *p_vcthScanLatchedDiff = (TH1D*) p_vcthScanLatched->Clone(); 
  p_vcthScanLatchedDiff->SetNameTitle("diffLatchedVcthScanDDL15", "Threshold scan (differential), latched, DDL15;Threshold [Vcth];Efficiency []");
  for ( int i = 0; i < p_vcthScanLatchedDiff->GetNbinsX()-1; i++ ) {
    p_vcthScanLatchedDiff->SetBinContent( i, p_vcthScanLatched->GetBinContent(i) - p_vcthScanLatched->GetBinContent(i+1) );
  }
  p_vcthScanLatchedDiff->Draw("");
  
  TCanvas *c2 = new TCanvas();
  TH1D *p_vcthScanSampled = cVcthHist->ProjectionX("projectionSampledVcthScanDDL15", 1, 1, "");
  p_vcthScanSampled->SetTitle("Threshold scan, sampled, DDL15;Threshold [Vcth];Efficiency []");
  p_vcthScanSampled->Draw("");

  TCanvas *c3 = new TCanvas();
  TH1D *p_vcthScanSampledDiff = (TH1D*) p_vcthScanLatched->Clone();
  p_vcthScanSampledDiff->SetNameTitle("diffSampledVcthScanDDL15", "Threshold scan (differential), sampled, DDL15;Threshold [Vcth];Efficiency []");
  for ( int i = 0; i < p_vcthScanSampledDiff->GetNbinsX()-1; i++ ) {
    p_vcthScanSampledDiff->SetBinContent( i, p_vcthScanSampled->GetBinContent(i) - p_vcthScanSampled->GetBinContent(i+1) );
  }
  p_vcthScanSampledDiff->Draw("");

  TCanvas *c4 = new TCanvas();
  TH1D *p_vcthScanOR = cVcthHist->ProjectionX("projectionOrVcthScanDDL15", 3, 3, "");
  p_vcthScanOR->SetTitle("Threshold scan, OR, DDL15;Threshold [Vcth];Efficiency []");
  p_vcthScanOR->Draw("");

  TCanvas *c5 = new TCanvas();
  TH1D *p_vcthScanORDiff = (TH1D*) p_vcthScanLatched->Clone();
  p_vcthScanORDiff->SetNameTitle("diffORVcthScanDDL15", "Threshold scan (differential), OR, DDL15;Threshold [Vcth];Efficiency []");
  for ( int i = 0; i < p_vcthScanORDiff->GetNbinsX()-1; i++ ) {
    p_vcthScanORDiff->SetBinContent( i, p_vcthScanOR->GetBinContent(i) - p_vcthScanOR->GetBinContent(i+1) );
  }
  p_vcthScanORDiff->Draw("");

  //TF1 *landau = new TF1("f", "landau", 0., 600.);
  TF1 *langausFitFunction = new TF1("f", Double_LanGaus, 50., 300., 4);  
  langausFitFunction->SetNpx(500);
  langausFitFunction->SetParameters(9., 144., 1., 6.);
  //langausFitFunction->FixParameter(0, 9.);
  //langausFitFunction->FixParameter(1, 144.);
  //langausFitFunction->FixParameter(3, 6.);
  langausFitFunction->SetParLimits(0, 0., 12.);
  langausFitFunction->SetParLimits(1., 120., 160.);
  langausFitFunction->SetParLimits(3, 0., 100.);
  langausFitFunction->SetParNames("width", "MPV", "Norm", "Noise");

  std::cout << "Latched: " << std::endl;
  p_vcthScanLatchedDiff->Fit(langausFitFunction, "R");
  std::cout << "Sampled: " << std::endl;
  p_vcthScanSampledDiff->Fit(langausFitFunction, "R");
  std::cout << "Or mode: " << std::endl;
  p_vcthScanORDiff->Fit(langausFitFunction, "R");

  c1->SaveAs("./images/latchedDiffFitLandauAddedNoiseWithLangaus.pdf");
  p_vcthScanLatchedDiff->SaveAs("./rootfiles/latchedDiffFitLandauAddedNoiseWithLangaus.root");
  c3->SaveAs("./images/sampledDiffFitLandauAddedNoiseWithLangaus.pdf");
  p_vcthScanSampledDiff->SaveAs("./rootfiles/sampledDiffFitLandauAddedNoiseWithLangaus.root");  
  c5->SaveAs("./images/orDiffFitLandauAddedNoiseWithLangaus.pdf");
  p_vcthScanORDiff->SaveAs("./rootfiles/orDiffFitLandauAddedNoiseWithLangaus.root");  


}
