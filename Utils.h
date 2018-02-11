#ifndef UTILS_HH
#define UTILS_HH

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

void labelAxis ( TAxis * ); 					// Axis labels

double Double_LanGaus ( double *, double * ); 	// Returns convolution of landau with gaussian, to be used in TF1 
double Function_RC ( double *, double * );		// Returns pulse shape

TF1 * Function_Landau ( ); 						// Landau distribution
TF1 * Function_LanGaus ( ); 					// Langaus distribution 
TF1 * Function_PulseShape( float );				// Pulseshape function

void labelAxis(TAxis * a) {
  a->SetBinLabel(1, "Sampled");
  a->SetBinLabel(2, "Latched");
  a->SetBinLabel(3, "OR");
  a->SetBinLabel(4, "Sampled_HIP");
  a->SetBinLabel(5, "OR_HIP");
}

TF1 * Function_Landau() {
  // Currently using a Landau in Vcth units, but could use fC and then scale the PreampShaper output
  TF1* f = new TF1("f", "landau", 0, 300);
  f->SetParameters(1, 144, 9);
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
  // Sigma of the gaussian is 6 (CBC3 noise is 2 Vcth units??? To be checked!!!)
  f->SetParameters(9, 144, 1, 6);
  f->SetParNames("Width","MP","Area","GSigma");
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

TF1 * Function_PulseShape(float charge) {
  // TF1* f = new TF1("f1", "gaus", 0, Nck*25);
  // f->SetParameters(charge, 20, 10);

  TF1* f = new TF1("f1", Function_RC, -25, (Nck-1)*25, 3);
  f->SetParameters(3, 4, 1);

  // Normalize so that peak matches "charge". This way if we draw "140" from the Landau distribution,
  // we will have a PulseShape that has a maximum value of 140
  double max = f->GetMaximum(0, 1000);
  f->SetParameters(3, 4, charge/max);

  return f;
}

#endif