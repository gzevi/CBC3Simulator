//SYSTEM
//#include <iostream>

//ROOT
#include "Math/WrappedTF1.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/GSLIntegrator.h"
#include "Math/PdfFuncMathCore.h"

const double cToleranceIntegration=1e-18;
    
// some useful functions (Landau/LanGauss convolution )
double Landau(double x , double* par)
{
    double xmpv_Landau  = 0.22278298;       // Landau maximum location
    
    double mpc = par[0];
    double sigmaLandau = par[1];
    
    double fLandau = (1/(sigmaLandau))*::ROOT::Math::landau_pdf( (x - mpc)/sigmaLandau - xmpv_Landau); 
    return fLandau;
}
double fLandauDistribution(double* x , double* par )
{
    return Landau(x[0],par);
}
double LanGauss(double* yValue , double* par)
{
    double xmpv_Landau  = 0.22278298;       // Landau maximum location
    
    double mpc = par[0];
    double sigmaLandau = par[1];
    double sigmaGaus = par[2];
    double x = par[3];

    double fLandau = (1/(sigmaLandau))*::ROOT::Math::landau_pdf( (yValue[0] - mpc)/sigmaLandau - xmpv_Landau); 
    double fGaus = (1.0/std::sqrt(2*sigmaGaus*sigmaGaus*TMath::Pi()))*TMath::Exp( -0.5*std::pow((x-yValue[0])/sigmaGaus,2.0) );
    return fLandau*fGaus;

}
double ConvoluteLanGaus(double x , double* par)
{
    double xmpv_Landau  = 0.22278298;       // Landau maximum location
    
    double mpc = par[0];
    double sigmaLandau = par[1];
    double sigmaGaus = par[2];
    double rIntegral=par[3];
    //int nIntegrationPoints=par[4];
    
    // numerically integrating 
    double xMin = x - rIntegral*sigmaGaus;//-1*rIntegral*mpc;//x - 100*sigmaGaus;
    double xMax = x + rIntegral*sigmaGaus;// 1*rIntegral*mpc;//x + 100*sigmaGaus;


    TF1 f("lanGauss",LanGauss , xMin , xMax,4);
    double params[4]={mpc,sigmaLandau,sigmaGaus,x};
    f.SetParameters(params);
    //return simpsonsRule(&f, xMin , xMax , nIntegrationPoints);

    ROOT::Math::WrappedTF1 wf1(f);
    // //using ROOT's built in integrator 
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
    // Set parameters of the integration
    ig.SetFunction(wf1);
    ig.SetRelTolerance(cToleranceIntegration);
    return ig.Integral(xMin,xMax);
    
}
double LanGausConvolution(double* x , double* par)
{
    return ConvoluteLanGaus(x[0],par);
}
double fLanGausIntegral(double x , double* par)
{
    double mpc = par[0];
    double sigmaLandau = par[1];
    double sigmaGaus = par[2];
    double rIntegral=par[3];
    //int nIntegrationPoints=par[4];

    double xMin = x;
    double xMax = 1*rIntegral*mpc;
    double convolutionParameters[5]={mpc,sigmaLandau,sigmaGaus,5};//{mpc,sigmaLandau,sigmaGaus,10,500};
    TF1 f("lanGaussIntegral",LanGausConvolution , xMin , xMax,4);
    f.SetParameters(convolutionParameters);
    ROOT::Math::WrappedTF1 wf1(f);
    
    // // Create the Integrator
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
    // Set parameters of the integration
    ig.SetFunction(wf1);
    ig.SetRelTolerance(cToleranceIntegration);
    return ig.Integral(xMin,xMax);
}
double mipEfficiency(double* x , double* par)
{
    double mpc = par[0];
    double sigmaLandau = par[1];
    double sigmaGaus = par[2];
    double rIntegral=par[3];
    double cEff=par[4];
    //double nIntegrationPoints=par[4];
    
    double parametersLanGaus[4]={mpc,sigmaLandau,sigmaGaus,rIntegral};
    return cEff*fLanGausIntegral(x[0],parametersLanGaus);
}
double detResponse(double* x , double* par)
{
    double mpc = par[0];
    double sigmaLandau = par[1];
    double sigmaGaus = par[2];
    double rIntegral=5;
    double cNorm=par[3];

    double parametersLanGaus[4]={mpc,sigmaLandau,sigmaGaus,rIntegral};
    return cNorm*ConvoluteLanGaus(x[0],parametersLanGaus);
    
}