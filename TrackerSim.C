#include <iostream>
#include <chrono>
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
#include "TProfile2D.h"
#include "TLine.h"

#include "LanGaus.cc"
#include "CBC3sim.C"
const int nCBCs_perModule = 16;
const int nStrips_perCBC = 127;

//CBC3 post ampl setttings 
const double cConversion_CBC3 = 133e-3; // in kElectrons
const double nStages_CBC3=3; // # amplification stages
const double tau_CBC3=4.; // in ns 
const double tRiseTime_CBC3=2; // in ns
// and some parameters for the default pulse shape (more realistic)
const double tau_CBC3improv =12.4;
const double r_CBC3improv=-12.0;
const double theta_CBC3improv=0.79;//TMath::PiOver2();
const double nTerms=3;


//some timing settings for the CBC3
const double cDeadTime_LatchedMode=2.7; // default is 2.7 in ns
const double cDeadTime_SampledMode=0.3; // default is 0.3 in ns

//some paramters that define "acceptable" efficiencies and fake rates
const double cMinAcceptableEfficiency = 0.9; 
const double cMaxAcceptableFakeRate = 0.0; 

typedef std::pair<bool,bool> HitDetectStates;
const bool ChargeSharing = false;
const double cSizeChargeSharingRegion=0.05*0.5; // default 0.01*0.5 ; 

// pretty color palette?!
void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

// try and figure out best settings for comparator
// from measured pulse shape using on-chip test pulse 
void fitPulseShape(TString pFileName="./examplePulseShape_MIP_CBC3.root", double pPedestal_DACunits=598, double pUnc_Perc=5)
{
	TString cOut;
	TString cHistName;
	TH1D* hMeasuredPulse=0;
	TFile* f = new TFile(pFileName,"READ");
	f->cd();
	cHistName.Form("pulseShape_MIP");
	f->GetObject (cHistName.Data() , hMeasuredPulse);
	hMeasuredPulse->SetDirectory(0);  
	f->Close();

	//clone histogram to remove pedestal 
	TH1D* hPulseShape_pSub = (TH1D*)hMeasuredPulse->Clone("pSubPulse");
	hPulseShape_pSub->Reset();
	double cMax = hMeasuredPulse->GetMaximum();
	for(int i = 1 ; i <= hPulseShape_pSub->GetNbinsX() ;i++)
	{
		double cBinContent = hMeasuredPulse->GetBinContent(i);
		if( cBinContent != 0)
		{
			hPulseShape_pSub->SetBinContent(i, (pPedestal_DACunits - cBinContent) );
			hPulseShape_pSub->SetBinError(i, (pUnc_Perc*(pPedestal_DACunits - cBinContent)/100.) );
		}
	}

	TString cFuncName;
	//first try and fit CRC
	cFuncName.Form("fPulse_CRC"); 
	double cNorm=cMax;
	double parsAmp[]={nStages_CBC3, tau_CBC3 , cNorm , 10.0 };
	TF1* f0 = new TF1(cFuncName.Data(), Function_RC, 0, 75.0, 4);
  	f0->SetParameters(parsAmp);
  	f0->SetLineColor(kBlue);
	f0->FixParameter(3,10.0);
	//hPulseShape_pSub->Fit(f0,"nr","",0,75);
	//f0->Draw("same");
	
	cFuncName.Form("fPostAmp_CBC3");
	double pulseShapePars[]={10.0 ,tau_CBC3improv,r_CBC3improv,theta_CBC3improv,nTerms, cMax};
	TF1 *f1 = new TF1(cFuncName.Data() , fPulseCBC3,  0, 75.0 , 6 );
  	f1->SetParameters(pulseShapePars);
  	f1->SetParNames("xOffset","tau","r","theta", "nTerms","k");
  	f1->FixParameter(4,nTerms);
  	//f1->FixParameter(0,10);
  	f1->SetParLimits(3,0,TMath::PiOver2());
	TCanvas* c = new TCanvas("c","c",350,350);
	c->cd(1);
	hPulseShape_pSub->Fit(f1,"nr+","",0,60);
	hPulseShape_pSub->SetStats(0);
	hPulseShape_pSub->DrawCopy("eHisto");
	f1->Draw("same");

	std::vector<double> aValues;aValues.clear();
	std::vector<double> bValues;bValues.clear();
	for( int i = 0 ; i < 2 ; i++ )
	{	
		double tauValue = f1->GetParameter(f1->GetParNumber("tau")) + std::pow(-1, (double)i)*f1->GetParError(f1->GetParNumber("tau"));
		for( int j = 0 ; j < 2 ; j++ )
		{
			double rValue = f1->GetParameter(f1->GetParNumber("r")) + std::pow(-1, (double)j)*f1->GetParError(f1->GetParNumber("r"));
			for( int k = 0 ; k < 2 ; k++ )
			{
				double thetaValue = f1->GetParameter(f1->GetParNumber("theta")) + std::pow(-1, (double)k)*f1->GetParError(f1->GetParNumber("theta"));
				aValues.push_back( tauValue + rValue*std::cos(thetaValue) );
				bValues.push_back( tauValue + rValue*std::sin(thetaValue) );
				cOut.Form("tau = %.2f , theta = %.2f , r = %.2f : a = %.1f , b = %.1f\n", tauValue, thetaValue, rValue , aValues[aValues.size()-1] , bValues[bValues.size()-1]);
				std::cout << cOut.Data();			
			}
		}
	}
	
	double aAvgValue=0;
	double bAvgValue=0;
	for(int i = 0 ; i < (int)aValues.size(); i++)	
	{
		aAvgValue += aValues[i]/aValues.size();
		bAvgValue += bValues[i]/aValues.size();
	}
	double aStdDev=0;
	double bStdDev=0;
	for(int i = 0 ; i < (int)aValues.size(); i++)	
	{
		aStdDev += std::pow(aValues[i]-aAvgValue,2.);
		bStdDev += std::pow(bValues[i]-bAvgValue,2.);
	}
	aStdDev = std::sqrt(aStdDev)/(aValues.size()-1);
	bStdDev = std::sqrt(bStdDev)/(bValues.size()-1);

	cOut.Form("a = %.2f ± %.2f\t b = %.2f ± %.2f\n" , aAvgValue , aStdDev,  bAvgValue , bStdDev );
	std::cout << cOut.Data();
}
// quick check of charge sharing 
void testChargeSharing( double pLandauMPV_kElectrons=25, double pLandauWidth_kElectrons=2, double pNoise_kElectrons=0.8, double pEnergyMin_kElectrons = 0.0 , double pEnergyMax_kElectrons = 100 )
{
	double cThreshold=10;
	double parsLandau[]={1.0, pLandauMPV_kElectrons,pLandauWidth_kElectrons };
	TF1* fLandau = Function_Landau(parsLandau, "fLandau",pEnergyMin_kElectrons, pEnergyMax_kElectrons);
	
	double parsChargeSharing[]={-0.5,0.5,cSizeChargeSharingRegion};
	TF1* fChargeSharing = Function_Comparator(parsChargeSharing, "fChargeSharing",-0.5, 0.5);
	fChargeSharing->Draw();

	// std::cout << fChargeSharing->Eval(-0.5) << " \t\t" << fChargeSharing->Eval(0.5) << std::endl;
	// TH1D* hEdeposited_DAC =  new TH1D("hEdeposited_DAC","Energy Deposited; Energy [DAC units]; Number of Entries", (int)(pEnergyMax_kElectrons/(0.1*cConversion_CBC3) - pEnergyMin_kElectrons/(cConversion_CBC3*0.1)) , pEnergyMin_kElectrons/cConversion_CBC3 , pEnergyMax_kElectrons/cConversion_CBC3 );
	// TH1D* hEdeposited2_DAC =  new TH1D("hEdeposited2_DAC","Energy Deposited; Energy [DAC units]; Number of Entries", (int)(pEnergyMax_kElectrons/(0.1*cConversion_CBC3) - pEnergyMin_kElectrons/(cConversion_CBC3*0.1)) , pEnergyMin_kElectrons/cConversion_CBC3 , pEnergyMax_kElectrons/cConversion_CBC3 );
	// TH2D* hEdepositedXpos_DAC =  new TH2D("hEdepositedXpos_DAC","Energy Deposited; Energy [DAC units]; Hit Position", (int)(pEnergyMax_kElectrons/(0.1*cConversion_CBC3) - pEnergyMin_kElectrons/(cConversion_CBC3*0.1)) , pEnergyMin_kElectrons/cConversion_CBC3 , pEnergyMax_kElectrons/cConversion_CBC3 , 1.0/(0.1*cSizeChargeSharingRegion), -0.5 , 0.5  );
	// TH1D* hHitPosition =  new TH1D("hHitPosition","Position Hits", 1.0/(0.1*cSizeChargeSharingRegion), -0.5 , 0.5 );
	// TH1D* hClusterSize1 =  new TH1D("hClustersSize1","Position 1 strip clusters", 0.50/(0.1*cSizeChargeSharingRegion),  0 , 0.5 );
	// TH1D* hClusterSize2 =  new TH1D("hClustersSize2","Position 2 strip clusters", 0.50/(0.1*cSizeChargeSharingRegion),  0 , 0.5 );

	// TRandom3 myDice;
	// for( int iEvent=0; iEvent < 10e6 ; iEvent++)
	// {
	// 	int cClusterSize=0;
	// 	double eDep = fLandau->GetRandom();
	// 	double eNoise = myDice.Gaus(0,pNoise_kElectrons);

	// 	// now assume hit is randomly distributed across a strip 
	// 	// so implant is @ 0 , -0.5 is strip n-1 , 0.5 is strip n 
	// 	double xPosHit = myDice.Uniform(-0.5,0.5);
	// 	hHitPosition->Fill(xPosHit);
	// 	// fraction of charge deposited in this strip 
	// 	double cFractionDepCharge = (cChargeSharing) ? fChargeSharing->Eval(xPosHit) : 1 ;
	// 	double eDep_Total = eDep*cFractionDepCharge + eNoise;
	// 	hEdeposited_DAC->Fill( (eDep_Total/cConversion_CBC3) );
	// 	// want to see what energy deposited looks like for hits close to the midpoint betweem two implants
	// 	hEdepositedXpos_DAC->Fill( (eDep_Total/cConversion_CBC3) , xPosHit);
	// 	cClusterSize += ((eDep_Total/cConversion_CBC3) > cThreshold ) ? 1 : 0 ;
	// 	eNoise = myDice.Gaus(0,pNoise_kElectrons);
	// 	double eDep_TotalSecond = eDep*(1-cFractionDepCharge) + eNoise;
	// 	cClusterSize += ((eDep_TotalSecond/cConversion_CBC3) > cThreshold ) ? 1 : 0 ;
	// 	if(iEvent%100000==0)
	// 		std::cout << "Event " << iEvent << " ----- " <<  cClusterSize << " strip cluster\n";

	// 	if(cClusterSize==2)
	// 		hClusterSize2->Fill(std::fabs(xPosHit));
	// 	if(cClusterSize==1)
	// 		hClusterSize1->Fill(std::fabs(xPosHit));

	// 	if( (0.5-std::fabs(xPosHit))<=0.25*cSizeChargeSharingRegion )
	// 	{
	// 		hEdeposited2_DAC->Fill( (eDep_TotalSecond/cConversion_CBC3) );
	// 	}
				
 			

	// 	// fill histograms for channel n-1 if fraction of energy deposited != 0 
	// 	if( xPosHit < 0 )
	// 	{
	// 		hEdeposited_DAC->Fill( (eDep_TotalSecond/cConversion_CBC3) );
	// 		if( (0.5-std::fabs(xPosHit))<=0.25*cSizeChargeSharingRegion )
	// 		{
	// 			hEdeposited2_DAC->Fill( (eDep_TotalSecond/cConversion_CBC3) );
	// 		}
	// 	}
	// 	// fill histograms for channel n+1 if fraction of energy deposited != 0 
	// 	else if(xPosHit > 0 )
	// 	{
	// 		hEdeposited_DAC->Fill( (eDep_TotalSecond/cConversion_CBC3) );
	// 		if( (0.5-std::fabs(xPosHit))<=0.1 )
	// 		{
	// 			hEdeposited2_DAC->Fill( (eDep_TotalSecond/cConversion_CBC3) );
	// 		}
	// 	}
	// }

	// TCanvas* c = new TCanvas("c","c",350*3,350);
	// c->Divide(3,1);
	// c->cd(1);
	// hEdepositedXpos_DAC->SetStats(0);
	// hEdepositedXpos_DAC->DrawCopy("colz");
	// c->cd(2);
	// hHitPosition->DrawCopy("eHisto");

	// c->cd(3);
	// hClusterSize1->SetStats(0);
	// hClusterSize1->SetLineColor(kBlue);
	// hClusterSize1->DrawCopy("eHisto");
	// hClusterSize2->SetLineColor(kRed);
	// hClusterSize2->DrawCopy("eHistoSAME");



}

//populate our toy tracker with hits 
//our toy tracker has pNmodules 
//runs for pNbunchCrossings
//and had an occupancy of pOccupancy
void populateHits(TString pFileName="./ToyMC_test.root", int pNmodules=10, int pNbunchCrossings=5, double pOccupancy= 1e-2, double pLandauMPV_kElectrons=25, double pLandauWidth_kElectrons=2, double pPedestal_kElectrons=1, double pNoise_kElectrons=0.8, double cChargeSharing=cSizeChargeSharingRegion, double pEnergyMin_kElectrons = 0.0 , double pEnergyMax_kElectrons = 100 )
{
	TString cOut;
	int nModules = pNmodules;
	int nTotalNumberOfChannels = nModules*nCBCs_perModule*nStrips_perCBC;
	double parsLandau[]={1.0, pLandauMPV_kElectrons,pLandauWidth_kElectrons };
	TF1* fLandau = Function_Landau(parsLandau, "fLandauDistribution",pEnergyMin_kElectrons, pEnergyMax_kElectrons);
	
	double parsChargeSharing[]={-0.5,0.5,cChargeSharing};
	TF1* fChargeSharing = Function_Comparator(parsChargeSharing, "fChargeSharing",-0.5, 0.5);

	double parsNoise[]={1.0, pPedestal_kElectrons, pNoise_kElectrons};
	TF1* fNoise = new TF1("fNoise","[0]*TMath::Gaus(x,[1],[2],true)", pEnergyMin_kElectrons , pEnergyMax_kElectrons);
	fNoise->SetParameters(parsNoise);

	int nBunchCrossings=pNbunchCrossings;
	TH2D* hHits = new TH2D("hHits", "Channels with a particle strike ; Bunch Crossing ID [a.u.]; Channel Number", nBunchCrossings , 0 , nBunchCrossings, nTotalNumberOfChannels , 0 , nTotalNumberOfChannels ); 
	TH2D* hClusters = new TH2D("hClusters", "Channels with a particle strike ; Bunch Crossing ID [a.u.]; Channel Number", nBunchCrossings , 0 , nBunchCrossings, nTotalNumberOfChannels*2.0 , 0 , nTotalNumberOfChannels ); 
	TH2D* hDepositedCharge_wNoise = new TH2D("hDepCharge_wNoise", "Deposited Charge; Bunch Crossing ID [a.u.]; Channel Number", nBunchCrossings , 0 , nBunchCrossings, nTotalNumberOfChannels , 0 , nTotalNumberOfChannels ); 
	TH2D* hPostAmplifier_wNoise = new TH2D("hPostAmplifier_wNoise", "Post Amplifier Output; Time [ns]; Channel Number", nBunchCrossings*25.0 , 0 , nBunchCrossings*25.0, nTotalNumberOfChannels , 0 , nTotalNumberOfChannels ); 
	TH1D* hEdeposited_DAC =  new TH1D("hEdeposited_DAC","Energy Deposited; Energy [DAC units]; Number of Entries", (int)(pEnergyMax_kElectrons/(1*cConversion_CBC3) - pEnergyMin_kElectrons/(cConversion_CBC3*0.1)) , pEnergyMin_kElectrons/cConversion_CBC3 , pEnergyMax_kElectrons/cConversion_CBC3 );
	TH1D* hEdeposited_SingleStripClusters_DAC =  new TH1D("hEdeposited_SingleStripClusters_DAC","Energy Deposited; Energy [DAC units]; Number of Entries", (int)(pEnergyMax_kElectrons/(1*cConversion_CBC3) - pEnergyMin_kElectrons/(cConversion_CBC3*0.1)) , pEnergyMin_kElectrons/cConversion_CBC3 , pEnergyMax_kElectrons/cConversion_CBC3 );
	TH1D* hEdeposited2_DAC =  new TH1D("hEdeposited2_DAC","Energy Deposited; Energy [DAC units]; Number of Entries", (int)(pEnergyMax_kElectrons/(1*cConversion_CBC3) - pEnergyMin_kElectrons/(cConversion_CBC3*0.1)) , pEnergyMin_kElectrons/cConversion_CBC3 , pEnergyMax_kElectrons/cConversion_CBC3 );
	TH2D* hDoubleHits =  new TH2D("hDoubleHits","Hits in consecutive BXs; Bunch Crossing ID [a.u.]; Channel Number", nBunchCrossings , 0 , nBunchCrossings, nTotalNumberOfChannels , 0 , nTotalNumberOfChannels ); 
	TH1D* hHitPosition =  new TH1D("hHitPosition","Position Hits", 1.0/(0.1*cSizeChargeSharingRegion), -0.5 , 0.5 );
	TH1D* hNoise_DAC =  new TH1D("hNoise_DACunits","Noise; Noise [DAC units]; Number of Entries", (int)(pEnergyMax_kElectrons/(0.1*cConversion_CBC3) - pEnergyMin_kElectrons/(cConversion_CBC3*0.1)) , pEnergyMin_kElectrons/cConversion_CBC3 , pEnergyMax_kElectrons/cConversion_CBC3 );
	TH1D* hClusterSize = new TH1D("hClusterSize","Cluster Size; Cluster Size; Number of Entries" , 100 , 0 , 100);
	TRandom3 myDice;
	//  populate channels in modules with hits in each BX 
	std::vector<int> cHits;cHits.clear();
	std::vector<int> cHitsPrevBx; cHitsPrevBx.clear();
	
 	for (int iBx =0 ; iBx < nBunchCrossings ; iBx++ ) 
 	{
 		int nChannels_Fired = myDice.Poisson(nTotalNumberOfChannels*pOccupancy);
 		cOut.Form("%d/%d channels fired in BX%d..... \n", nChannels_Fired, nTotalNumberOfChannels, iBx);
 		std::cout << cOut.Data();
 		cHits.clear();
 		for( int iChannel = 0; iChannel < nChannels_Fired ; iChannel++ )
 		{
 			// make sure I'm not choosing the same channel twice 
 			int channelNumber = myDice.Integer(nTotalNumberOfChannels);//(int)myDice.Uniform(0, (double)nModules*nCBCs_perModule*nStrips_perCBC );
 			auto it = find (cHits.begin(), cHits.end(), channelNumber);
 			do
 			{
 				channelNumber = myDice.Integer(nTotalNumberOfChannels);//(int)myDice.Uniform(0, (double)nModules*nCBCs_perModule*nStrips_perCBC );
 				it = find (cHits.begin(), cHits.end(), channelNumber);
 			}while( it != cHits.end() );
 			cHits.push_back(channelNumber);
				
 			// check if this channel also fired in the previous bx 
 			if( iBx > 0 ) // only makes sense if BX > 0 
 			{
	 			it = find (cHitsPrevBx.begin(), cHitsPrevBx.end(), channelNumber);
	 			if( it != cHitsPrevBx.end() )
	 			{
	 				hDoubleHits->Fill(iBx-1,channelNumber);
	 			}
 			}	

 			double eDep = fLandau->GetRandom();
	 		double eNoise = fNoise->GetRandom();
 			hNoise_DAC->Fill( eNoise/cConversion_CBC3);
 			int cClusterSize=1;

 			// now assume hit is randomly distributed across a strip 
 			// so implant is @ 0 , -0.5 is strip n-1 , 0.5 is strip n 
 			double xPosHit = myDice.Uniform(-0.5,0.5);
 			int nSecondChannel = (xPosHit < 0) ? channelNumber -1 : channelNumber + 1;
 			hHitPosition->Fill(xPosHit);
 			// fraction of charge deposited in this strip 
 			double cFractionDepCharge = (ChargeSharing) ? fChargeSharing->Eval(xPosHit) : 1 ;
 			double eDep_Total = eDep*cFractionDepCharge + eNoise;
 			// noise in second strip
 			eNoise = fNoise->GetRandom();//std::fabs(myDice.Gaus(0,pNoise_kElectrons));
 			double eDep_TotalSecond = eDep*(1-cFractionDepCharge) + eNoise;
 			
 			
 			// fill histograms for this channel  
			hHits->Fill(iBx,channelNumber);
			hDepositedCharge_wNoise->Fill(iBx, channelNumber, eDep_Total/cConversion_CBC3  );
			hEdeposited_DAC->Fill( (eDep_Total/cConversion_CBC3) );
			double cClusterCenter=channelNumber;
			if( nSecondChannel >=0 && nSecondChannel < nTotalNumberOfChannels && (1-cFractionDepCharge)>0.01 )
			{
				cClusterSize += 1; 
	 			cClusterCenter = channelNumber + (nSecondChannel-channelNumber)/2.0;
	 			hNoise_DAC->Fill( eNoise/cConversion_CBC3);
 				//hHits->Fill(iBx,nSecondChannel);
				hDepositedCharge_wNoise->Fill(iBx, nSecondChannel, eDep_TotalSecond/cConversion_CBC3  );
				hEdeposited_DAC->Fill( (eDep_TotalSecond/cConversion_CBC3) );
			}
			hClusterSize->Fill(cClusterSize);
			hClusters->Fill(iBx,cClusterCenter);
			if( cClusterSize == 1 )
				hEdeposited_SingleStripClusters_DAC->Fill((eDep_Total/cConversion_CBC3) );

 			// want to see what energy deposited looks like for hits close to the midpoint betweem two implants
			if( std::fabs(0.5-std::fabs(xPosHit))  <= 0.25*cSizeChargeSharingRegion ) 
 				hEdeposited2_DAC->Fill( (eDep_TotalSecond/cConversion_CBC3) );
				
 			//double parsAmp[]={nStages_CBC3, tau_CBC3 , 1.0 , iBx*25.0 };
			double parsAmp[]={iBx*25.0 ,tau_CBC3improv,r_CBC3improv,theta_CBC3improv,nTerms};

			TString cFuncName;
	  		cFuncName.Form("fPostAmp_CBC3_Ch%d_Bx%d",channelNumber, iBx);
	  		TF1* fSignal = Function_PulseShapeCBC3((eDep_Total/cConversion_CBC3), parsAmp, cFuncName.Data() ,0.0 , nBunchCrossings*25);
	  		//TF1* fSignal = Function_PulseShape((eDep_Total/cConversion_CBC3), parsAmp, cFuncName.Data() ,0.0 , nBunchCrossings*25);
	  		cFuncName.Form("fPostAmp_CBC3_Ch%d_Bx%d",nSecondChannel, iBx);
	  		TF1* fSignalSecond = 0;
	  		if( ChargeSharing ) 
	  			fSignalSecond = Function_PulseShapeCBC3((eDep_TotalSecond/cConversion_CBC3), parsAmp, cFuncName.Data() ,0.0 , nBunchCrossings*25);
	  		//TF1* fSignalSecond = Function_PulseShape((eDep_TotalSecond/cConversion_CBC3), parsAmp, cFuncName.Data() ,0.0 , nBunchCrossings*25);
	  		for(int iTimeBin = 0 ; iTimeBin < hPostAmplifier_wNoise->GetYaxis()->GetNbins() ;iTimeBin++)
	 	 	{
	 		 	double xValue = hPostAmplifier_wNoise->GetXaxis()->GetBinCenter(iTimeBin) - 0.5*hPostAmplifier_wNoise->GetXaxis()->GetBinWidth(iTimeBin) ;
	 		 	double yValue_thisBX =  fSignal->Eval(xValue);
	 		 	hPostAmplifier_wNoise->Fill(xValue, channelNumber, yValue_thisBX );
	 			if( nSecondChannel >=0 && nSecondChannel < nTotalNumberOfChannels && (1-cFractionDepCharge)>0.01  )
				{
					if( nSecondChannel == 20 && iTimeBin == 0 )
					{
						std::cout << "Channel " << nSecondChannel << "; fDep = " << cFractionDepCharge <<  " noise " << eNoise << " Edep = " << eDep_TotalSecond/cConversion_CBC3 << " value @ 0 = " << fSignalSecond->Eval(xValue) << " - compared to " << fSignal->Eval(xValue) <<  std::endl;
					} 
					hPostAmplifier_wNoise->Fill(xValue, nSecondChannel, fSignalSecond->Eval(xValue) );
				}
	 		}

	 		cOut.Form("..... finished populating hits .... [%d] Bx%d\n", iChannel, iBx );
	 		if( iChannel%250 == 0 )
	 			std::cout << cOut.Data();


		}
		cHitsPrevBx.clear();
		for( int iChannel=0;iChannel<(int)cHits.size();iChannel++)
		{
			cHitsPrevBx.push_back(cHits[iChannel]);
		}
	}
	
	TFile* f = new TFile(pFileName,"RECREATE");
	f->cd();
	hHits->Write ("hHits_perCh_perBx", TObject::kOverwrite);
	hClusters->Write ("hClusters_perCh_perBx", TObject::kOverwrite);
	hDepositedCharge_wNoise->Write ("hDepositedCharge_perCh_perBx", TObject::kOverwrite);
	hPostAmplifier_wNoise->Write ("hPostAmplifierOutput_perCh_perBx", TObject::kOverwrite);
	hEdeposited_DAC->Write("hEdeposited_DACunits", TObject::kOverwrite);
	hEdeposited2_DAC->Write("hEdeposited2_DACunits", TObject::kOverwrite);
	hEdeposited_SingleStripClusters_DAC->Write("hEdeposited_SingleStripClusters_DACunits", TObject::kOverwrite);
	hDoubleHits->Write("hDoubleHits", TObject::kOverwrite);
	hHitPosition->Write("hHitPositions", TObject::kOverwrite);
	hNoise_DAC->Write("hNoise_DACunits", TObject::kOverwrite);
	hClusterSize->Write("hClusterSize", TObject::kOverwrite);
	f->Close();
	
}



//for each hit figure out what the comparator output looks like 
//for thresholds between pThreshold_DAC_min -- pThreshold_DAC_max
void populateComparator(TString pFileName="./ToyMC_test.root", int pThreshold_DAC_min = 0 , int pThreshold_DAC_max = 400, int pThreshold_DAC_step = 5)
{
	int nCBCs_perModule = 16;
	int nStrips_perCBC = 127;

	TString cCanvasName, cCanvasTitle;
	TString cHistName, cHistTitle;
	TString cOut;
	
	TH2D* hHits=0;
	TH2D* hDepositedCharge_wNoise=0;
	TH2D* hPostAmplifier_wNoise=0;
	TH1D* hEdeposited_DAC=0;
	TH2D* hDoubleHits=0;

	TFile* cFile = new TFile(pFileName,"READ");
	cHistName.Form("hHits_perCh_perBx");
	cFile->GetObject (cHistName.Data() , hHits);
	hHits->SetDirectory(0);  

	cHistName.Form("hDepositedCharge_perCh_perBx");
	cFile->GetObject (cHistName.Data() , hDepositedCharge_wNoise);
	hDepositedCharge_wNoise->SetDirectory(0);  

	cHistName.Form("hPostAmplifierOutput_perCh_perBx");
	cFile->GetObject (cHistName.Data() , hPostAmplifier_wNoise);
	hPostAmplifier_wNoise->SetDirectory(0);  

	cHistName.Form("hEdeposited_DACunits");
	cFile->GetObject (cHistName.Data() , hEdeposited_DAC);
	hEdeposited_DAC->SetDirectory(0);  

	cHistName.Form("hDoubleHits");
	cFile->GetObject (cHistName.Data() , hDoubleHits);
	hDoubleHits->SetDirectory(0);  
	cFile->Close();
	
	int nTotalNumberOfChannels = hHits->GetYaxis()->GetNbins();
	int nBunchCrossings = hHits->GetXaxis()->GetNbins();
	for( int pThreshold_DAC = pThreshold_DAC_min ; pThreshold_DAC < pThreshold_DAC_max  ; pThreshold_DAC +=pThreshold_DAC_step )
	{
		if( pThreshold_DAC == 0 )
			continue;

		cHistName.Form("hComp_Thresold%d_DACunits" , pThreshold_DAC );
		cHistTitle.Form("Comparator State : Th = %d; Time [ns]; Channel Number [a.u]" ,pThreshold_DAC );
		TH2D* hComparator = new TH2D( cHistName.Data() , cHistTitle.Data(), nBunchCrossings*25.0 , 0 , nBunchCrossings*25.0, nTotalNumberOfChannels , 0 , nTotalNumberOfChannels );

		cOut.Form("Filling comparator output histogram for Threshold = %d\n", pThreshold_DAC );
		std::cout << cOut.Data();
		int iModule=0;
		for( int iChannel = 0; iChannel < nTotalNumberOfChannels ; iChannel++ )
		{
			cOut.Form("......... channel %d/%d ....... module%d \n " , iChannel,nTotalNumberOfChannels, iModule);
			if( iChannel%(nCBCs_perModule*nStrips_perCBC) == 0 )
			{
				if( iModule%10 == 0 )
					std::cout << cOut.Data();
				iModule++;
			}

			int iChannelBin = hPostAmplifier_wNoise->GetYaxis()->FindBin(iChannel);
		 	TH1D* hSamplingPostAmp = hPostAmplifier_wNoise->ProjectionX(Form("pPostAmp_Ch%d",iChannel), iChannelBin, iChannelBin);
		 	// if postAmplifier has fired at all then fill comparator output for all times in the simulation 
		 	if( hSamplingPostAmp->GetEntries() > 0 )
		 	{
		 		TString cTitle;
		 		cTitle.Form("Channel%d", iChannel);
		 		hSamplingPostAmp->SetTitle(cTitle.Data());
		 		hSamplingPostAmp->GetXaxis()->SetTitleOffset(1.4);
		 		

		 		for(int iTimeBin =0; iTimeBin < hComparator->GetXaxis()->GetNbins() ;iTimeBin++)
		 		{
		 			double xValue = hComparator->GetXaxis()->GetBinCenter(iTimeBin) - 0.5*hComparator->GetXaxis()->GetBinWidth(iTimeBin) ;
		 			int iChannelBin = hComparator->GetYaxis()->FindBin(iChannel);
		 			int iBin = hComparator->FindBin(xValue,iChannel);
		 			double cSignal = ( hPostAmplifier_wNoise->GetBinContent(iBin) > pThreshold_DAC ) ? 1 : 0 ;
		 			hComparator->Fill( xValue, iChannel, cSignal );
		 		}
		 	}
		}

		cFile = new TFile(pFileName,"UPDATE");
		cFile->cd();
		//TString cDirName;
		//cDirName.Form("Threshold_%d")
		//TDirectory *cDirectory = top->mkdir(cDirName.Data());
		//cDirectory->cd();			
		hComparator->Write(cHistName.Data(), TObject::kOverwrite);
		cFile->Close();
		
	

	}

}

// a few useful functions which tells me when comparator transitions happen.....
double findFirstTimeAboveThreshold(TH1D* hComparator , double pThreshold=0.5, double xMin=0, double xMax=100)
{
	int iBin0 = hComparator->GetXaxis()->FindBin(xMin);
	int iBin1 = hComparator->GetXaxis()->FindBin(xMax);
	double xValue=-1;
	for( int iBin = iBin0 ; iBin <= iBin1 ; iBin++)
	{
		double yValue = hComparator->GetBinContent(iBin);
		if( yValue > pThreshold )
		{
			xValue = hComparator->GetXaxis()->GetBinCenter(iBin) - 0.5*hComparator->GetXaxis()->GetBinWidth(iBin);
			break;
		}
	}
	return xValue;
}
double findLastTimeAboveThreshold(TH1D* hComparator , double pThreshold=0.5, double xMin=0, double xMax=100)
{
	int iBin0 = hComparator->GetXaxis()->FindBin(xMin);
	int iBin1 = hComparator->GetXaxis()->FindBin(xMax);
	double xValue=-1;
	for( int iBin = iBin1 ; iBin >= iBin0 ; iBin--)
	{
		double yValue = hComparator->GetBinContent(iBin);
		double x = hComparator->GetXaxis()->GetBinCenter(iBin) - 0.5*hComparator->GetXaxis()->GetBinWidth(iBin);
		if( yValue > pThreshold )
		{
			xValue = hComparator->GetXaxis()->GetBinCenter(iBin) - 0.5*hComparator->GetXaxis()->GetBinWidth(iBin);
			break;
		}
	}
	return xValue;
}
double findLastTimeBelowThreshold(TH1D* hComparator , double pThreshold=0.5, double xMin=0, double xMax=100)
{
	int iBin0 = hComparator->GetXaxis()->FindBin(xMin);
	int iBin1 = hComparator->GetXaxis()->FindBin(xMax);
	double xValue=-1;
	for( int iBin = iBin1 ; iBin >= iBin0 ; iBin--)
	{
		double yValue = hComparator->GetBinContent(iBin);
		if( yValue < pThreshold )
		{
			xValue = hComparator->GetXaxis()->GetBinCenter(iBin) - 0.5*hComparator->GetXaxis()->GetBinWidth(iBin);
			break;
		}
	}
	return xValue;
}

// returns the state of the latched/sampled hit detect logic for a given BX and delay [ns] of the sampling clock
HitDetectStates checkState(TH1D* hComparator, int pBx=0, int pDLL=0 ,double pDeadTime_Latched=cDeadTime_LatchedMode, double pDeadTime_Sampled=cDeadTime_SampledMode)
{
	double tEdge_DLLshiftedClk =  pBx*25 + pDLL; 
	
	TString cOut;
	//check for transition of comparator in the previous Bx of the delayed sampling clock
	bool cLatchedState = false;  
	double tStartSearch = (pBx-1)*25 + pDLL + pDeadTime_Latched ;// since we're blind for the first pDeadTime_Latched ns of a clock cycle
	double tStopSearch = (pBx-1)*25 + pDLL +  24; // to get to the end of the clock cycle 
	double tFirst_comparatorOn=findFirstTimeAboveThreshold(hComparator, 0.5 , tStartSearch, tStopSearch ); 
	double tLast_comparatorOff=findLastTimeBelowThreshold(hComparator, 0.5 , tStartSearch, tStopSearch ); 
	if( tLast_comparatorOff < tFirst_comparatorOn ) // since I'm only looking for a rising edge ... 
		cLatchedState = true;

	// check if comparator was high on the rising edge of delayed sampling clock 
	bool cSampledState = false; 
	tStartSearch = pBx*25 + pDLL + pDeadTime_Sampled ;// since we're blind for the first pDeadTime_Sampled ns of a clock cycle
	tStopSearch =  pBx*25 +24; // to get to the end of the clock cycle 
	tFirst_comparatorOn = findFirstTimeAboveThreshold(hComparator, 0.5 ,  tStartSearch, tStopSearch ); 
	double tLast_comparatorOn=findLastTimeAboveThreshold(hComparator, 0.5 ,  tStartSearch, tStopSearch ); 
	cOut.Form("Comparator on between : %.1f ns --- %.1f ns in this clock cycle. Edge of DLL shifted clock @ %.1f ns \n" , tFirst_comparatorOn , tLast_comparatorOn, tEdge_DLLshiftedClk);
	//std::cout << cOut.Data();
	// check if my sampling time is within the period where the comparator is on 
	if( tFirst_comparatorOn <= tEdge_DLLshiftedClk && tEdge_DLLshiftedClk <= tLast_comparatorOn ) 
		cSampledState = !( tFirst_comparatorOn < 0 || tLast_comparatorOn < 0); // because the above condition is true if tFirst_comparatorOn==tLast_comparatorOn==-1 

	HitDetectStates cStates;
	cStates.first = cSampledState;
	cStates.second = cLatchedState;
	return cStates;
}
void measureRates(TString pFileName="./ToyMC_CBC3_test.root",  int pThreshold_DAC_min = 0 , int pThreshold_DAC_max = 400, int pThreshold_DAC_step = 10 , int pDLL_min=0, int pDLL_max=25, int pDLL_step=1)
{
	TString cCanvasName, cCanvasTitle;
	TString cHistName, cHistTitle;
	TString cOut;
	
	TH2D* hHits=0;
	TH2D* hDepositedCharge_wNoise=0;
	TH2D* hPostAmplifier_wNoise=0;
	TH1D* hEdeposited_DAC=0;
	TH2D* hDoubleHits=0;

	TString cFileName=pFileName;
	TFile* cFile = new TFile(cFileName,"READ");
	cHistName.Form("hHits_perCh_perBx");
	cFile->GetObject (cHistName.Data() , hHits);
	hHits->SetDirectory(0);  

	cHistName.Form("hDepositedCharge_perCh_perBx");
	cFile->GetObject (cHistName.Data() , hDepositedCharge_wNoise);
	hDepositedCharge_wNoise->SetDirectory(0);  

	cHistName.Form("hPostAmplifierOutput_perCh_perBx");
	cFile->GetObject (cHistName.Data() , hPostAmplifier_wNoise);
	hPostAmplifier_wNoise->SetDirectory(0);  

	cHistName.Form("hEdeposited_DACunits");
	cFile->GetObject (cHistName.Data() , hEdeposited_DAC);
	hEdeposited_DAC->SetDirectory(0);  

	cHistName.Form("hDoubleHits");
	cFile->GetObject (cHistName.Data() , hDoubleHits);
	hDoubleHits->SetDirectory(0);  

	cFile->Close();
	
	int nTotalNumberOfChannels = hHits->GetYaxis()->GetNbins();
	int nBunchCrossings = hHits->GetXaxis()->GetNbins();

	cHistName.Form("hEff_SampledMode");
	cHistTitle.Form("Efficiency in Sampled mode; Time [ns]; Threshold [DAC]");
	TProfile2D* hEfficiency_SampledMode = new TProfile2D(cHistName.Data(),cHistTitle.Data(),  nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	cHistName.Form("hFakes_SampledMode");
	cHistTitle.Form("Fake Rate in Sampled mode; Time [ns]; Threshold [DAC]");
	TProfile2D* hFakeRate_SampledMode = new TProfile2D(cHistName.Data(),cHistTitle.Data(),  nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	cHistName.Form("hEff_LatchedMode");
	cHistTitle.Form("Efficiency in Latched mode; Time [ns]; Threshold [DAC]");
	TProfile2D* hEfficiency_LatchedMode = new TProfile2D(cHistName.Data(),cHistTitle.Data(),  nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	cHistName.Form("hFakes_LatchedMode");
	cHistTitle.Form("Fake Rate in Latched mode; Time [ns]; Threshold [DAC]");
	TProfile2D* hFakeRate_LatchedMode = new TProfile2D(cHistName.Data(),cHistTitle.Data(),  nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	cHistName.Form("hTimeWindow_LatchedMode");
	cHistTitle.Form("Time window for which Latched mode gives 0 perc. fakes; Bunch Crossing ID [a.u.]; Threshold;");
	TH2D* hTimeWindow_zeroFakes_Latched = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings, 0 , nBunchCrossings , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	cHistName.Form("hTimeWindow_LatchedMode_Avg");
	cHistTitle.Form("Time window for which Latched mode gives 0 perc. fakes; Threshold; Window Size [ns]");
	TProfile* hTimeWindow_zeroFakes_LatchedAvg = new TProfile(cHistName.Data(),cHistTitle.Data(), (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	

	cHistName.Form("hTimeWindow_SampledMode");
	cHistTitle.Form("Time window for which Sampled mode gives 0 perc. fakes; Bunch Crossing ID [a.u.]; Threshold;");
	TH2D* hTimeWindow_zeroFakes_Sampled = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings, 0 , nBunchCrossings , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	cHistName.Form("hTimeWindow_SampledMode_Avg");
	cHistTitle.Form("Time window for which Sampled mode gives 0 perc. fakes; Threshold; Window Size [ns]");
	TProfile* hTimeWindow_zeroFakes_SampledAvg = new TProfile(cHistName.Data(),cHistTitle.Data(), (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	cHistName.Form("hSamplingWindow_Eff_SampledMode");
	cHistTitle.Form("Time window for which Sampled mode gives > 90 perc. efficiency; Time [ns]; Threshold;");
	TH2D* hSamplingWindow_Eff_Sampled = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	cHistName.Form("hSamplingWindow_Faked_SampledMode");
	cHistTitle.Form("Time window for which Sampled mode gives 0 perc. fakes; Time [ns]; Threshold;");
	TH2D* hSamplingWindow_Fakes_Sampled = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	
	cHistName.Form("hSamplingWindow_SampledMode");
	cHistTitle.Form("Time window for which Sampled mode gives 0 perc. fakes AND > 90 perc. efficiency; Time [ns]; Threshold;");
	TH2D* hSamplingWindow_Sampled = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	cHistName.Form("hSamplingWindow_Eff_LatchedMode");
	cHistTitle.Form("Time window for which Latched mode gives > 90 perc. efficiency; Time [ns]; Threshold;");
	TH2D* hSamplingWindow_Eff_Latched = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	cHistName.Form("hSamplingWindow_Faked_LatchedMode");
	cHistTitle.Form("Time window for which Latched mode gives 0 perc. fakes; Time [ns]; Threshold;");
	TH2D* hSamplingWindow_Fakes_Latched = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	
	cHistName.Form("hSamplingWindow_LatchedMode");
	cHistTitle.Form("Time window for which Latched mode gives 0 perc. fakes AND > 90 perc. efficiency; Time [ns]; Threshold;");
	TH2D* hSamplingWindow_Latched = new TH2D(cHistName.Data(),cHistTitle.Data(), nBunchCrossings*25 , 0 , nBunchCrossings*25 , (pThreshold_DAC_max-pThreshold_DAC_min)/pThreshold_DAC_step , pThreshold_DAC_min , pThreshold_DAC_max );
	

	// first loop over threshold settings 
	for( int pThreshold_DAC = pThreshold_DAC_min ; pThreshold_DAC < pThreshold_DAC_max  ; pThreshold_DAC +=pThreshold_DAC_step )
	{
		if( pThreshold_DAC == 0 )
			continue;

		TH2D* hComparator=0;
		cHistName.Form("hComp_Thresold%d_DACunits" , pThreshold_DAC );
		cFile = new TFile(cFileName,"READ");
		cFile->cd();
		cFile->GetObject (cHistName.Data() , hComparator);
		if( hComparator == NULL)
		{
			cFile->Close();
			continue;
		}	
		hComparator->SetDirectory(0);  
		cFile->Close();

		cOut.Form(" ..... Threshold set to %d DAC units ..... \n" , pThreshold_DAC );
		std::cout << cOut.Data();
		
		cFile = new TFile(cFileName,"UPDATE");
		cFile->cd();
		for( int pDLL = pDLL_min ; pDLL < pDLL_max ; pDLL+=pDLL_step )
		{
			cHistName.Form("hSampledMode_DLL%d_Th%d", pDLL , pThreshold_DAC );
			cHistTitle.Form("Sampled Mode : Th = %d , DLL = %d; Bunch Crossing ID [a.u]; Channel Number [a.u]" ,pThreshold_DAC, pDLL );
			TH2D* hSampled = new TH2D( cHistName.Data() , cHistTitle.Data(), nBunchCrossings , 0 , nBunchCrossings, nTotalNumberOfChannels , 0 , nTotalNumberOfChannels );
			cHistName.Form("hLatchedMode_DLL%d_Th%d", pDLL , pThreshold_DAC );
			cHistTitle.Form("Latched Mode : Th = %d , DLL = %d; Bunch Crossing ID [a.u]; Channel Number [a.u]" ,pThreshold_DAC, pDLL );
			TH2D* hLatched  = new TH2D( cHistName.Data() , cHistTitle.Data(), nBunchCrossings , 0 , nBunchCrossings, nTotalNumberOfChannels , 0 , nTotalNumberOfChannels );
			

			cOut.Form("............... DLL set to %d ns .................\n", pDLL );
			if( pDLL%5 == 0 )
				std::cout << cOut.Data();

			for( int iChannel = 0; iChannel < nTotalNumberOfChannels ; iChannel++ )
			{
				int iChannelBin = hComparator->GetYaxis()->FindBin(iChannel);
				TH1D* hHits_thisChannel = hHits->ProjectionX(Form("hHits_Ch%d", iChannel) , iChannelBin, iChannelBin);
				//dont bother if this channel was never hit...
				if( hHits_thisChannel->GetEntries() == 0 )
					continue;

				TH1D* hComparator_thisChannel = hComparator->ProjectionX(Form("pComparator_Ch%d",iChannel), iChannelBin, iChannelBin);
		 		for(int iBx =0;iBx < nBunchCrossings ; iBx++)
		 		{
		 			int iBin = hHits->FindBin(iBx,iChannel);
		 			bool isHit = ( hHits->GetBinContent(iBin) > 0) ;
		 			
		 			// figure out what the different hit detection logic circuits say for this setting
		 			HitDetectStates cHitDetectStates = checkState(hComparator_thisChannel, iBx, pDLL );
		 			bool cSampledState = cHitDetectStates.first ; 
		 			bool cLatchedState = cHitDetectStates.second; 
		 			
		 			// if there was a hit in this channel in this Bx 
					if( isHit )
					{
						// fill profiles counting fakes with 0 - since this is a real hit 
						hFakeRate_LatchedMode->Fill( iBx*25 + pDLL, pThreshold_DAC, 0);
						hFakeRate_SampledMode->Fill( iBx*25 + pDLL, pThreshold_DAC, 0);
						
						// first for the latched mode 
						// fill profile with state of latched mode 
						hEfficiency_LatchedMode->Fill( iBx*25 + pDLL, pThreshold_DAC, (cLatchedState) ? 1 : 0 );
						// if the latched mode thinks there's a hit here 
						if( cLatchedState ) 
						{
							// then for this threshold and DLL setting fill the histogram which shows the state of the digital logic in this BX 
							hLatched->Fill(iBx, iChannel );
						}
						
						// now for the sampled mode 
						// fillprofile counting hits with the state of the sampled mode 
						hEfficiency_SampledMode->Fill( iBx*25 + pDLL, pThreshold_DAC,  (cSampledState) ? 1 : 0 );
						// if the sampled mode thinks there's a hit here 
						if( cSampledState ) 
						{
							// then for this threshold and DLL setting fill the histogram which shows the state of the digital logic in this BX 
							hSampled->Fill(iBx,iChannel);
						}
						
					}
					// if there was no hit in this Bx
					else
					{
						// if latched mode shows a hit ... then this is a fake  - and fill profile counting fakes with 1 
						if( cLatchedState )
						{
							hFakeRate_LatchedMode->Fill( iBx*25 + pDLL, pThreshold_DAC,  1);
						}

						// if sampled mode shows a hit ... then this is a fake  - and fill profile counting fakes with 1 
						if( cSampledState )
						{
							if( iBx == 0 && pDLL == 0 )
							{
								cOut.Form("Fake hit found in channel %d in Bx0 \n", iChannel);
								std::cout << cOut.Data();
							}
							hFakeRate_SampledMode->Fill( iBx*25 + pDLL, pThreshold_DAC, 1);
						}
					}
		 		}
			}
			cHistName.Form("hSampledMode_DLL%d_Th%d", pDLL , pThreshold_DAC );
			hSampled->Write(cHistName.Data(), TObject::kOverwrite);
			cHistName.Form("hLatchedMode_DLL%d_Th%d", pDLL , pThreshold_DAC );
			hLatched->Write(cHistName.Data(), TObject::kOverwrite);
		}
		cFile->Close();


		int iBin = hFakeRate_LatchedMode->GetYaxis()->FindBin(pThreshold_DAC);
		TH1D* hFakesLatched = (TH1D*)hFakeRate_LatchedMode->ProjectionX("pFakesLatched",iBin,iBin);
		TH1D* hFakesSampled = (TH1D*)hFakeRate_SampledMode->ProjectionX("pFakesSampled",iBin,iBin);

		iBin = hEfficiency_LatchedMode->GetYaxis()->FindBin(pThreshold_DAC);
		TH1D* hEffSampled = (TH1D*)hEfficiency_SampledMode->ProjectionX("pEffSampled",iBin,iBin);
		TH1D* hEffLatched = (TH1D*)hEfficiency_LatchedMode->ProjectionX("pEffLatched",iBin,iBin);
		for( int iBx=1 ; iBx < nBunchCrossings ; iBx++ )
		{
			double t_Latched = findLastTimeAboveThreshold(hFakesLatched , cMaxAcceptableFakeRate, iBx*25,iBx*25+24.0);
			double t_Latched_End =(iBx+1)*25;
			double tAccFakRate_Latched = t_Latched;
			double tAccFakRate_Latched_End = (iBx+1)*25-1;
			
			double t_Sampled = findLastTimeAboveThreshold(hFakesSampled , cMaxAcceptableFakeRate, iBx*25,iBx*25+24.0);
			double t_Sampled_End =(iBx+1)*25;
			// for(int tAboveThreshold = (int)t_Sampled ; tAboveThreshold <= (int)t_Sampled_End ; tAboveThreshold++)
			// {
			// 	hSamplingWindow_Fakes_Sampled->Fill(tAboveThreshold, pThreshold_DAC);
			// }
			double tAccFakRate_Sampled = t_Sampled;
			double tAccFakRate_Sampled_End = (iBx+1)*25-1;
			
			double twindow_Latched = (t_Latched>0) ? std::fabs((iBx+1)*25-t_Latched) : 25 ; 
			double twindow_Sampled = (t_Sampled>0) ? std::fabs((iBx+1)*25-t_Sampled) : 25 ;
			hTimeWindow_zeroFakes_Latched->Fill(iBx, pThreshold_DAC, twindow_Latched );
			hTimeWindow_zeroFakes_LatchedAvg->Fill(pThreshold_DAC, twindow_Latched );
			hTimeWindow_zeroFakes_Sampled->Fill(iBx, pThreshold_DAC, twindow_Sampled );
			hTimeWindow_zeroFakes_SampledAvg->Fill(pThreshold_DAC, twindow_Sampled );
			
			t_Sampled = findFirstTimeAboveThreshold(hEffSampled , cMinAcceptableEfficiency, iBx*25,iBx*25+24.0);
			t_Sampled_End = findLastTimeAboveThreshold(hEffSampled , cMinAcceptableEfficiency, iBx*25,iBx*25+24.0);
			cOut.Form("............... Bx%d \n", iBx);
			//std::cout << cOut.Data();
			cOut.Form("............................................ in sampled mode : %.1f -- %.1f ns [eff] , below acceptable fake rate between %.1f -- %.1f ns [fakes] \n" ,  t_Sampled , t_Sampled_End , tAccFakRate_Sampled ,tAccFakRate_Sampled_End);
			//std::cout << cOut.Data();
			cOut.Form("............................................ in latched mode : %.1f -- %.1f ns [eff] , below acceptable fake rate between %.1f -- %.1f ns [fakes] \n" ,  t_Latched , t_Latched_End , tAccFakRate_Latched ,tAccFakRate_Latched_End);
			//std::cout << cOut.Data();
			//cOut.Form("............... ................ size of acceptable window = %.3f ns[ latched mode ] %.3f ns [sampled mode]\n" , twindow_Latched , twindow_Sampled);
			//std::cout << cOut.Data();
			
			// //std::cout << cOut.Data();
			// for(int tAboveThreshold = (int)t_Sampled ; tAboveThreshold <= (int)t_Sampled_End ; tAboveThreshold++)
			// {
			// 	hSamplingWindow_Eff_Sampled->Fill(tAboveThreshold, pThreshold_DAC);
			// 	//cOut.Form("%d ns , ", tAboveThreshold );
			// 	//std::cout << cOut.Data();
			// }

			for( int tWithinBx = iBx*25 ; tWithinBx < (iBx+1)*25 ; tWithinBx++)
			{
				bool withinSamplWindow_Fakes = (tWithinBx >= tAccFakRate_Sampled && tWithinBx <= tAccFakRate_Sampled_End);
				bool withinSamplWindow_Eff = (tWithinBx >= t_Sampled && tWithinBx <= t_Sampled_End);

				if( withinSamplWindow_Eff )
					hSamplingWindow_Eff_Sampled->Fill(tWithinBx, pThreshold_DAC);
				if( withinSamplWindow_Fakes )
					hSamplingWindow_Fakes_Sampled->Fill(tWithinBx, pThreshold_DAC);

				if( withinSamplWindow_Fakes && withinSamplWindow_Eff)
				{
					hSamplingWindow_Sampled->Fill(tWithinBx, pThreshold_DAC);
				}

				withinSamplWindow_Fakes = (tWithinBx >= tAccFakRate_Latched && tWithinBx <= tAccFakRate_Latched_End);
				withinSamplWindow_Eff = (tWithinBx >= t_Latched && tWithinBx <= t_Latched_End);

				if( withinSamplWindow_Eff )
					hSamplingWindow_Eff_Latched->Fill(tWithinBx, pThreshold_DAC);
				if( withinSamplWindow_Fakes )
					hSamplingWindow_Fakes_Latched->Fill(tWithinBx, pThreshold_DAC);

				if( withinSamplWindow_Fakes && withinSamplWindow_Eff)
				{
					hSamplingWindow_Latched->Fill(tWithinBx, pThreshold_DAC);
				}

			}
			//std::cout << "\n";
			// cOut.Form("............... Bx%d - above acceptable efficiency between %.1f -- %.1f ns [sampled mode] , %.1f --- %.1f ns [latched mode]\n" , iBx, t_Sampled , t_Sampled_End, t_Latched , t_Latched_End );
			// std::cout << cOut.Data();
		}
	}
	cFile = new TFile(cFileName,"UPDATE");
	cFile->cd();
	cHistName.Form("hEfficiency_SampledMode_ns");
	hEfficiency_SampledMode->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hFakeRate_SampledMode_ns");
	hFakeRate_SampledMode->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hEfficiency_LatchedMode_ns");
	hEfficiency_LatchedMode->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hFakeRate_LatchedMode_ns");
	hFakeRate_LatchedMode->Write(cHistName.Data(), TObject::kOverwrite);

	cHistName.Form("hTimeWindow_SampledMode");
	hTimeWindow_zeroFakes_Sampled->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hTimeWindow_LatchedMode");
	hTimeWindow_zeroFakes_Latched->Write(cHistName.Data(), TObject::kOverwrite);

	cHistName.Form("hAvgTimeWindow_SampledMode");
	hTimeWindow_zeroFakes_SampledAvg->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hAvgTimeWindow_LatchedMode");
	hTimeWindow_zeroFakes_LatchedAvg->Write(cHistName.Data(), TObject::kOverwrite);

	cHistName.Form("hSamplingWindow_Eff_SampledMode");
	hSamplingWindow_Eff_Sampled->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hSamplingWindow_Fakes_SampledMode");
	hSamplingWindow_Fakes_Sampled->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hSamplingWindow_SampledMode");
	hSamplingWindow_Sampled->Write(cHistName.Data(), TObject::kOverwrite);

	cHistName.Form("hSamplingWindow_Eff_LatchedMode");
	hSamplingWindow_Eff_Latched->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hSamplingWindow_Fakes_LatchedMode");
	hSamplingWindow_Fakes_Latched->Write(cHistName.Data(), TObject::kOverwrite);
	cHistName.Form("hSamplingWindow_LatchedMode");
	hSamplingWindow_Latched->Write(cHistName.Data(), TObject::kOverwrite);

	
	cFile->Close();


}
void TrackerSimulation ( int pNmodules=20 , double pOccupancy=1e-2, int pNbunchCrossings=10, int pThreshold_DAC_min=0, int pThreshold_DAC_max=300, int pThreshold_DAC_step=5)
{ 

	TString pFileName_Modifier= (ChargeSharing) ? "wChargeSharing" : "noChargeSharing";
	TString pFileName; 
	pFileName.Form("%dmodules_%s.root", pNmodules, pFileName_Modifier.Data()); 

	TString cOut;
	int nModules=pNmodules;
	
	auto start = std::chrono::system_clock::now();
    //testChargeSharing();

    // populate hits for pNmodules modules [ only need to do this once]
    double cMPcharge_DACunits = 143.652;//
    double cLandauWidth_DACunits = 7.09689;//
    double cNoise_DACunits = 17.4744;//
    double cPedestal_kElectrons = 2;//
	//populateHits(pFileName, nModules,pNbunchCrossings,pOccupancy,cMPcharge_DACunits*cConversion_CBC3,cLandauWidth_DACunits*cConversion_CBC3,cPedestal_kElectrons,cNoise_DACunits*cConversion_CBC3, cSizeChargeSharingRegion);
	
	// then fill comparator outputs [only need to do this once]
	//populateComparator(pFileName, pThreshold_DAC_min ,pThreshold_DAC_max,pThreshold_DAC_step);
	// figure out how long this took
	auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t (end);
    cOut.Form ("It took %.3f mins to populate hits/comparators for this tracker with %d modules\n", elapsed_seconds.count() / 60.0 ,nModules);
    std::cout << cOut.Data();
    start = std::chrono::system_clock::now();
	measureRates(pFileName,  pThreshold_DAC_min , pThreshold_DAC_max , pThreshold_DAC_step );
	end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    end_time = std::chrono::system_clock::to_time_t (end);
    cOut.Form ("It took %.3f mins to measure hit/fake rates for this tracker with %d modules\n", elapsed_seconds.count() / 60.0 ,nModules);
    std::cout << cOut.Data();
}
void showEffect_chargeSharing(TString pFileName0="./20modules_noChargeSharing.root",TString pFileName1="./20modules_wChargeSharing.root" , int pExampleChannel=120, int pExampleThreshold=20 ,  int pExampleDLL = 10 )
{
	//set_plot_style();
	//gStyle->SetPalette(1);
	TString cHistName;

	TFile* f_wChargeSharing = new TFile(pFileName1,"READ");
	f_wChargeSharing->cd();

	TH2D* hHits_wChargESharing=0;
	cHistName.Form("hHits_perCh_perBx");
	f_wChargeSharing->GetObject (cHistName.Data() , hHits_wChargESharing);
	hHits_wChargESharing->SetDirectory(0);  

	TH2D* hChargeDeposited_wChargeSharing=0;
	cHistName.Form("hDepositedCharge_perCh_perBx");
	f_wChargeSharing->GetObject (cHistName.Data() , hChargeDeposited_wChargeSharing);
	hChargeDeposited_wChargeSharing->SetDirectory(0);  

	TH2D* hPostAmplifierOutput_wChargeSharing=0;
	cHistName.Form("hPostAmplifierOutput_perCh_perBx");
	f_wChargeSharing->GetObject (cHistName.Data() , hPostAmplifierOutput_wChargeSharing);
	hPostAmplifierOutput_wChargeSharing->SetDirectory(0);  

	TH2D* hComparatorOutput_wChargeSharing=0;
	cHistName.Form("hComp_Thresold%d_DACunits" , pExampleThreshold );
	f_wChargeSharing->GetObject (cHistName.Data() , hComparatorOutput_wChargeSharing);
	hComparatorOutput_wChargeSharing->SetDirectory(0);  

	TH2D* hComparatorOutput_Sampled_wChargeSharing=0;
	cHistName.Form("hSampledMode_DLL%d_Th%d", pExampleDLL , pExampleThreshold );
	f_wChargeSharing->GetObject (cHistName.Data() , hComparatorOutput_Sampled_wChargeSharing);
	hComparatorOutput_Sampled_wChargeSharing->SetDirectory(0);  

	TH2D* hComparatorOutput_Latched_wChargeSharing=0;
	cHistName.Form("hLatchedMode_DLL%d_Th%d", pExampleDLL , pExampleThreshold );
	f_wChargeSharing->GetObject (cHistName.Data() , hComparatorOutput_Latched_wChargeSharing);
	hComparatorOutput_Latched_wChargeSharing->SetDirectory(0);  

		
	TH2D* hEff_Sampled_wChargeSharing=0;
	cHistName.Form("hEfficiency_SampledMode_ns");
	f_wChargeSharing->GetObject (cHistName.Data() , hEff_Sampled_wChargeSharing);
	hEff_Sampled_wChargeSharing->SetDirectory(0);  

	cHistName.Form("hFakeRate_SampledMode_ns");
	TH2D* hFakes_Sampled_wChargeSharing=0;
	f_wChargeSharing->GetObject (cHistName.Data() , hFakes_Sampled_wChargeSharing);
	hFakes_Sampled_wChargeSharing->SetDirectory(0);  
		
	cHistName.Form("hAvgTimeWindow_SampledMode");
	TProfile* hTimeWindow_SampledMode_wChargeSharing=0;
	f_wChargeSharing->GetObject (cHistName.Data() , hTimeWindow_SampledMode_wChargeSharing);
	hTimeWindow_SampledMode_wChargeSharing->SetDirectory(0);  

	cHistName.Form("hTimeWindow_SampledMode");
	TH2D* hTimeWindow_SampledMode_wChargeSharing2D=0;
	f_wChargeSharing->GetObject (cHistName.Data() , hTimeWindow_SampledMode_wChargeSharing2D);
	hTimeWindow_SampledMode_wChargeSharing2D->SetDirectory(0);  
	f_wChargeSharing->Close();

	TFile* f_woutChargeSharing = new TFile(pFileName0,"READ");
	f_woutChargeSharing->cd();

	TH2D* hEff_Sampled_woutChargeSharing=0;
	cHistName.Form("hEfficiency_SampledMode_ns");
	f_woutChargeSharing->GetObject (cHistName.Data() , hEff_Sampled_woutChargeSharing);
	hEff_Sampled_woutChargeSharing->SetDirectory(0);  

	cHistName.Form("hFakeRate_SampledMode_ns");
	TH2D* hFakes_Sampled_woutChargeSharing=0;
	f_woutChargeSharing->GetObject (cHistName.Data() , hFakes_Sampled_woutChargeSharing);
	hFakes_Sampled_woutChargeSharing->SetDirectory(0);  
	
	cHistName.Form("hAvgTimeWindow_SampledMode");
	TProfile* hTimeWindow_SampledMode_woutChargeSharing=0;
	f_woutChargeSharing->GetObject (cHistName.Data() , hTimeWindow_SampledMode_woutChargeSharing);
	hTimeWindow_SampledMode_woutChargeSharing->SetDirectory(0);  
	
	cHistName.Form("hTimeWindow_SampledMode");
	TH2D* hTimeWindow_SampledMode_woutChargeSharing2D=0;
	f_woutChargeSharing->GetObject (cHistName.Data() , hTimeWindow_SampledMode_woutChargeSharing2D);
	hTimeWindow_SampledMode_woutChargeSharing2D->SetDirectory(0);  
	f_woutChargeSharing->Close();
	
	double cTime = 40; 
	int iBin = hEff_Sampled_wChargeSharing->GetXaxis()->FindBin(cTime);
	TH1D* hThresholdScan_woutChargeSharing = (TH1D*)hEff_Sampled_woutChargeSharing->ProjectionY("pScanwoutChargeSharing", iBin,iBin);
	TH1D* hThresholdScan_wChargeSharing = (TH1D*)hEff_Sampled_wChargeSharing->ProjectionY("pScanwChargeSharing", iBin,iBin);

	// TCanvas* cGeneratedHits = new TCanvas("cGeneratedHits","Generated hits",350*4,350*2);
	// cGeneratedHits->Divide(4,2);
	// cGeneratedHits->cd(1);
	// hHits_wChargESharing->SetStats(0);
	// hHits_wChargESharing->GetYaxis()->SetRangeUser(0,1e3);
	// hHits_wChargESharing->GetXaxis()->SetRangeUser(0,2);
	// hHits_wChargESharing->DrawCopy("colz");
	// cGeneratedHits->cd(2);
	// hChargeDeposited_wChargeSharing->SetStats(0);
	// hChargeDeposited_wChargeSharing->GetYaxis()->SetRangeUser(0,1e3);
	// hChargeDeposited_wChargeSharing->GetXaxis()->SetRangeUser(0,2);
	// hChargeDeposited_wChargeSharing->GetZaxis()->SetRangeUser(0,500);
	// hChargeDeposited_wChargeSharing->DrawCopy("colz");
	// cGeneratedHits->cd(3);
	// hPostAmplifierOutput_wChargeSharing->SetStats(0);
	// hPostAmplifierOutput_wChargeSharing->GetYaxis()->SetRangeUser(0,1e3);
	// hPostAmplifierOutput_wChargeSharing->GetXaxis()->SetRangeUser(0,2*25);
	// hPostAmplifierOutput_wChargeSharing->GetZaxis()->SetRangeUser(0,1000);
	// hPostAmplifierOutput_wChargeSharing->DrawCopy("colz");
	// cGeneratedHits->cd(4);
	// hComparatorOutput_wChargeSharing->SetStats(0);
	// hComparatorOutput_wChargeSharing->GetYaxis()->SetRangeUser(0,1e3);
	// hComparatorOutput_wChargeSharing->GetXaxis()->SetRangeUser(0,2*25);
	// hComparatorOutput_wChargeSharing->DrawCopy("colz");

	// TCanvas* cComparatorResponse = new TCanvas("cComparatorResponse","Comparator logic",350*4,350*2);
	// int cChannelBin = hChargeDeposited_wChargeSharing->GetYaxis()->FindBin(pExampleChannel);
	// TH1D* hPostAmpOut_wChargeSharing = (TH1D*)hPostAmplifierOutput_wChargeSharing->ProjectionX(Form("paOutput_wChargeSharing"), cChannelBin , cChannelBin);
	// TH1D* hCompOut_wChargeSharing = (TH1D*)hComparatorOutput_wChargeSharing->ProjectionX(Form("compOutput_wChargeSharing"), cChannelBin , cChannelBin);
	// TH1D* hCompSmapled_wChargeSharing = (TH1D*)hComparatorOutput_Sampled_wChargeSharing->ProjectionX(Form("sampledOutput_wChargeSharing"), cChannelBin , cChannelBin);
	// TH1D* hCompLatched_wChargeSharing = (TH1D*)hComparatorOutput_Latched_wChargeSharing->ProjectionX(Form("latchedOutput_wChargeSharing"), cChannelBin , cChannelBin);
	// cComparatorResponse->Divide(4,2);
	// cComparatorResponse->cd(1);
	// hPostAmpOut_wChargeSharing->SetStats(0);
	// hPostAmpOut_wChargeSharing->DrawCopy("hist");
	// cComparatorResponse->cd(2);
	// hCompOut_wChargeSharing->SetStats(0);
	// hCompOut_wChargeSharing->DrawCopy("hist");
	// cComparatorResponse->cd(3);
	// hCompSmapled_wChargeSharing->SetStats(0);
	// hCompSmapled_wChargeSharing->GetXaxis()->SetRangeUser(0,2);
	// hCompSmapled_wChargeSharing->DrawCopy("hist");
	// cComparatorResponse->cd(4);
	// hCompLatched_wChargeSharing->SetStats(0);
	// hCompLatched_wChargeSharing->GetXaxis()->SetRangeUser(0,2);
	// hCompLatched_wChargeSharing->DrawCopy("hist");
	

	TCanvas* c1 = new TCanvas("c1","c1",350*2,350*2);
	c1->Divide(2,2);
	c1->cd(1);
	c1->cd(1)->SetLogz(0);
	hEff_Sampled_wChargeSharing->SetTitle("P(detecting a hit) in sampled mode [w/ simple charge sharing model]");
	hEff_Sampled_wChargeSharing->SetStats(0);
	hEff_Sampled_wChargeSharing->GetXaxis()->SetRangeUser(25,50);
	hEff_Sampled_wChargeSharing->GetYaxis()->SetRangeUser(5,200);
	hEff_Sampled_wChargeSharing->GetZaxis()->SetRangeUser(0.95,1);
	hEff_Sampled_wChargeSharing->DrawCopy("colz");
	c1->cd(2);
	c1->cd(2)->SetLogz(0);
	hFakes_Sampled_wChargeSharing->SetTitle("P(Fake) in sampled mode [w/ simple charge sharing model]");
	hFakes_Sampled_wChargeSharing->SetStats(0);
	hFakes_Sampled_wChargeSharing->GetXaxis()->SetRangeUser(25,50);
	hFakes_Sampled_wChargeSharing->GetYaxis()->SetRangeUser(5,200);
	hEff_Sampled_wChargeSharing->GetZaxis()->SetRangeUser(0.9,1);
	hFakes_Sampled_wChargeSharing->DrawCopy("colz");
	c1->cd(3);
	c1->cd(3)->SetLogz(0);
	hEff_Sampled_woutChargeSharing->SetTitle("P(detecting a hit) in sampled mode");
	hEff_Sampled_woutChargeSharing->SetStats(0);
	hEff_Sampled_woutChargeSharing->GetXaxis()->SetRangeUser(25,50);
	hEff_Sampled_woutChargeSharing->GetYaxis()->SetRangeUser(5,200);
	hEff_Sampled_woutChargeSharing->GetZaxis()->SetRangeUser(0.95,1);
	hEff_Sampled_woutChargeSharing->DrawCopy("colz");
	c1->cd(4);
	c1->cd(4)->SetLogz(1);
	hFakes_Sampled_woutChargeSharing->SetTitle("P(Fake) in sampled mode");
	hFakes_Sampled_woutChargeSharing->SetStats(0);
	hFakes_Sampled_woutChargeSharing->GetXaxis()->SetRangeUser(25,50);
	hFakes_Sampled_woutChargeSharing->GetYaxis()->SetRangeUser(5,200);
	hFakes_Sampled_woutChargeSharing->GetZaxis()->SetRangeUser(1e-3,1);
	hFakes_Sampled_woutChargeSharing->DrawCopy("colz");


	TCanvas* c = new TCanvas("c","c",350*2,350);
	c->Divide(2,1);
	c->cd(1);
	hTimeWindow_SampledMode_woutChargeSharing2D->GetXaxis()->SetRangeUser(1,10);
	hTimeWindow_SampledMode_woutChargeSharing2D->SetStats(0);
	hTimeWindow_SampledMode_woutChargeSharing2D->GetZaxis()->SetRangeUser(10,25);
	hTimeWindow_SampledMode_woutChargeSharing2D->DrawCopy("colz");
	c->cd(2);
	hTimeWindow_SampledMode_wChargeSharing2D->GetXaxis()->SetRangeUser(1,10);
	hTimeWindow_SampledMode_wChargeSharing2D->SetStats(0);
	hTimeWindow_SampledMode_wChargeSharing2D->GetZaxis()->SetRangeUser(10,25);
	hTimeWindow_SampledMode_wChargeSharing2D->DrawCopy("colz");
	
	TCanvas* c3 = new TCanvas("c3","c3",350*2,350);
	c3->Divide(2,1);
	c3->cd(1);
	hThresholdScan_woutChargeSharing->GetYaxis()->SetTitle("Hit finding Efficiency in sampled mode");
	hThresholdScan_woutChargeSharing->GetXaxis()->SetTitleOffset(1.2);
	hThresholdScan_woutChargeSharing->SetStats(0);
	hThresholdScan_woutChargeSharing->GetXaxis()->SetRangeUser(0,200);
	hThresholdScan_woutChargeSharing->SetLineColor(kBlue);
	hThresholdScan_woutChargeSharing->DrawCopy("eHisto");
	hThresholdScan_wChargeSharing->GetXaxis()->SetRangeUser(0,200);
	hThresholdScan_wChargeSharing->SetLineColor(kRed);
	hThresholdScan_wChargeSharing->DrawCopy("eHistoSAME");
	c3->cd(2);
	hTimeWindow_SampledMode_woutChargeSharing->GetXaxis()->SetTitleOffset(1.2);
	hTimeWindow_SampledMode_woutChargeSharing->SetStats(0);
	hTimeWindow_SampledMode_woutChargeSharing->GetXaxis()->SetRangeUser(0,200);
	hTimeWindow_SampledMode_woutChargeSharing->SetLineColor(kBlue);
	hTimeWindow_SampledMode_woutChargeSharing->DrawCopy("ehisto");
	hTimeWindow_SampledMode_wChargeSharing->SetLineColor(kRed);
	hTimeWindow_SampledMode_wChargeSharing->DrawCopy("ehistoSAME");
	TLine*  lLimit = new TLine(0, 10 , 200 , 10);
	lLimit->SetLineColor(kBlack);
	lLimit->Draw("same");
	
}
void TrackerSim(TString pFileName0="./Test_40modules_wChargeSharing.root",TString pFileName1="./Test_40modules_noChargeSharing.root" ,  double pLandauMPV_kElectrons=25, double pLandauWidth_kElectrons=2, double pNoise_kElectrons=0.8)
{
	// int cNmodules=10;
	// int cNbunchCrossings=10;
			
	// std::vector<double> cSizeChargeSharingRegion={0.02}; //0.01
	// std::vector<double> cOccupancy = {5e-2};
	// for( auto pOccupancy: cOccupancy)
	// {
	// 	TString cFileName;
	// 	cFileName.Form("./Test_%dmodules_%dBXs_%dx1em3_PercOcc.root",cNmodules, cNbunchCrossings,  (int)(pOccupancy*1e3) );
	// 	std::cout << cFileName.Data() << "\n";
	// 	TrackerSimulation(cFileName.Data(), pOccupancy , 0.01, cNmodules, cNbunchCrossings);
	// }

	// TString cHistName;

	// TFile* f_wChargeSharing = new TFile(pFileName0,"READ");
	// f_wChargeSharing->cd();
	// TH2D* hEff_Sampled_wChargeSharing=0;
	// cHistName.Form("hEfficiency_SampledMode_ns");
	// f_wChargeSharing->GetObject (cHistName.Data() , hEff_Sampled_wChargeSharing);
	// hEff_Sampled_wChargeSharing->SetDirectory(0);  

	// cHistName.Form("hFakeRate_SampledMode_ns");
	// TH2D* hFakes_Sampled_wChargeSharing=0;
	// f_wChargeSharing->GetObject (cHistName.Data() , hFakes_Sampled_wChargeSharing);
	// hFakes_Sampled_wChargeSharing->SetDirectory(0);  


	// cHistName.Form("hEfficiency_LatchedMode_ns");
	// TH2D* hEff_Latched_wChargeSharing=0;
	// f_wChargeSharing->GetObject (cHistName.Data() , hEff_Latched_wChargeSharing);
	// hEff_Latched_wChargeSharing->SetDirectory(0);  

	// cHistName.Form("hFakeRate_LatchedMode_ns");
	// TH2D* hFakes_Latched_wChargeSharing=0;
	// f_wChargeSharing->GetObject (cHistName.Data() , hFakes_Latched_wChargeSharing);
	// hFakes_Latched_wChargeSharing->SetDirectory(0);  
	

	// f_wChargeSharing->Close();

	// TH1D* cProfile_Sampled = new TH1D("cTimeSampled","Sampled;Threshold [DAC]; Time within BX with 0 fakes.", 150/5 , 0, 150);
	// TH1D* cProfile_Latched = new TH1D("cTimeLatched","Latched;Threshold [DAC]; Time within BX with 0 fakes.", 150/5 , 0, 150);
	// TCanvas* cCompareThresholds = new TCanvas("cCompareThresholds","Compare Thresholds",350*6,350*5);
 //    cCompareThresholds->Divide(6,5);
 //    int cPadCounter=0;
	// for( int cThreshold=0; cThreshold <= 150; cThreshold+=5)
	// {
	// 	int iBin = hEff_Sampled_wChargeSharing->GetYaxis()->FindBin((double)cThreshold);
	// 	TH1D* hEff_wChargeSharing = (TH1D*)hEff_Sampled_wChargeSharing->ProjectionX("pEff_wChargeSharing",iBin,iBin);
	// 	TH1D* hFakes_wChargeSharing = (TH1D*)hFakes_Sampled_wChargeSharing->ProjectionX("pFakes_wChargeSharing",iBin,iBin);
	// 	TH1D* hEffLatched_wChargeSharing = (TH1D*)hEff_Latched_wChargeSharing->ProjectionX("pEff_Latched_wChargeSharing",iBin,iBin);
	// 	TH1D* hFakesLatched_wChargeSharing = (TH1D*)hFakes_Latched_wChargeSharing->ProjectionX("pFakes_Latched_wChargeSharing",iBin,iBin);
		
	// 	double xFakes_Sampled = findLastTimeAboveThreshold(hFakes_wChargeSharing , 0.0, 25,49.0);
	// 	double xFakes_Latched = findLastTimeAboveThreshold(hFakesLatched_wChargeSharing , 0.0, 25,49.0);
	// 	if( cThreshold%5 == 0 )
	// 	{
	// 		std::cout << "At a threshold of " << cThreshold << " .... last time w/ fakes > 0 perc. === " << xFakes_Sampled << " and " << xFakes_Latched << " ns.\n";
	// 		cCompareThresholds->cd(cPadCounter+1);
	// 		hEff_wChargeSharing->SetTitle(Form("Threshold = %d DAC units [%.1f sigma above pedestal]",cThreshold,cThreshold/(pNoise_kElectrons/cConversion_CBC3)));
	// 		hEff_wChargeSharing->SetStats(0);
	// 	    hEff_wChargeSharing->SetFillColor(kGreen);
	// 	    hEff_wChargeSharing->SetFillStyle(3004);
	// 	    hEff_wChargeSharing->GetXaxis()->SetRangeUser(25,50);
	// 	    hEff_wChargeSharing->DrawCopy("eHisto");
	// 	    hFakes_wChargeSharing->SetFillColor(kRed);
	// 	    hFakes_wChargeSharing->SetFillStyle(3004);
	// 	    hFakes_wChargeSharing->DrawCopy("eHistoSAME");
	// 		cPadCounter++;
	// 	}
	// 	cProfile_Sampled->Fill(cThreshold, 50.0 - xFakes_Sampled);
	// 	cProfile_Latched->Fill(cThreshold, 50.0 - xFakes_Latched);
	// }
	// // TFile* f_noChargeSharing = new TFile(pFileName1,"READ");
	// // f_noChargeSharing->cd();
	// // TH2D* hEff_Sampled_noChargeSharing=0;
	// // cHistName.Form("hFakeRate_SampledMode_ns");
	// // f_noChargeSharing->GetObject (cHistName.Data() , hEff_Sampled_noChargeSharing);
	// // hEff_Sampled_noChargeSharing->SetDirectory(0);  
	// // f_noChargeSharing->Close();


 //    TCanvas* c = new TCanvas("c","c",350*2,350*1);
 //    c->Divide(2,1);
 //    c->cd(1);
 //    cProfile_Sampled->SetStats(0);
 //    cProfile_Sampled->SetLineColor(kRed); cProfile_Sampled->DrawCopy("histo");
 //    cProfile_Latched->DrawCopy("histoSAME");
 //    // hEff_Sampled_wChargeSharing->SetStats(0);
 //    // hEff_Sampled_wChargeSharing->DrawCopy("colz");
 //    // c->cd(2);
 //    // hFakes_Sampled_wChargeSharing->SetStats(0);
 //    // hFakes_Sampled_wChargeSharing->DrawCopy("colz");
 //    c->cd(1);
    // hEff_wChargeSharing->SetStats(0);
    // hEff_wChargeSharing->SetFillColor(kGreen);
    // hEff_wChargeSharing->SetFillStyle(3004);
    // hEff_wChargeSharing->GetXaxis()->SetRangeUser(25,50);
    // hEff_wChargeSharing->DrawCopy("eHisto");
    // hFakes_wChargeSharing->SetFillColor(kRed);
    // hFakes_wChargeSharing->SetFillStyle(3004);
    // hFakes_wChargeSharing->DrawCopy("eHistoSAME");

    // c->cd(2);
    // hEffLatched_wChargeSharing->SetStats(0);
    // hEffLatched_wChargeSharing->SetFillColor(kGreen);
    // hEffLatched_wChargeSharing->SetFillStyle(3004);
    // hEffLatched_wChargeSharing->GetXaxis()->SetRangeUser(25,50);
    // hEffLatched_wChargeSharing->DrawCopy("eHisto");
    // hFakesLatched_wChargeSharing->SetFillColor(kRed);
    // hFakesLatched_wChargeSharing->SetFillStyle(3004);
    // hFakesLatched_wChargeSharing->DrawCopy("eHistoSAME");

    // c->cd(2);
    // hEff_Sampled_wChargeSharing->SetStats(0);
    // hEff_Sampled_wChargeSharing->DrawCopy("colz");
    //c->cd(1);
    //hEff_Sampled_noChargeSharing->SetStats(0);
    //hEff_Sampled_noChargeSharing->DrawCopy("colz");

}