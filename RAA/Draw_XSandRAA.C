//!!!! TGraphMultiErrors needs ROOTv20+, for example with cmsenv in ~/miniAOD/CMSSW_11_2_0_pre10

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphMultiErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

void Draw_XSandRAA(){
  gInterpreter->GenerateDictionary("vector<vector<vector<float> > >", "vector");

  vector<vector<vector<double> > > *x_LW_2ndstep; //[col][method][bin]
  TFile *infile3 = new TFile("../twoSteps/pTBiases.root","READ");
  infile3->GetObject("x_LW_correctedpTspectrum_2ndStep", x_LW_2ndstep);

  vector<vector<double> > *InvAccEff_pp;
  vector<vector<double> > *InvAccEff_PbPb;
  vector<float> *AEcorr12_pp;
  vector<float> *AEcorr12_PbPb;

  TFile *infile2 = new TFile("../twoSteps/AccEffFrom2ndStepToys.root","READ");
  infile2->GetObject("InvAccEffFromCorrMC_LinearisedCorrelationFactor_pp_2ndStep", AEcorr12_pp);
  infile2->GetObject("InvAccEffFromCorrMC_LinearisedCorrelationFactor_PbPb_2ndStep", AEcorr12_PbPb);
  infile2->GetObject("InvAccEffFromCorrMC_withSystErr_pp_2ndStep", InvAccEff_pp);
  infile2->GetObject("InvAccEffFromCorrMC_withSystErr_PbPb_2ndStep", InvAccEff_PbPb);

  vector<vector<vector<float> > > *Yields_postfit_pp;
  vector<vector<vector<float> > > *Yields_postfit_PbPb;
  vector<float> *y_nom_pp;
  vector<vector<float> > *y_fitErr_pp;
  vector<float> *TnPrelErr_pp;
  vector<float> *y_metafitRelErrLo_pp;
  vector<float> *y_metafitRelErrHi_pp;
  vector<float> *r1r2Corr_pp;
  vector<float> *metafitErrCorr_pp;
  vector<vector<float> > *rsig_relerr_pp;
  vector<float> *y_nom_PbPb;
  vector<vector<float> > *y_fitErr_PbPb;
  vector<float> *TnPrelErr_PbPb;
  vector<float> *y_metafitRelErrLo_PbPb;
  vector<float> *y_metafitRelErrHi_PbPb;
  vector<float> *r1r2Corr_PbPb;
  vector<float> *metafitErrCorr_PbPb;
  vector<vector<float> > *rsig_relerr_PbPb;
  vector<float> *y_metafitRelErrLo_RAA;
  vector<float> *y_metafitRelErrHi_RAA;
  vector<float> *metafitErrCorr_RAA;
  vector<float> *r1r2Corr_RAA;
  vector<float> *AcceffSystCorr_RAA;
  
  TFile *infile = new TFile("../AccEffCorr/corrected_yields_3rdStep.root","UPDATE");
  infile->GetObject("FinalCorrectedYield_pp", y_nom_pp);
  infile->GetObject("FinalCorrectedYield_fitError_pp", y_fitErr_pp);
  infile->GetObject("TagAndProbe_relError_pp", TnPrelErr_pp);
  infile->GetObject("r1r2Correlation_pp", r1r2Corr_pp);
  infile->GetObject("FinalCorrectedYield_PbPb", y_nom_PbPb);
  infile->GetObject("FinalCorrectedYield_fitError_PbPb", y_fitErr_PbPb);
  infile->GetObject("TagAndProbe_relError_PbPb", TnPrelErr_PbPb);
  infile->GetObject("r1r2Correlation_PbPb", r1r2Corr_PbPb);
  infile->GetObject("r1r2Correlation_RAA", r1r2Corr_RAA);
  infile->GetObject("AcceffSyst_Correlation_RAA", AcceffSystCorr_RAA);

  TFile *infile4 = new TFile("../AccEffCorr/corrected_yields_2ndStep.root","READ");
  infile4->GetObject("Yields_postfit_pp", Yields_postfit_pp);
  infile4->GetObject("Yields_postfit_PbPb", Yields_postfit_PbPb);
  infile4->GetObject("CorrectedYields_MetafitRelSystErrorLo_pp", y_metafitRelErrLo_pp); //or "CorrectedYield_MetafitRelSystErrorMaxDeviation_pp"
  infile4->GetObject("CorrectedYields_MetafitRelSystErrorHi_pp", y_metafitRelErrHi_pp); //or "CorrectedYield_MetafitRelSystErrorMaxDeviation_pp"
  infile4->GetObject("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixpp", metafitErrCorr_pp);
  infile4->GetObject("rsig_relerr_pp", rsig_relerr_pp);
  infile4->GetObject("CorrectedYields_MetafitRelSystErrorLo_PbPb", y_metafitRelErrLo_PbPb);
  infile4->GetObject("CorrectedYields_MetafitRelSystErrorHi_PbPb", y_metafitRelErrHi_PbPb);
  infile4->GetObject("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixPbPb", metafitErrCorr_PbPb);
  infile4->GetObject("rsig_relerr_PbPb", rsig_relerr_PbPb);
  infile4->GetObject("RAA_MetafitRelSystErrorLo_", y_metafitRelErrLo_RAA);
  infile4->GetObject("RAA_MetafitRelSystErrorHi_", y_metafitRelErrHi_RAA);
  infile4->GetObject("RAA_MetafitSyst_LinearizedCorrelationMatrix", metafitErrCorr_RAA);

  vector<vector<double> > nomiAndErr;

  float corrTnPerr = 0.7;

  const int nbins = _NanaBins;
  const int nYregions = 2;//hard-coded
  int nPts[nYregions+1] = {2,1,1};
  int nGr = 3; //pp,PbPb,RAA
  const int nMaxPts = 2;//hard-coded
  int nbins_test = 0;
  for(int b=1;b<=nYregions;b++) nbins_test += nPts[b];
  if(nbins_test!=nbins) cout<<"!!!!!!!!!!!! PROBLEM with number of analysis bins and/or number of Y regions ! Expected segfault"<<endl;

  double xErr[nGr][3][nYregions+1][nMaxPts], x[nGr][nYregions+1][nMaxPts];
  double y_metafitErr[nGr][3][nbins];
  double y_BcTauErr[nGr][3][nbins];
  double yErrPart[nGr][3][nbins];
  double yfitErr[nGr][3][nYregions+1][nMaxPts];
  double yacceffErr[nGr][3][nbins];
  double yTnPErr[nGr][nbins];
  double y[nGr][nYregions+1][nMaxPts], yErr[nGr][3][nYregions+1][nMaxPts];
  vector<vector<float> > corrtot(3, vector<float>((int)(nbins*(nbins-1)/2) , 0));

  for(int b=0;b<=nYregions;b++){ //b==0 is the graph with all points, without Y discrimination
    for(int m=0;m<nPts[b];m++){

      int trueb = (b==0)?m:(b-1);
      //x values //take LW prescription for non-centered x
      x[0][b][m] = (*x_LW_2ndstep)[0][1][trueb]; //take method 1: quadratic in log-log, fix m
      x[1][b][m] = (*x_LW_2ndstep)[1][1][trueb]; //take method 1: quadratic in log-log, fix m
      x[2][b][m] = (x[0][b][m]+x[1][b][m])/2; //average of pp and PbPb bin centers
      for(int igr=0;igr<3;igr++){
	xErr[igr][1][b][m] =  x[igr][b][m] - _BcPtmin[trueb+1];
	xErr[igr][2][b][m] = -x[igr][b][m] + _BcPtmax[trueb+1];      
      }

      //y values
      float normpp = L_pp * (_BcPtmax[trueb+1]-_BcPtmin[trueb+1]) * (_BcYmax[trueb+1]-_BcYmin[trueb+1]);
      float normPbPb = NMB_PbPb * TAA_090 * (_BcPtmax[trueb+1]-_BcPtmin[trueb+1]) * (_BcYmax[trueb+1]-_BcYmin[trueb+1]); //Leq_PbPb if we choose the lumi option
      y[0][b][m] = (*y_nom_pp)[trueb] / normpp;//(*Yields_postfit_pp)[0][trueb+1][0] * (*InvAccEff_pp)[trueb][0] / normpp; //
      y[1][b][m] = (*y_nom_PbPb)[trueb] / normPbPb;//(*Yields_postfit_PbPb)[0][trueb+1][0] * (*InvAccEff_PbPb)[trueb][0] / normPbPb;
      y[2][b][m] = y[1][b][m]/y[0][b][m];

      //metafit error
      y_metafitErr[0][1][trueb] = y[0][b][m] * (*y_metafitRelErrLo_pp)[trueb];
      y_metafitErr[0][2][trueb] = y[0][b][m] * (*y_metafitRelErrHi_pp)[trueb];
      y_metafitErr[1][1][trueb] = y[1][b][m] * (*y_metafitRelErrLo_PbPb)[trueb];
      y_metafitErr[1][2][trueb] = y[1][b][m] * (*y_metafitRelErrHi_PbPb)[trueb];
      y_metafitErr[2][1][trueb] = y[2][b][m] * (*y_metafitRelErrLo_RAA)[trueb];
      y_metafitErr[2][2][trueb] = y[2][b][m] * (*y_metafitRelErrHi_RAA)[trueb];
      for(int igr=0;igr<3;igr++)
	y_metafitErr[igr][0][trueb] = (y_metafitErr[igr][1][trueb] + y_metafitErr[igr][2][trueb])/2;

      //fit error
      for(int pm=0;pm<3;pm++){
	yfitErr[0][pm][b][m] = y[0][b][m] * (*rsig_relerr_pp)[trueb+1][pm];//(*y_fitErr_pp)[trueb][pm] / normpp;
	yfitErr[1][pm][b][m] = y[1][b][m] * (*rsig_relerr_PbPb)[trueb+1][pm];//(*y_fitErr_PbPb)[trueb][pm] / normPbPb;
	yfitErr[2][pm][b][m] = y[2][b][m] * sqrt(pow(yfitErr[0][pm][b][m]/y[0][b][m],2) + pow(yfitErr[1][pm][b][m]/y[1][b][m],2) );}

      //HERE put loop on pm for new AE error
      //acceff error
      for(int pm=0;pm<3;pm++){
	yacceffErr[0][pm][trueb] = y[0][b][m] * (*InvAccEff_pp)[trueb][pm+1] / (*InvAccEff_pp)[trueb][0];
	yacceffErr[1][pm][trueb] = y[1][b][m] * (*InvAccEff_PbPb)[trueb][pm+1] / (*InvAccEff_PbPb)[trueb][0];
	yacceffErr[2][pm][trueb] = y[2][b][m] * sqrt(pow(yacceffErr[0][pm][trueb]/y[0][b][m],2) + pow(yacceffErr[1][pm][trueb]/y[1][b][m],2));
      }

      //TnP error
      yTnPErr[0][trueb] = y[0][b][m] * (*TnPrelErr_pp)[trueb];
      yTnPErr[1][trueb] = y[1][b][m] * (*TnPrelErr_PbPb)[trueb];
      yTnPErr[2][trueb] = y[2][b][m] * sqrt(pow((*TnPrelErr_pp)[trueb],2) + pow((*TnPrelErr_PbPb)[trueb],2));

      //Bc to tau channel systematic
      for(int igr=0;igr<2;igr++){
	y_BcTauErr[igr][1][trueb] = _XS_BcTauRelSystLo * y[igr][b][m];
	y_BcTauErr[igr][2][trueb] = 0;
	y_BcTauErr[igr][0][trueb] = (y_BcTauErr[igr][1][trueb] + y_BcTauErr[igr][2][trueb])/2;
      }
      for(int pm=0;pm<3;pm++)
	y_BcTauErr[2][pm][trueb] = _RAA_BcTauRelSyst * y[2][b][m];

      //partial error, without metafit error
      for(int pm=0;pm<3;pm++){
	for(int igr=0;igr<2;igr++){
	  yErrPart[igr][pm][trueb] = sqrt(pow(yfitErr[igr][pm][b][m] ,2) + pow(yacceffErr[igr][pm][trueb] ,2) + pow(y_BcTauErr[igr][pm][trueb],2) + pow(yTnPErr[igr][trueb],2) );}
	yErrPart[2][pm][trueb] = y[2][b][m] * sqrt(pow(yErrPart[0][pm][trueb]/y[0][b][m],2) + pow(yErrPart[1][pm][trueb]/y[1][b][m],2) + _RAA_BcTauRelSyst ); //without metafit error	
      }

      //total error
      for(int pm=0;pm<3;pm++){
        for(int igr=0;igr<3;igr++){
	  yErr[igr][pm][b][m] = sqrt(pow(y_metafitErr[igr][pm][trueb] ,2) + pow(yErrPart[igr][pm][trueb] ,2));}
      }

      if(b==0){
	cout<<"point "<<m<<": x_LW pp, PbPb, RAA = "<<x[0][b][m]<<" "<<x[1][b][m]<<" "<<x[2][b][m]<<endl;
	cout<<"point "<<m<<": y_pp, y_PbPb, y_RAA = "<<y[0][b][m]<<" "<<y[1][b][m]<<" "<<y[2][b][m]<<endl;
	cout<<"point "<<m<<": low  errors on y_pp, y_PbPb, y_RAA = "<<yErr[0][1][b][m]<<" "<<yErr[1][1][b][m]<<" "<<yErr[2][1][b][m]<<endl;
	cout<<"point "<<m<<": high errors on y_pp, y_PbPb, y_RAA = "<<yErr[0][2][b][m]<<" "<<yErr[1][2][b][m]<<" "<<yErr[2][2][b][m]<<endl;
	cout<<"fit     rel errors on pp, PbPb, RAA = "<<yfitErr[0][0][b][m]/y[0][b][m]<<" "<<yfitErr[1][0][b][m]/y[1][b][m]<<" "<<yfitErr[2][0][b][m]/y[2][b][m]<<endl;
	cout<<"metafit rel errors on pp, PbPb, RAA = "<<y_metafitErr[0][0][trueb]/y[0][b][m]<<" "<<y_metafitErr[1][0][trueb]/y[1][b][m]<<" "<<y_metafitErr[2][0][trueb]/y[2][b][m]<<endl;
	cout<<"acceff  rel errors on pp, PbPb, RAA = "<<yacceffErr[0][0][trueb]/y[0][b][m]<<" "<<yacceffErr[1][0][trueb]/y[1][b][m]<<" "<<yacceffErr[2][0][trueb]/y[2][b][m]<<endl;
	cout<<"TnP  rel errors on pp, PbPb, RAA = "<<yTnPErr[0][trueb]/y[0][b][m]<<" "<<yTnPErr[1][trueb]/y[1][b][m]<<" "<<yTnPErr[2][trueb]/y[2][b][m]<<endl;
      }

    }
  }

  corrtot[0][0] = ( (*r1r2Corr_pp)[0] * yfitErr[0][0][0][0] * yfitErr[0][0][0][1]  
		 +(*metafitErrCorr_pp)[0] * y_metafitErr[0][0][0] * y_metafitErr[0][0][1] 
		    +(*AEcorr12_pp)[0] * yacceffErr[0][0][0] * yacceffErr[0][0][1]
		 +_corr_BcTauSyst * y_BcTauErr[0][0][0] * y_BcTauErr[0][0][1]
		 +corrTnPerr * yTnPErr[0][0] * yTnPErr[0][1]
		 )/( yErr[0][0][0][0]*yErr[0][0][0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  corrtot[1][0] = ( (*r1r2Corr_PbPb)[0] * yfitErr[1][0][0][0] * yfitErr[1][0][0][1]  
		 +(*metafitErrCorr_PbPb)[0] * y_metafitErr[1][0][0] * y_metafitErr[1][0][1]
		    +(*AEcorr12_PbPb)[0] * yacceffErr[1][0][0] * yacceffErr[1][0][1]
		 +_corr_BcTauSyst * y_BcTauErr[1][0][0] * y_BcTauErr[1][0][1]
		 +corrTnPerr * yTnPErr[0][0] * yTnPErr[0][1]
		 )/( yErr[1][0][0][0]*yErr[1][0][0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  corrtot[2][0] = ( (*r1r2Corr_RAA)[0] * yfitErr[2][0][0][0] * yfitErr[2][0][0][1]  
		 +(*metafitErrCorr_RAA)[0] * y_metafitErr[2][0][0] * y_metafitErr[2][0][1] 
		 +(*AcceffSystCorr_RAA)[0] * yacceffErr[2][0][0] * yacceffErr[2][0][1]
		 +_corr_BcTauSyst * y_BcTauErr[2][0][0] * y_BcTauErr[2][0][1]
		 +corrTnPerr * yTnPErr[0][0] * yTnPErr[0][1]
		 )/( yErr[2][0][0][0]*yErr[2][0][0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  cout<<"correlation factor between 2 pT bins for pp, from fit, metafit, acceff = "<<(*r1r2Corr_pp)[0]<<" "<<(*metafitErrCorr_pp)[0]<<" "<<(*AEcorr12_pp)[0]<<endl;
  cout<<"correlation factor between 2 pT bins for PbPb, from fit, metafit, acceff = "<<(*r1r2Corr_PbPb)[0]<<" "<<(*metafitErrCorr_PbPb)[0]<<" "<<(*AEcorr12_PbPb)[0]<<endl;
  cout<<"correlation factor between 2 pT bins for RAA, from fit, metafit, acceff = "<<(*r1r2Corr_RAA)[0]<<" "<<(*metafitErrCorr_RAA)[0]<<" "<<(*AcceffSystCorr_RAA)[0]<<endl;
  cout<<"correlation factor between 2 pT bins, for pp, PbPb, RAA = "<<corrtot[0][0]<<" "<<corrtot[1][0]<<" "<<corrtot[2][0]<<endl;

  //Record final results
  vector<vector<vector<float> > > XSRAA(3, vector<vector<float> >(4, vector<float>(nbins, 0))); //dimension1: pp,PbPb,RAA. dimension2: value, errlo, errhi, errsym. dimension3: pt bins
  vector<vector<vector<float> > > metafitRelerr(3, vector<vector<float> >(3, vector<float>(nbins, 0)));//same, without value for dimension2
  for(int b=0;b<nbins;b++){
    for(int p=0;p<3;p++){
      XSRAA[p][0][b] = y[p][0][b];
      XSRAA[p][1][b] = yErr[p][1][0][b];
      XSRAA[p][2][b] = yErr[p][2][0][b];
      XSRAA[p][3][b] = yErr[p][0][0][b];
    }
    metafitRelerr[0][1][b] = (*y_metafitRelErrLo_pp)[b];
    metafitRelerr[0][2][b] = (*y_metafitRelErrHi_pp)[b];
    metafitRelerr[0][0][b] = (metafitRelerr[0][1][b]+metafitRelerr[0][2][b])/2;
    metafitRelerr[1][1][b] = (*y_metafitRelErrLo_PbPb)[b];
    metafitRelerr[1][2][b] = (*y_metafitRelErrHi_PbPb)[b];
    metafitRelerr[1][0][b] = (metafitRelerr[1][1][b]+metafitRelerr[1][2][b])/2;
    metafitRelerr[2][1][b] = (*y_metafitRelErrLo_RAA)[b];
    metafitRelerr[2][2][b] = (*y_metafitRelErrHi_RAA)[b];
    metafitRelerr[2][0][b] = (metafitRelerr[2][1][b]+metafitRelerr[2][2][b])/2;
  }
  infile->WriteObject(&XSRAA,"XSRAA_final");
  infile->WriteObject(&metafitRelerr,"MetaFit_RelErr");
  infile->WriteObject(&corrtot,"LinearizedCorrelationMatrix_total");

  //graph styles
  gStyle->SetEndErrorSize(15);//size in pixels
  vector<vector<Style_t> > Mstyle = {{20,20,89} , {21,21,90} , {20,20,89}}; //89,90: empty markers, line width 4 (go to 107,108 for width 5)
  vector<Color_t> Mcol = {kBlue+2  , kSpring-7 , kAzure+4};
  vector<Color_t> Lcol = {kAzure-2 , kSpring-6 , kAzure+5};
  vector<Color_t> Fcol = {kCyan-6  , kSpring+9 , kCyan};
  vector<Width_t> Lwidth = {3,3,4};

  vector<vector<TGraphMultiErrors*> > gr(nGr);
  for(int igr=0;igr<3;igr++){
    for(int b=0;b<=nYregions;b++){ 
      gr[igr].push_back( new TGraphMultiErrors(nPts[b], x[igr][b],y[igr][b], xErr[igr][1][b], xErr[igr][2][b],yErr[igr][1][b], yErr[igr][2][b]) );
      gr[igr][b]->AddYError(nPts[b], yfitErr[igr][1][b], yfitErr[igr][2][b]);
      gr[igr][b]->AddYError(nPts[b], yErr[igr][1][b], yErr[igr][2][b]);
      gr[igr][b]->AddYError(nPts[b], yfitErr[igr][1][b], yfitErr[igr][2][b]);

      gr[igr][b]->SetMarkerColor(Mcol[igr]);
      gr[igr][b]->SetMarkerSize(5);
      gr[igr][b]->SetMarkerStyle(Mstyle[igr][b]);
      gr[igr][b]->SetLineColor(Lcol[igr]);
      gr[igr][b]->SetLineWidth(Lwidth[igr]);
      gr[igr][b]->SetFillColor(Fcol[igr]); 
      gr[igr][b]->SetFillStyle(1001);
      for(int p=0;p<4;p++){
	gr[igr][b]->GetAttLine(p)->SetLineColor(Lcol[igr]);
	gr[igr][b]->GetAttLine(p)->SetLineWidth(Lwidth[igr]);
      }

      gr[igr][b]->GetAttFill(0)->SetFillStyle(3001); //full error, draw boxes
      gr[igr][b]->GetAttFill(0)->SetFillColor(Fcol[igr]); 
      gr[igr][b]->GetAttFill(1)->SetFillStyle(1001); //fit error only, draw boxes
      gr[igr][b]->GetAttFill(1)->SetFillColor(Fcol[igr]); 
      gr[igr][b]->GetAttLine(2)->SetLineStyle(2);
      gr[igr][b]->GetAttLine(3)->SetLineStyle(1);
      
    }
  }

  //Draw XS and RAA
  TCanvas *c1 = new TCanvas("c1","XS",2000,2000);
  c1->SetLeftMargin(0.14);
  c1->SetTopMargin(0.11);
  gr[0][0]->SetTitle("Cross-sections");
  gr[0][0]->GetYaxis()->SetRangeUser(0.5*y[1][0][1],2*y[1][0][0]);
  gr[0][0]->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);
  gr[0][0]->GetYaxis()->SetTitleOffset(1.4);
  gr[0][0]->GetXaxis()->SetTitle("p_{T}^{vis}(B_{c}) [GeV]");
  gr[0][0]->GetYaxis()->SetTitle("BF  #times  #left(#frac{d#sigma_{pp}}{dp_{T}^{vis} dy^{vis}}    or    #frac{1}{N_{MB} T_{PbPb}} #frac{dN_{PbPb}^{corr}}{dp_{T}^{vis} dy^{vis}}#right)  [pb/GeV]");
  gr[0][0]->SetMarkerColor(kWhite);
  gr[0][0]->SetLineColor(kWhite);
  gr[1][0]->Draw("APX");
  gr[1][1]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
  gr[1][2]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
  gr[0][1]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
  gr[0][2]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"

  vector<TLegend*> leg;
  for(int b=0;b<nYregions;b++){
    leg.push_back(new TLegend((b==0)?0.48:0.725,0.62,(b==0)?0.67:0.89,0.765));
    leg[b]->AddEntry(gr[0][b+1],(b==0)?"pp":"pp");
    leg[b]->AddEntry(gr[1][b+1],(b==0)?"PbPb":"PbPb");
    leg[b]->SetTextSize(0.039);
    leg[b]->SetBorderSize(0);
    leg[b]->Draw("same");
  }
  TLatex legEntry;
  legEntry.SetNDC();
  legEntry.SetTextFont(42);
  legEntry.SetTextSize(0.035);
  legEntry.DrawLatex(0.45,0.77,"1.3<|Y^{vis}|<2.3");
  legEntry.DrawLatex(0.71,0.77,"0<|Y^{vis}|<2.3");
  legEntry.DrawLatex(0.65,0.49,Form("#rho_{1-2}^{pp} = %.2f",corrtot[0][0]));
  legEntry.DrawLatex(0.65,0.42,Form("#rho_{1-2}^{PbPb} = %.2f",corrtot[1][0]));
  legEntry.DrawLatex(0.55,0.84,"Centrality 0-90%");

  gPad->SetLogy();

  //DRAW CMS Preliminary                                                                                                                                                                                                 
  TLatex CMStag;
  CMStag.SetNDC();
  CMStag.SetTextFont(42);
  CMStag.SetTextSize(0.035);
  CMStag.DrawLatex(0.12,0.9,"#font[61]{CMS} #font[52]{Preliminary}");
  CMStag.DrawLatex(0.4,0.9, Form("pp (%.0f pb^{-1}), PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));

  c1->SaveAs("CrossSections.pdf");

  TCanvas *c2 = new TCanvas("c2","RAA",2000,2000);
  //c2->SetLeftMargin(0.14);
  c2->SetTopMargin(0.11);
  gr[2][0]->SetTitle("R_{PbPb}");
  gr[2][0]->SetMarkerColor(kWhite);
  gr[2][0]->SetLineColor(kWhite);
  gr[2][0]->GetYaxis()->SetTitleOffset(1.2);
  gr[2][0]->GetYaxis()->SetRangeUser(0,2.5);
  gr[2][0]->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);
  gr[2][0]->GetXaxis()->SetTitle("p_{T}^{vis}(B_{c}) [GeV]");
  gr[2][0]->GetYaxis()->SetTitle("R_{PbPb}(B_{c})");
  gr[2][0]->Draw("AP");
  gr[2][1]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
  gr[2][2]->Draw("PX s same; 2; 2; P s=0; P s=0");

  TLegend *leg2 = new TLegend(0.5,0.64,0.8,0.8);
  leg2->AddEntry(gr[2][1], "1.3<|Y^{vis}|<2.3");
  leg2->AddEntry(gr[2][2], "0<|Y^{vis}|<2.3");
  leg2->SetTextSize(0.039);
  leg2->SetBorderSize(0);
  leg2->Draw("same");

  legEntry.DrawLatex(0.65,0.53,Form("#rho_{1-2} = %.2f",corrtot[2][0]));
  legEntry.DrawLatex(0.5,0.83,"Centrality 0-90%");

  CMStag.DrawLatex(0.1,0.9,"#font[61]{CMS} #font[52]{Preliminary}");
  CMStag.DrawLatex(0.4,0.9, Form("pp (%.0f pb^{-1}), PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));

  TLine *one = new TLine(_BcPtmin[0],1,_BcPtmax[0],1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->Draw("same");

  c2->SaveAs("RAA.pdf");

  ///////////////////////////////////////
  // Draw all systematics
  ///////////////////////////////////////
  vector<TH1F*> S(6);
  vector<TString> name = {"fit","metafit","acceff","TnP","BcTau","total"};
  vector<TString> prettyname = {"fit","fit method","Acc#times Eff","Tag-and-probe","B_{c}#rightarrow J/#psi #tau","total"};
  vector<Color_t> col = {kRed,kBlue,kGreen+1,kOrange-6,kMagenta-4, kBlack};
  
  for(int i=0;i<6;i++)
    S[i] = new TH1F(name[i],name[i],6,0,6);

  S[0]->SetBinContent(1,yfitErr[0][0][0][0]/y[0][0][0]);//take average of hi and lo errors
  S[0]->SetBinContent(2,yfitErr[0][0][0][1]/y[0][0][1]);
  S[0]->SetBinContent(3,yfitErr[1][0][0][0]/y[1][0][0]);
  S[0]->SetBinContent(4,yfitErr[1][0][0][1]/y[1][0][1]);
  S[0]->SetBinContent(5,yfitErr[2][0][0][0]/y[2][0][0]);
  S[0]->SetBinContent(6,yfitErr[2][0][0][1]/y[2][0][1]);

  S[1]->SetBinContent(1,y_metafitErr[0][0][0]/y[0][0][0]);//take average of hi and lo errors
  S[1]->SetBinContent(2,y_metafitErr[0][0][1]/y[0][0][1]);
  S[1]->SetBinContent(3,y_metafitErr[1][0][0]/y[1][0][0]);
  S[1]->SetBinContent(4,y_metafitErr[1][0][1]/y[1][0][1]);
  S[1]->SetBinContent(5,y_metafitErr[2][0][0]/y[2][0][0]);
  S[1]->SetBinContent(6,y_metafitErr[2][0][1]/y[2][0][1]);

  S[2]->SetBinContent(1,yacceffErr[0][0][0]/y[0][0][0]);
  S[2]->SetBinContent(2,yacceffErr[0][0][1]/y[0][0][1]);
  S[2]->SetBinContent(3,yacceffErr[1][0][0]/y[1][0][0]);
  S[2]->SetBinContent(4,yacceffErr[1][0][1]/y[1][0][1]);
  S[2]->SetBinContent(5,yacceffErr[2][0][0]/y[2][0][0]);
  S[2]->SetBinContent(6,yacceffErr[2][0][1]/y[2][0][1]);

  S[3]->SetBinContent(1,yTnPErr[0][0]/y[0][0][0]);
  S[3]->SetBinContent(2,yTnPErr[0][1]/y[0][0][1]);
  S[3]->SetBinContent(3,yTnPErr[1][0]/y[1][0][0]);
  S[3]->SetBinContent(4,yTnPErr[1][1]/y[1][0][1]);
  S[3]->SetBinContent(5,yTnPErr[2][0]/y[2][0][0]);
  S[3]->SetBinContent(6,yTnPErr[2][1]/y[2][0][1]);

  S[4]->SetBinContent(1,y_BcTauErr[0][1][0]/y[0][0][0]); //take the lower error here
  S[4]->SetBinContent(2,y_BcTauErr[0][1][1]/y[0][0][1]);
  S[4]->SetBinContent(3,y_BcTauErr[1][1][0]/y[1][0][0]);
  S[4]->SetBinContent(4,y_BcTauErr[1][1][1]/y[1][0][1]);
  S[4]->SetBinContent(5,y_BcTauErr[2][1][0]/y[2][0][0]);
  S[4]->SetBinContent(6,y_BcTauErr[2][1][1]/y[2][0][1]);

  S[5]->SetBinContent(1,yErr[0][0][0][0]/y[0][0][0]);//take average of hi and lo errors
  S[5]->SetBinContent(2,yErr[0][0][0][1]/y[0][0][1]);
  S[5]->SetBinContent(3,yErr[1][0][0][0]/y[1][0][0]);
  S[5]->SetBinContent(4,yErr[1][0][0][1]/y[1][0][1]);
  S[5]->SetBinContent(5,yErr[2][0][0][0]/y[2][0][0]);
  S[5]->SetBinContent(6,yErr[2][0][0][1]/y[2][0][1]);

  for(int i=0;i<6;i++){
    S[i]->SetLineColor(col[i]);
    S[i]->SetLineWidth(2);
    S[i]->GetXaxis()->SetBinLabel(1,"pp bin1");
    S[i]->GetXaxis()->SetBinLabel(2,"pp bin2");
    S[i]->GetXaxis()->SetBinLabel(3,"PbPb bin1");
    S[i]->GetXaxis()->SetBinLabel(4,"PbPb bin2");
    S[i]->GetXaxis()->SetBinLabel(5,"R_{PbPb} bin1");
    S[i]->GetXaxis()->SetBinLabel(6,"R_{PbPb} bin2");
    S[i]->GetXaxis()->SetLabelSize(0.055);
    S[i]->GetYaxis()->SetTitle("relative error");
    S[i]->GetYaxis()->SetTitleSize(0.055);
    S[i]->GetYaxis()->SetLabelSize(0.04);
    S[i]->GetYaxis()->SetTitleOffset(1.);
    S[i]->GetYaxis()->SetRangeUser(0,1.15*S[5]->GetMaximum());
  }
  
  gStyle->SetOptStat(0);

  TCanvas *c3 = new TCanvas("c3","Syst",2000,1500);
  c3->SetBottomMargin(0.1);
  c3->SetTopMargin(0.04);
  c3->SetRightMargin(0.04);
  c3->SetLeftMargin(0.12);
  S[0]->SetTitle("");
  for(int i=0;i<6;i++)
    S[i]->Draw((TString)((i==0)?"hist":"histsame"));

  TLegend *leg3 = new TLegend(0.15,0.5,0.37,0.95);
  for(int i=0;i<6;i++)
    leg3->AddEntry(S[i], prettyname[i]);
  leg3->SetHeader("Uncertainty sources");
  leg3->SetTextSize(0.043);
  leg3->SetBorderSize(0);
  leg3->Draw("same");

  c3->SaveAs("SystematicsSummary.pdf");

  infile->Close();      
  infile2->Close();      
  infile3->Close();      
}
