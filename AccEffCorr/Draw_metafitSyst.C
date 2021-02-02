#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"
#include "../helpers/SystematicsList.h"

void Draw_metafitSyst(bool firstStep=true){
  
  bool averageNomError = false; //for the new nominal (average over some fit methods), take the average of the fit errors too
  bool doAverageNominal = true;
  bool RMSbyBlocks = true;
  Color_t cols[4] = {1,634,810,882};//kBlack,kRed+2,kOrange+10,kViolet+2//,kBlue+2 602;

  vector<vector<vector<float> > > *Yields_postfit[2][_nMetafitSyst];
  vector<vector<float> > *rsig_relerr[3][_nMetafitSyst];
  vector<vector<float> > *eff_oneBinned[3];
  vector<float> *acc_oneBinned[3];

  vector<vector<double> > *InvAccEff_pp;
  vector<vector<double> > *InvAccEff_PbPb;
  TFile *infile4 = new TFile("../twoSteps/AccEffFrom2ndStepToys.root","READ");
  infile4->GetObject("InvAccEffFrom2ndStep_withSystErr_pp", InvAccEff_pp);
  infile4->GetObject("InvAccEffFrom2ndStep_withSystErr_PbPb", InvAccEff_PbPb);

  vector<float> y_nom[3];
  vector<float> y_systErrLo[3];
  vector<vector<float> > y_systErrLo_maxDev[3];
  vector<float> y_systErrHi[3];
  vector<vector<float> > y_systErrHi_maxDev[3];
  vector<float> y_systCorr[3];
  //vector<float> y_fitErr(nbins);

  TFile *infile = new TFile("corrected_yields.root","READ");

  const int nbins = _NanaBins;
  vector<TGraphAsymmErrors*> gSyst[3];
  TGraphAsymmErrors* gSystRatio[3][_nMetafitSyst+1][nbins];
  //vector<vector<TMultiGraph*> > gSystRatio_col;//[3][nbins];
  TMultiGraph* gSystRatio_all[3][nbins];
  double xErrLo[_nMetafitSyst+1][nbins], xErrHi[_nMetafitSyst+1][nbins], x[_nMetafitSyst+1][nbins], zero[nbins], newy[1];
  double yErrLo[3][_nMetafitSyst+1][nbins], yErrHi[3][_nMetafitSyst+1][nbins], yErrLo_ratio[3][_nMetafitSyst+1][nbins], yErrHi_ratio[3][_nMetafitSyst+1][nbins];
  double y_oneBinned[3][_nMetafitSyst+1][nbins], y_ratio[3][_nMetafitSyst+1][nbins]; //oneBinned is actually corrected by the 2-steps Acc x Eff, if firstStep==false

  //DRAW CMS Preliminary 
  TLatex CMStag;
  CMStag.SetNDC();
  CMStag.SetTextFont(42);
  CMStag.SetTextSize(0.035);

  for(int ispp=0;ispp<3;ispp++){
    cout<<"\n    ***********    "<<((ispp==1)?"pp":((ispp==0)?"PbPb":"RAA"))<<"    *************"<<endl;

    //gSystRatio_col.push_back(vector<TMultiGraph*>());

    TFile *infile3, *infile2;
    if(ispp<2){
      infile3 = new TFile("../acceptance/acceptanceMap.root","READ");
      infile3->GetObject("acceptance_oneBinned", acc_oneBinned[ispp]);

      infile2 = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
      infile2->GetObject("efficiency_oneBinned"+(TString)((ispp==1)?"_pp":"_PbPb"), eff_oneBinned[ispp]);
    }
  
    float nAvNominal[nbins];
    float ynom[nbins];
    float relerrLoNom[nbins];
    float relerrHiNom[nbins];
    for(int b=1;b<=nbins;b++){
      nAvNominal[b-1] = 0; ynom[b-1] = 0;
      relerrLoNom[b-1] = 0; relerrHiNom[b-1] = 0;}
    
    //x (pt) and error
    for(int b=1;b<=nbins;b++){
      zero[b-1] = 0;      
      gSystRatio_all[ispp][b-1] = new TMultiGraph();
      for(int isys=0;isys<=_nMetafitSyst;isys++){
	float xErrBase = (-_BcPtmin[b]+_BcPtmax[b])/2;
	float xBase = (_BcPtmin[b]+_BcPtmax[b])/2;
	float xcor = ((b==1)?1.:0.65)*xErrBase*(isys-(_nMetafitSyst+1)/2)/(_nMetafitSyst+1);
	x[isys][b-1] = xBase + xcor;
	xErrLo[isys][b-1] = xErrBase + xcor;
	xErrHi[isys][b-1] = xErrBase - xcor;
      }
    }

    //calculate new nominal (average over some methods)
    for(int isys=0;isys<_nMetafitSyst;isys++){
      if(ispp<2){
	infile->GetObject("Yields_postfit"+(TString)((ispp==1)?"_pp":"_PbPb")+systName[isys]+systExtName[isys], Yields_postfit[ispp][isys]);
	infile->GetObject("rsig_relerr"+(TString)((ispp==1)?"_pp":"_PbPb")+systName[isys]+systExtName[isys], rsig_relerr[ispp][isys]);

	for(int b=1;b<=nbins;b++){
	  //cout<<"analysis bin #"<<b<<endl;
	  //      cout<<"one-binned acceptance, efficiency = "<<(*acc_oneBinned)[b]<<" "<<(*eff_oneBinned)[b][0]<<endl;
	  float invAE;
	  if(firstStep) invAE = 1 / (((*acc_oneBinned[ispp])[b] * (*eff_oneBinned[ispp])[b][0]));
	  else invAE = (*(ispp?InvAccEff_pp:InvAccEff_PbPb))[b-1][0];
	  y_oneBinned[ispp][isys][b-1] = (*Yields_postfit[ispp][isys])[0][b][0] * invAE;
	  yErrLo[ispp][isys][b-1] = (*rsig_relerr[ispp][isys])[b][1] * y_oneBinned[ispp][isys][b-1];
	  yErrHi[ispp][isys][b-1] = (*rsig_relerr[ispp][isys])[b][2] * y_oneBinned[ispp][isys][b-1];

	  if(usedInNominalAverage[isys]>0.01){
	    ynom[b-1] += usedInNominalAverage[isys] * y_oneBinned[ispp][isys][b-1]; 
	    relerrLoNom[b-1] += usedInNominalAverage[isys] * (*rsig_relerr[ispp][isys])[b][1];
	    relerrHiNom[b-1] += usedInNominalAverage[isys] * (*rsig_relerr[ispp][isys])[b][2];
	    nAvNominal[b-1] += usedInNominalAverage[isys];
          }
	}
      }
    }
    
    for(int b=1;b<=nbins;b++){
      if(ispp<2){
	y_oneBinned[ispp][_nMetafitSyst][b-1] = ynom[b-1]/nAvNominal[b-1];
	if(!doAverageNominal) y_oneBinned[ispp][_nMetafitSyst][b-1] = y_oneBinned[ispp][0][b-1];
	yErrLo[ispp][_nMetafitSyst][b-1] = -y_oneBinned[ispp][_nMetafitSyst][b-1] * relerrLoNom[b]/nAvNominal[b-1];
	yErrHi[ispp][_nMetafitSyst][b-1] = y_oneBinned[ispp][_nMetafitSyst][b-1] * relerrHiNom[b]/nAvNominal[b-1];
	if(!averageNomError){
	  yErrLo[ispp][_nMetafitSyst][b-1] = y_oneBinned[ispp][_nMetafitSyst][b-1] * yErrLo[ispp][0][b-1]/y_oneBinned[ispp][0][b-1]; //keep the relative error from old nominal fit
	  yErrHi[ispp][_nMetafitSyst][b-1] = y_oneBinned[ispp][_nMetafitSyst][b-1] * yErrHi[ispp][0][b-1]/y_oneBinned[ispp][0][b-1]; //keep the relative error from old nominal fit
	}
      }
    }

    for(int isys=0;isys<=_nMetafitSyst;isys++){
      cout<<"  *** Building graph for systematic: "<<systPrettyName[isys]<<endl;

      //ratio to new nominal (average)
      for(int b=1;b<=nbins;b++){
	if(ispp<2){
	  y_ratio[ispp][isys][b-1] = y_oneBinned[ispp][isys][b-1] / y_oneBinned[ispp][_nMetafitSyst][b-1];
	  yErrLo_ratio[ispp][isys][b-1] = yErrLo[ispp][isys][b-1] / y_oneBinned[ispp][_nMetafitSyst][b-1];
	  yErrHi_ratio[ispp][isys][b-1] = yErrHi[ispp][isys][b-1] / y_oneBinned[ispp][_nMetafitSyst][b-1];
	}
	else{
	  y_oneBinned[ispp][isys][b-1] = y_oneBinned[0][isys][b-1] / y_oneBinned[1][isys][b-1];	
	  y_ratio[ispp][isys][b-1] = y_ratio[0][isys][b-1] / y_ratio[1][isys][b-1]; // = y_oneBinned[2][isys][b-1]/y_oneBinned[2][_nMetafitSyst][b-1] -- if the average (idx _nMetafitSyst) was already calculated
	  yErrLo_ratio[ispp][isys][b-1] = y_ratio[ispp][isys][b-1] * sqrt(pow(yErrLo[0][isys][b-1]/y_oneBinned[0][isys][b-1],2) + pow(yErrLo[1][isys][b-1]/y_oneBinned[1][isys][b-1],2));
	  yErrHi_ratio[ispp][isys][b-1] = y_ratio[ispp][isys][b-1] * sqrt(pow(yErrHi[0][isys][b-1]/y_oneBinned[0][isys][b-1],2) + pow(yErrHi[1][isys][b-1]/y_oneBinned[1][isys][b-1],2));
	  yErrLo[ispp][isys][b-1] = (y_oneBinned[ispp][isys][b-1] / y_ratio[ispp][isys][b-1]) * yErrLo_ratio[ispp][isys][b-1];
	  yErrHi[ispp][isys][b-1] = (y_oneBinned[ispp][isys][b-1] / y_ratio[ispp][isys][b-1]) * yErrHi_ratio[ispp][isys][b-1];
	}
	cout<<"b, y, yErr (+,-) = "<<b<<" "<<y_oneBinned[ispp][isys][b-1]<<" +"<<yErrHi[ispp][isys][b-1]<<" -"<<yErrLo[ispp][isys][b-1]<<endl;

	double ynew[] = { ynew[0] = isys+0.5+1 };
	if(isys==_nMetafitSyst) ynew[0] = 0.5;

	//make graphs
	gSystRatio[ispp][isys][b-1] = new TGraphAsymmErrors(1, &(y_ratio[ispp][isys][b-1]), ynew, &(yErrLo_ratio[ispp][isys][b-1]), &(yErrHi_ratio[ispp][isys][b-1]), zero, zero);

	Color_t col;
	if(isys==_nMetafitSyst) col = cols[0];
	else if(usedInNominalAverage[isys]) col = cols[1];
	else if(usedInRMS[isys]) col = cols[2];
	else col = cols[3];

	gSystRatio[ispp][isys][b-1]->SetMarkerColor(col);
	gSystRatio[ispp][isys][b-1]->SetLineColor(col);
	gSystRatio[ispp][isys][b-1]->SetMarkerSize(3);
	gSystRatio[ispp][isys][b-1]->SetLineWidth(3);
	gSystRatio[ispp][isys][b-1]->SetMarkerStyle(systStyle[_nMetafitSyst]);

	gSystRatio_all[ispp][b-1]->Add(gSystRatio[ispp][isys][b-1]);
      }

      gSyst[ispp].push_back( new TGraphAsymmErrors(nbins, x[isys], y_oneBinned[ispp][isys], xErrLo[isys], xErrHi[isys], yErrLo[ispp][isys], yErrHi[ispp][isys]) );

      gSyst[ispp][isys]->SetMarkerSize(4);
      gSyst[ispp][isys]->SetMarkerStyle(systStyle[isys]);
      gSyst[ispp][isys]->SetMarkerColor(systCol[isys]);
      gSyst[ispp][isys]->GetYaxis()->SetRangeUser(0.7*y_oneBinned[ispp][0][1],(ispp?1.4:2.)*y_oneBinned[ispp][0][0]);
      gSyst[ispp][isys]->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);//+3*yErr[0]);
      gSyst[ispp][isys]->SetTitle("Metafit systematic variations of "+(TString)((ispp==2)?"R_{PbPb}":"corrected yields"));
      gSyst[ispp][isys]->GetYaxis()->SetTitle((ispp==2)?"R_{PbPb}":"corrected yield");
      gSyst[ispp][isys]->GetXaxis()->SetTitle("p_{T}(#mu#mu#mu) [GeV]");
    }

    //Draw corrected yields for pp and PbPb, for each systematic
    if(ispp<2){
      TCanvas *c1 = new TCanvas("c1"+(TString)ispp,"c1",2000,2000);

      TLegend *leg = new TLegend(0.37,0.4,0.9,0.89);
      for(int isys=0;isys<=_nMetafitSyst;isys++){
    	leg->AddEntry(gSyst[ispp][isys],systPrettyName[isys],"lp");
      }
      leg->SetTextSize(0.023);
      leg->SetBorderSize(0);

      gSyst[ispp][0]->Draw("AP");
      for(int isys=0;isys<=_nMetafitSyst;isys++){
    	gSyst[ispp][isys]->Draw((isys==_nMetafitSyst)?"Psame":"Psame");
      }
      leg->Draw("same");
      if(!ispp) gPad->SetLogy();

      CMStag.DrawLatex(0.19,0.85,"#font[61]{CMS "+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":""))+"}");
      CMStag.DrawLatex(0.19,0.80,"#font[52]{Preliminary}");

      c1->SetLeftMargin(0.15);
      c1->SaveAs("figs/correctedYields_metafitSyst"+(TString)(firstStep?"_1stStep":"_2ndStep")+(TString)((ispp==1)?"_pp":((ispp==0)?"_PbPb":""))+".pdf");
      c1->SaveAs("figs/correctedYields_metafitSyst"+(TString)(firstStep?"_1stStep":"_2ndStep")+(TString)((ispp==1)?"_pp":((ispp==0)?"_PbPb":""))+".png");
    }

    //******************************************************
    //Record the (nominal result +) systematic error
    //******************************************************

    for(int b=0;b<nbins;b++){
      y_nom[ispp].push_back( y_oneBinned[ispp][_nMetafitSyst][b] );
      y_systErrLo[ispp].push_back( 0 ); //rms of the deviations to nominal y
      y_systErrHi[ispp].push_back( 0 );
      y_systErrLo_maxDev[ispp].push_back( vector<float>(_nRMSblocks+1,0) );//max deviation from nominal
      y_systErrHi_maxDev[ispp].push_back( vector<float>(_nRMSblocks+1,0) );
      float nVarUsed = 0;
      int iblock = 1;
      for(int isys=0;isys<_nMetafitSyst;isys++){ //include old nominal (0) in rms
	cout<<"b, isys, name, yRatioToNominal, yerrlo, yerrhi = "<<b<<" "<<isys<<" "<<systPrettyName[isys]<<" "<<y_ratio[ispp][isys][b]<<" -"<<yErrLo_ratio[ispp][isys][b]<<" +"<<yErrHi_ratio[ispp][isys][b]<<endl;
	if(RMSblock[isys]>0 || (doAverageNominal && isys==0)){
	  if(RMSblock[isys]>iblock) iblock += 1; //summing the deviation to the next block
	  y_systErrLo_maxDev[ispp][b][0] = max(y_systErrLo[ispp][b] , (float)max(0.,-(y_ratio[ispp][isys][b] - 1)) ); //for maximum deviation, take asymmetric error (ErrLo counts only downwards variations)
	  y_systErrHi_maxDev[ispp][b][0] = max(y_systErrHi[ispp][b] , (float)max(0., (y_ratio[ispp][isys][b] - 1)) );
	  y_systErrLo_maxDev[ispp][b][iblock] = min((double)y_systErrLo_maxDev[ispp][b][iblock] , y_ratio[ispp][isys][b] - 1 );
	  y_systErrHi_maxDev[ispp][b][iblock] = max((double)y_systErrHi_maxDev[ispp][b][iblock] , y_ratio[ispp][isys][b] - 1 );
	  if(!RMSbyBlocks){
	    y_systErrLo[ispp][b] += pow(y_ratio[ispp][isys][b] - 1 ,2); //symmetric systematic error here
	    y_systErrHi[ispp][b] += pow(y_ratio[ispp][isys][b] - 1 ,2);
	  }
	  else if(RMSblock[isys+1]>iblock || RMSblock[isys+1]==0){
	    y_systErrLo[ispp][b] += pow(max(fabs(y_systErrLo_maxDev[ispp][b][iblock]),y_systErrHi_maxDev[ispp][b][iblock]) ,2); //symmetric systematic error here
	    y_systErrHi[ispp][b] += pow(max(fabs(y_systErrLo_maxDev[ispp][b][iblock]),y_systErrHi_maxDev[ispp][b][iblock]) ,2);
	  }
	  nVarUsed += 1;
	}
	else{
	  bool varLowerThanNom = y_ratio[ispp][isys][b] < 1;
	  float yerrVar = varLowerThanNom?(yErrHi_ratio[ispp][isys][b]):(yErrLo_ratio[ispp][isys][b]);
	  float yerrNom = varLowerThanNom?(yErrLo_ratio[ispp][_nMetafitSyst][b]):(yErrHi_ratio[ispp][_nMetafitSyst][b]);
	  cout<<"For this variation (not counted in systematic, we check (yratio-1) , sqrt(sigmayratio^2 + sigmayratio_nominal^2), ratio of the two = "<<y_ratio[ispp][isys][b]-1<<" "<<
	  sqrt(pow(yerrVar,2)+pow(yerrNom,2))<<" "<<
	  (y_ratio[ispp][isys][b]-1)/sqrt(pow(yerrVar,2)+pow(yerrNom,2))<<endl;
	}
      }
      y_systErrLo[ispp][b] /= ((RMSbyBlocks?_nRMSblocks:nVarUsed) - (RMSbyBlocks?0:1));
      y_systErrLo[ispp][b] = sqrt( y_systErrLo[ispp][b] );
      y_systErrHi[ispp][b] /= ((RMSbyBlocks?_nRMSblocks:nVarUsed) - (RMSbyBlocks?0:1));
      y_systErrHi[ispp][b] = sqrt( y_systErrHi[ispp][b] );
      cout<<"******   b, y relative systematic error = "<<b<<" "<<(y_systErrLo[ispp][b] + y_systErrHi[ispp][b])/2<<endl;
      //y_fitErr[b] = (*rsig_relerr)[b+1] * y_nom[ispp][b];
    }

    int idx = 0;
    for(int b=0;b<nbins;b++){
      for(int b2=b+1;b2<nbins;b2++){
	y_systCorr[ispp].push_back(0);
	float nVarUsed = 0;
	if(RMSbyBlocks){
	  for(int iblo=0;iblo<_nRMSblocks;iblo++){
	    float err1 = (fabs(y_systErrLo_maxDev[ispp][b][iblo])>y_systErrHi_maxDev[ispp][b][iblo])?( y_systErrLo_maxDev[ispp][b][iblo] ):( y_systErrHi_maxDev[ispp][b][iblo] );
	    float err2 = (fabs(y_systErrLo_maxDev[ispp][b2][iblo])>y_systErrHi_maxDev[ispp][b2][iblo])?( y_systErrLo_maxDev[ispp][b2][iblo] ):( y_systErrHi_maxDev[ispp][b2][iblo] );
	    y_systCorr[ispp][idx] += err1*err2;
	  }
	} else{
	  for(int isys=0;isys<_nMetafitSyst;isys++){
	    if(usedInRMS[isys] || (doAverageNominal && isys==0)){
	      y_systCorr[ispp][idx] += (y_ratio[ispp][isys][b] -1) * (y_ratio[ispp][isys][b2] -1);
	      nVarUsed += 1;
	    }
	  }
	}
	
	y_systCorr[ispp][idx] /= ((RMSbyBlocks?_nRMSblocks:nVarUsed) - (RMSbyBlocks?0:1));
	y_systCorr[ispp][idx] /= ((y_systErrLo[ispp][b] + y_systErrHi[ispp][b])/2);
	y_systCorr[ispp][idx] /= ((y_systErrLo[ispp][b2] + y_systErrHi[ispp][b2])/2);
	cout<<"correlation factor of metafit syst between bins "<<b<<" and "<<b2<<" = "<<y_systCorr[ispp][idx]<<endl;
	idx += 1;
      }
    }

    TFile * outf = new TFile("corrected_yields.root","UPDATE");
    if(ispp<2) outf->WriteObject(&y_nom[ispp],"FinalCorrectedYield"+(TString)(firstStep?"_1stStep":"_2ndStep")+(TString)((ispp==1)?"_pp":"_PbPb"));
    outf->WriteObject(&y_systErrLo[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitRelSystErrorLo_"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    outf->WriteObject(&y_systErrHi[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitRelSystErrorHi_"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    outf->WriteObject(&y_systErrLo_maxDev[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitRelSystErrorLoMaxDeviation"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    outf->WriteObject(&y_systErrHi_maxDev[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitRelSystErrorHiMaxDeviation"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    outf->WriteObject(&y_systCorr[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitSyst_LinearizedCorrelationMatrix"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    //  outf->WriteObject(&y_fitErr,"FinalCorrectedYield_fitError"+(TString)((ispp==1)?"_pp":"_PbPb"));
    outf->Close();
  }

  //******************************************************
  //Draw the ratio of variations to nominal, for each bin, in pp, PbPb, and RAA
  //******************************************************

  TCanvas *c2 = new TCanvas("c2","c2",3300,1100);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.15);
  c2->Divide(7,1,0,0);

  TLatex binname;
  binname.SetTextFont(42);
  binname.SetTextSize(0.08);
  binname.SetTextAlign(21);
  binname.SetNDC();

  for(int ispp=0;ispp<3;ispp++){
    for(int b=0;b<nbins;b++){
      c2->cd(2+nbins*ispp+b);
      gPad->SetLeftMargin(0.);
      gPad->SetRightMargin((2+nbins*ispp+b==7)?0.01:0.);
      gPad->SetTopMargin(0.04);
      gPad->SetBottomMargin(0.15);

      //multigraph for a given bin and ispp value, gathering all systs
      gSystRatio_all[ispp][b]->Draw("AP0");

      gSystRatio_all[ispp][b]->SetTitle("");//Ratio of metafit systematics to nominal");
      (gSystRatio_all[ispp][b]->GetHistogram())->GetYaxis()->SetRangeUser(0,_nMetafitSyst+1);
      (gSystRatio_all[ispp][b]->GetHistogram())->GetXaxis()->SetRangeUser((ispp==1)?0.85:0.5,(ispp==1)?1.15:1.5);
      (gSystRatio_all[ispp][b]->GetHistogram())->GetYaxis()->SetTitle("");
      if(ispp==2 && b==2) (gSystRatio_all[ispp][b]->GetHistogram())->GetXaxis()->SetTitle("variation / nominal");//(ispp==2)?"R_{PbPb}^{systematic}/R_{PbPb}^{nominal}":"yield^{systematic}/yield^{nominal}");
      (gSystRatio_all[ispp][b]->GetHistogram())->GetXaxis()->SetLabelSize(0.065);
      (gSystRatio_all[ispp][b]->GetHistogram())->GetXaxis()->SetTitleSize(0.09);
      (gSystRatio_all[ispp][b]->GetHistogram())->GetXaxis()->SetTitleOffset(0.35);
      (gSystRatio_all[ispp][b]->GetHistogram())->GetXaxis()->SetLabelOffset(-0.02);

      (gSystRatio_all[ispp][b]->GetHistogram())->GetYaxis()->SetTickLength(0);
      (gSystRatio_all[ispp][b]->GetHistogram())->GetXaxis()->SetNdivisions(507);

      gSystRatio_all[ispp][b]->Draw("AP0");
      //gSystRatio_all[ispp][b]->Draw("Psame");

      TLine *one = new TLine(1,0,1,_nMetafitSyst+1);
      one->SetLineWidth(2);
      one->SetLineStyle(2);
      one->Draw("same");	

      if(b==0 && ispp==0){
	CMStag.SetTextSize(0.07);
	CMStag.DrawLatex(0.05,0.97,"#font[61]{CMS}#font[52]{ Preliminary}");
      }

      bool lastPad = b==nbins-1 && ispp==2;
      binname.DrawLatex(lastPad?0.25:0.5,lastPad?0.04:0.05,"#font[61]{#splitline{"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"PbPb/pp"))+"}{p_{T} bin "+(TString)to_string(b+1)+"}}");
    }

  }

  //DRAW CMS Preliminary 
  TLatex systname;
  systname.SetTextFont(62);
  systname.SetTextAlign(32);

  c2->cd(1);
  gPad->SetLeftMargin(0);
  gPad->SetRightMargin(0.);
  float dy = (0.96-0.15)/(_nMetafitSyst+1);
  for(int isys=0;isys<=_nMetafitSyst;isys++){
    float y = 0.15 + (isys+1.5)*dy;
    if(isys==_nMetafitSyst) y = 0.15+0.5*dy;

    Color_t col;
    if(isys==_nMetafitSyst) col = cols[0];
    else if(usedInNominalAverage[isys]) col = cols[1];
    else if(usedInRMS[isys]) col = cols[2];
    else col = cols[3];

    if(isys==6 || isys==7) systname.SetTextSize(0.046);
    else if(isys==3) systname.SetTextSize(0.05);
    else systname.SetTextSize(0.054);
    float correctWeird = 0;
    if(isys==3 || isys==6 || isys==7) correctWeird = 0.035;
    systname.DrawLatex(0.98+correctWeird,y,"#color["+(TString)to_string((int)col)+"]{"+systPrettyName[isys]+"}");
  }

  //horizontal lines to guide the eye
  for ( int isys : { 1,3,5,9,11 } ){
    TLine *hor = new TLine(0,0.15 + isys*dy ,1, 0.15 + isys*dy);
    hor->SetLineWidth(2);
    hor->SetNDC(true);
    hor->SetLineStyle(3);
    for (int p=1;p<=7;p++){
      c2->cd(p);
      hor->Draw("same");	
    }
  }

  c2->SaveAs("figs/correctedYields_and_RAA_metafitSyst_ratioToNominal.pdf");
  c2->SaveAs("figs/correctedYields_and_RAA_metafitSyst_ratioToNominal.png");

}
