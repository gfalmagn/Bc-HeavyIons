#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TDirectory.h"
#include "../helpers/Cuts_BDT.h" //this first
#include "../helpers/Cuts.h"

float Maxi3(float f1, float f2, float f3){
  return max(f1,max(f2,f3));
}

bool needsRegul(TH1F* h, bool loose=true){
  //check if the histo passes the 'not enough stats' conditions
  bool needsReg = false, looseCondition = true;
  bool emptybin = false, seenNonEmpty = false;

  for(int b=1;b<=h->GetNbinsX();b++){ //if some bins have a lot of stats, no regularization
    if(h->GetBinContent(b)>0 && h->GetBinError(b)/h->GetBinContent(b)<0.3) return false;
  }
  
  if(loose){
    for(int b=1;b<=h->GetNbinsX();b++){
      if(h->GetBinContent(b)>0 && h->GetBinError(b)/h->GetBinContent(b)<0.5) {
	looseCondition = false; break;}
    }
  }

  //empty bin surrounded by non-empty bins?
  for(int b=2;b<h->GetNbinsX();b++){
    if(emptybin && 
       h->GetBinContent(b)>0) 
      needsReg = true; //order of the 3 conditions does matter
    if(seenNonEmpty && 
       (h->GetBinContent(b)<=0 || (loose && h->GetBinError(b) > 0.95*h->GetBinContent(b))))
      emptybin = true;
    if(h->GetBinContent(b)>0) 
      seenNonEmpty = true;
  }

  return needsReg || (loose && looseCondition);
}

void floatingAverage3bins(TH1F* h){
  TH1F* hcopy = (TH1F*)h->Clone("hcopy");

  const int n = h->GetNbinsX();
  //special cases for bin 1 and n
  h->SetBinContent(1, (2*hcopy->GetBinContent(1)+hcopy->GetBinContent(2))/3 );
  h->SetBinError(1, max(hcopy->GetBinError(1),hcopy->GetBinError(2)));//sqrt( pow(2*hcopy->GetBinError(1),2) + pow(hcopy->GetBinError(2),2) )/3 );// sqrt( ...)/3
  h->SetBinContent(n, (2*hcopy->GetBinContent(n)+hcopy->GetBinContent(n-1))/3 );
  h->SetBinError(n, max(hcopy->GetBinError(n),hcopy->GetBinError(n-1)));//sqrt( pow(2*hcopy->GetBinError(n),2) + pow(hcopy->GetBinError(n-1),2) )/3 );// sqrt( ...)/3
  //average over 3 bins
  for(int b=2;b<n;b++){
    h->SetBinContent(b, ( hcopy->GetBinContent(b-1)+hcopy->GetBinContent(b)+hcopy->GetBinContent(b+1) )/3 );
    h->SetBinError(b, Maxi3(hcopy->GetBinError(b-1),hcopy->GetBinError(b),hcopy->GetBinError(b+1)));//sqrt( pow(hcopy->GetBinError(b-1),2) + pow(hcopy->GetBinError(b),2) + pow(hcopy->GetBinError(b+1),2) )/3 );
  }

  hcopy->Delete();
}

int MakePositive(TH1F* h, bool regularize=false, int forceReg=-1){

  //for histos with too few stats, make a 3-bins floating average
  int reg = 0;
  if(regularize && forceReg!=0){
    if((needsRegul(h,true) && forceReg==-1) || forceReg>=1){
      floatingAverage3bins(h);
      reg = 1;
    }
    // if((reg==1 && needsRegul(h,false) && forceReg==-1) || forceReg>=2) { //only if 1st regularization was done
    //   floatingAverage3bins(h); //if there are still empty bins surrounded by non-empty bins, regularize again
    //   reg = 2;
    // }
    // if((reg==2 && needsRegul(h,false) && forceReg==-1) || forceReg>=3) { //only if 2nd regularization was done
    //   floatingAverage3bins(h); //if there are still empty bins surrounded by non-empty bins, regularize again
    //   reg = 3;
    // }
    if(forceReg==-1 && reg>=1) cout<<"Regularized "<<h->GetName()<<" "<<reg<<" time(s)"<<endl;
  }

  //forbid bins with negative content
  bool nonEmpty = false;
  for(int b=1;b<=h->GetNbinsX();b++){
    if(h->GetBinContent(b)>0) nonEmpty = true;
    if(h->GetBinContent(b)<0) { //move negative contents to other bins
      int eps = (b>1 && h->GetBinContent(b-1)>=fabs(h->GetBinContent(b)))?-1:1;
      h->SetBinContent(b+eps,h->GetBinContent(b+eps)+h->GetBinContent(b));
      h->SetBinError(b+eps,sqrt(pow(h->GetBinContent(b+eps),2)+pow(h->GetBinContent(b),2)));
      h->SetBinContent(b,0);
      h->SetBinError(b,0);
    }
  }
  for(int b=1;b<=h->GetNbinsX();b++){ //avoid errors reaching negative values
    if(h->GetBinContent(b) - h->GetBinError(b) <0.01*h->GetBinContent(b)) {
      h->SetBinError(b, 0.99*h->GetBinContent(b));
    }
  }
  if(!nonEmpty){
    h->SetBinContent(h->FindBin(4.8),1e-5);
    h->SetBinError(h->FindBin(4.8),1e-5);
  }
  h->SetMinimum(0);

  return reg;
}

void application(vector<float> BDTcut, bool ispp, bool BDTuncorrFromM, int kinBin, bool scaleSystBDTintegrated=false, bool regulLowStatShapes=false, bool addAccEff=true, int massBinning=0, int varyBDTbin=0){

  cout<<"***************\n   Make histograms for options _"<<(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(scaleSystBDTintegrated?"_scaleSystBDTintegrated":"")+(TString)(regulLowStatShapes?"_regulLowStatShapes":"")+(TString)((massBinning!=0)?("_MbinsVar"+(TString)to_string(massBinning)):"")+(TString)((varyBDTbin!=0)?("_BDTbins"+(TString)((varyBDTbin==1)?"Up":"Down")):"")<<" , kinematic bin #"<<kinBin<<endl;

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  auto acceffFile = TFile::Open(addAccEff?"../efficiency/AcceptanceEfficiencyMap.root":"","READ");
  TH2Poly* hp_acceff = (TH2Poly*)(addAccEff?(  acceffFile->Get("hp_acceff_"+(TString)(ispp?"pp":"PbPb"))  ):(new TH2Poly()));
  hp_acceff->SetDirectory(0);
  TH2Poly* hpcoarse_inBDT23 = (TH2Poly*)(addAccEff?(  acceffFile->Get("hpcoarse_inBDT23_"+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)((kinBin==0)?"_integratePtBins":""))  ):(new TH2Poly()));
  hpcoarse_inBDT23->SetDirectory(0);
  TH2Poly* hpcoarse_inBDT3 = (TH2Poly*)(addAccEff?(  acceffFile->Get("hpcoarse_inBDT3_"+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)((kinBin==0)?"_integratePtBins":""))  ):(new TH2Poly()));
  hpcoarse_inBDT3->SetDirectory(0);
  acceffFile->Close();

  auto fullFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root");

  bool applyBDTweights = false; //recent change: don't apply these weights in pp

  const int nCuts = BDTcut.size()-1;
  int ntrees = 8;

  //Initialization of variables
  float Bc_CorrM[ntrees];
  float Bc_ctauSignif[ntrees];
  float Bc_ctauSignif3D[ntrees], Bc_log_ctauSignif3D[ntrees];
  float Bc_cosAlpha[ntrees], Bc_Alpha[ntrees];
  float Bc_cosAlpha3D[ntrees], Bc_Alpha3D[ntrees];
  float Bc_VtxProb[ntrees], Bc_log1min_VtxProb[ntrees];
  float QQ_VtxProb[ntrees], QQ_log1min_VtxProb[ntrees];
  float QQ_dca[ntrees], QQ_log_dca[ntrees];
  float dR_sum[ntrees];
  float dR_jpsiOverMuW[ntrees];
  float MuonDxySignif_sum[ntrees], log_MuonDxySignif_sum[ntrees];
  float MuonDzSignif_sum[ntrees];
  float muonsGlbInLooseAcc[ntrees];
  float muonsInTightAcc[ntrees];
  bool muW_isGlb[ntrees];
  bool muW_inLooseAcc[ntrees];
  bool muW_inTightAcc[ntrees];
  bool mumi_isGlb[ntrees];
  bool mumi_inLooseAcc[ntrees];
  bool mumi_inTightAcc[ntrees];
  bool mupl_isGlb[ntrees];
  bool mupl_inLooseAcc[ntrees];
  bool mupl_inTightAcc[ntrees];
  float QQmuW_ptImbal[ntrees];
  float Bc_M[ntrees];
  float Bc_Pt[ntrees];
  float Bc_Y[ntrees];
  float QQ_M[ntrees];
  float QQ2_M[ntrees];
  float QQ3_M[ntrees];
  int bkgType[ntrees];
  bool muW_isJpsiBro[ntrees];
  int muW_trueId[ntrees];
  float weight[ntrees];
  float weightJpsiChoice[ntrees];
  int flipJpsi[ntrees];
  float BDTv[ntrees];
  float BDTprob[ntrees];
  UInt_t eventNb[ntrees];

  //what fitting strategy
  int nbin = 32; //for other histos than Bc_M
  int nbinM = 25;
  vector<float> CRbinwRatio;
  for(int k=1;k<=3;k++) CRbinwRatio.push_back( ((_mBcMax-_mBcMin)/_nbinMSR(ispp)[k]) * (_nbinMCR(ispp)[k]/(_mMax-_mBcMax)) );

  //*******************************************
  //process names, pretty names of histos (designating actual content)
  vector<TString> procName = ispp?( 
				   (vector<TString>){"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig", //0-4
				       "NonPromptJpsi","bToJpsi","promptJpsi","JpsiMC_uncorr","JpsiMC","JpsiMC_flipJSameSideUp","JpsiMC_flipJSameSideDown","JpsiMC_wPromptMCUp","JpsiMC_wPromptMCDown", //5-13
				       "flipJpsi","flipJpsi_flipJSameSideUp","flipJpsi_flipJSameSideDown","flipJpsi_wPromptMCUp","flipJpsi_wPromptMCDown", //14-18
				       "FakeJpsi_JpsiSBUp","FakeJpsi_JpsiSBDown"}//19-20
				    ):(//PbPb
				    (vector<TString>){"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig", //0-4
					"NonPromptJpsi","bToJpsi","promptJpsi","JpsiMC_uncorr","NPJpsi","NPJpsi_FlipJorMCUp","NPJpsi_FlipJorMCDown","NPJpsi_UncorrNPJUp","NPJpsi_UncorrNPJDown", //5-13
					"flipJpsi","flipJpsi_FlipJorMCUp","flipJpsi_FlipJorMCDown","flipJpsi_UncorrNPJUp","flipJpsi_UncorrNPJDown", //14-18
					"FakeJpsi_JpsiSBUp","FakeJpsi_JpsiSBDown"}//19-20
      );
  int JMCcontent[5] = { (ispp?1:0), (ispp?1:0), (ispp?1:0), (ispp?2:1), 0}; //{nominal,nominal+flipJSameSideUp,nominal+flipJSameSideDown,systUp,systDown} //points to index of JMCname
  TString JMCname[6] = {"MC b#to J/#psi","MC non-prompt J/#psi","full J/#psi MC","MC non-prompt J/#psi, enhanced bToJpsi","MC non-prompt J/#psi, reduced bToJpsi","MC J/#psi prompt + loosely correlated"}; //0 = only bToJpsi, 1 = NonPrompt Jpsi, 2 = prompt+nonprompt Jpsi, 3 = NonPrompt Jpsi with 3x bToJpsi, 4 = NonPrompt Jpsi with 0.33x bToJpsi, 5 = prompt + !bToJpsi
  vector<TString> prettyName = ispp?(
				     (vector<TString>){"WrongSign","J/#psi sidebands","High mass control","signal region","MC signal expectation",
					 JMCname[1],JMCname[0],"MC prompt J/#psi",JMCname[5], 
					 JMCname[JMCcontent[0]], JMCname[JMCcontent[1]]+" flipJSameSideUp", JMCname[JMCcontent[2]]+" flipJSameSideDown", JMCname[JMCcontent[3]], JMCname[JMCcontent[4]],
					 "flipped J/#psi","flipped J/#psi close-#phi side","flipped J/#psi opposite-#phi side","flipped J/#psi wPromptMCUp","flipped J/#psi wPromptMCDown",
					 "J/#psi upper sidebands","J/#psi lower sideband"}
				     ):(//PbPb
				     (vector<TString>){"WrongSign","J/#psi sidebands","High mass control","signal region","MC signal expectation",
					 JMCname[1],JMCname[0],"MC prompt J/#psi",JMCname[5], 
					 JMCname[JMCcontent[0]], JMCname[JMCcontent[1]]+" FlipJorMCUp", JMCname[JMCcontent[2]]+" FlipJorMCDown", JMCname[JMCcontent[3]], JMCname[JMCcontent[4]],
					 "flipped J/#psi",JMCname[5],"flipped J/#psi","flipped J/#psi UncorrNPJUp","flipped J/#psi UncorrNPJDown",
					 "J/#psi upper sidebands","J/#psi lower sideband"}
					);


  //trees and hists
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","flipJpsi"};
  int nProc = procName.size();
  vector<TH1F*> h_BDT, h_BDTprob;
  vector<vector<TH1F*> > h_bdt(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_BcM(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_MeanInvAccEff(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_AccEffCorr(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_MeanInvAccEff_BDT23(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_AccEffCorr_BDT23(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_MeanInvAccEff_BDT3(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_AccEffCorr_BDT3(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_AccEffWeights(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_CorrYieldsVsAccEffWeights(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_BcPt(nProc, vector<TH1F*>(nCuts,new TH1F( "Bc_Pt", "Bc transverse momentum", 50, 0, 30 )));
  vector<vector<TH1F*> > h_QQM(nProc, vector<TH1F*>(nCuts,new TH1F( "QQ_M", "Jpsi mass", nbinM, 2,4 )));

  //grab tree
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }

  //logarithmic binning
  const int nlogbins = 30;
  float xmin=2, xmax=2500;
  float lxmin=TMath::Log(xmin), lxmax=TMath::Log(xmax);
  float logbins[nlogbins+1];
  for(int l=0;l<=nlogbins;l++){
    logbins[l] = TMath::Exp(lxmin + l * (lxmax-lxmin)/((float)nlogbins));
  }

  //Prepare output for each dataset (signal and backgrounds)
  for(int i=0; i<nProc; i++){
    h_BDT.push_back( new TH1F( "BDT_"+procName[i], "BDT_"+prettyName[i], nbin, -0.9, 0.9 ) );
    for(int k=0;k<nCuts;k++){
      //mass binning can change with systematics
      int kForMbin = (massBinning==1)?(nCuts-1):((massBinning==2)?0:k);
      int nmbins =  _nbinM(ispp)[kForMbin]; 
      float mbins[nmbins+1];
      for(int mb=0;mb<=nmbins;mb++) 
	mbins[mb] = _Mbinning(ispp,kForMbin)[mb];

      //Initialise histos
      h_bdt[i][k] = new TH1F( "BDT_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "BDT "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"]", 50, -0.9,0.9 );
      h_BcM[i][k] = new TH1F( "BcM_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "B_{c} mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"]", nmbins, mbins );
      h_MeanInvAccEff[i][k] = new TH1F( "MeanInvAccEff_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), 
					"Mean 1/#alpha#times#varepsilon "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"]", nmbins, mbins );
      h_AccEffWeights[i][k] = new TH1F( "AccEffWeights_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), 
					"1/#alpha#times#varepsilon weights "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"]", nlogbins,logbins );
      h_MeanInvAccEff_BDT23[i][k] = new TH1F( "MeanInvAccEff_BDT23_"+procName[i]+"_BDT"+(TString)(to_string(k+1)),
					      "Mean 1/#alpha#times#varepsilon "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"], with efficiency of being in BDT bins 2-3", nmbins, mbins );
      h_MeanInvAccEff_BDT3[i][k] = new TH1F( "MeanInvAccEff_BDT3_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), 
					     "Mean 1/#alpha#times#varepsilon "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"], with efficiency of being in BDT bin 3", nmbins, mbins );
      h_CorrYieldsVsAccEffWeights[i][k] = new TH1F( "CorrYieldsVsAccEffWeights_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), 
						    "Corrected yields VS 1/#alpha#times#varepsilon weights "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"]", nlogbins,logbins );
      h_BcPt[i][k] = new TH1F( "BcPt_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "B_{c} p_{T} "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"]", 50, 0, 30 );
      h_QQM[i][k] = new TH1F( "QQM_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "J/#psi mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,kinBin,BDTuncorrFromM)[k+1])+"]", nbinM, 2.6,3.6 );
    }
  }

  //*******************************************
  //Fetch the BDT weights to be applied to flipJpsi
  TFile* flipJWfile = new TFile("../BDT/weightingCR/flipJpsiWeights_fromBDTinCR_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  vector<TH1F*> flipJBdtWeight;
  if(applyBDTweights){
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinCR_flipJpsi"+(TString)(ispp?"0":"3")+"_JpsiMC"+(TString)(to_string(JMCcontent[0]))+"_AddMCtoFlipJ"));
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinCR_flipJpsi1_JpsiMC"+(TString)(to_string(JMCcontent[1]))+"_AddMCtoFlipJ"));//flipJpsiSameSide
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinCR_flipJpsi2_JpsiMC"+(TString)(to_string(JMCcontent[2]))+"_AddMCtoFlipJ"));//flipJpsiOppSide
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinCR_flipJpsi"+(TString)(ispp?"0":"3")+"_JpsiMC"+(TString)(to_string(JMCcontent[3]))+"_AddMCtoFlipJ"));//wPromptMCUp
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinCR_flipJpsi"+(TString)(ispp?"0":"3")+"_JpsiMC"+(TString)(to_string(JMCcontent[4]))+"_AddMCtoFlipJ"));//wPromptMCDown

    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinSR_flipJpsi"+(TString)(ispp?"0":"3")+"_JpsiMC"+(TString)(to_string(JMCcontent[0]))+"_AddMCtoFlipJ"));
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinSR_flipJpsi1_JpsiMC"+(TString)(to_string(JMCcontent[1]))+"_AddMCtoFlipJ"));//flipJpsiSameSide
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinSR_flipJpsi2_JpsiMC"+(TString)(to_string(JMCcontent[2]))+"_AddMCtoFlipJ"));//flipJpsiOppSide
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinSR_flipJpsi"+(TString)(ispp?"0":"3")+"_JpsiMC"+(TString)(to_string(JMCcontent[3]))+"_AddMCtoFlipJ"));//wPromptMCUp
    flipJBdtWeight.push_back((TH1F*) flipJWfile->Get("flipJpsiWeights_fromBDTinSR_flipJpsi"+(TString)(ispp?"0":"3")+"_JpsiMC"+(TString)(to_string(JMCcontent[4]))+"_AddMCtoFlipJ"));//wPromptMCDown
  }

  //*******************************************
  //Fetch BDT correction = f(M) to be subtracted from BDT, to uncorrelate it from mass
  TFile* f_BDTuncorrel = new TFile("BDTuncorrFromM_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  TH1F* h_correctBDT = BDTuncorrFromM?( (TH1F*) f_BDTuncorrel->Get("avBDTvsM_bkg_KinBin"+(TString)to_string(kinBin)) ):NULL; //make the BDT of the expected background (postfit) uncorrelated with mass
  if(kinBin==0 && BDTuncorrFromM) {
    for(int b=2;b<=_NanaBins;b++)
      h_correctBDT->Add((TH1F*) f_BDTuncorrel->Get("avBDTvsM_bkg_KinBin"+(TString)to_string(b)));
    h_correctBDT->Scale(1/_NanaBins); //average of correction functions of all bins 
  }

  //*******************************************
  //For the signal region and the tmva output, fill the histos for BDT and Bc_M
  for(int iT=0; iT<(int)T.size(); iT++){
    if(BDTuncorrFromM==false && scaleSystBDTintegrated==false && regulLowStatShapes==false && massBinning==0 && varyBDTbin==0) //nominal
      std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    T[iT]->SetBranchAddress("BDT", &BDTv[iT]);
    T[iT]->SetBranchAddress("Bc_M", &Bc_M[iT]);
    T[iT]->SetBranchAddress("Bc_Pt", &Bc_Pt[iT]);
    T[iT]->SetBranchAddress("Bc_Y", &Bc_Y[iT]);
    T[iT]->SetBranchAddress("QQ_M", &QQ_M[iT]);
    T[iT]->SetBranchAddress("muW_isJpsiBro", &muW_isJpsiBro[iT]);
    if(iT==7) T[iT]->SetBranchAddress("flipJpsi", &flipJpsi[iT]);
    T[iT]->SetBranchAddress("weight", &weight[iT]);
    //    T[iT]->SetBranchAddress("eventNb", &eventNb[iT]);

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){//T[iT]->GetEntries() 
      T[iT]->GetEntry(j);

      //Keep events from the wanted analysis bin
      if(!inFidCuts(kinBin,Bc_Pt[iT],Bc_Y[iT])) continue;
      if(weight[iT]==0) continue;

      vector<int> iProc; //processes in which we want to include this event
      iProc.push_back(iT);
      
      //*******************************************
      //Which processes this event should go in ?
      if(iT==1){
	if(QQ_M[iT] > m_Jpsi) iProc.push_back(19);
	else iProc.push_back(20);
      }
      else if(iT==5){

	for(int JMCvar=0;JMCvar<5;JMCvar++){ //test JpsiMC and 4 variations, to see if they include nonprompt Jpsi
	  if(JMCcontent[JMCvar]!=0 && JMCcontent[JMCvar]!=5) iProc.push_back(9+JMCvar);}

	if(muW_isJpsiBro[iT]){ //if true b->Jpsi, add in corresponding process
	  iProc.push_back(6); 
	  for(int JMCvar=0;JMCvar<5;JMCvar++){ //test JpsiMC and 4 variations, to see if they are bToJpsiOnly (JMCcontent==0)
	    if(JMCcontent[JMCvar]==0) iProc.push_back(9+JMCvar);}	  
	}
	else {
	  iProc.push_back(8); //prompt + loosely correlated
	  if(!ispp) iProc.push_back(15); //for variation of flipJpsi shape in PbPb
	  for(int JMCvar=0;JMCvar<5;JMCvar++){ //test JpsiMC and 4 variations, to see if they exclude bToJpsi (JMCcontent==5)
	    if(JMCcontent[JMCvar]==5) iProc.push_back(9+JMCvar);} 
	}

      }
      else if(iT==6){

	iProc[0] = 7; //index for prompt Jpsi process
	iProc.push_back(8); //prompt + loosely correlated
	if(!ispp) iProc.push_back(15); //for variation of flipJpsi shape in PbPb
	for(int JMCvar=0;JMCvar<5;JMCvar++){ //test JpsiMC and 4 variations, to see if they include prompt Jpsi
	  if(JMCcontent[JMCvar]==2 || JMCcontent[JMCvar]==5) iProc.push_back(9+JMCvar);}
	//if(!ispp) {//OLD
	//  iProc.push_back(14); //in PbPb, PromptMC is the equivalent of flipJpsi in pp 
	//  iProc.push_back(17); iProc.push_back(18); //change BDT weights when the nominal Jpsi MC changes
	//}

      }
      else if(iT==7){
	iProc.erase(iProc.begin());
	//if(ispp) {
	iProc.push_back(14); if(!ispp) iProc.push_back(16);
	iProc.push_back(17); iProc.push_back(18);
	//} //change BDT weights when the nominal Jpsi MC changes
	if(ispp){
	  if(flipJpsi[iT]==1 || flipJpsi[iT]==6) iProc.push_back(16); //flipJpsi opposite-side
	  else iProc.push_back(15); //flipJpsi same-side
	}
      }

      //*******************************************
      //loop on processes that include this given event
      for(int l=0; l<iProc.size();l++){
	int iproc = iProc[l];
	//weight
	float w = weight[iT];
	float bdtcorr = BDTuncorrFromM?( h_correctBDT->GetBinContent(h_correctBDT->FindBin(Bc_M[iT])) ):0;

	// if(iT==1){
	//   if(iproc==19) w *= (ispp?2.415:2.028); //hard-coded correction of norm of lower and upper Jpsi sidebands //2020/04/27 //needed only if the histos are not scaled by the nominal process (done later on)
	//   if(iproc==20) w *= (ispp?1.706:1.973);
	// }

	if(iT>=5){ //weights of Jpsi MC or flipJpsi events
	  if(ispp && iproc==15) w *= 7/5.;// *(ispp?1.:0.4);
	  if(ispp && iproc==16) w *= 7/2.;// *(ispp?1.:0.4);	
	  //BDT weights for flipJpsi
	  for(int flipMeth=0; flipMeth<5;flipMeth++){
	    if(iproc==(14+flipMeth) && applyBDTweights){
	      if(BDTv[iT]<(_withTM?-0.2:-0.35) && Bc_M[iT]<5.6 && _withTM) w *= flipJBdtWeight[5+flipMeth]->GetBinContent(flipJBdtWeight[5+flipMeth]->FindBin(BDTv[iT])); //weights from SR
	      else                                                         w *= flipJBdtWeight[0+flipMeth]->GetBinContent(flipJBdtWeight[0+flipMeth]->FindBin(BDTv[iT])); //from CR
	    }
	  }
	  //BDT weights for JpsiMC (uncorrelated part)
	  if(!muW_isJpsiBro[iT]){
	    for(int flipMeth=0; flipMeth<5;flipMeth++){
	      if(iproc==(9+flipMeth) && applyBDTweights){
		if(BDTv[iT]<(_withTM?-0.2:-0.35) && Bc_M[iT]<5.6 && _withTM) w *= flipJBdtWeight[5+flipMeth]->GetBinContent(flipJBdtWeight[5+flipMeth]->FindBin(BDTv[iT]));
		else                                                         w *= flipJBdtWeight[0+flipMeth]->GetBinContent(flipJBdtWeight[0+flipMeth]->FindBin(BDTv[iT]));
	      }
	    }
	  }
	  // //variations of the non-prompt Jpsi MC in PbPb //OLD
	  // if(!ispp && muW_isJpsiBro[iT]){
	  //   if(iproc==12) w *= 3;
	  //   if(iproc==13) w *= 0.33;
	  // }	
	}

	//	if(iproc==9 || iproc==12 || iproc==13)
	//cout<<"iproc, iT, w, w change, muW_isJpsiBro = "<<iproc<<" "<<iT<<" "<<w<<" "<<w/weight[iT]<<" "<<muW_isJpsiBro[iT]<<endl;
	//Fill histos
	h_BDT[iproc]->Fill( BDTv[iT] , w);

	for(int k=0;k<nCuts;k++){
	  if(BDTv[iT]-bdtcorr > BDTcut[k] && BDTv[iT]-bdtcorr < BDTcut[k+1]){
	    int kForMbin = (massBinning==1)?(nCuts-1):((massBinning==2)?0:k);
	    if(Bc_M[iT]>_mBcMax) w *= CRbinwRatio[kForMbin]; //compensate for larger bins in the high mass CR
	    h_bdt[iproc][k]->Fill(BDTv[iT]-bdtcorr, w ); 
	    h_BcM[iproc][k]->Fill(Bc_M[iT], w ); 
	    h_BcPt[iproc][k]->Fill(Bc_Pt[iT], w ); 
	    h_QQM[iproc][k]->Fill(QQ_M[iT], w );

	    if(addAccEff){
	      float acceff = hp_acceff->GetBinContent( hp_acceff->FindBin( fabs(Bc_Y[iT]) , Bc_Pt[iT] ));
	      float acceffW = (acceff==0)?300:(1/acceff);
	      h_MeanInvAccEff[iproc][k]->Fill(Bc_M[iT], w * acceffW);
	      h_AccEffWeights[iproc][k]->Fill(acceffW, w);
	      h_CorrYieldsVsAccEffWeights[iproc][k]->Fill(acceffW, w * acceffW);
	      float effBDT23 = hpcoarse_inBDT23->GetBinContent( hpcoarse_inBDT23->FindBin( fabs(Bc_Y[iT]) , Bc_Pt[iT] ));
	      float effBDT3 = hpcoarse_inBDT3->GetBinContent( hpcoarse_inBDT3->FindBin( fabs(Bc_Y[iT]) , Bc_Pt[iT] ));
	      h_MeanInvAccEff_BDT23[iproc][k]->Fill(Bc_M[iT], w * acceffW / ((effBDT23==0)?1.5:effBDT23));
	      h_MeanInvAccEff_BDT3[iproc][k]->Fill(Bc_M[iT], w * acceffW / ((effBDT3==0)?3.:effBDT3));
	    }

	  }
	}
      }//END loop on processes
      
    }
    //END event loop
    
  }
  //END loop on trees

  //Getting integrals of processes, integrated over BDT bins
  vector<float> integ(nProc,0);
  if(scaleSystBDTintegrated){
    for(int i=0; i<nProc; i++){
      for(int k=0;k<nCuts;k++){
	integ[i] += h_BcM[i][k]->Integral();}
    }
  }

  //********************************************************
  //Recording mass histos intended for combine
  //********************************************************
  TFile *f = new TFile("InputForCombine_"
		       +(TString)(ispp?"pp":"PbPb")
		       +(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")
		       +(TString)(scaleSystBDTintegrated?"_scaleSystBDTintegrated":"")
		       +(TString)(regulLowStatShapes?"_regulLowStatShapes":"")
		       +(TString)((massBinning!=0)?("_MbinsVar"+(TString)to_string(massBinning)):"")
		       +(TString)((varyBDTbin!=0)?("_BDTbins"+(TString)((varyBDTbin==1)?"Up":"Down")):"")
		       +".root", (kinBin<2)?"recreate":"update");  
  TDirectory *dir1[nCuts];
  TDirectory *dir2[nProc];
  for(int k=0;k<nCuts;k++){
    
    f->cd();
    dir1[k] = f->mkdir("BDT"+(TString)(to_string(k+1))+(TString)((kinBin>0)?("Kin"+to_string(kinBin)):""));
    for(int i=0; i<nProc; i++){
      dir1[k]->cd();
      dir2[k] = dir1[k]->mkdir(procName[i]);
      dir2[k]->cd();
      
      MakePositive(h_bdt[i][k]);
      int didReg = MakePositive(h_BcM[i][k], regulLowStatShapes && i!=3 && i!=4);
      MakePositive(h_BcPt[i][k]);
      MakePositive(h_QQM[i][k]);

      int scaleWithProc = -1;
      if(i>9 && i<=13) scaleWithProc = 9;
      if(i>14 && i<=18) scaleWithProc = 14;
      if(i>=19 && i<=20) scaleWithProc = 1;
      
      //if(scaleWithProc!=-1) cout<<"scaling histo of proc "<<i<<" with "<<integ[scaleWithProc]<<"/"<<integ[i]<<" = "<<integ[scaleWithProc]/integ[i]<<endl;
      float scale = (scaleWithProc==-1)?1:(h_BcM[scaleWithProc][k]->Integral() / h_BcM[i][k]->Integral());
      if(scaleSystBDTintegrated)
	scale = (scaleWithProc==-1)?1:(integ[scaleWithProc]/integ[i]);
      h_BcM[i][k]->Scale(scale);

      h_bdt[i][k]->Write("BDTv");
      h_BcM[i][k]->Write("BcM");
      h_BcPt[i][k]->Write("BcPt");
      h_QQM[i][k]->Write("JpsiM");

      if(addAccEff) {
	MakePositive(h_MeanInvAccEff[i][k], (i!=3 && i!=4), didReg); //need same regularization as simple Bc_M
	MakePositive(h_MeanInvAccEff_BDT23[i][k], (i!=3 && i!=4), didReg);
	MakePositive(h_MeanInvAccEff_BDT3[i][k], (i!=3 && i!=4), didReg);
	h_MeanInvAccEff[i][k]->Scale(scale);
	h_MeanInvAccEff_BDT23[i][k]->Scale(scale);
	h_MeanInvAccEff_BDT3[i][k]->Scale(scale);

	h_MeanInvAccEff[i][k]->Write("BcM_AccEffWeighted");
	h_MeanInvAccEff_BDT23[i][k]->Write("BcM_AccEffWeighted_BDTeff23");
	h_MeanInvAccEff_BDT3[i][k]->Write("BcM_AccEffWeighted_BDTeff3");
	h_AccEffWeights[i][k]->Write("AccEffWeights");
	h_CorrYieldsVsAccEffWeights[i][k]->Write("CorrYieldsVsAccEffWeights");
	
	h_AccEffCorr[i][k] = (TH1F*)h_MeanInvAccEff[i][k]->Clone("AccEffCorr_"+procName[i]+"_BDT"+(TString)(to_string(k+1)));
	h_AccEffCorr[i][k]->Divide(h_BcM[i][k]);
	for(int B=0;B<h_BcM[i][k]->GetNbinsX();B++){
	  if(h_BcM[i][k]->GetBinContent(B) <=0) {
	    h_AccEffCorr[i][k]->SetBinError(B, 0.); 
	    h_AccEffCorr[i][k]->SetBinContent(B, 1.); }
	}
	h_AccEffCorr[i][k]->Write("BcM_AccEffCorr");
	  
	h_AccEffCorr_BDT23[i][k] = (TH1F*)h_MeanInvAccEff_BDT23[i][k]->Clone("AccEffCorr_BDTeff23_"+procName[i]+"_BDT"+(TString)(to_string(k+1)));
	h_AccEffCorr_BDT23[i][k]->Divide(h_BcM[i][k]);
	for(int B=0;B<h_BcM[i][k]->GetNbinsX();B++){
	  if(h_BcM[i][k]->GetBinContent(B) <=0) {
	    h_AccEffCorr_BDT23[i][k]->SetBinError(B, 0.); 
	    h_AccEffCorr_BDT23[i][k]->SetBinContent(B, 1.); }
	}
	h_AccEffCorr_BDT23[i][k]->Write("BcM_AccEffCorr_BDTeff23");
	  
	h_AccEffCorr_BDT3[i][k] = (TH1F*)h_MeanInvAccEff_BDT3[i][k]->Clone("AccEffCorr_BDTeff3_"+procName[i]+"_BDT"+(TString)(to_string(k+1)));
	h_AccEffCorr_BDT3[i][k]->Divide(h_BcM[i][k]);
	for(int B=0;B<h_BcM[i][k]->GetNbinsX();B++){
	  if(h_BcM[i][k]->GetBinContent(B) <=0) {
	    h_AccEffCorr_BDT3[i][k]->SetBinError(B, 0.); 
	    h_AccEffCorr_BDT3[i][k]->SetBinContent(B, 1.); }
	}
	h_AccEffCorr_BDT3[i][k]->Write("BcM_AccEffCorr_BDTeff3");
	  
      }

    }
  }
  f->Close();

}

void HistsForCombine(bool ispp = true){

  bool BDTuncorrFromM=false;
  bool scaleSystBDTintegrated=false;
  bool regulLowStatShapes=false;
  bool addAccEff = true;
  int massBinning = 0;
  int BDTbinning = 0;
  //  application( _BDTcuts(ispp,0,BDTuncorrFromM) , ispp, BDTuncorrFromM, 0); //integrated bin
  for(int b=1;b<=_NanaBins;b++){
    BDTuncorrFromM=false;
    application( _BDTcuts(ispp,b,BDTuncorrFromM, 0) , ispp, BDTuncorrFromM, b, scaleSystBDTintegrated, regulLowStatShapes, addAccEff, massBinning, BDTbinning); //run bin1 before bin2
    application( _BDTcuts(ispp,b,BDTuncorrFromM,-1) , ispp, BDTuncorrFromM, b, scaleSystBDTintegrated, regulLowStatShapes, addAccEff, massBinning, -1); //BDT binning
    application( _BDTcuts(ispp,b,BDTuncorrFromM, 1) , ispp, BDTuncorrFromM, b, scaleSystBDTintegrated, regulLowStatShapes, addAccEff, massBinning, +1);
    application( _BDTcuts(ispp,b,BDTuncorrFromM, 0) , ispp, BDTuncorrFromM, b, scaleSystBDTintegrated, regulLowStatShapes, addAccEff, 1,           BDTbinning); //mass binning
    application( _BDTcuts(ispp,b,BDTuncorrFromM, 0) , ispp, BDTuncorrFromM, b, scaleSystBDTintegrated, regulLowStatShapes, addAccEff, 2,           BDTbinning);
    application( _BDTcuts(ispp,b,BDTuncorrFromM, 0) , ispp, BDTuncorrFromM, b, true,                   true,               addAccEff, massBinning, BDTbinning); //scale systematic shape variations with BDT-integrated yields
    application( _BDTcuts(ispp,b,BDTuncorrFromM, 0) , ispp, BDTuncorrFromM, b, scaleSystBDTintegrated, true,               addAccEff, massBinning, BDTbinning); //regularise low-stats shapes 
    BDTuncorrFromM=true;
    application( _BDTcuts(ispp,b,BDTuncorrFromM, 0) , ispp, BDTuncorrFromM, b, scaleSystBDTintegrated, regulLowStatShapes, addAccEff, massBinning, BDTbinning); //uncorrelate BDT variable from M
  }
}
