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
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

void MakePositive(TH1F* h){
  for(int b=0;b<=h->GetNbinsX();b++){
    if(h->GetBinContent(b)<0) h->SetBinContent(b,0);
  }
  h->SetMinimum(0);
  return;
}

void application(vector<float> BDTcut, bool ispp, bool BDTuncorrFromM, int kinBin){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  auto fullFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root");

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
  int nbin = 38; //for other histos than Bc_M
  int nbinM = 30;
  float CRbinwRatio = ((m_Bc-3.3)/_nbinMSR(ispp)) * (_nbinMCR(ispp)/1.);

  //*******************************************
  //process names, pretty names of histos (designating actual content)
  vector<TString> procName = ispp?( 
				   (vector<TString>){"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig", //0-4
				       "NonPromptJpsi","bToJpsi","promptJpsi","JpsiMC_uncorr","JpsiMC","JpsiMC_flipJSameSideUp","JpsiMC_flipJSameSideDown","JpsiMC_wPromptMCUp","JpsiMC_wPromptMCDown", //5-13
				       "flipJpsi","flipJpsi_flipJSameSideUp","flipJpsi_flipJSameSideDown","flipJpsi_wPromptMCUp","flipJpsi_wPromptMCDown", //14-18
				       "FakeJpsi_JpsiSBUp","FakeJpsi_JpsiSBDown"}//19-20
				    ):(//PbPb
				    (vector<TString>){"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig", //0-4
					"NonPromptJpsi","bToJpsi","promptJpsi","JpsiMC_uncorr","NPJpsi","NPJpsi_PromptOrFlipJUp","NPJpsi_PromptOrFlipJDown","NPJpsi_bJpsiFracUp","NPJpsi_bJpsiFracDown", //5-13
					"PromptJpsi","PromptJpsi_PromptOrFlipJUp","PromptJpsi_PromptOrFlipJDown","PromptJpsi_bJpsiFracUp","PromptJpsi_bJpsiFracDown", //14-18
					"FakeJpsi_JpsiSBUp","FakeJpsi_JpsiSBDown"}//19-20
      );
  int JMCcontent[5] = { 1, 1, 1, (ispp?2:3), (ispp?0:4)}; //{nominal,nominal+flipJSameSideUp,nominal+flipJSameSideDown,systUp,systDown} //points to index of JMCname
  TString JMCname[6] = {"MC b#to J/#psi","MC non-prompt J/#psi","full J/#psi MC","MC non-prompt J/#psi, enhanced bToJpsi","MC non-prompt J/#psi, reduced bToJpsi","MC J/#psi prompt + loosely correlated"}; //0 = only bToJpsi, 1 = NonPrompt Jpsi, 2 = prompt+nonprompt Jpsi, 3 = NonPrompt Jpsi with 3x bToJpsi, 4 = NonPrompt Jpsi with 0.33x bToJpsi, 5 = prompt + !bToJpsi
  vector<TString> prettyName = ispp?(
				     (vector<TString>){"WrongSign","J/#psi sidebands","High mass control","signal region","MC signal expectation",
					 JMCname[1],JMCname[0],"MC prompt J/#psi",JMCname[5], 
					 JMCname[JMCcontent[0]], JMCname[JMCcontent[1]]+" flipJSameSideUp", JMCname[JMCcontent[2]]+" flipJSameSideDown", JMCname[JMCcontent[3]], JMCname[JMCcontent[4]],
					 "flipped J/#psi","flipped J/#psi same-#eta side","flipped J/#psi opposite-#eta side","flipped J/#psi wPromptMCUp","flipped J/#psi wPromptMCDown",
					 "J/#psi upper sidebands","J/#psi lower sideband"}
				     ):(//PbPb
				     (vector<TString>){"WrongSign","J/#psi sidebands","High mass control","signal region","MC signal expectation",
					 JMCname[1],JMCname[0],"MC prompt J/#psi",JMCname[5], 
					 JMCname[JMCcontent[0]], JMCname[JMCcontent[1]]+" PromptOrFlipJUp", JMCname[JMCcontent[2]]+" PromptOrFlipJDown", JMCname[JMCcontent[3]], JMCname[JMCcontent[4]],
					 "Prompt J/#psi","Prompt J/#psi same-#eta side","Prompt J/#psi opposite-#eta side","Prompt J/#psi bJpsiFracUp","Prompt J/#psi bJpsiFracDown",
					 "J/#psi upper sidebands","J/#psi lower sideband"}
					);


  //trees and hists
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","flipJpsi"};
  int nProc = procName.size();
  vector<TH1F*> h_BDT, h_BDTprob;
  vector<vector<TH1F*> > h_bdt(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_BcM(nProc, vector<TH1F*>(nCuts)); 
  vector<vector<TH1F*> > h_BcPt(nProc, vector<TH1F*>(nCuts,new TH1F( "Bc_Pt", "Bc transverse momentum", 50, 0, 30 )));
  vector<vector<TH1F*> > h_QQM(nProc, vector<TH1F*>(nCuts,new TH1F( "QQ_M", "Jpsi mass", nbinM, 2,4 )));

  //grab tree
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }

  //Prepare output for each dataset (signal and backgrounds)
  for(int i=0; i<nProc; i++){
    h_BDT.push_back( new TH1F( "BDT_"+procName[i], "BDT_"+prettyName[i], nbin, -0.5, 0.5 ) );

    for(int k=0;k<nCuts;k++){
      h_bdt[i][k] = new TH1F( "BDT_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "BDT "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k+1])+"]", 50, -0.5,0.5 );
      h_BcM[i][k] = new TH1F( "BcM_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "B_{c} mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k+1])+"]", _nbinM(ispp), _Mbinning(ispp) );
      h_BcPt[i][k] = new TH1F( "BcPt_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "B_{c} p_{T} "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k+1])+"]", 50, 0, 30 );
      h_QQM[i][k] = new TH1F( "QQM_"+procName[i]+"_BDT"+(TString)(to_string(k+1)), "J/#psi mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k])+","+(TString)(_BDTcut_s(ispp,BDTuncorrFromM)[k+1])+"]", nbinM, 2.6,3.6 );
    }
  }

  //*******************************************
  //Fetch the BDT weights to be applied to flipJpsi
  TFile* flipJWfile = new TFile("../BDT/weightingCR/flipJpsiWeights_fromBDTinCR_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  vector<TH1F*> flipJBdtWeight;
  if(_withTM){
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
  TH1F* h_correctBDT = BDTuncorrFromM?( (TH1F*) f_BDTuncorrel->Get("avBDTvsM_bkg") ):NULL; //make the BDT of the expected background (postfit) uncorrelated with mass

  //*******************************************
  //For the signal region and the tmva output, fill the histos for BDT and Bc_M
  for(int iT=0; iT<(int)T.size(); iT++){
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
      if(kinBin>-1 && !( Bc_Pt[iT]>_BcPtmin[kinBin] && Bc_Pt[iT]<_BcPtmax[kinBin] && fabs(Bc_Y[iT])>_BcYmin[kinBin] && fabs(Bc_Y[iT])<_BcYmax[kinBin])) continue;
      if(kinBin==-1){//integrated bin
	bool inFidCuts = Bc_Pt[iT]>_BcPtmin[0] && Bc_Pt[iT]<_BcPtmax[0] && fabs(Bc_Y[iT])>_BcYmin[0] && fabs(Bc_Y[iT])<_BcYmax[0];
	inFidCuts = inFidCuts || (Bc_Pt[iT]>_BcPtmin[1] && Bc_Pt[iT]<_BcPtmax[1] && fabs(Bc_Y[iT])>_BcYmin[1] && fabs(Bc_Y[iT])<_BcYmax[1]);
	if(!inFidCuts) continue;
      }

      vector<int> iProc;
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
	  for(int JMCvar=0;JMCvar<5;JMCvar++){ //test JpsiMC and 4 variations, to see if they exclude bToJpsi (JMCcontent==5)
	    if(JMCcontent[JMCvar]==5) iProc.push_back(9+JMCvar);} 
	}

      }
      else if(iT==6){

	iProc[0] = 7; //index for prompt Jpsi process
	iProc.push_back(8); //prompt + loosely correlated
	for(int JMCvar=0;JMCvar<5;JMCvar++){ //test JpsiMC and 4 variations, to see if they include prompt Jpsi
	  if(JMCcontent[JMCvar]==2 || JMCcontent[JMCvar]==5) iProc.push_back(9+JMCvar);}
	if(!ispp) {
	  iProc.push_back(14); //in PbPb, PromptMC is the equivalent of flipJpsi in pp
	  iProc.push_back(17); iProc.push_back(18); //change BDT weights when the nominal Jpsi MC changes
	}

      }
      else if(iT==7){
	iProc.erase(iProc.begin());
	if(ispp) {iProc.push_back(14);
	  iProc.push_back(17); iProc.push_back(18);} //change BDT weights when the nominal Jpsi MC changes
	if(flipJpsi[iT]>4) iProc.push_back(15); //flipJpsi same-side
	else iProc.push_back(16); //flipJpsi opposite-side
      }

      //*******************************************
      //loop on processes that include this given event
      for(int l=0; l<iProc.size();l++){
	int iproc = iProc[l];
	//weight
	float w = weight[iT];
	if(Bc_M[iT]>m_Bc) w *= CRbinwRatio; //compensate for larger bins in the high mass CR
	float bdtcorr = BDTuncorrFromM?( h_correctBDT->GetBinContent(h_correctBDT->FindBin(Bc_M[iT])) ):0;

	if(iT==1){
	  if(iproc==19) w *= (ispp?2.415:2.028); //hard-coded correction of norm of lower and upper Jpsi sidebands //2020/04/27
	  if(iproc==20) w *= (ispp?1.706:1.973);
	}

	if(iT>=5){ //weights of Jpsi MC or flipJpsi events
	  if(iproc==15) w *= 7/3. *(ispp?1.:0.4);	
	  if(iproc==16) w *= 7/4. *(ispp?1.:0.4);	
	  //BDT weights for flipJpsi
	  for(int flipMeth=0; flipMeth<5;flipMeth++){
	    if(iproc==(14+flipMeth) && _withTM){
	      if(BDTv[iT]<(_withTM?-0.2:-0.3) && Bc_M[iT]<5.5) w *= flipJBdtWeight[5+flipMeth]->GetBinContent(flipJBdtWeight[5+flipMeth]->FindBin(BDTv[iT])); //weights from SR
	      else                              w *= flipJBdtWeight[0+flipMeth]->GetBinContent(flipJBdtWeight[0+flipMeth]->FindBin(BDTv[iT])); //from CR
	    }
	  }
	  //BDT weights for JpsiMC (uncorrelated part)
	  if(!muW_isJpsiBro[iT]){
	    for(int flipMeth=0; flipMeth<5;flipMeth++){
	      if(iproc==(9+flipMeth) && _withTM){
		if(BDTv[iT]<(_withTM?-0.2:-0.3) && Bc_M[iT]<5.5) w *= flipJBdtWeight[5+flipMeth]->GetBinContent(flipJBdtWeight[5+flipMeth]->FindBin(BDTv[iT]));
		else                              w *= flipJBdtWeight[0+flipMeth]->GetBinContent(flipJBdtWeight[0+flipMeth]->FindBin(BDTv[iT]));
	      }
	    }
	  }
	  //variations of the non-prompt Jpsi MC in PbPb
	  if(!ispp && !muW_isJpsiBro[iT]){
	    if(iproc==12) w *= 3;
	    if(iproc==13) w *= 0.33;
	  }	
	}

	//Fill histos
	h_BDT[iproc]->Fill( BDTv[iT] , w);

	for(int k=0;k<nCuts;k++){
	  if(BDTv[iT]-bdtcorr > BDTcut[k] && BDTv[iT]-bdtcorr < BDTcut[k+1]){
	    h_bdt[iproc][k]->Fill(BDTv[iT]-bdtcorr, w ); 
	    h_BcM[iproc][k]->Fill(Bc_M[iT], w ); 
	    h_BcPt[iproc][k]->Fill(Bc_Pt[iT], w ); 
	    h_QQM[iproc][k]->Fill(QQ_M[iT], w );
	  }
	}
      }
      
    }
    //END event loop
    
  }
  //END loop on trees

  //********************************************************
  //Recording mass histos intended for combine
  //********************************************************
  TFile *f = new TFile("InputForCombine_"+(TString)(BDTuncorrFromM?"BDTuncorrFromM_":"")+(TString)(ispp?"pp":"PbPb")+".root", (kinBin<1)?"recreate":"update");  
  TDirectory *dir1[nCuts];
  TDirectory *dir2[nProc];
  for(int k=0;k<nCuts;k++){
    
    f->cd();
    dir1[k] = f->mkdir("BDT"+(TString)(to_string(k+1))+(TString)((kinBin>-1)?("Kin"+to_string(kinBin)):""));
    for(int i=0; i<nProc; i++){
      dir1[k]->cd();
      dir2[k] = dir1[k]->mkdir(procName[i]);
      dir2[k]->cd();
      
      MakePositive(h_bdt[i][k]);
      MakePositive(h_BcM[i][k]);
      MakePositive(h_BcPt[i][k]);
      MakePositive(h_QQM[i][k]);

      h_bdt[i][k]->Write("BDTv");
      h_BcM[i][k]->Write("BcM");
      h_BcPt[i][k]->Write("BcPt");
      h_QQM[i][k]->Write("JpsiM");
    }
  }
  f->Close();

}

void HistsForCombine(bool ispp = true, bool BDTuncorrFromM=false){

  //  application(BDTuncorrFromM?(_corrBDTcuts(ispp)):(_BDTcuts(ispp)), ispp, BDTuncorrFromM, -1); //integrated bin
  application(BDTuncorrFromM?(_corrBDTcuts(ispp)):(_BDTcuts(ispp)), ispp, BDTuncorrFromM, 0); //run bin0 before bin1
  application(BDTuncorrFromM?(_corrBDTcuts(ispp)):(_BDTcuts(ispp)), ispp, BDTuncorrFromM, 1);

}
