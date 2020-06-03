#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
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
#include "TCanvas.h"
#include "TPad.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

void UncorrelateBDTfromM(bool ispp=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  int ntrees = 9;
  bool ignoreBin2 = false;
  bool ignoreBin1 = ignoreBin2 || false;
  bool useFlipJpsi = ispp;
  bool flipJpsiSameSide = false; // whether to keep only events with flipJpsi angle on same |eta| side
  bool bToJpsiOnly = false;//ispp;

  vector<float> BDTcut = _BDTcuts(ispp); //vector of BDT cut values put into array
  //  if(ignore1stBin) BDTcut.erase(0);
  //int nchan = BDTcut.size() -1;
  float BDTcut_l[_nChan(ispp)+2]; 
  BDTcut_l[0] = -1;
  for(int k=0;k<=_nChan(ispp);k++){
    BDTcut_l[k+1] = BDTcut[k];
  }
  
  vector<vector<TH1F*> > h_BcM(ntrees, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_postfit(ntrees, vector<TH1F*>(_nChan(ispp)+1 ));
  vector<TH1F* > h_BcM_bkg(_nChan(ispp)+1);
  
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi"};
  TString prettyName[] = {"WrongSign","J/#psi sidebands","High mass control","signal region","MC signal",
  			  "MC NonPromptJpsi","MC PromptJpsi","dimuon+track (misID)","flipped J/#psi"};
  
  TString procNameDef[] = {"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig",(ispp?"JpsiMC":"NPJpsi"),"PromptJpsi","dimuTrk","flipJpsi"};
  vector<vector<TString>> procName{{"Wrongsign"},{"FakeJpsi"},{"TrueJpsi"},{"data_obs"},{"BcSig"},
				  {"NonPromptJpsi","bToJpsi"}, {"PromptJpsi", "JpsiMC_uncorr"},
				  {"dimuTrk"}, {"flipJpsi","flipJpsiSameSide" }};
  vector<int> systIdx{0,0,0,0,0,
      (int)bToJpsiOnly,//nonprompt
      0,//prompt
      0,//dimuTrk
      (int)flipJpsiSameSide,//flipJpsi
      };
  
  bool usedForFit[] = {false,true,false,true,true,true,//nonprompt
		     !ispp,//prompt
		      false,//dimutrk
		      useFlipJpsi//flipJpsi
  };
  int nPiled = 0; for(int j=0;j<ntrees;j++) nPiled+= usedForFit[j];
  nPiled -= 1;
  int idxPiled[] = {1,5,8,4};//{1,5,8,4};
  if(!ispp) idxPiled[2] = 6; //PromptJpsi instead of flipJpsi

  //********************************************************  
  //Extract BcM HISTOGRAMS used as combine input
  //********************************************************
  vector<vector<TH1F*> > h_BcM_prefit(ntrees, vector<TH1F*>(_nChan(ispp)+1));
  auto histFile = TFile::Open("InputForCombine_"+(TString)(ispp?"pp":"PbPb")+".root");
  for(int i=0; i<ntrees; i++){
    if(!usedForFit[i]) continue;
    for(int k=0;k<=_nChan(ispp);k++){
      int kk = (k==0)?1:k;
      h_BcM_prefit[i][k] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/"+procName[i][systIdx[i]]+"/BcM");
      if(k>1) {
	h_BcM_prefit[i][0]->Add(h_BcM_prefit[i][k]);}
    }
  }

  //********************************************************
  //Extract (POSTFIT NORMALISATIONS & SHAPES) from combine output
  //********************************************************
  TString normFileName = "./CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(bToJpsiOnly?"bToJpsi":"NonPromptJpsi")+(TString)(useFlipJpsi?(flipJpsiSameSide?"_flipJpsiSameSide":"_flipJpsi"):(ispp?"":"_PromptJpsi"))+(TString)(ignoreBin2?"_ignoreBin1-2":(ignoreBin1?"_ignoreBin1":""))+(TString)(ispp?"_pp":"_PbPb")+".root";
  cout<<"Extract normalisations from file "<<normFileName<<endl;
  auto normFile = TFile::Open(normFileName,"READ");
  
  RooArgSet *Yields = (RooArgSet*)normFile->Get("norm_fit_s");
  for(int i=0; i<ntrees; i++){
    if(!usedForFit[i] || i==3) continue;
    cout<<"Recovering postfit BcM shape for process "<<procNameDef[i]<<endl;
    for(int k=0;k<=_nChan(ispp);k++){

      //Clone trimuon mass histos prefit into postfit. Processes not used in fitting will not change histos postfit
      h_BcM_postfit[i][k] = (TH1F*) normFile->Get("shapes_fit_s/BDT"+(TString)(to_string((k==0)?1:k))+"/"+procNameDef[i]);
      h_BcM_postfit[i][k]->SetDirectory(0);

      //inclusive on BDT bins
      if(k>1) h_BcM_postfit[i][0]->Add(h_BcM_postfit[i][k]);

      //total background
      if(k>0 && i!=3 && i!=4){
	if(h_BcM_bkg[k]==NULL || h_BcM_bkg[k]->GetEntries()==0) h_BcM_bkg[k] = (TH1F*) h_BcM_postfit[i][k]->Clone();
	else h_BcM_bkg[k]->Add(h_BcM_postfit[i][k]);}
    }
  }

  //total background, inclusive on BDT bins
  h_BcM_bkg[0] = (TH1F*) h_BcM_bkg[1]->Clone();
  for(int k=2;k<=_nChan(ispp);k++){
    h_BcM_bkg[0]->Add(h_BcM_bkg[k]);
  }
  for(int k=0;k<=_nChan(ispp);k++) h_BcM_bkg[k]->SetDirectory(0);

  normFile->Close();

  //fraction of bkg i in mass bin m: h_BcM_postfit[i][0]->GetBinContent(m) / h_BcM_bkg[0]->GetBinContent(m);

  //********************************************************
  //Get BDT vs M histograms
  //********************************************************
  auto fullFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root");

  //grab tree 
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }

  //variables for branches
  float Bc_M[ntrees];
  float BDT[ntrees];
  float weight[ntrees];
  int flipJpsi[ntrees];

  //init histograms
  double mbins[_nbinM(ispp)+1];
  for(int bx=0;bx<=_nbinM(ispp);bx++){ 
    mbins[bx] = _Mbinning(ispp)[bx];
  }
  double bdtlo = _withTM?-0.47:-1, bdthi = _withTM?(ispp?0.31:0.39):1;
  int nbdt = _withTM?(ispp?50:40):35;
  TH2F* h_BDTM_bkg = new TH2F("Bc_BDTM_bkg", "BDT vs M(B_{c}) background", _nbinM(ispp), mbins, nbdt,bdtlo,bdthi);
  TH2F* h_BDTM_sig = new TH2F("Bc_BDTM_sig", "BDT vs M(B_{c}) signal", _nbinM(ispp), mbins, nbdt,bdtlo,bdthi);
  TH2F* h_BDTM_data = new TH2F("Bc_BDTM_data", "BDT vs M(B_{c}) data", _nbinM(ispp), mbins, nbdt,bdtlo,bdthi);
  float nbkg = 0, ndata = 0, nsig = 0;
  float CRbinwRatio = ((m_Bc-3.3)/_nbinMSR(ispp)) * (_nbinMCR(ispp)/1.);

  //tree and event loop
  for(int iT=0; iT<(int)T.size(); iT++){
    if(!usedForFit[iT]) continue;
    std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    T[iT]->SetBranchAddress("BDT", &BDT[iT]);
    T[iT]->SetBranchAddress("Bc_M", &Bc_M[iT]);
    if(iT==7) T[iT]->SetBranchAddress("flipJpsi", &flipJpsi[iT]);
    T[iT]->SetBranchAddress("weight", &weight[iT]);

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){
      T[iT]->GetEntry(j);
      float w = weight[iT];
      if(Bc_M[iT]>m_Bc) w *= CRbinwRatio;

      if(iT==3){
	h_BDTM_data->Fill(Bc_M[iT], BDT[iT], w);
	ndata += w;}
      else if(iT==4){
	h_BDTM_sig->Fill(Bc_M[iT], BDT[iT], w);
	nsig += w;}
      else{
	int mbin = h_BcM_prefit[iT][0]->FindBin(Bc_M[iT]);
	float yieldCorr = (h_BcM_prefit[iT][0]->GetBinContent(mbin) == 0)?0:( h_BcM_postfit[iT][0]->GetBinContent(mbin) / h_BcM_prefit[iT][0]->GetBinContent(mbin) );

	h_BDTM_bkg->Fill(Bc_M[iT], BDT[iT], w*yieldCorr);
	nbkg += w*yieldCorr;
      }
      
    }//end event loop

  }//end tree loop

  cout<<"total bkg = "<<nbkg<<endl;
  cout<<"signal = "<<nsig<<endl;
  cout<<"data = "<<ndata<<endl;

  //suppress bins with <0 content, for drawing purposes
  for(int b=0;b<=h_BDTM_bkg->GetNbinsX();b++){
    for(int by=0;by<=h_BDTM_bkg->GetNbinsY();by++){
      if(h_BDTM_bkg->GetBinContent(b,by)<0) h_BDTM_bkg->SetBinContent(b,by, 0 );
      if(h_BDTM_sig->GetBinContent(b,by)<0) h_BDTM_sig->SetBinContent(b,by, 0 );
    }
  }

  //********************************************************
  //DRAW BDT vs M 2D
  //********************************************************
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kLightTemperature);
  gStyle->SetNumberContours(40);
  TCanvas* c1 = new TCanvas("c1","c1",3900,1300);
  c1->Divide(3,1);
  c1->cd(1);
  h_BDTM_bkg->GetXaxis()->SetTitle("M(B_{c})");
  h_BDTM_bkg->GetYaxis()->SetTitle("BDT");
  h_BDTM_bkg->Draw("COLZ0");
  c1->cd(2);
  h_BDTM_sig->GetXaxis()->SetTitle("M(B_{c})");
  h_BDTM_sig->GetYaxis()->SetTitle("BDT");
  h_BDTM_sig->Draw("COLZ0");
  c1->cd(3);
  h_BDTM_data->GetXaxis()->SetTitle("M(B_{c})");
  h_BDTM_data->GetYaxis()->SetTitle("BDT");
  h_BDTM_data->Draw("COLZ0");
  
  c1->SaveAs("figs_UncorrBDT/BDTvsM_2D_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  //********************************************************
  //DRAW average BDT vs M
  //********************************************************
  TH1F *avBDTvsM_bkg = (TH1F*)h_BDTM_bkg->ProfileX("avBDTvsM_bkg",1,-1,""); //default error is error on the mean, option g is for weighted average errors 
  TH1F *avBDTvsM_sig = (TH1F*)h_BDTM_sig->ProfileX("avBDTvsM_sig",1,-1,""); 
  TH1F *avBDTvsM_data = (TH1F*)h_BDTM_data->ProfileX("avBDTvsM_data",1,-1,"");
  TCanvas* c2 = new TCanvas("c2","c2",3900,1300);
  c2->Divide(3,1);
  c2->cd(1);
  avBDTvsM_bkg->GetYaxis()->SetTitle("mean BDT value");
  if(!ispp && _withTM) avBDTvsM_bkg->GetYaxis()->SetRangeUser(-0.35, -0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  avBDTvsM_bkg->Draw();
  c2->cd(2);
  avBDTvsM_sig->GetYaxis()->SetTitle("mean BDT value");
  avBDTvsM_sig->GetYaxis()->SetRangeUser(avBDTvsM_sig->GetMinimum()-0.02, avBDTvsM_sig->GetMaximum()+0.02);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  avBDTvsM_sig->Draw();
  c2->cd(3);
  avBDTvsM_data->GetYaxis()->SetTitle("mean BDT value");
  if(!ispp && _withTM) avBDTvsM_data->GetYaxis()->SetRangeUser(-0.35, -0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  avBDTvsM_data->Draw();

  c2->SaveAs("figs_UncorrBDT/BDTvsM_profile_"+(TString)(ispp?"pp":"PbPb")+".pdf");


  //********************************************************
  //Get BDT vs M histograms
  //********************************************************
  TH2F* h_corrBDT_bkg = new TH2F("Bc_corrBDT_bkg", "corrected BDT vs M(B_{c}) background", _nbinM(ispp), mbins, nbdt,_withTM?(ispp?-0.31:-0.26):-0.8,_withTM?(ispp?0.44:0.63):(ispp?1:1.4));
  TH2F* h_corrBDT_sig = new TH2F("Bc_corrBDT_sig", "corrected BDT vs M(B_{c}) signal", _nbinM(ispp), mbins, nbdt,    _withTM?(ispp?-0.31:-0.26):-0.8,_withTM?(ispp?0.44:0.63):(ispp?1:1.4));
  TH2F* h_corrBDT_data = new TH2F("Bc_corrBDT_data", "corrected BDT vs M(B_{c}) data", _nbinM(ispp), mbins, nbdt,    _withTM?(ispp?-0.31:-0.26):-0.8,_withTM?(ispp?0.44:0.63):(ispp?1:1.4));

  //tree and event loop
  for(int iT=0; iT<(int)T.size(); iT++){
    if(!usedForFit[iT]) continue;
    std::cout << "--- Processing again: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){
      T[iT]->GetEntry(j);
      float w = weight[iT];
      if(Bc_M[iT]>m_Bc) w *= CRbinwRatio;

      int mbin = h_BcM_prefit[iT][0]->FindBin(Bc_M[iT]);
      float bdtcorr = avBDTvsM_bkg->GetBinContent(mbin); //correct with the BDT average of expected background
      if(iT==3)
	h_corrBDT_data->Fill(Bc_M[iT], BDT[iT] - bdtcorr, w);
      else if(iT==4)
	h_corrBDT_sig->Fill(Bc_M[iT], BDT[iT] - bdtcorr, w);
      else{
	float yieldCorr = (h_BcM_prefit[iT][0]->GetBinContent(mbin) == 0)?0:( h_BcM_postfit[iT][0]->GetBinContent(mbin) / h_BcM_prefit[iT][0]->GetBinContent(mbin) );

	h_corrBDT_bkg->Fill(Bc_M[iT], BDT[iT] - bdtcorr, w*yieldCorr);
      }
      
    }//end event loop

  }//end tree loop


  //suppress bins with <0 content, for drawing purposes
  for(int b=0;b<=h_corrBDT_bkg->GetNbinsX();b++){
    for(int by=0;by<=h_corrBDT_bkg->GetNbinsY();by++){
      if(h_corrBDT_bkg->GetBinContent(b,by)<0) h_corrBDT_bkg->SetBinContent(b,by, 0 );
      if(h_corrBDT_sig->GetBinContent(b,by)<0) h_corrBDT_sig->SetBinContent(b,by, 0 );
    }
  }


  //********************************************************
  //DRAW corrected BDT vs M 2D
  //********************************************************
  TCanvas* c3 = new TCanvas("c3","c3",3900,1300);
  c3->Divide(3,1);
  c3->cd(1);
  h_corrBDT_bkg->GetXaxis()->SetTitle("M(B_{c})");
  h_corrBDT_bkg->GetYaxis()->SetTitle("BDT corrected");
  h_corrBDT_bkg->Draw("COLZ0");
  c3->cd(2);
  h_corrBDT_sig->GetXaxis()->SetTitle("M(B_{c})");
  h_corrBDT_sig->GetYaxis()->SetTitle("BDT corrected");
  h_corrBDT_sig->Draw("COLZ0");
  c3->cd(3);
  h_corrBDT_data->GetXaxis()->SetTitle("M(B_{c})");
  h_corrBDT_data->GetYaxis()->SetTitle("BDT corrected");
  h_corrBDT_data->Draw("COLZ0");
  
  c3->SaveAs("figs_UncorrBDT/correctedBDTvsM_2D_"+(TString)(ispp?"pp":"PbPb")+".pdf");


  //********************************************************
  //DRAW average corrected BDT vs M
  //********************************************************
  TH1F *avcorrBDTvsM_bkg = (TH1F*)h_corrBDT_bkg->ProfileX("avcorrBDTvsM_bkg",1,-1,""); //default error is error on the mean, option g is for weighted average errors 
  TH1F *avcorrBDTvsM_sig = (TH1F*)h_corrBDT_sig->ProfileX("avcorrBDTvsM_sig",1,-1,""); 
  TH1F *avcorrBDTvsM_data = (TH1F*)h_corrBDT_data->ProfileX("avcorrBDTvsM_data",1,-1,"");
  TCanvas* c4 = new TCanvas("c4","c4",3900,1300);
  c4->Divide(3,1);
  c4->cd(1);
  avcorrBDTvsM_bkg->GetYaxis()->SetTitle("mean corrected BDT");
  avcorrBDTvsM_bkg->GetYaxis()->SetRangeUser(avcorrBDTvsM_bkg->GetMinimum()-0.02, avcorrBDTvsM_bkg->GetMaximum()+0.02);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  avcorrBDTvsM_bkg->Draw();
  c4->cd(2);
  avcorrBDTvsM_sig->GetYaxis()->SetTitle("mean corrected BDT");
  if(!ispp && _withTM) avcorrBDTvsM_sig->GetYaxis()->SetRangeUser(0.1,0.5);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  avcorrBDTvsM_sig->Draw();
  c4->cd(3);
  avcorrBDTvsM_data->GetYaxis()->SetTitle("mean corrected BDT");
  avcorrBDTvsM_data->GetYaxis()->SetRangeUser(avcorrBDTvsM_data->GetMinimum()-0.02, avcorrBDTvsM_data->GetMaximum()+0.02);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  avcorrBDTvsM_data->Draw();

  c4->SaveAs("figs_UncorrBDT/correctedBDTvsM_profile_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  //********************************************************
  //DRAW corrected BDT
  //********************************************************
  TH1F *corrBDT_bkg = (TH1F*)h_corrBDT_bkg->ProjectionY("corrBDT_bkg");
  TH1F *corrBDT_sig = (TH1F*)h_corrBDT_sig->ProjectionY("corrBDT_sig");
  TH1F *corrBDT_data = (TH1F*)h_corrBDT_data->ProjectionY("corrBDT_data");
  TCanvas* c5 = new TCanvas("c5","c5",3900,1300);
  c5->Divide(3,1);
  c5->cd(1);
  corrBDT_bkg->GetXaxis()->SetTitle("corrected BDT");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  corrBDT_bkg->Draw();
  c5->cd(2);
  corrBDT_sig->GetXaxis()->SetTitle("corrected BDT");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  corrBDT_sig->Draw();
  c5->cd(3);
  corrBDT_data->GetXaxis()->SetTitle("corrected BDT");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  corrBDT_data->Draw();

  c5->SaveAs("figs_UncorrBDT/correctedBDT_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  cout<<"poportion of signal MC with corrBDT<-0.1 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(-0.1)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<-0.07 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(-0.07)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<-0.05 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(-0.05)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0. = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.05 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.05)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.1 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.1)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.15 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.15)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.2 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.2)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.22 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.22)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.23 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.23)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.25 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.25)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.3 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.3)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.35 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.35)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.4 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.4)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.45 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.45)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.5 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.5)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.55 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.55)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.65 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.65)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.7 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.7)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.75 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.75)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.8 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.8)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.85 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.85)) / nsig <<endl;
  cout<<"poportion of signal MC with corrBDT<0.9 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.9)) / nsig <<endl;

  TFile *outf = new TFile("BDTuncorrFromM_"+(TString)(ispp?"pp":"PbPb")+".root","RECREATE"); 
  avBDTvsM_bkg->Write();
  outf->Close();
}
