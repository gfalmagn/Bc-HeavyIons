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
#include "TROOT.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"
#include "../helpers/SgMuonAcceptanceCuts.h"

void MakePositive(TH2F* h){
  for(int i=1;i<=h->GetNbinsX(); i++){
    for(int j=1;j<=h->GetNbinsX(); j++){
      if(h->GetBinContent(i,j)<=0) h->SetBinContent(i,j,0);
    }
  }
  h->SetMinimum(0);
}

void visibleVsGenKinematics(bool ispp=true, bool secondStep=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  bool wFidCuts = false;
  bool w2ndstepPtWeight = true;

  //Grab the variations of the pT bias of MC, from first step r1 and r2
  TFile *BiasFile = TFile::Open("../twoSteps/pTBiases.root","READ");
  TH1F* bias3muPt = (TH1F*)BiasFile->Get("pTbias_"+(TString)(ispp?"pp":"PbPb")+"_var"+(TString)to_string(1)+(TString)(secondStep?"_2ndStep":"")); //variation 1 is the nominal for: p_T^{n+m log(p_T)} , fixed m
  
  //Get tree branches
  auto fullFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  int ntrees = 9;

  TTree* Ts = (TTree*)fullFile->Get("signal_MC");
  
  TH1F* Rec3muPt = new TH1F("Rec3muPt",";p_{T} [GeV];N_{sig}",70,0,50);
  TH1F* GenBcPt = new TH1F("GenBcPt",";p_{T} [GeV];N_{sig}",70,0,50);
  TH2F* PtGenVisRatio_vsPt = new TH2F("PtGenVisRatio_vsPt",";p_{T}(trimuon) [GeV];(gen B_{c} p_{T}) / (reco trimuon p_{T})",25,4.5,35,81,0.5,2.3);
  TH2F* PtVisGenRatio_vsM = new TH2F("PtVisGenRatio_vsM",";m(trimuon) [GeV];(reco trimuon p_{T}) / (gen B_{c} p_{T})",27,_mBcMin,_mBcMax,50,0.3,1.3);

  float sgenBc_Pt;
  float sgen3mu_Pt;
  float sBc_M;
  float sQQ_M;
  float sQQ_dca;
  float sBc_Pt;
  float sBc_Y;
  float sBc_VtxProb;
  float sweight;
  float sBDT;
  float smuW_eta;
  float smupl_eta;
  float smumi_eta;
  float smuW_pt;
  float smupl_pt;
  float smumi_pt;
  Ts->SetBranchAddress("QQ_M", &sQQ_M);
  Ts->SetBranchAddress("QQ_dca", &sQQ_dca);
  Ts->SetBranchAddress("Bc_M", &sBc_M);
  Ts->SetBranchAddress("genBc_Pt", &sgenBc_Pt);
  Ts->SetBranchAddress("gen3mu_Pt", &sgen3mu_Pt);
  Ts->SetBranchAddress("Bc_Pt", &sBc_Pt);
  Ts->SetBranchAddress("Bc_Y", &sBc_Y);
  Ts->SetBranchAddress("Bc_VtxProb", &sBc_VtxProb);
  Ts->SetBranchAddress("BDT", &sBDT);
  Ts->SetBranchAddress("weight", &sweight);
  Ts->SetBranchAddress("muW_eta", &smuW_eta);
  Ts->SetBranchAddress("mupl_eta", &smupl_eta);
  Ts->SetBranchAddress("mumi_eta", &smumi_eta);
  Ts->SetBranchAddress("muW_Pt", &smuW_pt);
  Ts->SetBranchAddress("mupl_Pt", &smupl_pt);
  Ts->SetBranchAddress("mumi_Pt", &smumi_pt);

  //*******************************************
  //Fill the histograms
  //signal MC
  for(int j=0; j<Ts->GetEntries(); j++){
    
    Ts->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, sBc_Pt, sBc_Y)) continue; //fiducial cuts

    if(w2ndstepPtWeight)
      sweight *= getBias(bias3muPt,sBc_Pt);

    Rec3muPt->Fill(sBc_Pt, sweight);
    GenBcPt->Fill(sgenBc_Pt, sweight);
    PtGenVisRatio_vsPt->Fill(sBc_Pt, sgenBc_Pt/sBc_Pt , sweight);
    PtVisGenRatio_vsM->Fill(sBc_M, sBc_Pt/sgenBc_Pt, sweight);
  }

  gStyle->SetOptStat(0);

  //Reco3muPt_and_GenBcPt
  TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
  c1->SetTopMargin(0.03);
  c1->SetRightMargin(0.05);
  Rec3muPt->SetLineWidth(3);
  GenBcPt->SetLineWidth(3);
  GenBcPt->SetLineColor(kRed);
  Rec3muPt->Draw("hist");
  GenBcPt->Draw("histsame");

  TLegend *leg1 = new TLegend(0.55,0.74,0.94,0.96);
  leg1->SetHeader(ispp?"pp":"PbPb");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(Rec3muPt, "reco (visible) p_{T}(#mu#mu#mu)");
  leg1->AddEntry(GenBcPt, "generated p_{T}(B_{c})");
  leg1->Draw("same");

  c1->SaveAs("figs/Reco3muPt_and_GenBcPt_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  //PtVisGenRatio_vsM
  TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
  c2->SetLeftMargin(0.12);
  MakePositive(PtVisGenRatio_vsM);
  PtVisGenRatio_vsM->GetYaxis()->SetRangeUser(0.45,1.15);
  PtVisGenRatio_vsM->Draw("COLZ");
  TProfile* PtVisGenRatio_mean = PtVisGenRatio_vsM->ProfileX();
  PtVisGenRatio_mean->SetLineColor(kRed+1);
  PtVisGenRatio_mean->SetLineWidth(3);
  PtVisGenRatio_mean->Draw("histsame][");

  c2->SaveAs("figs/PtVis_To_PtGen_RatioVsM_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  //PtVisGenRatio_vsPt
  TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
  c3->SetLeftMargin(0.12);
  MakePositive(PtGenVisRatio_vsPt);
  PtGenVisRatio_vsPt->GetYaxis()->SetRangeUser(0.8,1.85);
  PtGenVisRatio_vsPt->Draw("COLZ");
  TProfile* PtGenVisRatio_mean = PtGenVisRatio_vsPt->ProfileX();
  PtGenVisRatio_mean->SetLineColor(kRed+1);
  PtGenVisRatio_mean->SetLineWidth(3);
  PtGenVisRatio_mean->Draw("histsame][");

  c3->SaveAs("figs/PtGen_To_PtVis_RatioVsPtVis_"+(TString)(ispp?"pp":"PbPb")+".pdf");

}
