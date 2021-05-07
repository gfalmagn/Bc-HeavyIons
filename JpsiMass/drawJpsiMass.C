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
#include "../helpers/SgMuonAcceptanceCuts.h"

void drawJpsiMass(bool ispp=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  auto fullFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  int ntrees = 9;

  TTree* T = (TTree*)fullFile->Get("bkgBCMASS");
  TTree* Ts = (TTree*)fullFile->Get("signal_MC");
  TTree* Td = (TTree*)fullFile->Get("sigRegion");
  TTree* Tj = (TTree*)fullFile->Get("flipJpsi");

  TH1F* hd_tightEta = new TH1F("hd_tightEta",";dimuon mass [GeV];N", (int)( (2*JpeakWidT+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLoT-JpeakLoBuf-JpeakWidT/2,m_Jpsi+JpeakHiT+JpeakHiBuf+JpeakWidT/2);
  TH1F* hd_looseEta = new TH1F("hd_looseEta",";dimuon mass [GeV];N",  (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);
  TH1F* hd_tightEta_ambig = new TH1F("hd_tightEta_ambig",";dimuon mass [GeV];N", (int)( (2*JpeakWidT+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLoT-JpeakLoBuf-JpeakWidT/2,m_Jpsi+JpeakHiT+JpeakHiBuf+JpeakWidT/2);
  TH1F* hd_looseEta_ambig = new TH1F("hd_looseEta_ambig",";dimuon mass [GeV];N",  (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);

  TH1F* hj_tightEta = new TH1F("hj_tightEta",";dimuon mass [GeV];N", (int)( (2*JpeakWidT+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLoT-JpeakLoBuf-JpeakWidT/2,m_Jpsi+JpeakHiT+JpeakHiBuf+JpeakWidT/2);
  TH1F* hj_looseEta = new TH1F("hj_looseEta",";dimuon mass [GeV];N",  (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);
  TH1F* hj_tightEta_ambig = new TH1F("hj_tightEta_ambig",";dimuon mass [GeV];N", (int)( (2*JpeakWidT+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLoT-JpeakLoBuf-JpeakWidT/2,m_Jpsi+JpeakHiT+JpeakHiBuf+JpeakWidT/2);
  TH1F* hj_looseEta_ambig = new TH1F("hj_looseEta_ambig",";dimuon mass [GeV];N",  (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);

  TH1F* hs_tightEta = new TH1F("hs_tightEta",";dimuon mass [GeV];N", (int)( (2*JpeakWidT+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLoT-JpeakLoBuf-JpeakWidT/2,m_Jpsi+JpeakHiT+JpeakHiBuf+JpeakWidT/2);
  TH1F* hs_looseEta = new TH1F("hs_looseEta",";dimuon mass [GeV];N",  (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);
  TH1F* hs_tightEta_ambig = new TH1F("hs_tightEta_ambig",";dimuon mass [GeV];N", (int)( (2*JpeakWidT+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLoT-JpeakLoBuf-JpeakWidT/2,m_Jpsi+JpeakHiT+JpeakHiBuf+JpeakWidT/2);
  TH1F* hs_looseEta_ambig = new TH1F("hs_looseEta_ambig",";dimuon mass [GeV];N",  (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);

  float Bc_M;
  float QQ_M;
  float Bc_Pt;
  float Bc_Y;
  float weight;
  float w_simple;
  float BDT;
  float muW_eta;
  float mupl_eta;
  float mumi_eta;
  T->SetBranchAddress("Bc_M", &Bc_M);
  T->SetBranchAddress("QQ_M", &QQ_M);
  T->SetBranchAddress("Bc_Pt", &Bc_Pt);
  T->SetBranchAddress("Bc_Y", &Bc_Y);
  T->SetBranchAddress("BDT", &BDT);
  T->SetBranchAddress("muW_eta", &muW_eta);
  T->SetBranchAddress("mupl_eta", &mupl_eta);
  T->SetBranchAddress("mumi_eta", &mumi_eta);
  T->SetBranchAddress("weight", &weight);
  T->SetBranchAddress("w_simple", &w_simple);

  float dBc_M;
  float dQQ_M;
  float dBc_Pt;
  float dBc_Y;
  float dweight;
  float dw_simple;
  float dBDT;
  float dmuW_eta;
  float dmupl_eta;
  float dmumi_eta;
  Td->SetBranchAddress("QQ_M", &dQQ_M);
  Td->SetBranchAddress("Bc_M", &dBc_M);
  Td->SetBranchAddress("Bc_Pt", &dBc_Pt);
  Td->SetBranchAddress("Bc_Y", &dBc_Y);
  Td->SetBranchAddress("BDT", &dBDT);
  Td->SetBranchAddress("weight", &dweight);
  Td->SetBranchAddress("w_simple", &dw_simple);
  Td->SetBranchAddress("muW_eta", &dmuW_eta);
  Td->SetBranchAddress("mupl_eta", &dmupl_eta);
  Td->SetBranchAddress("mumi_eta", &dmumi_eta);

  float jBc_M;
  float jQQ_M;
  float jBc_Pt;
  float jBc_Y;
  float jweight;
  float jw_simple;
  float jBDT;
  float jmuW_eta;
  float jmupl_eta;
  float jmumi_eta;
  Tj->SetBranchAddress("QQ_M", &jQQ_M);
  Tj->SetBranchAddress("Bc_M", &jBc_M);
  Tj->SetBranchAddress("Bc_Pt", &jBc_Pt);
  Tj->SetBranchAddress("Bc_Y", &jBc_Y);
  Tj->SetBranchAddress("BDT", &jBDT);
  Tj->SetBranchAddress("weight", &jweight);
  Tj->SetBranchAddress("w_simple", &jw_simple);
  Tj->SetBranchAddress("muW_eta", &jmuW_eta);
  Tj->SetBranchAddress("mupl_eta", &jmupl_eta);
  Tj->SetBranchAddress("mumi_eta", &jmumi_eta);

  float sBc_M;
  float sQQ_M;
  float sBc_Pt;
  float sBc_Y;
  float sweight;
  float sw_simple;
  float sBDT;
  float smuW_eta;
  float smupl_eta;
  float smumi_eta;
  Ts->SetBranchAddress("QQ_M", &sQQ_M);
  Ts->SetBranchAddress("Bc_M", &sBc_M);
  Ts->SetBranchAddress("Bc_Pt", &sBc_Pt);
  Ts->SetBranchAddress("Bc_Y", &sBc_Y);
  Ts->SetBranchAddress("BDT", &sBDT);
  Ts->SetBranchAddress("weight", &sweight);
  Ts->SetBranchAddress("w_simple", &sw_simple);
  Ts->SetBranchAddress("muW_eta", &smuW_eta);
  Ts->SetBranchAddress("mupl_eta", &smupl_eta);
  Ts->SetBranchAddress("mumi_eta", &smumi_eta);
 
  bool wFidCuts = false;
  //*******************************************
  //Fill the histograms
  //Jpsi sidebands
  for(int j=0; j<T->GetEntries(); j++){
    
    T->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, Bc_Pt, Bc_Y)) continue; //fiducial cuts
    
    if(Bc_M>6.2 || Bc_M<3.5) continue; //mass SR for dimuon mass plot
    if(max(fabs(mupl_eta), max(fabs(mumi_eta),fabs(muW_eta))) <1.5){
      hd_tightEta->Fill(QQ_M,weight);
      if(weight/w_simple!=1)
	hd_tightEta_ambig->Fill(QQ_M,weight);
    }
    else{
      hd_looseEta->Fill(QQ_M,weight);
      if(weight/w_simple!=1)
	hd_looseEta_ambig->Fill(QQ_M,weight);
    }
  }

  //data SR
  for(int j=0; j<Td->GetEntries(); j++){
    
    Td->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, dBc_Pt, dBc_Y)) continue; //fiducial cuts
    
    if(dBc_M>6.2 || dBc_M<3.5) continue; //mass SR for dimuon mass plot
    if(max(fabs(dmupl_eta), max(fabs(dmumi_eta),fabs(dmuW_eta))) <1.5){
      hd_tightEta->Fill(dQQ_M,dweight);
      if(dw_simple!=0 && dweight!=1 && dweight!=4)
	hd_tightEta_ambig->Fill(dQQ_M,dweight);
    }
    else{
      hd_looseEta->Fill(dQQ_M,dweight);
      if(dw_simple!=0 && dweight!=1 && dweight!=4)
	hd_looseEta_ambig->Fill(dQQ_M,dweight);
    }
  }

  //flipped Jpsi
  for(int j=0; j<Tj->GetEntries(); j++){
    
    Tj->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, jBc_Pt, jBc_Y)) continue; //fiducial cuts
    
    if(jBc_M>6.2 || jBc_M<3.5) continue; //mass SR for dimuon mass plot
    if(max(fabs(jmupl_eta), max(fabs(jmumi_eta),fabs(jmuW_eta))) <1.5){
      hj_tightEta->Fill(jQQ_M,fabs(jw_simple));
      if(jweight/jw_simple!=1)
  	hj_tightEta_ambig->Fill(jQQ_M,fabs(jw_simple));
    }
    else{
      hj_looseEta->Fill(jQQ_M,fabs(jw_simple));
      if(jweight/jw_simple!=1)
  	hj_looseEta_ambig->Fill(jQQ_M,fabs(jw_simple));
    }
  }

  //signal MC
  for(int j=0; j<Ts->GetEntries(); j++){
    
    Ts->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, sBc_Pt, sBc_Y)) continue; //fiducial cuts
    
    if(sBc_M>6.2 || sBc_M<3.5) continue; //mass SR for dimuon mass plot
    if(max(fabs(smupl_eta), max(fabs(smumi_eta),fabs(smuW_eta))) <1.5){
      hs_tightEta->Fill(sQQ_M,fabs(sweight));
      if(sweight/sw_simple!=1)
  	hs_tightEta_ambig->Fill(sQQ_M,fabs(sweight));
    }
    else{
      hs_looseEta->Fill(sQQ_M,fabs(sweight));
      if(sweight/sw_simple!=1)
  	hs_looseEta_ambig->Fill(sQQ_M,fabs(sweight));
    }
  }



  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","cdata",2200,1100);
  c1->Divide(2,1);

  c1->cd(1);
  hd_looseEta_ambig->SetLineColor(kRed+2);
  hd_looseEta->SetLineWidth(3);
  hd_looseEta_ambig->SetLineWidth(3);
  TH1F* hd_looseEta_ambig_sc = (TH1F*)hd_looseEta_ambig->Clone("hd_looseEta_ambig");
  hd_looseEta_ambig_sc->Scale(hd_looseEta->Integral()/hd_looseEta_ambig->Integral());
  hd_looseEta_ambig_sc->SetLineColor(kOrange+8);
  hd_looseEta_ambig_sc->SetMinimum(0);
  hd_looseEta_ambig_sc->Draw("");
  hd_looseEta_ambig->Draw("same");
  hd_looseEta->Draw("same");

  TLegend *legd1 = new TLegend(0.6,0.7,0.9,0.9);
  legd1->SetBorderSize(0);
  legd1->SetTextSize(0.035);
  legd1->SetHeader("loose range (high |#eta|)");
  legd1->AddEntry(hd_looseEta,"all");
  legd1->AddEntry(hd_looseEta_ambig,"ambiguous J/#psi choice");
  legd1->AddEntry(hd_looseEta_ambig_sc,"ambiguous, scaled");
  legd1->Draw("same");

  c1->cd(2);
  hd_tightEta_ambig->SetLineColor(kRed+2);
  hd_tightEta->SetLineWidth(3);
  hd_tightEta_ambig->SetLineWidth(3);
  TH1F* hd_tightEta_ambig_sc = (TH1F*)hd_tightEta_ambig->Clone("hd_tightEta_ambig");
  hd_tightEta_ambig_sc->Scale(hd_tightEta->Integral()/hd_tightEta_ambig->Integral());
  hd_tightEta_ambig_sc->SetLineColor(kOrange+8);
  hd_tightEta_ambig_sc->SetMinimum(0);
  hd_tightEta_ambig_sc->Draw("");
  hd_tightEta_ambig->Draw("same");
  hd_tightEta->Draw("same");

  TLegend *legd2 = new TLegend(0.6,0.7,0.9,0.9);
  legd2->SetBorderSize(0);
  legd2->SetTextSize(0.035);
  legd2->SetHeader("tight range (central |#eta|)");
  legd2->AddEntry(hd_tightEta,"all");
  legd2->AddEntry(hd_tightEta_ambig,"ambiguous J/#psi choice");
  legd2->AddEntry(hd_tightEta_ambig_sc,"ambiguous, scaled");
  legd2->Draw("same");

  c1->SaveAs("dimuonMass_ambiguousJpsiChoice_data_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  TCanvas *c3 = new TCanvas("c3","csigMC",2200,1100);
  c3->Divide(2,1);

  c3->cd(1);
  hs_looseEta_ambig->SetLineColor(kRed+2);
  hs_looseEta->SetLineWidth(3);
  hs_looseEta_ambig->SetLineWidth(3);
  TH1F* hs_looseEta_ambig_sc = (TH1F*)hs_looseEta_ambig->Clone("hs_looseEta_ambig");
  hs_looseEta_ambig_sc->Scale(hs_looseEta->Integral()/hs_looseEta_ambig->Integral());
  hs_looseEta_ambig_sc->SetLineColor(kOrange+8);
  hs_looseEta_ambig_sc->SetMinimum(0);
  hs_looseEta_ambig_sc->Draw("");
  hs_looseEta_ambig->Draw("same");
  hs_looseEta->Draw("same");
  legd2->Draw("same");

  c3->cd(2);
  hs_tightEta_ambig->SetLineColor(kRed+2);
  hs_tightEta->SetLineWidth(3);
  hs_tightEta_ambig->SetLineWidth(3);
  TH1F* hs_tightEta_ambig_sc = (TH1F*)hs_tightEta_ambig->Clone("hs_tightEta_ambig");
  hs_tightEta_ambig_sc->Scale(hs_tightEta->Integral()/hs_tightEta_ambig->Integral());
  hs_tightEta_ambig_sc->SetLineColor(kOrange+8);
  hs_tightEta_ambig_sc->SetMinimum(0);
  hs_tightEta_ambig_sc->Draw("");
  hs_tightEta_ambig->Draw("same");
  hs_tightEta->Draw("same");
  legd2->Draw("same");

  c3->SaveAs("dimuonMass_ambiguousJpsiChoice_signalMC_"+(TString)(ispp?"pp":"PbPb")+".pdf");


  TCanvas *c2 = new TCanvas("c2","cflipJ",2200,1100);
  c2->Divide(2,1);

  c2->cd(1);
  hj_looseEta_ambig->SetLineColor(kRed+2);
  hj_looseEta->SetLineWidth(3);
  hj_looseEta_ambig->SetLineWidth(3);
  TH1F* hj_looseEta_ambig_sc = (TH1F*)hj_looseEta_ambig->Clone("hj_looseEta_ambig");
  hj_looseEta_ambig_sc->Scale(hj_looseEta->Integral()/hj_looseEta_ambig->Integral());
  hj_looseEta_ambig_sc->SetLineColor(kOrange+8);
  hj_looseEta_ambig_sc->SetMinimum(0);
  hj_looseEta_ambig_sc->Draw("");
  hj_looseEta_ambig->Draw("same");
  hj_looseEta->Draw("same");
  legd1->Draw("same");

  c2->cd(2);
  hj_tightEta_ambig->SetLineColor(kRed+2);
  hj_tightEta->SetLineWidth(3);
  hj_tightEta_ambig->SetLineWidth(3);
  TH1F* hj_tightEta_ambig_sc = (TH1F*)hj_tightEta_ambig->Clone("hj_tightEta_ambig");
  hj_tightEta_ambig_sc->Scale(hj_tightEta->Integral()/hj_tightEta_ambig->Integral());
  hj_tightEta_ambig_sc->SetLineColor(kOrange+8);
  hj_tightEta_ambig_sc->SetMinimum(0);
  hj_tightEta_ambig_sc->Draw("");
  hj_tightEta_ambig->Draw("same");
  hj_tightEta->Draw("same");
  legd2->Draw("same");

  c2->SaveAs("dimuonMass_ambiguousJpsiChoice_flipJpsi_"+(TString)(ispp?"pp":"PbPb")+".pdf");



}
