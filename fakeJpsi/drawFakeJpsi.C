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

void drawFakeJpsi(bool ispp=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  auto fullFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  int ntrees = 9;

  TTree* T = (TTree*)fullFile->Get("bkgBCMASS");
  TTree* Td = (TTree*)fullFile->Get("sigRegion");
  TTree* Tnp = (TTree*)fullFile->Get("bToJpsi_MC");

  TH1F* h_bToJ = new TH1F("h_bToJ",";trimuon mass [GeV];N", 28,3.3,5.75);
  TH1F* h_lowSB = new TH1F("h_lowSB",";trimuon mass [GeV];N", 20,3.3,7.3);
  TH1F* h_hiSB = new TH1F("h_hiSB",";trimuon mass [GeV];N", 20,3.3,7.3);
  TH1F* h_tightEta = new TH1F("h_tightEta",";dimuon mass [GeV];N", (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);
  TH1F* h_looseEta = new TH1F("h_looseEta",";dimuon mass [GeV];N",  (int)( (2*JpeakWid+JpeakLoBuf+JpeakHiBuf+1e-5)/0.01), m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2,m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2);

  float Bc_M;
  float QQ_M;
  float Bc_Pt;
  float Bc_Y;
  float weight;
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

  // float dBc_M;
  // float dBc_Pt;
  // float dBc_Y;
  // float dweight;
  // float dBDT;
  // T->SetBranchAddress("dBc_M", &dBc_M);
  // T->SetBranchAddress("dBc_Pt", &dBc_Pt);
  // T->SetBranchAddress("dBc_Y", &dBc_Y);
  // T->SetBranchAddress("dBDT", &dBDT);
  // T->SetBranchAddress("dweight", &dweight);

  float npBc_M;
  float npBc_Pt;
  float npBc_Y;
  float npweight;
  float npBDT;
  bool muW_isJpsiBro;
  Tnp->SetBranchAddress("Bc_M", &npBc_M);
  Tnp->SetBranchAddress("Bc_Pt", &npBc_Pt);
  Tnp->SetBranchAddress("Bc_Y", &npBc_Y);
  Tnp->SetBranchAddress("BDT", &npBDT);
  Tnp->SetBranchAddress("weight", &npweight);
  Tnp->SetBranchAddress("muW_isJpsiBro", &muW_isJpsiBro);

  float dBc_M;
  float dQQ_M;
  float dBc_Pt;
  float dBc_Y;
  float dweight;
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
  Td->SetBranchAddress("muW_eta", &dmuW_eta);
  Td->SetBranchAddress("mupl_eta", &dmupl_eta);
  Td->SetBranchAddress("mumi_eta", &dmumi_eta);

  bool wFidCuts = false;
  //*******************************************
  //Fill the histograms
  //Jpsi sidebands
  for(int j=0; j<T->GetEntries(); j++){
    
    T->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, Bc_Pt, Bc_Y)) continue; //fiducial cuts
    
    if(QQ_M<3.1)
      h_lowSB->Fill(Bc_M,weight);        
    else
      h_hiSB->Fill(Bc_M,weight);        

    if(Bc_M>6.2 || Bc_M<3.5) continue; //mass SR for dimuon mass plot
    if(max(fabs(mupl_eta), max(fabs(mumi_eta),fabs(muW_eta))) <1.5)
      h_tightEta->Fill(QQ_M,weight);
    else
      h_looseEta->Fill(QQ_M,weight);
  }

  //data SR
  for(int j=0; j<Td->GetEntries(); j++){
    
    Td->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, dBc_Pt, dBc_Y)) continue; //fiducial cuts

    if(dBc_M>6.2 || dBc_M<3.5) continue; //mass SR for dimuon mass plot
    if(max(fabs(dmupl_eta), max(fabs(dmumi_eta),fabs(dmuW_eta))) <1.5)
      h_tightEta->Fill(dQQ_M,dweight);
    else
      h_looseEta->Fill(dQQ_M,dweight);
  }

  //non-prompt MC
  for(int j=0; j<Tnp->GetEntries(); j++){
    
    Tnp->GetEntry(j);
    if(wFidCuts) 
      if(!inFidCuts(0, npBc_Pt, npBc_Y)) continue; //fiducial cuts

    if(muW_isJpsiBro)
      h_bToJ->Fill(npBc_M,npweight);    
  }


  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
  h_bToJ->SetLineWidth(3);
  h_bToJ->Draw("Ehist");
  c1->SetTopMargin(0.03);
  c1->SaveAs("figs/BcM_trueBtoJpsi_"+(TString)(ispp?"pp":"PbPb")+".pdf");



  TCanvas *c2 = new TCanvas("c2","c2",3000,1500);
  c2->Divide(2,1);
  h_tightEta->SetLineWidth(2);
  h_looseEta->SetLineWidth(2);

  c2->cd(1)->SetTopMargin(0.03);
  h_tightEta->Draw("E");
  TLatex t1(.7,.8,"low |#eta|");  
  t1.SetNDC(kTRUE);
  t1.SetTextSize(0.05);
  t1.Draw("same");

  c2->cd(2)->SetTopMargin(0.03);
  h_looseEta->Draw("E");
  TLatex t2(.7,.8,"high |#eta|");  
  t2.SetNDC(kTRUE);
  t2.SetTextSize(0.05);
  t2.Draw("same");

  c2->SaveAs("figs/JpsiMass_"+(TString)(ispp?"pp":"PbPb")+".pdf");



  gStyle->SetOptStat("emr");
  gStyle->SetStatY(0.97); 
  gStyle->SetStatH(0.2);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.3);

  TCanvas *c3 = new TCanvas("c3","c3",3000,1500);
  c3->Divide(2,1);
  c3->SetTopMargin(0.03);
  h_lowSB->SetLineWidth(3);
  h_hiSB->SetLineWidth(3);

  c3->cd(1)->SetTopMargin(0.03);
  h_lowSB->Draw("histE");
  TLatex t3(.5,.7,"lower sideband");  
  t3.SetNDC(kTRUE);
  t3.SetTextSize(0.05);
  t3.Draw("same");

  c3->cd(2)->SetTopMargin(0.03);
  h_hiSB->Draw("histE");
  TLatex t4(.5,.7,"upper sideband");  
  t4.SetNDC(kTRUE);
  t4.SetTextSize(0.05);
  t4.Draw("same");

  c3->SaveAs("figs/BcM_separateJpsiSidebands_"+(TString)(ispp?"pp":"PbPb")+".pdf");

}
