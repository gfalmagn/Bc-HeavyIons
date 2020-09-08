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

void drawFlipJpsi(bool ispp=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  auto fullFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  int ntrees = 9;

  TTree* T = (TTree*)fullFile->Get("flipJpsi");
  TTree* Td = (TTree*)fullFile->Get("sigRegion");

  TH1F* h_yields = new TH1F("h_yields",(TString)(ispp?"pp":"PbPb")+" yields;J/#psi rotation angle;counts", 9,0,9);
  TH1F* h_oppEta = new TH1F("h_oppEta",";trimuon mass [GeV];normalised to 1", ispp?20:16,3,7.3);
  TH1F* h_samEta = new TH1F("h_samEta",";trimuon mass [GeV];normalised to 1", ispp?20:16,3,7.3);
  TH1F* h_oppPhi = new TH1F("h_oppPhi",";trimuon mass [GeV];normalised to 1", ispp?20:16,3,7.3);
  TH1F* h_otherPhi = new TH1F("h_otherPhi",";trimuon mass [GeV];normalised to 1", ispp?20:16,3,7.3);
  TH1F* h_BDT = new TH1F("h_BDT",";BDT;normalised to 1", ispp?20:16,-0.6,0.6);
  TH1F* h_dBDT = new TH1F("h_dBDT",";BDT;normalised to 1", ispp?20:16,-0.6,0.6);

  float Bc_M;
  float Bc_Pt;
  float Bc_Y;
  int flipJpsi;
  float weight;
  float BDT;
  T->SetBranchAddress("Bc_M", &Bc_M);
  T->SetBranchAddress("Bc_Pt", &Bc_Pt);
  T->SetBranchAddress("Bc_Y", &Bc_Y);
  T->SetBranchAddress("BDT", &BDT);
  T->SetBranchAddress("weight", &weight);
  T->SetBranchAddress("flipJpsi", &flipJpsi);

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

  bool wFidCuts = false;
  bool highMass = false;
  //*******************************************
  //Fill the histograms
  for(int j=0; j<T->GetEntries(); j++){
    
    T->GetEntry(j);
    if(wFidCuts) 
      if(!( Bc_Pt>_BcPtmin[0] && Bc_Pt<_BcPtmax[0] && fabs(Bc_Y)<_BcYmax[0] && (Bc_Pt>_BcPtmin[1] || fabs(Bc_Y)>_BcYmin[2]) )) continue; //fiducial cuts
    if(highMass)
      if(Bc_M<6.2) continue;
    
    h_yields->Fill(flipJpsi+1e-3,weight*7);

    if(flipJpsi<5) h_oppEta->Fill(Bc_M,weight*7);
    else h_samEta->Fill(Bc_M,weight*7);
    if(flipJpsi==1 || flipJpsi==6) h_oppPhi->Fill(Bc_M,weight*7);
    else h_otherPhi->Fill(Bc_M,weight*7);

    h_BDT->Fill(BDT,weight);        
  }

  // //data
  // for(int j=0; j<Td->GetEntries(); j++){
    
  //   Td->GetEntry(j);
  //   if(wFidCuts) 
  //     if(!( dBc_Pt>_BcPtmin[0] && dBc_Pt<_BcPtmax[0] && fabs(dBc_Y)<_BcYmax[0] && (dBc_Pt>_BcPtmin[1] || fabs(dBc_Y)>_BcYmin[2]) )) continue; //fiducial cuts
  //   if(highMass)
  //     if(dBc_M<6.2) continue;
    
  //   h_dBDT->Fill(dBDT,dweight);
  // }

  gStyle->SetOptStat(0);

  h_yields->SetTitleOffset(0.2);
  h_yields->SetLineWidth(4);
  h_yields->GetXaxis()->SetBinLabel(1,"");
  h_yields->GetXaxis()->SetBinLabel(2,"#eta #rightarrow #minus#eta, #Delta#phi = #pi");
  h_yields->GetXaxis()->SetBinLabel(3,"#eta #rightarrow #minus#eta, #Delta#phi = #frac{#pi}{2}");
  h_yields->GetXaxis()->SetBinLabel(4,"#eta #rightarrow #minus#eta, #Delta#phi = 0");
  h_yields->GetXaxis()->SetBinLabel(5,"#eta #rightarrow #minus#eta, #Delta#phi = #minus#frac{#pi}{2}");
  h_yields->GetXaxis()->SetBinLabel(6,"#Delta#eta = 0, #Delta#phi = #frac{#pi}{2}");
  h_yields->GetXaxis()->SetBinLabel(7,"#Delta#eta = 0, #Delta#phi = #pi");
  h_yields->GetXaxis()->SetBinLabel(8,"#Delta#eta = 0, #Delta#phi = #minus#frac{#pi}{2}");
  h_yields->GetXaxis()->SetBinLabel(9,"");
  h_yields->GetXaxis()->SetTitleOffset(-0.7);
  h_yields->GetXaxis()->SetTitleSize(0.046);
  h_yields->GetXaxis()->SetLabelSize(0.06);
  h_yields->LabelsOption("v");//vertical bin labels

  TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
  gPad->SetBottomMargin(0.28);
  gPad->SetTopMargin(0.06);
  gPad->SetLeftMargin(0.125);
  h_yields->Draw("E");

  c1->SaveAs("flipJpsiYields_"+(TString)(highMass?"highMassRegion_":"")+(TString)(ispp?"pp":"PbPb")+".pdf");

  TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
  h_oppEta->SetLineWidth(3);
  h_oppEta->Scale(1/h_oppEta->Integral());
  h_samEta->SetLineWidth(3);
  h_samEta->Scale(1/h_samEta->Integral());
  h_oppEta->SetLineColor(kRed);

  h_oppEta->Draw("hist");
  h_samEta->Draw("histsame");

  TLegend *leg2 = new TLegend(0.6,0.76,0.9,0.96);
  leg2->AddEntry(h_oppEta,"opposite #eta flipping");
  leg2->AddEntry(h_samEta,"same #eta flipping");
  leg2->SetTextSize(0.04);
  leg2->Draw();

  c2->SetLeftMargin(0.12);
  c2->SetTopMargin(0.04);
  c2->SaveAs("BcM_flipJpsi_"+(TString)(highMass?"highMassRegion_":"")+"OppositeVsSameEta_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
  h_oppPhi->SetLineWidth(3);
  h_oppPhi->Scale(1/h_oppPhi->Integral());
  h_otherPhi->SetLineWidth(3);
  h_otherPhi->Scale(1/h_otherPhi->Integral());
  h_oppPhi->SetLineColor(kRed);

  h_oppPhi->GetYaxis()->SetRangeUser(0,1.05*max(h_otherPhi->GetMaximum(),h_oppPhi->GetMaximum()));
  h_oppPhi->Draw("hist");
  h_otherPhi->Draw("histsame");

  TLegend *leg3 = new TLegend(0.6,0.76,0.9,0.96);
  leg3->AddEntry(h_oppPhi,"opposite #phi flipping");
  leg3->AddEntry(h_otherPhi,"other #phi flippings");
  leg3->SetTextSize(0.04);
  leg3->Draw();

  c3->SetLeftMargin(0.12);
  c3->SetTopMargin(0.04);
  c3->SaveAs("BcM_flipJpsi_"+(TString)(highMass?"highMassRegion_":"")+"OppositeVsOtherPhi_"+(TString)(ispp?"pp":"PbPb")+".pdf");


  // TCanvas *c4 = new TCanvas("c4","c4",1500,1500);
  // h_BDT->SetLineWidth(3);
  // h_BDT->Scale(1/h_BDT->Integral());
  // h_dBDT->SetLineWidth(3);
  // h_dBDT->Scale(1/h_dBDT->Integral());
  // h_BDT->SetLineColor(kRed);

  // h_BDT->GetYaxis()->SetRangeUser(0,1.05*max(h_dBDT->GetMaximum(),h_BDT->GetMaximum()));
  // h_BDT->Draw("hist");
  // h_dBDT->Draw("histsame");

  // TLegend *leg4 = new TLegend(0.6,0.76,0.9,0.96);
  // leg4->AddEntry(h_BDT,"flipJpsi");
  // leg4->AddEntry(h_dBDT,"data");
  // leg4->SetTextSize(0.04);
  // leg4->Draw();

  // c4->SetLeftMargin(0.12);
  // c4->SetTopMargin(0.04);
  // c4->SaveAs("BDT_flipJpsiVsData_"+(TString)(highMass?"highMassRegion_":"")+(TString)(ispp?"pp":"PbPb")+".pdf");


}