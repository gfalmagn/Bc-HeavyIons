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
#include "../BDT/Definitions.h"

void plotMorphedShapes(bool ispp=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  int ntrees = 9;

  vector<float> BDTcut = _BDTcuts(ispp); //vector of BDT cut values put into array
  bool ignore1stBin = BDTcut.size()>4;
  //  if(ignore1stBin) BDTcut.erase(0);
  //int nchan = BDTcut.size() -1;
  float BDTcut_l[_nChan(ispp)+2]; 
  BDTcut_l[0] = -1;
  for(int k=0;k<=_nChan(ispp);k++){
    BDTcut_l[k+1] = BDTcut[k];
  }

  //********************************************************
  //Extract BcM HISTOGRAMS used as combine input
  //********************************************************
  vector<vector<TH1F*> > flipJ_FlipJsyst,flipJ_JMCsyst,JMC_FlipJsyst,JMC_JMCsyst,Jpsi_FlipJsyst,Jpsi_JMCsyst,FakeJ_syst;

  auto histFile = TFile::Open("InputForCombine_"+(TString)(ispp?"pp":"PbPb")+".root");
  TString JMCsName[] = {"","_wPromptMCUp","_wPromptMCDown"};
  TString flipJsName[] = {"","_flipJSameSideUp","_flipJSameSideDown"};
  TString JSBsName[] = {"","_JpsiSBUp","_JpsiSBDown"};

  for(int k=0;k<=_nChan(ispp);k++){
    flipJ_FlipJsyst.push_back(vector<TH1F*>(3));
    JMC_FlipJsyst.push_back(vector<TH1F*>(3));
    flipJ_JMCsyst.push_back(vector<TH1F*>(3));
    JMC_JMCsyst.push_back(vector<TH1F*>(3));
    Jpsi_FlipJsyst.push_back(vector<TH1F*>(3));
    Jpsi_JMCsyst.push_back(vector<TH1F*>(3));
    FakeJ_syst.push_back(vector<TH1F*>(3));

    int kk = (k==0)?1:k;

    for(int s=0;s<3;s++){//nominal,systUp,systDown
      flipJ_FlipJsyst[k][s] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/flipJpsi"+flipJsName[s]+"/BcM");
      JMC_FlipJsyst[k][s] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/JpsiMC"+flipJsName[s]+"/BcM");
      flipJ_JMCsyst[k][s] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/flipJpsi"+JMCsName[s]+"/BcM");
      JMC_JMCsyst[k][s] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/JpsiMC"+JMCsName[s]+"/BcM");

      Jpsi_FlipJsyst[k][s] = (TH1F*)flipJ_FlipJsyst[k][s]->Clone(); Jpsi_FlipJsyst[k][s]->Add(JMC_FlipJsyst[k][s]);
      Jpsi_JMCsyst[k][s] = (TH1F*)flipJ_JMCsyst[k][s]->Clone(); Jpsi_JMCsyst[k][s]->Add(JMC_JMCsyst[k][s],2.);
      FakeJ_syst[k][s] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/FakeJpsi"+JSBsName[s]+"/BcM");

      if(k>1){
	flipJ_FlipJsyst[0][s]->Add(flipJ_FlipJsyst[k][s]);
	JMC_FlipJsyst[0][s]->Add(JMC_FlipJsyst[k][s]);
	flipJ_JMCsyst[0][s]->Add(flipJ_JMCsyst[k][s]);
	JMC_JMCsyst[0][s]->Add(JMC_JMCsyst[k][s]);

	Jpsi_FlipJsyst[0][s]->Add(Jpsi_FlipJsyst[k][s]);
	Jpsi_JMCsyst[0][s]->Add(Jpsi_JMCsyst[k][s]);
	FakeJ_syst[0][s]->Add(FakeJ_syst[k][s]);
      }
    }
  }

  gStyle->SetOptStat(0);
  // TCanvas* c1 = new TCanvas("c1", "c1",1500, 1500);
  // TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  // leg->SetTextSize(0.05);
  // leg->SetBorderSize(0);
  // for(int k=0;k<=_nChan(ispp);k++){
  //   c1->Divide(2,2);
  //   for(int i=0;i<3;i++){
  //     flipJ_FlipJsyst[k][i]->SetLineWidth(2);
  //     flipJ_FlipJsyst[k][i]->SetTitle("flipJpsi_systFlipJpsi");
  //     flipJ_FlipJsyst[k][i]->Scale(1/flipJ_FlipJsyst[k][i]->Integral());
  //     flipJ_FlipJsyst[k][i]->SetLineColor((i==0)?kRed:((i==1)?(kGreen+1):kBlue));
  //     c1->cd(1);
  //     flipJ_FlipJsyst[k][i]->Draw((i==0)?"hist":"histsame");
  //     leg->AddEntry(flipJ_FlipJsyst[k][i],(i==0)?"nominal":((i==1)?"systUp":"systDown"));
  //     leg->DrawClone("same");
      
  //     JMC_FlipJsyst[k][i]->SetLineWidth(2);
  //     JMC_FlipJsyst[k][i]->SetTitle("J/#psi MC_systFlipJpsi");
  //     JMC_FlipJsyst[k][i]->Scale(1/JMC_FlipJsyst[k][i]->Integral());
  //     JMC_FlipJsyst[k][i]->SetLineColor((i==0)?kRed:((i==1)?(kGreen+1):kBlue));
  //     //JMC_FlipJsyst[k][i]->AddEntry((i==0)?"nominal":((i==1)?"systUp":"systDown"));
  //     c1->cd(2);
  //     JMC_FlipJsyst[k][i]->Draw((i==0)?"hist":"histsame");

  //     flipJ_JMCsyst[k][i]->SetLineWidth(2);
  //     flipJ_JMCsyst[k][i]->SetTitle("flipJpsi_systJpsiMC");
  //     flipJ_JMCsyst[k][i]->Scale(1/flipJ_JMCsyst[k][i]->Integral());
  //     flipJ_JMCsyst[k][i]->SetLineColor((i==0)?kRed:((i==1)?(kGreen+1):kBlue));
  //     //flipJ_JMCsyst[k][i]->AddEntry((i==0)?"nominal":((i==1)?"systUp":"systDown"));
  //     c1->cd(3);
  //     flipJ_JMCsyst[k][i]->Draw((i==0)?"hist":"histsame");

  //     JMC_JMCsyst[k][i]->SetLineWidth(2);
  //     JMC_JMCsyst[k][i]->SetTitle("J/#psi MC_systJpsiMC");
  //     JMC_JMCsyst[k][i]->Scale(1/JMC_JMCsyst[k][i]->Integral());
  //     JMC_JMCsyst[k][i]->SetLineColor((i==0)?kRed:((i==1)?(kGreen+1):kBlue));
  //     //JMC_JMCsyst[k][i]->AddEntry((i==0)?"nominal":((i==1)?"systUp":"systDown"));
  //     c1->cd(4);
  //     JMC_JMCsyst[k][i]->Draw((i==0)?"hist":"histsame");      
      
  //   }

  //   c1->SaveAs("figs_morphing/morphedShapes_BDT"+(TString)((k==0)?"all":to_string(k))+".pdf");
  //   c1->Clear();
  //   leg->Clear();
  // }


  TCanvas* c2 = new TCanvas("c2", "c2",3000, 1500);
  TLegend* leg2 = new TLegend(0.55,0.6,0.9,0.9);
  leg2->SetTextSize(0.065);
  leg2->SetBorderSize(0);
  gStyle->SetTitleFontSize(0.11);
  gStyle->SetTitleOffset(-0.05);
  for(int k=0;k<=_nChan(ispp);k++){
    c2->Divide(3,1);
    for(int i=0;i<3;i++){
      Jpsi_FlipJsyst[k][i]->SetLineWidth(2);
      Jpsi_FlipJsyst[k][i]->GetXaxis()->SetTitle("M(B_{c}) [GeV]");
      Jpsi_FlipJsyst[k][i]->SetTitle("TrueJpsi_systFlipJpsi");
      Jpsi_FlipJsyst[k][i]->Scale(1/Jpsi_FlipJsyst[k][i]->Integral());
      Jpsi_FlipJsyst[k][i]->SetLineColor((i==0)?kRed:((i==1)?(kGreen+1):kBlue));
      c2->cd(1);
      Jpsi_FlipJsyst[k][i]->Draw((i==0)?"hist":"histsame");
      leg2->AddEntry(Jpsi_FlipJsyst[k][i],(i==0)?"nominal":((i==1)?"systUp":"systDown"));
      leg2->DrawClone("same");
      
      Jpsi_JMCsyst[k][i]->SetLineWidth(2);
      Jpsi_JMCsyst[k][i]->SetTitle("TrueJpsi_systJpsiMC");
      Jpsi_JMCsyst[k][i]->GetXaxis()->SetTitle("M(B_{c}) [GeV]");
      Jpsi_JMCsyst[k][i]->Scale(1/Jpsi_JMCsyst[k][i]->Integral());
      Jpsi_JMCsyst[k][i]->SetLineColor((i==0)?kRed:((i==1)?(kGreen+1):kBlue));
      //Jpsi_JMCsyst[k][i]->AddEntry((i==0)?"nominal":((i==1)?"systUp":"systDown"));
      c2->cd(2);
      Jpsi_JMCsyst[k][i]->Draw((i==0)?"hist":"histsame");
      
      FakeJ_syst[k][i]->SetLineWidth(2);
      FakeJ_syst[k][i]->SetTitle("FakeJpsi_systSB");
      FakeJ_syst[k][i]->GetXaxis()->SetTitle("M(B_{c}) [GeV]");
      FakeJ_syst[k][i]->Scale(1/FakeJ_syst[k][i]->Integral());
      FakeJ_syst[k][i]->SetLineColor((i==0)?kRed:((i==1)?(kGreen+1):kBlue));
      //FakeJ_syst[k][i]->AddEntry((i==0)?"nominal":((i==1)?"systUp":"systDown"));
      c2->cd(3);
      FakeJ_syst[k][i]->Draw((i==0)?"hist":"histsame");
      
    }

    c2->SaveAs("figs_morphing/morphedShapes_mergeTrueJpsi_BDT"+(TString)((k==0)?"all":to_string(k))+".pdf");
    c2->Clear();
    leg2->Clear();
  }

}

