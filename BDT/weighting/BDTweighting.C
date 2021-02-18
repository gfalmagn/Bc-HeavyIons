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
#include "TSystem.h"
#include "TROOT.h"
#include "TLatex.h"
#include "../../helpers/Cuts_BDT.h"
#include "../../helpers/Cuts.h"

void BDTweight(bool ispp=true, bool step2=true, bool step3=false){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  TString metafitSyst = "";
  TString systExt = "";
  bool useFutureFit = true; //if true, needs to have run previously a fit without BDT weights or with weights from the 1st-step fit

  //definitions of sig and bkg processes
  vector<vector<TString>> procName{{"data_obs"},{"BcSig"},{"FakeJpsi","FakeJpsi_JpsiSBUp","FakeJpsi_JpsiSBDown"},
							    ispp?(vector<TString>{"JpsiMC","JpsiMC_wPromptMCUp","JpsiMC_wPromptMCDown"}):(vector<TString>{"NPJpsi","NPJpsi_UncorrNPJUp","NPJpsi_UncorrNPJDown"}),
							      ispp?(vector<TString>{"flipJpsi","flipJpsi_flipJSameSideUp","flipJpsi_flipJSameSideDown"}):(vector<TString>{"flipJpsi","flipJpsi_FlipJorMCUp","flipJpsi_FlipJorMCDown"})};
  int nproc = procName.size();

  //file with distros
  auto histFile = TFile::Open("../../templateFit/InputForCombine_"+(TString)(ispp?"pp":"PbPb")+(TString)(step2?"_2ndStep":"")+metafitSyst+".root");

  //file of fit output
  TString normFileName = "../../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_2bins"+(TString)((step2 && useFutureFit)?"_2ndStep":"")+metafitSyst+systExt+".root"; //!!!! HERE change to 2nd step
  auto normFile = TFile::Open(normFileName);

  //Grab NP for shape morphings
  RooArgList fittedPars = ((RooFitResult*)normFile->Get("fit_s"))->floatParsFinal();
  float np_FakeJpsiSB = ((RooRealVar*)fittedPars.find("JpsiSB"))->getValV();
  float np_JpsiMC = ((RooRealVar*)fittedPars.find(ispp?"wPromptMC":"UncorrNPJ"))->getValV();
  float np_flipJpsi = ((RooRealVar*)fittedPars.find(ispp?"flipJSameSide":"FlipJorMC"))->getValV();
  vector<float> shapeMorph = {0,0,np_FakeJpsiSB,np_JpsiMC,np_flipJpsi};
  cout<<"np_FakeJpsiSB np_JpsiMC np_flipJpsi = "<<np_FakeJpsiSB <<" "<<np_JpsiMC <<" "<<np_flipJpsi <<endl;

  //Grab BDT distributions  
  vector<vector<TH1F*> > h_bdtSum(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<vector<TH1F*> > > h_bdt(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));
  vector<vector<vector<TH1F*> > > h_bdtUp(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));
  vector<vector<vector<TH1F*> > > h_bdtDown(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));
  vector<TH1F*> h_BDTSum(_NanaBins+1);
  vector<TH1F*> h_BDTData(_NanaBins+1);
  vector<TH1F*> h_BDTratio(_NanaBins+1);

  for(int i=0; i<nproc; i++){
    for(int b=1;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){
	h_bdt[i][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procName[i][0]+"/BDTv");
	h_bdtUp[i][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procName[i][(i>1)?1:0]+"/BDTv");
	h_bdtDown[i][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procName[i][(i>1)?2:0]+"/BDTv");

	//deal with shape morphing, assuming linear interpolation
	if(shapeMorph[i]>0){
	  h_bdt[i][b][k]->Scale(1-shapeMorph[i]); //bdtnew = (1-np)*bdtold + np*bdtup;
	  h_bdt[i][b][k]->Add(h_bdtUp[i][b][k] , shapeMorph[i]);
	}
	else if(shapeMorph[i]<0){
	  h_bdt[i][b][k]->Scale(1+shapeMorph[i]); //bdtnew = (1+np)*bdtold - np*bdtdown;
	  h_bdt[i][b][k]->Add(h_bdtDown[i][b][k] , -shapeMorph[i]);
	}

	h_bdt[i][b][k]->Rebin(ispp?3:4);
      }
    }
  }
  
  //Grab postfit yields
  RooArgSet *Yields = (RooArgSet*)normFile->Get("norm_fit_s");
  vector<vector<vector<float> > > yields(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1)));

  for(int i=0; i<nproc; i++){
    for(int b=1;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){
	yields[i][b][k] = Yields->getRealValue("BDT"+(TString)(to_string(k))+"Kin"+(TString)(to_string(b))+"/"+procName[i][0]);
	if(i==0) yields[i][b][k] = h_bdt[i][b][k]->Integral(); //for data_obs
	//update the integral of BDT distro to the postfit yield
	else h_bdt[i][b][k]->Scale( yields[i][b][k] / h_bdt[i][b][k]->Integral() );
	//	cout<<"ispp iproc b k yield = "<<ispp<<" "<<i<<" "<<b<<" "<<k<<" "<<yields[i][b][k]<<endl;

	//sum of postfit templates
	if(i==1) h_bdtSum[b][k] = (TH1F*)h_bdt[i][b][k]->Clone("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/sumTemplates/BDTv");
	else if(i>1) h_bdtSum[b][k]->Add(h_bdt[i][b][k]);

	//drawing stuff
	h_bdt[i][b][k]->SetLineWidth(3);
	h_bdt[i][b][k]->GetXaxis()->SetRangeUser(_BDTcuts(ispp,b,step2)[0]-0.05,_BDTcuts(ispp,b,step2)[_nChan(ispp)]+0.05);
	h_bdt[i][b][k]->GetXaxis()->SetTitleSize(0.06);
	h_bdt[i][b][k]->GetYaxis()->SetTitleSize(0.06);
	if(i==1){
	  h_bdtSum[b][k]->SetLineWidth(3);
	  h_bdtSum[b][k]->GetXaxis()->SetRangeUser(_BDTcuts(ispp,b,step2)[0]-0.05,_BDTcuts(ispp,b,step2)[_nChan(ispp)]+0.05);
	  h_bdtSum[b][k]->GetXaxis()->SetTitleSize(0.06);
	  h_bdtSum[b][k]->GetYaxis()->SetTitleSize(0.06);
	}
      }
    }
  }

  //Summing over BDT bins
  for(int b=1;b<=_NanaBins;b++){      
    for(int k=1;k<=_nChan(ispp);k++){
      if(k==1) {
	h_BDTData[b] = (TH1F*)h_bdt[0][b][k]->Clone("Kin"+(TString)to_string(b)+"/data/BDTv");
	h_BDTSum[b] = (TH1F*)h_bdtSum[b][k]->Clone("Kin"+(TString)to_string(b)+"/sumTemplates/BDTv");
      }
      else {
	h_BDTData[b]->Add(h_bdt[0][b][k]);
	h_BDTSum[b]->Add(h_bdtSum[b][k]);
      }
    }
  }

  gStyle->SetOptStat(0);
  for(int b=1;b<=_NanaBins;b++){
    
    //plot all separate BDT distros
    TCanvas* c1 = new TCanvas("allSourcesBin"+(TString)(ispp?"_pp":"_PbPb")+(TString)to_string(b),"allSourcesBin"+(TString)(ispp?"_pp":"_PbPb")+(TString)to_string(b),3000,1000);
    c1->Divide(3,2);

    c1->cd(1);
    h_BDTData[b]->SetTitle("data;BDT;postfit yield;");
    h_BDTData[b]->Draw("");
      
    c1->cd(2);
    h_BDTSum[b]->SetTitle("summed templates;BDT;postfit yield;");
    h_BDTSum[b]->Draw("");

    for(int k=2;k<=_nChan(ispp);k++){
      for(int i=1; i<nproc; i++){      
	//sum over BDT bins
	if(k>1) h_bdt[i][b][1]->Add(h_bdt[i][b][k]);
      }
    }

    c1->cd(3);
    h_bdt[1][b][1]->SetTitle("signal MC;BDT;postfit yield;");
    h_bdt[1][b][1]->Draw("");

    c1->cd(4);
    h_bdt[2][b][1]->SetTitle("Fake J/#psi;BDT;postfit yield;");
    h_bdt[2][b][1]->Draw("");

    c1->cd(5);
    h_bdt[3][b][1]->SetTitle("B decays;BDT;postfit yield;");
    h_bdt[3][b][1]->Draw("");

    c1->cd(6);
    h_bdt[4][b][1]->SetTitle("J/#psi-#mu combinatorics;BDT;postfit yield;");
    h_bdt[4][b][1]->Draw("");

    c1->SaveAs("figs/BDTdistr_allSources_"+(TString)(step2?(step3?"3rdStep_":(useFutureFit?"2ndStepUseFinalFit_":"2ndStep_")):"")+(TString)(ispp?"pp":"PbPb")+"_kinBin_"+(TString)to_string(b)+".pdf");

    //plot comparison data/summed templates
    TCanvas* c2 = new TCanvas("dataPostfitComp"+(TString)(ispp?"_pp":"_PbPb")+(TString)to_string(b),"dataPostfitComp"+(TString)(ispp?"_pp":"_PbPb")+(TString)to_string(b),1500,1500);
    c2->Divide(1,2);
    c2->cd(1)->SetPad(0.,0.5,1.,1.);//upper pad
    c2->cd(1)->SetLeftMargin(0.13);
    c2->cd(1)->SetRightMargin(0.04);
    c2->cd(1)->SetTopMargin(0.04);
    c2->cd(1)->SetBottomMargin(0.);

    c2->cd(2)->SetPad(0.,0.,1.,0.5);//lower pad
    c2->cd(2)->SetLeftMargin(0.13);
    c2->cd(2)->SetRightMargin(0.04);
    c2->cd(2)->SetTopMargin(0.);
    c2->cd(2)->SetBottomMargin(0.14);
    gPad->SetTickx(1);

    c2->cd(1);
    h_BDTData[b]->SetTitle(";BDT;postfit yield;");
    h_BDTData[b]->Draw("");
    h_BDTSum[b]->SetLineColor(kRed);
    h_BDTSum[b]->Draw("same");

    TLegend* leg = new TLegend((ispp && b==2)?0.15:0.6,0.7,(ispp && b==2)?0.45:0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->AddEntry(h_BDTData[b], "data","l");
    leg->AddEntry(h_BDTSum[b], "sum of postfit templates","l");
    leg->Draw("same");
    
    c2->cd(2);
    h_BDTratio[b] = (TH1F*)h_BDTData[b]->Clone("BDT_RatioDataToSummedTemplates");
    h_BDTratio[b]->Divide(h_BDTSum[b]);
    h_BDTratio[b]->GetYaxis()->SetTitle("data/(sum postfit templates)");
    h_BDTratio[b]->GetYaxis()->SetRangeUser(0,2);
    h_BDTratio[b]->SetLineColor(kOrange+4);
    h_BDTratio[b]->Draw("");

    TLine one = TLine();
    one.SetLineStyle(7);
    one.SetLineColor(kGray+2);
    one.DrawLine(_BDTcuts(ispp,b,step2)[0]-0.05, 1 ,_BDTcuts(ispp,b,step2)[_nChan(ispp)]+0.05, 1);

    c2->SaveAs("figs/BDTdistr_dataComparisonToSummedTemplates_"+(TString)(step2?(step3?"3rdStep_":(useFutureFit?"2ndStepUseFinalFit_":"2ndStep_")):"")+(TString)(ispp?"pp":"PbPb")+"_kinBin"+(TString)to_string(b)+".pdf");

    //regularize weights histo
    for(int B=1;B<=h_BDTratio[b]->GetNbinsX();B++){
      if(fabs(h_BDTratio[b]->GetBinError(B)/h_BDTratio[b]->GetBinContent(B))>0.4 || 
	 (fabs(h_BDTratio[b]->GetBinError(B)/h_BDTratio[b]->GetBinContent(B))>0.2 && (h_BDTratio[b]->GetBinContent(B)<0.3 || h_BDTratio[b]->GetBinContent(B)>2)))
	h_BDTratio[b]->SetBinContent(B,1);
    }

  }

  TFile* weightF = new TFile("BDTdistrWeights_"+(TString)(ispp?"pp":"PbPb")+".root","UPDATE");
  for(int b=1;b<=_NanaBins;b++)
    h_BDTratio[b]->Write("BDT_RatioDataToSummedTemplates_"+(TString)(step2?(step3?"3rdStep_":(useFutureFit?"2ndStepUseFinalFit_":"2ndStep_")):"")+(TString)(ispp?"pp":"PbPb")+"_kinBin"+(TString)to_string(b));
  weightF->Close();

}

void BDTweighting(bool step2=true, bool step3=false){

  cout<<endl<<"BDT comparison for pp"<<endl<<endl;
  BDTweight(true,step2,step3);
  cout<<endl<<"BDT comparison for PbPb"<<endl<<endl;
  BDTweight(false,step2,step3);
}

