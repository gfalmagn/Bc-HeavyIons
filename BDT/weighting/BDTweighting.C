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

void BDTweight(bool ispp=true, bool step2=true, bool step3=false, bool inCentBins=false){

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
  TString normFileName = "../../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+(TString)(inCentBins?"_centBins":"_2bins")+(TString)((step2 && useFutureFit)?"_2ndStep":"")+metafitSyst+systExt+".root";
  auto normFile = TFile::Open(normFileName);
  TString normFileNameInteg = "../../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_integrated"+(TString)((step2 && useFutureFit)?"_2ndStep":"")+metafitSyst+systExt+".root";
  auto normFileInteg = TFile::Open(normFileNameInteg);

  //Grab NP for shape morphings
  vector<float> np_FakeJpsiSB, np_JpsiMC, np_flipJpsi;
  vector<vector<float> > shapeMorph;
  for(int i=0;i<2;i++){
    RooArgList fittedPars = ((RooFitResult*)((i==0)?normFileInteg:normFile)->Get("fit_s"))->floatParsFinal();
    np_FakeJpsiSB.push_back( ((RooRealVar*)fittedPars.find("JpsiSB"))->getValV() );
    np_JpsiMC.push_back( ispp?(((RooRealVar*)fittedPars.find("wPromptMC"))->getValV()):0. ); //ispp?"wPromptMC":"UncorrNPJ"
    np_flipJpsi.push_back( ((RooRealVar*)fittedPars.find(ispp?"flipJSameSide":"FlipJorMC"))->getValV() );
    shapeMorph.push_back( vector<float>({0,0,np_FakeJpsiSB[i],np_JpsiMC[i],np_flipJpsi[i]}) );
    cout<<(TString)((i==0)?"integrated fit: ":"2pT bins fit: ")+" np_FakeJpsiSB np_JpsiMC np_flipJpsi = "<<np_FakeJpsiSB[i] <<" "<<np_JpsiMC[i] <<" "<<np_flipJpsi[i] <<endl;
  }

  //Grab BDT distributions  
  vector<vector<TH1F*> > h_bdtSum(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<vector<TH1F*> > > h_bdt(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));
  vector<vector<vector<TH1F*> > > h_bdtUp(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));
  vector<vector<vector<TH1F*> > > h_bdtDown(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));
  vector<TH1F*> h_BDTSum(_NanaBins+1);
  vector<TH1F*> h_BDTData(_NanaBins+1);
  vector<TH1F*> h_BDTratio(_NanaBins+1);

  for(int i=0; i<nproc; i++){
    for(int b=0;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){
	h_bdt[i][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":((TString)(inCentBins?"Cent":"Kin")+(TString)to_string(b)))+"/"+procName[i][0]+"/BDTv");
	h_bdtUp[i][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":((TString)(inCentBins?"Cent":"Kin")+(TString)to_string(b)))+"/"+procName[i][(i>1)?1:0]+"/BDTv");
	h_bdtDown[i][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":((TString)(inCentBins?"Cent":"Kin")+(TString)to_string(b)))+"/"+procName[i][(i>1)?2:0]+"/BDTv");

	//deal with shape morphing, assuming linear interpolation
	if(shapeMorph[(int)(b>0)][i]>0){
	  h_bdt[i][b][k]->Scale(1-shapeMorph[(int)(b>0)][i]); //bdtnew = (1-np)*bdtold + np*bdtup;
	  h_bdt[i][b][k]->Add(h_bdtUp[i][b][k] , shapeMorph[(int)(b>0)][i]);
	}
	else if(shapeMorph[(int)(b>0)][i]<0){
	  h_bdt[i][b][k]->Scale(1+shapeMorph[(int)(b>0)][i]); //bdtnew = (1+np)*bdtold - np*bdtdown;
	  h_bdt[i][b][k]->Add(h_bdtDown[i][b][k] , -shapeMorph[(int)(b>0)][i]);
	}

	h_bdt[i][b][k]->Rebin(ispp?3:4);
      }
    }
  }
  
  //Grab postfit yields
  RooArgSet *Yields = (RooArgSet*)normFile->Get("norm_fit_s");
  RooArgSet *Yields_integ = (RooArgSet*)normFileInteg->Get("norm_fit_s");
  vector<vector<vector<float> > > yields(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1)));

  for(int i=0; i<nproc; i++){
    for(int b=0;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){
	yields[i][b][k] = ((b==0)?Yields_integ:Yields)->getRealValue("BDT"+(TString)to_string(k)+(TString)((b==0)?"":((TString)(inCentBins?"Cent":"Kin")+(TString)to_string(b)))+"/"+procName[i][0]);
	if(i==0) yields[i][b][k] = h_bdt[i][b][k]->Integral(); //for data_obs
	//update the integral of BDT distro to the postfit yield
	else h_bdt[i][b][k]->Scale( yields[i][b][k] / h_bdt[i][b][k]->Integral() );
	//	cout<<"ispp iproc b k yield = "<<ispp<<" "<<i<<" "<<b<<" "<<k<<" "<<yields[i][b][k]<<endl;

	//sum of postfit templates
	if(i==1) h_bdtSum[b][k] = (TH1F*)h_bdt[i][b][k]->Clone("BDT"+(TString)to_string(k)+(TString)((b==0)?"":((TString)(inCentBins?"Cent":"Kin")+(TString)to_string(b)))+"/sumTemplates/BDTv");
	else if(i>1) h_bdtSum[b][k]->Add(h_bdt[i][b][k]);

	//drawing stuff
	h_bdt[i][b][k]->SetLineWidth(3);
	h_bdt[i][b][k]->GetXaxis()->SetRangeUser(_BDTcuts(ispp,inCentBins?0:b , inCentBins?b:0,step2)[0]-0.05,_BDTcuts(ispp,inCentBins?0:b , inCentBins?b:0,step2)[_nChan(ispp)]+0.05);
	h_bdt[i][b][k]->GetXaxis()->SetTitleSize(0.075);
	h_bdt[i][b][k]->GetYaxis()->SetTitleSize(0.075);
	h_bdt[i][b][k]->GetXaxis()->SetTitleOffset(0.85);
	h_bdt[i][b][k]->GetYaxis()->SetTitleOffset(0.7);
	h_bdt[i][b][k]->GetXaxis()->SetLabelSize(0.061);
	h_bdt[i][b][k]->GetYaxis()->SetLabelSize(0.061);
	if(i==1){
	  h_bdtSum[b][k]->SetLineWidth(3);
	  h_bdtSum[b][k]->GetXaxis()->SetRangeUser(_BDTcuts(ispp,inCentBins?0:b , inCentBins?b:0,step2)[0]-0.05,_BDTcuts(ispp,inCentBins?0:b , inCentBins?b:0,step2)[_nChan(ispp)]+0.05);
 	  h_bdtSum[b][k]->GetXaxis()->SetTitleSize(0.075);
	  h_bdtSum[b][k]->GetYaxis()->SetTitleSize(0.075);
	  h_bdtSum[b][k]->GetXaxis()->SetTitleOffset(0.85);
	  h_bdtSum[b][k]->GetYaxis()->SetTitleOffset(0.7);
 	  h_bdtSum[b][k]->GetXaxis()->SetLabelSize(0.061);
	  h_bdtSum[b][k]->GetYaxis()->SetLabelSize(0.061);
	}
      }
    }
  }

  //Summing over BDT bins
  for(int b=0;b<=_NanaBins;b++){      
    for(int k=1;k<=_nChan(ispp);k++){
      if(k==1) {
	h_BDTData[b] = (TH1F*)h_bdt[0][b][k]->Clone((TString)((b==0)?"":((TString)(inCentBins?"Cent":"Kin")+(TString)to_string(b)))+"/data/BDTv");
	h_BDTSum[b] = (TH1F*)h_bdtSum[b][k]->Clone((TString)((b==0)?"":((TString)(inCentBins?"Cent":"Kin")+(TString)to_string(b)))+"/sumTemplates/BDTv");
      }
      else {
	h_BDTData[b]->Add(h_bdt[0][b][k]);
	h_BDTSum[b]->Add(h_bdtSum[b][k]);
      }
    }
  }

  gStyle->SetOptStat(0);
  for(int b=0;b<=_NanaBins;b++){
    
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

    c1->SaveAs("figs/BDTdistr_allSources_"+(TString)(step2?(step3?"3rdStep_":(useFutureFit?"2ndStepUseFinalFit_":"2ndStep_")):"")+(TString)(ispp?"pp":"PbPb")+(TString)(inCentBins?"_centBin":"_kinBin")+(TString)to_string(b)+".pdf");

    //plot comparison data/summed templates
    TCanvas* c2 = new TCanvas("dataPostfitComp"+(TString)(ispp?"_pp":"_PbPb")+(TString)to_string(b),"dataPostfitComp"+(TString)(ispp?"_pp":"_PbPb")+(TString)to_string(b),1500,1500);
    c2->Divide(1,2);
    c2->cd(1)->SetPad(0.,0.5,1.,1.);//upper pad
    c2->cd(1)->SetLeftMargin(0.13);
    c2->cd(1)->SetRightMargin(0.04);
    c2->cd(1)->SetTopMargin(0.08);
    c2->cd(1)->SetBottomMargin(0.);

    c2->cd(2)->SetPad(0.,0.,1.,0.5);//lower pad
    c2->cd(2)->SetLeftMargin(0.13);
    c2->cd(2)->SetRightMargin(0.04);
    c2->cd(2)->SetTopMargin(0.);
    c2->cd(2)->SetBottomMargin(0.14);
    gPad->SetTickx(1);

    c2->cd(1);
    h_BDTData[b]->SetTitle(";BDT;postfit yield;");
    h_BDTData[b]->GetYaxis()->SetRangeUser(0.01, 1.2*h_BDTData[b]->GetMaximum());
    h_BDTData[b]->Draw("");
    h_BDTSum[b]->SetLineColor(kRed);
    h_BDTSum[b]->Draw("same");

    TLatex CMStag;
    CMStag.SetNDC();
    CMStag.SetTextFont(42);
    CMStag.SetTextSize(0.065);
    CMStag.DrawLatex(0.75,0.62,"#font[61]{CMS}");
    //CMStag.DrawLatex(0.75,0.54,"#font[52]{Preliminary}");
    TLatex lumitag;
    lumitag.SetNDC();
    lumitag.SetTextFont(42);
    lumitag.SetTextSize(0.063);
    lumitag.DrawLatex(ispp?0.66:0.62,0.935, Form(ispp?"pp (%.0f pb^{-1}), 5.02 TeV":"PbPb (%.2f nb^{-1}), 5.02 TeV",ispp?L_pp:(L_PbPb*1e3)));

    TLegend* leg = new TLegend(0.64,0.72,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->AddEntry(h_BDTData[b], "data","l");
    leg->AddEntry(h_BDTSum[b], "#Sigma postfit templates","l");
    leg->Draw("same");
    
    c2->cd(2);
    h_BDTratio[b] = (TH1F*)h_BDTData[b]->Clone("BDT_RatioDataToSummedTemplates");
    h_BDTratio[b]->Divide(h_BDTSum[b]);
    h_BDTratio[b]->GetYaxis()->SetTitle("data / (#Sigma postfit templates)");
    h_BDTratio[b]->GetYaxis()->SetRangeUser(0,1.99);
    h_BDTratio[b]->SetLineColor(kOrange+4);
    h_BDTratio[b]->Draw("");

    TLine one = TLine();
    one.SetLineStyle(7);
    one.SetLineColor(kGray+2);
    one.DrawLine(_BDTcuts(ispp,inCentBins?0:b , inCentBins?b:0,step2)[0]-((b==1)?0.15:0.05), 1 ,_BDTcuts(ispp,inCentBins?0:b , inCentBins?b:0,step2)[_nChan(ispp)]+(ispp?0.1:0.06), 1);

    c2->SaveAs("figs/BDTdistr_dataComparisonToSummedTemplates_"+(TString)(step2?(step3?"3rdStep_":(useFutureFit?"2ndStepUseFinalFit_":"2ndStep_")):"")+(TString)(ispp?"pp":"PbPb")+(TString)(inCentBins?"_centBin":"_kinBin")+(TString)to_string(b)+".pdf");

    //regularize weights histo
    for(int B=1;B<=h_BDTratio[b]->GetNbinsX();B++){
      if(fabs(h_BDTratio[b]->GetBinError(B)/h_BDTratio[b]->GetBinContent(B))>0.4 || 
	 (fabs(h_BDTratio[b]->GetBinError(B)/h_BDTratio[b]->GetBinContent(B))>0.2 && (h_BDTratio[b]->GetBinContent(B)<0.3 || h_BDTratio[b]->GetBinContent(B)>2)))
	h_BDTratio[b]->SetBinContent(B,1);
    }

  }

  TFile* weightF = new TFile("BDTdistrWeights_"+(TString)(ispp?"pp":"PbPb")+".root","UPDATE");
  for(int b=0;b<=_NanaBins;b++)
    h_BDTratio[b]->Write("BDT_RatioDataToSummedTemplates_"+(TString)(step2?(step3?"3rdStep_":(useFutureFit?"2ndStepUseFinalFit_":"2ndStep_")):"")+(TString)(ispp?"pp":"PbPb")+(TString)(inCentBins?"_centBin":"_kinBin")+(TString)to_string(b));
  weightF->Close();

}

void BDTweighting(bool step2=true, bool step3=false){

  cout<<endl<<"BDT comparison for pp"<<endl<<endl;
  BDTweight(true,step2,step3,false);
  cout<<endl<<"BDT comparison for PbPb"<<endl<<endl;
  BDTweight(false,step2,step3,false);
  cout<<endl<<"BDT comparison for PbPb centrality bins"<<endl<<endl;
  BDTweight(false,step2,step3,true);
}

