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
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"

void Uncorrelate(bool ispp=true, int kinBin=1, int centBin=0,bool secondStep=false){

  cout<<"\n   ************ Uncorrelate BDT from mass, ispp, kinbin, centbin = "<<ispp<<" "<<kinBin<<" "<<centBin<<endl;

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  int ntrees = 9;

  vector<vector<TH1F*> > h_BcM(ntrees, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_postfit(ntrees, vector<TH1F*>(_nChan(ispp)+1 ));
  vector<TH1F* > h_BcM_bkg(_nChan(ispp)+1);
  
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi"};
  TString prettyName[] = {"WrongSign","J/#psi sidebands","High mass control","signal region","MC signal",
  			  "MC NonPromptJpsi","MC PromptJpsi","dimuon+track (misID)","flipped J/#psi"};
  
  TString procNameDef[] = {"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig",(ispp?"JpsiMC":"NPJpsi"),"PromptJpsi","dimuTrk","flipJpsi"};
  vector<vector<TString>> procName{{"Wrongsign"},{"FakeJpsi"},{"TrueJpsi"},{"data_obs"},{"BcSig"},
				  {ispp?"JpsiMC":"NPJpsi","bToJpsi"}, {"PromptJpsi", "JpsiMC_uncorr"},
				  {"dimuTrk"}, {"flipJpsi","flipJpsiSameSide" }};
  vector<int> systIdx{0,0,0,0,0,
      0,//nonprompt
      0,//prompt
      0,//dimuTrk
      0,//flipJpsi
      };
  
  bool usedForFit[] = {false,true,false,true,true,true,//nonprompt
		       false,//!ispp,//prompt
		      false,//dimutrk
		      true//flipJpsi
  };
  int nPiled = 0; for(int j=0;j<ntrees;j++) nPiled+= usedForFit[j];
  nPiled -= 1;
  int idxPiled[] = {1,5,8,4};//{1,5,8,4};
  if(!ispp) idxPiled[2] = 6; //PromptJpsi instead of flipJpsi

  //********************************************************  
  //Extract BcM HISTOGRAMS used as combine input
  //********************************************************
  vector<vector<TH1F*> > h_BcM_prefit(ntrees, vector<TH1F*>(_nChan(ispp)+1));
  auto histFile = TFile::Open("../templateFit/InputForCombine_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+".root");
  for(int i=0; i<ntrees; i++){
    if(!usedForFit[i]) continue;
    for(int k=0;k<=_nChan(ispp);k++){
      int kk = (k==0)?1:k;
      TString s_bin = (centBin==0)?((kinBin==0)?"":("Kin"+(TString)to_string(kinBin))):("Cent"+(TString)to_string(centBin));
      h_BcM_prefit[i][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(kk)+s_bin+"/"+procName[i][systIdx[i]]+"/BcM");
      
      if(k>1) {
	AddTH1(h_BcM_prefit[i][0],h_BcM_prefit[i][k]);}
    }
  }

  //********************************************************
  //Extract (POSTFIT NORMALISATIONS & SHAPES) from combine output
  //********************************************************
  TString normFileName = "../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics"
    +(TString)(ispp?"_pp":"_PbPb")+(TString)((centBin==0)?((kinBin==0)?"_integrated":"_2bins"):"_centBins")+(TString)(secondStep?"_2ndStep":"")+".root";
  cout<<"Extract normalisations from file "<<normFileName<<endl;
  auto normFile = TFile::Open(normFileName,"READ");
  
  RooArgSet *Yields = (RooArgSet*)normFile->Get("norm_fit_s");
  for(int i=0; i<ntrees; i++){
    if(!usedForFit[i] || i==3) continue;
    cout<<"Recovering postfit BcM shape for process "<<procNameDef[i]<<endl;
    for(int k=0;k<=_nChan(ispp);k++){

      //Clone trimuon mass histos prefit into postfit. Processes not used in fitting will not change histos postfit
      TString s_bin = (centBin==0)?((kinBin==0)?"":("Kin"+(TString)to_string(kinBin))):("Cent"+(TString)to_string(centBin));
      h_BcM_postfit[i][k] = (TH1F*) normFile->Get("shapes_fit_s/BDT"+(TString)(to_string((k==0)?1:k))+s_bin+"/"+procNameDef[i]);
      h_BcM_postfit[i][k]->SetDirectory(0);

      //inclusive on BDT bins
      if(k>1) AddTH1(h_BcM_postfit[i][0],h_BcM_postfit[i][k]);

      //total background
      if(k>0 && i!=3 && i!=4){
	if(h_BcM_bkg[k]==NULL || h_BcM_bkg[k]->GetEntries()==0) h_BcM_bkg[k] = (TH1F*) h_BcM_postfit[i][k]->Clone();
	else h_BcM_bkg[k]->Add(h_BcM_postfit[i][k]);}
    }
  }

  //total background, inclusive on BDT bins
  h_BcM_bkg[0] = (TH1F*) h_BcM_bkg[1]->Clone();
  for(int k=2;k<=_nChan(ispp);k++){
    AddTH1(h_BcM_bkg[0],h_BcM_bkg[k]);
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
  float Bc_Pt[ntrees];
  float Bc_Y[ntrees];
  float BDT[ntrees];
  float weight[ntrees];
  int flipJpsi[ntrees];
  int hiBin[ntrees];

  //init histograms
  double mbins[_nbinM(ispp)[0]+1];
  for(int bx=0;bx<=_nbinM(ispp)[0];bx++){ 
    mbins[bx] = _Mbinning(ispp,0)[bx];
  }
  vector<float> BDTcuts = _BDTcuts(ispp, kinBin, centBin, secondStep);
  double bdtlo = BDTcuts[0]-(ispp?0.2:0.25) , bdthi = BDTcuts[BDTcuts.size()-1]+0.03;
  int nbdt = 60;
  TH2F* h_BDTM_bkg = new TH2F("Bc_BDTM_bkg", "BDT vs M(B_{c}) background", _nbinM(ispp)[0], mbins, nbdt,bdtlo,bdthi);
  TH2F* h_BDTM_sig = new TH2F("Bc_BDTM_sig", "BDT vs M(B_{c}) signal", _nbinM(ispp)[0], mbins, nbdt,bdtlo,bdthi);
  TH2F* h_BDTM_data = new TH2F("Bc_BDTM_data", "BDT vs M(B_{c}) data", _nbinM(ispp)[0], mbins, nbdt,bdtlo,bdthi);
  TH1F* h_rmsBDT_bkg = new TH1F("Bc_rmsBDT_bkg", "RMS(BDT) vs M(B_{c}) background;M(trimuon);RMS(BDT)", _nbinM(ispp)[0], mbins);
  TH1F* h_rmsBDT_data = new TH1F("Bc_rmsBDT_data", "RMS(BDT) vs M(B_{c}) data;M(trimuon);RMS(BDT)", _nbinM(ispp)[0], mbins);
  float nbkg = 0, ndata = 0, nsig = 0;
  float CRbinwRatio = ((_mBcMax-_mBcMin)/_nbinMSR(ispp)[0]) * (_nbinMCR(ispp)[0]/(_mMax-_mBcMax));

  //tree and event loop
  for(int iT=0; iT<(int)T.size(); iT++){
    if(!usedForFit[iT]) continue;
    std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    T[iT]->SetBranchAddress(secondStep?"BDT2":"BDT", &BDT[iT]);
    T[iT]->SetBranchAddress("Bc_M", &Bc_M[iT]);
    T[iT]->SetBranchAddress("Bc_Y", &Bc_Y[iT]);
    T[iT]->SetBranchAddress("Bc_Pt", &Bc_Pt[iT]);
    if(iT==7) T[iT]->SetBranchAddress("flipJpsi", &flipJpsi[iT]);
    T[iT]->SetBranchAddress((secondStep && iT==4)?"weight2":"weight", &weight[iT]);
    if(!ispp) T[iT]->SetBranchAddress("Centrality", &hiBin[iT]);

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){
      T[iT]->GetEntry(j);
      if(!inFidCuts(kinBin,Bc_Pt[iT],Bc_Y[iT])) continue;
      if(!ispp && ((float)hiBin[iT] < 2*_Centmin[centBin] || (float)hiBin[iT] >= 2*_Centmax[centBin])) continue;
      
      float w = weight[iT];
      if(Bc_M[iT]>_mBcMax) w *= CRbinwRatio;

      if(iT==3){
	h_BDTM_data->Fill(Bc_M[iT], BDT[iT], w);
	ndata += w;}
      else if(iT==4){
	h_BDTM_sig->Fill(Bc_M[iT], BDT[iT], w);
	nsig += w;}
      else{
	
	//find BDT bin to get postfit correction
	int k0 = 1;
	for(int k=2;k<=_nChan(ispp);k++){
	  if(BDT[iT]>_BDTcuts(ispp,kinBin,centBin,secondStep)[k-1] && BDT[iT]<_BDTcuts(ispp,kinBin,centBin,secondStep)[k]) //1st step here because refers to 1st step fits
	    k0 = k;
	}

	int mbin = h_BcM_prefit[iT][k0]->FindBin(Bc_M[iT]);
	//correct mass distribution with postfit histograms
	float yieldCorr = (h_BcM_prefit[iT][k0]->GetBinContent(mbin) == 0)?1.:( h_BcM_postfit[iT][k0]->GetBinContent(mbin) / h_BcM_prefit[iT][k0]->GetBinContent(mbin) );

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

  TString s_Bin = (centBin==0)?((kinBin==0)?"_integrated":("_KinBin"+(TString)to_string(kinBin))):("_CentBin"+(TString)to_string(centBin));

  //********************************************************
  //DRAW BDT vs M 2D
  //********************************************************
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kLightTemperature);
  gStyle->SetNumberContours(40);
  TCanvas* c1 = new TCanvas("c1","c1",3900,1300);
  c1->Divide(3,1);
  c1->cd(1);
  h_BDTM_bkg->GetXaxis()->SetTitle("M(trimuon)");
  h_BDTM_bkg->GetYaxis()->SetTitle("BDT");
  h_BDTM_bkg->GetYaxis()->SetTitleOffset(0.8);
  h_BDTM_bkg->Draw("COLZ0");
  c1->cd(2);
  h_BDTM_sig->GetXaxis()->SetTitle("M(trimuon)");
  h_BDTM_sig->GetYaxis()->SetTitle("BDT");
  h_BDTM_sig->GetYaxis()->SetTitleOffset(0.8);
  h_BDTM_sig->Draw("COLZ0");
  c1->cd(3);
  h_BDTM_data->GetXaxis()->SetTitle("M(trimuon)");
  h_BDTM_data->GetYaxis()->SetTitle("BDT");
  h_BDTM_data->GetYaxis()->SetTitleOffset(0.8);
  h_BDTM_data->Draw("COLZ0");
  
  c1->SaveAs("figs_UncorrBDT/BDTvsM_2D_"+(TString)(ispp?"pp":"PbPb")+s_Bin+(TString)(secondStep?"_2ndStep":"")+".pdf");

  //********************************************************
  //DRAW average BDT vs M
  //********************************************************
  TH1F *avBDTvsM_bkg = (TH1F*)h_BDTM_bkg->ProfileX("avBDTvsM_bkg_KinBin"+(TString)to_string(kinBin),1,-1,""); //default error is error on the mean, option g is for weighted average errors 
  TH1F *avBDTvsM_sig = (TH1F*)h_BDTM_sig->ProfileX("avBDTvsM_sig_KinBin"+(TString)to_string(kinBin),1,-1,""); 
  TH1F *avBDTvsM_data = (TH1F*)h_BDTM_data->ProfileX("avBDTvsM_data_KinBin"+(TString)to_string(kinBin),1,-1,"");
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

  c2->SaveAs("figs_UncorrBDT/BDTvsM_profile_"+(TString)(ispp?"pp":"PbPb")+s_Bin+(TString)(secondStep?"_2ndStep":"")+".pdf");

  //********************************************************
  //Get RMS of BDT for each mass bin
  //********************************************************
  for(int b=1;b<=h_rmsBDT_bkg->GetNbinsX();b++){
    TH1D* htmp_bkg = (TH1D*)h_BDTM_bkg->ProjectionY("h_BDTForThisM_bkg",b,b,"");
    h_rmsBDT_bkg->SetBinContent( b, htmp_bkg->GetRMS() );
    htmp_bkg->Delete();

    TH1D* htmp_data = (TH1D*)h_BDTM_data->ProjectionY("h_BDTForThisM_data",b,b,"");
    h_rmsBDT_data->SetBinContent( b, htmp_data->GetRMS() );
    htmp_data->Delete();
  }

  TCanvas* c9 = new TCanvas("c9","c9",2600,1300);
  c9->Divide(2,1);
  c9->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  h_rmsBDT_bkg->Draw();
  c9->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  h_rmsBDT_data->Draw();

  c9->SaveAs("figs_UncorrBDT/RMSBDTvsM_"+(TString)(ispp?"pp":"PbPb")+s_Bin+(TString)(secondStep?"_2ndStep":"")+".pdf");

  //********************************************************
  //Get corrected-BDT vs M histograms
  //********************************************************
  TH2F* h_corrBDT_bkg = new TH2F("Bc_corrBDT_bkg", "corrected BDT vs M(B_{c}) background", _nbinM(ispp)[0], mbins, nbdt,-2.8,ispp?4.2:5.3);
  TH2F* h_corrBDT_sig = new TH2F("Bc_corrBDT_sig", "corrected BDT vs M(B_{c}) signal", _nbinM(ispp)[0], mbins, nbdt,    -2.8,ispp?4.2:5.3);
  TH2F* h_corrBDT_data = new TH2F("Bc_corrBDT_data", "corrected BDT vs M(B_{c}) data", _nbinM(ispp)[0], mbins, nbdt,    -2.8,ispp?4.2:5.3);

  //tree and event loop
  for(int iT=0; iT<(int)T.size(); iT++){
    if(!usedForFit[iT]) continue;
    std::cout << "--- Processing again: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){
      T[iT]->GetEntry(j);
      if(!inFidCuts(kinBin,Bc_Pt[iT],Bc_Y[iT])) continue;
      if(!ispp && ((float)hiBin[iT] < 2*_Centmin[centBin] || (float)hiBin[iT] >= 2*_Centmax[centBin])) continue;

      float w = weight[iT];
      if(Bc_M[iT]>_mBcMax) w *= CRbinwRatio;

      //find BDT bin to get postfit correction
      int k0 = 1;
      for(int k=2;k<=_nChan(ispp);k++){
	if(BDT[iT]>_BDTcuts(ispp,kinBin,centBin,false)[k-1] && BDT[iT]<_BDTcuts(ispp,kinBin,centBin,false)[k]) //1st step here because refers to 1st step fits
	  k0 = k;
      }

      float bdtcorr = avBDTvsM_bkg->Interpolate(Bc_M[iT]); //correct with the BDT average of expected background
      float bdtrms = h_rmsBDT_bkg->Interpolate(Bc_M[iT]); //divide by the rms(BDT) at this mass bin
      if(iT==3)
	h_corrBDT_data->Fill(Bc_M[iT], (BDT[iT] - bdtcorr)/bdtrms, w);
      else if(iT==4)
	h_corrBDT_sig->Fill(Bc_M[iT], (BDT[iT] - bdtcorr)/bdtrms, w);
      else{
	int mbin = h_BcM_prefit[iT][k0]->FindBin(Bc_M[iT]);
	float yieldCorr = (h_BcM_prefit[iT][k0]->GetBinContent(mbin) == 0)?0:( h_BcM_postfit[iT][k0]->GetBinContent(mbin) / h_BcM_prefit[iT][k0]->GetBinContent(mbin) );

	h_corrBDT_bkg->Fill(Bc_M[iT], (BDT[iT] - bdtcorr)/bdtrms, w*yieldCorr);
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
  
  c3->SaveAs("figs_UncorrBDT/correctedBDTvsM_2D_"+(TString)(ispp?"pp":"PbPb")+s_Bin+(TString)(secondStep?"_2ndStep":"")+".pdf");


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

  c4->SaveAs("figs_UncorrBDT/correctedBDTvsM_profile_"+(TString)(ispp?"pp":"PbPb")+s_Bin+(TString)(secondStep?"_2ndStep":"")+".pdf");

  //********************************************************
  //DRAW corrected BDT
  //********************************************************
  TH1F *corrBDT_bkg = (TH1F*)h_corrBDT_bkg->ProjectionY("corrBDT_bkg");
  TH1F *corrBDT_sig = (TH1F*)h_corrBDT_sig->ProjectionY("corrBDT_sig");
  TH1F *corrBDT_data = (TH1F*)h_corrBDT_data->ProjectionY("corrBDT_data");
  TCanvas* c5 = new TCanvas("c5","c5",3900,1300);
  c5->Divide(3,1);
  c5->cd(1);
  corrBDT_bkg->SetTitle("corrected BDT background;corrected BDT;N");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  corrBDT_bkg->Draw();
  c5->cd(2);
  corrBDT_sig->SetTitle("corrected BDT signal;corrected BDT;N");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  corrBDT_sig->Draw();
  c5->cd(3);
  corrBDT_data->SetTitle("corrected BDT data;corrected BDT;N");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  corrBDT_data->Draw();

  c5->SaveAs("figs_UncorrBDT/correctedBDT_"+(TString)(ispp?"pp":"PbPb")+s_Bin+(TString)(secondStep?"_2ndStep":"")+".pdf");

  // cout<<"poportion of signal MC with corrBDT<-0.1 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(-0.1)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<-0.07 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(-0.07)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<-0.05 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(-0.05)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0. = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.05 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.05)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.1 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.1)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.15 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.15)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.2 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.2)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.22 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.22)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.23 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.23)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.25 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.25)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.3 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.3)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.35 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.35)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.4 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.4)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.45 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.45)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.5 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.5)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.55 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.55)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.65 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.65)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.7 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.7)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.75 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.75)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.8 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.8)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.85 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.85)) / nsig <<endl;
  // cout<<"poportion of signal MC with corrBDT<0.9 = "<<corrBDT_sig->Integral(0,corrBDT_sig->FindBin(0.9)) / nsig <<endl;

  TFile *outf = new TFile("BDTuncorrFromM_"+(TString)(ispp?"pp":"PbPb")+".root",(kinBin==0 && centBin==0 && !secondStep)?"RECREATE":"UPDATE"); 
  avBDTvsM_bkg->Write("avBDTvsM_bkg"+s_Bin+(TString)(secondStep?"_2ndStep":""));
  h_rmsBDT_bkg->Write("rmsBDTvsM_bkg"+s_Bin+(TString)(secondStep?"_2ndStep":""));
  outf->Close();
}

void UncorrelateBDTfromM(bool secondStep=false){
  Uncorrelate(true,0,0,secondStep);
  Uncorrelate(true,1,0,secondStep);
  Uncorrelate(true,2,0,secondStep);
  Uncorrelate(false,0,0,secondStep);
  Uncorrelate(false,1,0,secondStep);
  Uncorrelate(false,2,0,secondStep);
  //centrality
  Uncorrelate(false,0,1,secondStep);
  Uncorrelate(false,0,2,secondStep);
}
