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

  //flipJ: 0=nominal, 1=flipSameSide, 2=flipOppSide, 3=PromptJpsi
  //JpsiMC: what to add to flipJpsi (or subtract from data, depending on AddMCtoFlipJ) after SB-subtraction from data? 0=bToJpsi, 1=NonPromptJpsi(-bToJpsi), 2=fullJpsiMC(-bToJpsi), 3=NonPromptJpsi(-bToJpsi x3), 4=NonPromptJpsi(-bToJpsi x0.33)
void BDTweight(bool ispp=true, int flipJ=0, int JpsiMC=0, float JpsiMCSF=1., float flipJpsiSF=1., float sigSF=1., int bin=0, bool step2=false){

  bool AddMCtoFlipJ = true; //whether to include the Jpsi MC (uncorrelated part) to the flipJpsi sample whose BDT distro is weighted
  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  //need to start from untouched saved file
  //gSystem->Exec("cp BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+"_copystep1.root");  
  auto fullFile = TFile::Open("../BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  //  fullFile->Cp("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+"_copystep1.root");
  int ntrees = 9;

  //Initialization of variables
  float Bc_CorrM[ntrees];
  float Bc_ctauSignif[ntrees];
  float Bc_ctauSignif3D[ntrees], Bc_log_ctauSignif3D[ntrees];
  float Bc_alpha[ntrees];
  float Bc_alpha3D[ntrees];
  float Bc_VtxProb[ntrees], Bc_log1min_VtxProb[ntrees];
  float QQ_VtxProb[ntrees], QQ_log1min_VtxProb[ntrees];
  float QQ_dca[ntrees], QQ_log_dca[ntrees];
  float dR_sum[ntrees];
  float dR_jpsiOverMuW[ntrees];
  float dR_jpsiMuW[ntrees];
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
  int flipJpsi[ntrees];
  float weight[ntrees];
  UInt_t eventNb[ntrees];
  float BDT[ntrees];
  float BDTprob[ntrees];
  float BDTrarity[ntrees];

  //*******************************************
  //Get trees
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi"};
  TString prettyName[] = {"WRONGSIGN","J/Psi sidebands","High mass control","signal region","MC signal expectation",
			  "MC NonPromptJpsi","MC PromptJpsi","dimuon+track (misID)","flipped J/Psi"};
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }
  vector<int> usedT{1,3,4,5,6,7,8}; //the trees that will be read
  vector<int> usedH{1,2,3,4,5,6,7,8}; //7 will be for the full MC, 2 for the SB-subtracted data

  //*******************************************
  //Initialize histos
  vector<vector<TH1F*>> h_bdt;
  for(int iT=0; iT<(int)T.size(); iT++){
    if(!count(usedH.begin(),usedH.end(),iT)) {
      h_bdt.push_back(vector<TH1F*>(2,new TH1F()));
    }
    else{
      TString idx_s = (TString)(to_string(iT)); 
      h_bdt.push_back((vector<TH1F*>){
	                              new TH1F("bdtSR"+idx_s,"bdtSR"+idx_s,(ispp?20:13),_withTM?-0.48:-1,_withTM?0.3:1.),
				      new TH1F("bdtCR"+idx_s,"bdtCR"+idx_s,(ispp?20:13),_withTM?-0.48:-1,_withTM?0.3:1.)
					});
    }
    h_bdt[iT][0]->SetDirectory(0);
    h_bdt[iT][1]->SetDirectory(0);
  }

  //*******************************************
  //Fill the BDT histos
  for(int iT=0; iT<(int)T.size(); iT++){
    if(!count(usedT.begin(),usedT.end(),iT) || iT==7) continue;

    std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    // T[iT]->SetBranchAddress("Bc_alpha", &Bc_alpha[iT]);
    // T[iT]->SetBranchAddress("Bc_alpha3D", &Bc_alpha3D[iT]);
    // T[iT]->SetBranchAddress("Bc_VtxProb", &Bc_VtxProb[iT]);
    // T[iT]->SetBranchAddress("QQ_VtxProb", &QQ_VtxProb[iT]);
    // T[iT]->SetBranchAddress("QQ_dca", &QQ_dca[iT]);
    // T[iT]->SetBranchAddress("muW_isGlb", &muW_isGlb[iT]);
    // T[iT]->SetBranchAddress("muW_inLooseAcc", &muW_inLooseAcc[iT]);
    // T[iT]->SetBranchAddress("muW_inTightAcc", &muW_inTightAcc[iT]);
    // T[iT]->SetBranchAddress("mupl_isGlb", &mupl_isGlb[iT]);
    // T[iT]->SetBranchAddress("mupl_inLooseAcc", &mupl_inLooseAcc[iT]);
    // T[iT]->SetBranchAddress("mupl_inTightAcc", &mupl_inTightAcc[iT]);
    // T[iT]->SetBranchAddress("mumi_isGlb", &mumi_isGlb[iT]);
    // T[iT]->SetBranchAddress("mumi_inLooseAcc", &mumi_inLooseAcc[iT]);
    // T[iT]->SetBranchAddress("mumi_inTightAcc", &mumi_inTightAcc[iT]);

    //T[iT]->SetBranchAddress("Bc_CorrM", &Bc_CorrM[iT]);
    // T[iT]->SetBranchAddress("QQmuW_ptImbal", &QQmuW_ptImbal[iT]);
    // T[iT]->SetBranchAddress("Bc_ctauSignif", &Bc_ctauSignif[iT]);
    // T[iT]->SetBranchAddress("Bc_ctauSignif3D", &Bc_ctauSignif3D[iT]);
    // T[iT]->SetBranchAddress("dR_sum", &dR_sum[iT]);
    // T[iT]->SetBranchAddress("dR_jpsiOverMuW", &dR_jpsiOverMuW[iT]);
    // T[iT]->SetBranchAddress("dR_jpsiMuW", &dR_jpsiMuW[iT]);
    // T[iT]->SetBranchAddress("MuonDxySignif_sum", &MuonDxySignif_sum[iT]);
    // T[iT]->SetBranchAddress("MuonDzSignif_sum", &MuonDzSignif_sum[iT]);
    T[iT]->SetBranchAddress("Bc_M", &Bc_M[iT]);
    T[iT]->SetBranchAddress("Bc_Pt", &Bc_Pt[iT]);
    T[iT]->SetBranchAddress("Bc_Y", &Bc_Y[iT]);
    //T[iT]->SetBranchAddress("QQ_M", &QQ_M[iT]);
    T[iT]->SetBranchAddress("muW_isJpsiBro", &muW_isJpsiBro[iT]);
    //    T[iT]->SetBranchAddress("muW_trueId", &muW_trueId[iT]);
    T[iT]->SetBranchAddress("bkgType", &bkgType[iT]);
    T[iT]->SetBranchAddress("weight", &weight[iT]);
    T[iT]->SetBranchAddress("BDT", &BDT[iT]);
    T[iT]->SetBranchAddress("eventNb", &eventNb[iT]);
    if(iT==8) T[iT]->SetBranchAddress("flipJpsi", &flipJpsi[iT]);

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){//T[iT]->GetEntries()

      T[iT]->GetEntry(j);
      if(bin>0 && !( Bc_Pt[iT]>_BcPtmin[bin] && Bc_Pt[iT]<_BcPtmax[bin] && fabs(Bc_Y[iT])>_BcYmin[bin] && fabs(Bc_Y[iT])<_BcYmax[bin])) continue;
      if(bin==0){//integrated bin
	bool inFidCuts = false;
	for(int b=1;b<=_NanaBins;b++)
	  inFidCuts = inFidCuts || (Bc_Pt[iT]>_BcPtmin[b] && Bc_Pt[iT]<_BcPtmax[b] && fabs(Bc_Y[iT])>_BcYmin[b] && fabs(Bc_Y[iT])<_BcYmax[b]);
	if(!inFidCuts) continue;
      }

      int ifill = iT;
      float w = weight[iT];

      if(iT==5)
	if(!muW_isJpsiBro[iT]) ifill=7;
      if(iT==5 || iT==6)
	w *= JpsiMCSF;
      if(iT==6 && flipJ==3) {
	ifill = 8; //assimilate promptJpsi to nominal flipJpsi
	w *= JpsiMCSF;}

      if(iT==8){
	w *= flipJpsiSF;
	if(flipJ==1){
	  if(flipJpsi[iT]>4) {w *= 7/3.;} else {continue;}
	}
	else if(flipJ==2){
	  if(flipJpsi[iT]<=4) {w *= 7/4.;} else {continue;}
	}
      }

      if(Bc_M[iT]>6.3) h_bdt[ifill][1]->Fill(BDT[iT], w);
      else h_bdt[ifill][0]->Fill(BDT[iT], w);
      
    }
    //END event loop
  }
  //END loop on trees
  fullFile->Close();

  //*******************************************
  //Add prompt and nonprompt MC, subtract Jpsi SB (and appropriate JpsiMC)
  for(int i=0; i<2; i++){
    //5:bToJpsi, 6:prompt, 7:nonprompt without bToJpsi
    delete h_bdt[2][i]; h_bdt[2][i] = (TH1F*)h_bdt[3][i]->Clone(); //clone data hist
    h_bdt[2][i]->Add(h_bdt[1][i],-1); //subtract Jpsi SB from data
    h_bdt[2][i]->Add(h_bdt[4][i],-sigSF); //subtract Bc MC from data

    h_bdt[2][i]->Add(h_bdt[5][i], -((JpsiMC==3)?3.:((JpsiMC==4)?0.33:1.)) ); //subtract pure bToJpsi from data
    if(AddMCtoFlipJ){
      if(JpsiMC>=1) h_bdt[8][i]->Add(h_bdt[7][i]); //add rest of nonprompt JpsiMC to flipJpsi
      if(JpsiMC==2) h_bdt[8][i]->Add(h_bdt[6][i]); //add prompt JpsiMC to flipJpsi
    }
    else{      
      if(JpsiMC>=1) h_bdt[2][i]->Add(h_bdt[7][i], -1 ); //subtract rest of nonprompt JpsiMC from data
      if(JpsiMC==2) h_bdt[2][i]->Add(h_bdt[6][i], -1 ); //subtract prompt JpsiMC from data
    }

    //Now 7=full Jpsi MC
    h_bdt[7][i]->Add(h_bdt[6][i]); //add prompt to full JpsiMC
    h_bdt[7][i]->Add(h_bdt[5][i]); //add bToJpsi to full JpsiMC

  }

  //*******************************************
  //Draw BDT distributions
  for(int t=0; t<ntrees; t++){
    if(!count(usedH.begin(),usedH.end(),t)) continue;
    for(int i=0; i<2; i++){
      h_bdt[t][i]->SetLineWidth(2);
      h_bdt[t][i]->SetMarkerSize(1.9);
      h_bdt[t][i]->SetMarkerStyle(20);
      h_bdt[t][i]->SetMarkerColor((i==0)?kBlue:kRed);
      h_bdt[t][i]->SetLineColor((i==0)?kBlue:kRed);

      h_bdt[t][i]->Scale(1/h_bdt[t][i]->Integral());      
    }

    for(int i=0; i<2; i++)  h_bdt[t][i]->GetYaxis()->SetRangeUser(0, 1.1*h_bdt[t][1]->GetMaximum());
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("CRvsSR","CRvsSR",1500,1500);
  c1->Divide(2,2);

  TLegend* leg = new TLegend(0.59,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

  TString title[4] = {"data", "data ("+(TString)(AddMCtoFlipJ?"Bc":"")+"MC&SB-subtracted)" , "J/#psi MC", (TString)((flipJ==3)?"uncorrelated J/#psi MC":("flipped J/#psi"+(TString)(AddMCtoFlipJ?" + J/#psi MC":""))) };
  int drawIdx[4] = {3,2,7,8};

  for(int p=0;p<4;p++){
    c1->cd(p+1);
    int t = drawIdx[p];
    h_bdt[t][0]->SetTitle(title[p]);
    h_bdt[t][0]->GetXaxis()->SetTitle("BDT");
    h_bdt[t][0]->Draw("PE");
    h_bdt[t][1]->Draw("PEsame");
    leg->Clear();
    leg->AddEntry(h_bdt[t][0], "signal region");
    leg->AddEntry(h_bdt[t][1], "high mass CR");
    leg->DrawClone("same");
  }

  c1->SaveAs("BDTdistributions_SRandCR_flipJpsi"+(TString)(to_string(flipJ))+"_JpsiMC"+(TString)(to_string(JpsiMC))+(TString)(AddMCtoFlipJ?"_AddMCtoFlipJ":"")+(TString)((bin>0)?("_kinBin"+(TString)to_string(bin)):"")+(TString)(ispp?"_pp":"_PbPb")+".pdf");

  //*******************************************
  //Get ratios of flipJpsi to data for weights, in CR and SR
  vector<TH1F*> dataOverFlipJ;
  for(int i=0; i<2; i++){
    dataOverFlipJ.push_back((TH1F*)h_bdt[2][i]->Clone());
    dataOverFlipJ[i]->Divide(h_bdt[8][i]); //data divided by flipJpsi
    for (int B=1;B<dataOverFlipJ[i]->GetNbinsX();B++){
      if(h_bdt[8][i]->GetBinContent(B)<=0 || (h_bdt[2][i]->GetBinContent(B)<=0 && h_bdt[2][i]->GetBinCenter(B)>0.) ) dataOverFlipJ[i]->SetBinContent(B,1);
      else if(fabs(h_bdt[2][i]->GetBinContent(B) / h_bdt[8][i]->GetBinContent(B)) > 10 ) dataOverFlipJ[i]->SetBinContent(B,10); //limit weight to x10
      else if(fabs(h_bdt[2][i]->GetBinContent(B) / h_bdt[8][i]->GetBinContent(B)) <= 0 ) dataOverFlipJ[i]->SetBinContent(B, step2?0.:1.); //limit weight to x10
    }
  }

  TCanvas* c2 = new TCanvas("CRweights","CR weights",2300,1500);
  c2->Divide(2,2);
  c2->cd(3)->SetPad(0.,0.,0.5,0.35);//lower pad for ratio
  c2->cd(1)->SetPad(0.,0.35,0.5,1.);//upper pad for BDT distros
  c2->cd(4)->SetPad(0.5,0.,1.,0.35);//lower pad for ratio
  c2->cd(2)->SetPad(0.5,0.35,1.,1.);//upper pad for BDT distros

  TLegend* leg3 = new TLegend(0.5,0.7,0.9,0.9);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.04);
	
  for(int i=0; i<2; i++){
    c2->cd(2-i);
    h_bdt[8][i]->GetXaxis()->SetTitle("BDT");
    h_bdt[8][i]->SetTitle("BDT in "+(TString)((i==0)?"SR":"CR"));
    h_bdt[(i==0)?8:2][i]->SetLineColor((i==0)?kRed:kBlue);
    h_bdt[(i==0)?8:2][i]->SetMarkerColor((i==0)?kRed:kBlue);
    h_bdt[8][i]->Draw("PE");
    h_bdt[2][i]->Draw("PEsame");
    leg3->Clear();
    leg3->AddEntry(h_bdt[2][i], (TString)((i==0)?"SR ":"CR ")+title[1]);
    leg3->AddEntry(h_bdt[8][i], (TString)((i==0)?"SR ":"CR ")+title[3]);
    leg3->DrawClone("same");

    c2->cd(4-i);
    gPad->SetGridy();
    gPad->SetLogy();
    dataOverFlipJ[i]->SetTitle("data / "+(TString)((flipJ==3)?"(uncorrelated J/#psi MC)":"(flipped-J/#psi + MC J/#psi)"));
    dataOverFlipJ[i]->SetTitleSize(0.1);
    dataOverFlipJ[i]->SetLineColor(kBlack);
    dataOverFlipJ[i]->SetMarkerColor(kBlack);
    dataOverFlipJ[i]->GetYaxis()->SetRangeUser(0.3,10.1);
    dataOverFlipJ[i]->GetYaxis()->SetTitle("data / "+(TString)((flipJ==3)?"(uncorr. J/#psi MC)":"(flipped+MC)-J/#psi"));
    dataOverFlipJ[i]->GetXaxis()->SetTitle("BDT");
    dataOverFlipJ[i]->Draw("PE");
  }  

  c2->cd(4);
  dataOverFlipJ[1]->SetMarkerColor(kRed+3);
  dataOverFlipJ[1]->SetLineColor(kRed+3);
  dataOverFlipJ[1]->Draw("PEsame");
  TLegend* leg2 = new TLegend(0.3,0.7,0.6,0.9);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.08);
  leg2->AddEntry(dataOverFlipJ[0], "SR");
  leg2->AddEntry(dataOverFlipJ[1], "CR");
  leg2->DrawClone("same");

  c2->SaveAs("BDTdistributions_ratioDataFlipJ_SRandCR_flipJpsi"+(TString)(to_string(flipJ))+"_JpsiMC"+(TString)(to_string(JpsiMC))+(TString)(AddMCtoFlipJ?"_AddMCtoFlipJ":"")+(TString)((bin>0)?("_kinBin"+(TString)to_string(bin)):"")+(TString)(ispp?"_pp":"_PbPb")+".pdf");


  TFile* weightF = new TFile("flipJpsiWeights_fromBDTinCR_"+(TString)(ispp?"pp":"PbPb")+".root","UPDATE");
  
  dataOverFlipJ[0]->Write("flipJpsiWeights_fromBDTinSR_flipJpsi"+(TString)(to_string(flipJ))+"_JpsiMC"+(TString)(to_string(JpsiMC))+(TString)(AddMCtoFlipJ?"_AddMCtoFlipJ":"")+(TString)((bin>0)?("_kinBin"+(TString)to_string(bin)):"") );
  dataOverFlipJ[1]->Write("flipJpsiWeights_fromBDTinCR_flipJpsi"+(TString)(to_string(flipJ))+"_JpsiMC"+(TString)(to_string(JpsiMC))+(TString)(AddMCtoFlipJ?"_AddMCtoFlipJ":"")+(TString)((bin>0)?("_kinBin"+(TString)to_string(bin)):"") );
  weightF->Close();

}

void BDTweighting(bool ispp=true, bool step2=false){

  float scaleJMC = ispp?1.8:1.5;
  float scaleFlipJ = 1.;
  vector<float> scaleSig(_NanaBins+1, ispp?1.:0.8);

  if(step2){
    bool BDTuncorrFromM = false;
    int skipBDTbins=0;
    bool ignoreBin2 = (skipBDTbins==2);
    bool ignoreBin1 = (skipBDTbins==1 || skipBDTbins==2);

    //  TString normFileName = "./CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(bToJpsiOnly?"bToJpsi":"NonPromptJpsi")+(TString)(useFlipJpsi?(flipJpsiSameSide?"_flipJpsiSameSide":"_flipJpsi"):(ispp?"":"_PromptJpsi"))+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ignoreBin2?"_ignoreBin1-2":(ignoreBin1?"_ignoreBin1":""))+(TString)(ispp?"_pp":"_PbPb")+"_2bins"+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
    TString normFileName = "./CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_NonPromptJpsi"+(TString)(ispp?"_flipJpsi":"_PromptJpsi")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ignoreBin2?"_ignoreBin1-2":(ignoreBin1?"_ignoreBin1":""))+(TString)(ispp?"_pp":"_PbPb")+"_2bins.root";
    auto normFile = TFile::Open(normFileName);
    RooArgList fittedPars = ((RooFitResult*)normFile->Get("fit_s"))->floatParsFinal();
    scaleFlipJ = ((RooRealVar*)fittedPars.find(ispp?"flipJSameSide":"PromptOrFlipJ"))->getValV();
    scaleJMC = ((RooRealVar*)fittedPars.find(ispp?"wPromptMC":"bJpsiFrac"))->getValV();
    scaleSig[1] = ((RooRealVar*)fittedPars.find("r1"))->getValV();
    scaleSig[2] = ((RooRealVar*)fittedPars.find("r2"))->getValV();
    scaleSig[0] = (scaleSig[1]+scaleSig[2])/2;
  }
  
  // for(int b=1;b<=_NanaBins;b++){
  for(int b=0;b<1;b++){
    //in PbPb: replace flipJpsi by PromptMC
    BDTweight(ispp, (ispp?0:3) , (ispp?0:3) , scaleJMC, scaleFlipJ, scaleSig[b], b, step2); 
    BDTweight(ispp, (ispp?0:3) , (ispp?2:4) , scaleJMC, scaleFlipJ, scaleSig[b], b, step2);
    BDTweight(ispp, (ispp?0:3) ,1, scaleJMC, scaleFlipJ, scaleSig[b], b, step2); //NonPromptMC - bToJpsi

    BDTweight(ispp, 1, 1, scaleJMC, scaleFlipJ, scaleSig[b], b, step2); //flipJpsi method variation: nonpromptMC-bToJpsi is nominal in pp and PbPb
    BDTweight(ispp, 2, 1, scaleJMC, scaleFlipJ, scaleSig[b], b, step2);
  }

}
