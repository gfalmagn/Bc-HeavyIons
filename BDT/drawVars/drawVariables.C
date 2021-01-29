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
#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/TMVAGui.h"
#include "../../helpers/Cuts_BDT.h"
#include "../../helpers/Cuts.h"


//Draw all samples for a given variable
void DrawVar(vector<TH1F*> h, int ntrees, TString varName, bool ispp, bool lowerCut, bool leftLeg, bool noROC=false){
  cout<<"Draw variable "<<varName<<endl;
  gStyle->SetOptStat(0);

  vector<int> tToDraw{4,1,5,(ispp?8:6)};
  vector<TString> legends{"signal MC","J/#psi sidebands","NonPrompt J/#psi MC",(TString)(ispp?"flipped J/#psi":"Prompt J/#psi MC")};
  if(!ispp){
    tToDraw.push_back(8);
    legends.push_back("flipped J/#psi");
  }
  Color_t cols[] = {kCyan, kMagenta+1, kGreen, kRed-2, kBlue, kOrange-7, kOrange, kGreen+4, kGreen+1};

  TCanvas* c = new TCanvas("c"+varName,"c",noROC?1500:3000,1500);
  if(!noROC){
    c->Divide(2,1);
    c->cd(1);
  }
  TLegend* leg = new TLegend(leftLeg?0.14:0.59,ispp?0.67:0.61,leftLeg?0.42:0.98,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  float maxy = 0;
  int j=0;
  for(const auto& iT:tToDraw){
    //hist[ih][iT]->SetFillColor(cols[iT]);
    h[iT]->SetLineColor(((iT==4)?kBlack:(cols[iT]+2)));
    h[iT]->SetLineWidth((iT==4)?3:2);
    h[iT]->Scale(1/h[iT]->Integral(0,h[iT]->GetNbinsX()+1));
    maxy = max(maxy, (float)h[iT]->GetMaximum());

    leg->AddEntry(h[iT],legends[j],"l");
    j++;
  }

  h[tToDraw[0]]->GetYaxis()->SetRangeUser(0,1.05*maxy);
  for(int i=0;i<tToDraw.size();i++) h[tToDraw[i]]->Draw((i==0)?"hist":"histsame");
  h[4]->Draw("histsame");
  
  leg->SetTextSize(0.035);
  leg->Draw("same");

  if(noROC){
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.022);
    c->SetBottomMargin(0.12);

    c->SaveAs(varName+(TString)(ispp?"_pp":"_PbPb")+".pdf");    
  } else{
    c->cd(1)->SetLeftMargin(0.12);
    c->cd(1)->SetRightMargin(0.022);
    c->cd(1)->SetBottomMargin(0.12);
  }

  if(noROC) return; 

  //ROC CURVES
  c->cd(2)->SetBottomMargin(0.12);

  TLegend* leg2 = new TLegend(0.12,0.14,0.52,ispp?0.38:0.44);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.038);
  leg2->SetHeader("ROC "+(TString)h[tToDraw[0]]->GetXaxis()->GetTitle());
  vector<TGraph*> ROC;

  const int npt = h[tToDraw[0]]->GetNbinsX()+3; //signalMC (iT==4) should be first of tToDraw
  double effSig[npt];

  //signal eff
  //cout<<"calc efficiency signal"<<endl;
  int start = (lowerCut?(npt-1):0);
  for(int b=start;(b<npt && b>=0);(lowerCut?(b--):(b++))){
    if(b==start) effSig[b] = 0;
    else effSig[b] = effSig[lowerCut?(b+1):(b-1)] + h[tToDraw[0]]->GetBinContent(lowerCut?b:(b-1));  //assuming histo is scaled!! //add new bin to sum of previous bins
    //cout<<"b effSig = "<<b<<" "<<effSig[b]<<endl;
  }

  //loop on backgrounds (rejection)
  int k=-1;
  for(const auto& iT:tToDraw){
    k++;
    if(k==0) continue;
    double rejBkg[npt];

    //cout<<"calc rej bkg #"<<iT<<endl;
    int start = (lowerCut?0:(npt-1));
    for(int b=start;(b<npt && b>=0);(lowerCut?(b++):(b--))){
      if(b==start) rejBkg[b] = 0;
      else rejBkg[b] = rejBkg[b+(lowerCut?-1:1)] + h[iT]->GetBinContent(lowerCut?(b-1):b);  //assuming histo is scaled!!
      //cout<<"b rejBkg = "<<b<<" "<<rejBkg[b]<<endl;
    }

    // for(int b=0;b<npt;b++)
    //  cout<<"b, effsig, rejbkg = "<<b<<" "<<effSig[b]<<" "<<rejBkg[b]<<endl;
    ROC.push_back(new TGraph(npt,effSig,rejBkg) );
  }    

  //DRAW
  // ROC[0]->Draw("AP");
  for(int i=0;i<ROC.size();i++){
    ROC[i]->SetMarkerSize(2);
    ROC[i]->SetMarkerStyle(20+i);
    ROC[i]->SetMarkerColor(cols[tToDraw[i+1]]+2);
    ROC[i]->SetLineColor(cols[tToDraw[i+1]]+2);
    ROC[i]->GetXaxis()->SetLimits(0,1);
    ROC[i]->GetHistogram()->SetMaximum(1.);
    ROC[i]->GetHistogram()->SetMinimum(0.);
    ROC[i]->SetTitle("ROC curve ("+(TString)(lowerCut?"lower cut":"upper cut")+");signal efficiency;background rejection");

    leg2->AddEntry(ROC[i],legends[i+1],"lp");
    ROC[i]->Draw((i==0)?"APL":"PLsame");
  }
  leg2->Draw("same");

  TLine *line = new TLine(0,1,1,0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray);
  line->Draw("same");

  c->SaveAs(varName+(TString)(lowerCut?"_lowerCut":"_upperCut")+(TString)(ispp?"_pp":"_PbPb")+".pdf");

}





void drawVariables(bool ispp=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  auto fullFile = TFile::Open("../BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","UPDATE");

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
  float dR_jpsi[ntrees];
  float dR_muWmi[ntrees];
  float dR_muWpl[ntrees];
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
  float QQ_M[ntrees];
  float QQ2_M[ntrees];
  float QQ3_M[ntrees];
  int bkgType[ntrees];
  bool muW_isJpsiBro[ntrees];
  int muW_trueId[ntrees];
  float weight[ntrees];
  UInt_t eventNb[ntrees];
  float BDT[ntrees];

  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi","CorrbToJpsi_MC","prefitBkg"};
  TString prettyName[] = {"WRONGSIGN","J/Psi sidebands","High mass control","signal region","MC signal expectation",
			  "MC NonPromptJpsi","MC PromptJpsi","dimuon+track (misID)","flipped J/Psi","correlated NonPrompt MC","prefit full background"};
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }

  //histograms initialisation
  vector<TH1F*> h_VProb;
  vector<TH1F*> h_QQVProb;
  vector<TH1F*> h_alpha;
  vector<TH1F*> h_alpha3D;
  vector<TH1F*> h_TauSignif;  
  vector<TH1F*> h_Tau3DSignif;
  vector<TH1F*> h_mcorr;
  vector<TH1F*> h_dca;
  vector<TH1F*> h_QQMuWImbal;
  vector<TH1F*> h_SumDeltaR;
  vector<TH1F*> h_dRJpsiOverMuW;
  vector<TH1F*> h_SumDxyMu;
  vector<TH1F*> h_SumDzMu;
  vector<TH1F*> h_BDT;
  vector<TH1F*> h_AllDimuDeltaR;

  int nbins = 41;
  vector<TString> varName = {"VtxProb","QQVtxProb","alpha","alpha3D","TauSignif","TauSignif3D","mcorr","dca","QQMuWImbal","SumDeltaR","dRJpsiOverMuW","SumDxyMu","SumDzMu","BDT","AllDimuDeltaR"};
  vector<bool> lowercut = {1,1,0,0,1,1,0,0,0,0,1,1,1,1,0};
  vector<bool> leftLeg = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0};
  for(int iT=0; iT<(int)ntrees+1; iT++){ //order is important
    h_VProb.push_back(new TH1F("h_"+varName[0]+"_"+(TString)to_string(iT),"Trimuon vertex probability;VtxProb(#mu#mu#mu);norm to 1",nbins,_vtxProb_cut,1));
    h_QQVProb.push_back(new TH1F("h_"+varName[1]+"_"+(TString)to_string(iT),"J/#psi vertex probability;VtxProb(J/#psi);norm to 1",nbins,_QQvtxProb_cut,1));
    h_alpha.push_back(new TH1F("h_"+varName[2]+"_"+(TString)to_string(iT),"2D pointing angle;#alpha_{2D} [rad];norm to 1",nbins,0,_alpha_cut(ispp)));
    h_alpha3D.push_back(new TH1F("h_"+varName[3]+"_"+(TString)to_string(iT),"3D pointing angle;#alpha_{3D} [rad];norm to 1",nbins,0,_alpha3D_cut(ispp)));
    h_TauSignif.push_back(new TH1F("h_"+varName[4]+"_"+(TString)to_string(iT),"2D lifetime significance;#tau/#sigma_{#tau} ;norm to 1",nbins,_ctauSignif_cut,15));
    h_Tau3DSignif.push_back(new TH1F("h_"+varName[5]+"_"+(TString)to_string(iT),"3D lifetime significance;#tau_{3D}/#sigma_{#tau_{3D}} ;norm to 1",nbins,_ctauSignif3D_cut,15));
    h_mcorr.push_back(new TH1F("h_"+varName[6]+"_"+(TString)to_string(iT),"Corrected mass;m_{corr} [GeV];norm to 1",nbins,3.5,_BcCorrM_cut(ispp)));
    h_dca.push_back(new TH1F("h_"+varName[8]+"_"+(TString)to_string(iT),"Dimuon distance of closest approach;dca(J/#psi) [mm];norm to 1",nbins,0,min(0.22,(double)_QQdca_cut)));
    h_QQMuWImbal.push_back(new TH1F("h_"+varName[9]+"_"+(TString)to_string(iT),"p_{T} imbalance between J/#psi and #mu_{W};#left|#frac{p_{T}(J/#psi)-p_{T}(#mu_{W})}{p_{T}(J/#psi)+p_{T}(#mu_{W})}#right|;norm to 1",nbins,0,1.2));
    h_SumDeltaR.push_back(new TH1F("h_"+varName[10]+"_"+(TString)to_string(iT),"Sum of the three dimuon #DeltaR's;sum(#DeltaR(#mu_{i}#mu_{j}));norm to 1",nbins,0,7));
    h_dRJpsiOverMuW.push_back(new TH1F("h_"+varName[11]+"_"+(TString)to_string(iT),"Ratio of J/#psi #DeltaR to other #DeltaR's;#DeltaR(J/#psi)/(#DeltaR_{2} + #DeltaR_{3});norm to 1",nbins,0.1,1));
    h_SumDxyMu.push_back(new TH1F("h_"+varName[12]+"_"+(TString)to_string(iT),"Sum of significances of muon 2D displacement to PV;sum(d_{xy}(#mu)/#sigma);norm to 1",nbins,0,40));
    h_SumDzMu.push_back(new TH1F("h_"+varName[13]+"_"+(TString)to_string(iT),"Sum of significances of muon z displacement to PV;sum(d_{z}(#mu)/#sigma);norm to 1",nbins,0,40));
    h_BDT.push_back(new TH1F("h_"+varName[14]+"_"+(TString)to_string(iT),"BDT;BDT;norm to 1",nbins,_BDTcuts(ispp)[0]-0.2,_BDTcuts(ispp)[_BDTcuts(ispp).size()-1]+0.02));
    h_AllDimuDeltaR.push_back(new TH1F("h_"+varName[10]+"_"+(TString)to_string(iT),"#DeltaR of all dimuon pairs of trimuon;#DeltaR(#mu_{i}#mu_{j});norm to 1",nbins,0,3.25));
  }
  int nSpeHist = 1;

  vector<vector<TH1F*> > hist{h_VProb,h_QQVProb,h_alpha,h_alpha3D,h_TauSignif,h_Tau3DSignif,h_mcorr,h_dca,h_QQMuWImbal,h_SumDeltaR,h_dRJpsiOverMuW,h_SumDxyMu,h_SumDzMu,h_BDT,h_AllDimuDeltaR};

  for(int iT=0; iT<(int)ntrees; iT++){
    if(iT==7) continue;
    std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    T[iT]->SetBranchAddress("Bc_alpha", &Bc_alpha[iT]);
    T[iT]->SetBranchAddress("Bc_alpha3D", &Bc_alpha3D[iT]);
    T[iT]->SetBranchAddress("Bc_VtxProb", &Bc_VtxProb[iT]);
    T[iT]->SetBranchAddress("QQ_VtxProb", &QQ_VtxProb[iT]);
    T[iT]->SetBranchAddress("QQ_dca", &QQ_dca[iT]);
    T[iT]->SetBranchAddress("muW_isGlb", &muW_isGlb[iT]);
    T[iT]->SetBranchAddress("muW_inLooseAcc", &muW_inLooseAcc[iT]);
    T[iT]->SetBranchAddress("muW_inTightAcc", &muW_inTightAcc[iT]);
    T[iT]->SetBranchAddress("mupl_isGlb", &mupl_isGlb[iT]);
    T[iT]->SetBranchAddress("mupl_inLooseAcc", &mupl_inLooseAcc[iT]);
    T[iT]->SetBranchAddress("mupl_inTightAcc", &mupl_inTightAcc[iT]);
    T[iT]->SetBranchAddress("mumi_isGlb", &mumi_isGlb[iT]);
    T[iT]->SetBranchAddress("mumi_inLooseAcc", &mumi_inLooseAcc[iT]);
    T[iT]->SetBranchAddress("mumi_inTightAcc", &mumi_inTightAcc[iT]);
    T[iT]->SetBranchAddress("Bc_CorrM_shiftedM", &Bc_CorrM[iT]);
    T[iT]->SetBranchAddress("QQmuW_ptImbal", &QQmuW_ptImbal[iT]);
    T[iT]->SetBranchAddress("Bc_ctauSignif", &Bc_ctauSignif[iT]);
    T[iT]->SetBranchAddress("Bc_ctauSignif3D", &Bc_ctauSignif3D[iT]);
    T[iT]->SetBranchAddress("dR_jpsi", &dR_jpsi[iT]);
    T[iT]->SetBranchAddress("dR_muWmi", &dR_muWmi[iT]);
    T[iT]->SetBranchAddress("dR_muWpl", &dR_muWpl[iT]);
    T[iT]->SetBranchAddress("dR_sum_shiftedM", &dR_sum[iT]);
    T[iT]->SetBranchAddress("dR_jpsiOverMuW_shiftedM", &dR_jpsiOverMuW[iT]);
    T[iT]->SetBranchAddress("dR_jpsiMuW_shiftedM", &dR_jpsiMuW[iT]);
    T[iT]->SetBranchAddress("MuonDxySignif_sum", &MuonDxySignif_sum[iT]);
    T[iT]->SetBranchAddress("MuonDzSignif_sum", &MuonDzSignif_sum[iT]);
    T[iT]->SetBranchAddress("Bc_M", &Bc_M[iT]);
    T[iT]->SetBranchAddress("Bc_Pt", &Bc_Pt[iT]);
    T[iT]->SetBranchAddress("QQ_M", &QQ_M[iT]);
    T[iT]->SetBranchAddress("muW_isJpsiBro", &muW_isJpsiBro[iT]);
    T[iT]->SetBranchAddress("muW_trueId", &muW_trueId[iT]);
    T[iT]->SetBranchAddress("bkgType", &bkgType[iT]);
    T[iT]->SetBranchAddress("weight", &weight[iT]);
    T[iT]->SetBranchAddress("eventNb", &eventNb[iT]);
    T[iT]->SetBranchAddress("BDT", &BDT[iT]);

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){//T[iT]->GetEntries()

      T[iT]->GetEntry(j);
      vector<float> v{Bc_VtxProb[iT],QQ_VtxProb[iT],Bc_alpha[iT],Bc_alpha3D[iT],Bc_ctauSignif[iT],Bc_ctauSignif3D[iT],Bc_CorrM[iT],QQ_dca[iT],QQmuW_ptImbal[iT],dR_sum[iT],dR_jpsiOverMuW[iT],MuonDxySignif_sum[iT],MuonDzSignif_sum[iT],BDT[iT]};

      int iT2 = (iT==5 && muW_isJpsiBro[iT])?ntrees:iT;
      for(int ih=0;ih<hist.size()-nSpeHist;ih++){
	hist[ih][iT2]->Fill(v[ih],weight[iT]);
      }

      //AllDimuDeltaR
      hist[hist.size()-1][iT2]->Fill( dR_jpsi[iT] ,weight[iT]);
      hist[hist.size()-1][iT2]->Fill( dR_muWmi[iT] ,weight[iT]);
      hist[hist.size()-1][iT2]->Fill( dR_muWpl[iT] ,weight[iT]);


      // h_QQVProb[iT]->Fill(QQ_VtxProb[iT],weight[iT]);
      // h_alpha[iT]->Fill(Bc_alpha[iT],weight[iT]);
      // h_alpha3D[iT]->Fill(Bc_alpha3D[iT],weight[iT]);
      // h_TauSignif[iT]->Fill(Bc_ctauSignif[iT],weight[iT]);
      // h_TauSignif3D[iT]->Fill(Bc_ctauSignif3D[iT],weight[iT]);
      // h_mcorr[iT]->Fill(Bc_CorrM[iT],weight[iT]);
      // h_dca[iT]->Fill(QQ_dca[iT],weight[iT]);
      // h_QQMuWImbal[iT]->Fill(QQmuW_ptImbal[iT],weight[iT]);
      // h_SumDeltaR[iT]->Fill(dR_sum[iT],weight[iT]);
      // h_dRJpsiOverMuW[iT]->Fill(dR_jpsiOverMuW[iT],weight[iT]);
      // h_SumDxyMu[iT]->Fill(MuonDxySignif_sum[iT],weight[iT]);
      // h_SumDzMu[iT]->Fill(MuonDzSignif_sum[iT],weight[iT]);
      // h_BDT[iT]->Fill(BDT[iT],weight[iT]);

      // if(iT==5 && muW_isJpsiBro[iT]){
      // 	h_VProb[ntrees]->Fill(Bc_VtxProb[iT],weight[iT]);
      // 	h_QQVProb[ntrees]->Fill(QQ_VtxProb[iT],weight[iT]);
      // 	h_alpha[ntrees]->Fill(Bc_alpha[iT],weight[iT]);
      // 	h_alpha3D[ntrees]->Fill(Bc_alpha3D[iT],weight[iT]);
      // 	h_TauSignif[ntrees]->Fill(Bc_ctauSignif[iT],weight[iT]);
      // 	h_TauSignif3D[ntrees]->Fill(Bc_ctauSignif3D[iT],weight[iT]);
      // 	h_mcorr[ntrees]->Fill(Bc_CorrM[iT],weight[iT]);
      // 	h_dca[ntrees]->Fill(QQ_dca[iT],weight[iT]);
      // 	h_QQMuWImbal[ntrees]->Fill(QQmuW_ptImbal[iT],weight[iT]);
      // 	h_SumDeltaR[ntrees]->Fill(dR_sum[iT],weight[iT]);
      // 	h_dRJpsiOverMuW[ntrees]->Fill(dR_jpsiOverMuW[iT],weight[iT]);
      // 	h_SumDxyMu[ntrees]->Fill(MuonDxySignif_sum[iT],weight[iT]);
      // 	h_SumDzMu[ntrees]->Fill(MuonDzSignif_sum[iT],weight[iT]);
      // 	h_BDT[ntrees]->Fill(QQ_BDT[iT],weight[iT]);
      // }

    }
    //END event loop

  } 
  //END loop on trees

  for(int ih=0;ih<hist.size()-nSpeHist;ih++){
    DrawVar(hist[ih],ntrees,varName[ih],ispp,lowercut[ih],leftLeg[ih]);
  }
  DrawVar(hist[hist.size()-1],ntrees,varName[hist.size()-1],ispp,0,leftLeg[hist.size()-1],true);

  delete fullFile;
  
}
