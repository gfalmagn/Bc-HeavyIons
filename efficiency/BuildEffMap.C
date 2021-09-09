#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH2Poly.h"
#include "TLine.h"
#include "TClonesArray.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include <TVirtualFitter.h>
#include "../PbPb18/Utilities/EVENTUTILS.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"
#include "../helpers/SgMuonAcceptanceCuts.h"
#include "../helpers/AccEff2DBinning.h"
#include "./TnP/tnp_weight_lowptPbPb.h"
#include "./TnP/tnp_weight_pp_efficiencies.h"

double trimuEffPartial(double t1, double l2, double h2, double t2, double l3, double h3){
  return t1*( (l2-h2)*h3 + (h2-t2)*l3 );
}

double asymTrimuEff(double l1, double l2, double l3, double h1, double h2, double h3, double t1=-1, double t2=-1, double t3=-1){
  //assumes l<h<t (factorising efficiencies)
  if(t1==-1) t1 = h1; //for pp case when there is no L2-L3 trigger distinction
  if(t2==-1) t2 = h2;
  if(t3==-1) t3 = h3;

  // if(l1<h1-1e-5 || l2<h2-1e-5 || l3<h3-1e-5 || h1<t1-1e-5 || h2<t2-1e-5 || h3<t3-1e-5){
  //   cout<<"A loose efficiency is smaller than a tight one !!!!"<<endl;
  //   cout<<"l1, h1, t1 = "<<l1<<" "<<h1<<" "<<t1<<endl;
  //   cout<<"l2, h2, t2 = "<<l2<<" "<<h2<<" "<<t2<<endl;
  //   cout<<"l3, h3, t3 = "<<l3<<" "<<h3<<" "<<t3<<endl;
  // }

  //only happens with systematic variations of SF
  if(l1<h1) l1=h1;
  if(l2<h2) l2=h2;
  if(l3<h3) l3=h3;
  if(l1<t1) l1=t1;
  if(l2<t2) l2=t2;
  if(l3<t3) l3=t3;
  if(h1<t1) h1=t1;
  if(h2<t2) h2=t2;
  if(h3<t3) h3=t3;

  double res = t1*t2*t3;
  res += trimuEffPartial(t1, l2, h2, t2, l3, h3); //permutation(1,2,3)
  res += trimuEffPartial(t2, l3, h3, t3, l1, h1); //permutation(2,3,1)
  res += trimuEffPartial(t3, l1, h1, t1, l2, h2); //permutation(3,1,2)
  return res;
}

void drawEffMap(TH2Poly* hpEff, TH2Poly* hpSel, TH2Poly* hpAcc, TLine* l1, TLine* l2, TLine* l3, TLine* l4, TLine* l5, TLine* l6, bool ispp, TString nameSuf, float zmin, float zmax){

  TCanvas *c2 = new TCanvas("c2","c2",3000,1500);
  c2->Divide(2,1);

  c2->cd(1);
  //  gPad->SetLogz();
  //hpSel->GetZaxis()->SetRangeUser(1,200);
  hpSel->GetXaxis()->SetTitle("|y^{#mu#mu#mu}(B_{c})|");
  hpSel->GetYaxis()->SetTitle("p_{T}^{#mu#mu#mu}(B_{c})");
  hpSel->SetTitle("Selected B_{c}'s");
  hpSel->Draw("COLZ0");
  l1->Draw("same"); l2->Draw("same"); l3->Draw("same"); l4->Draw("same"); l5->Draw("same"); l6->Draw("same");

  c2->cd(2)->SetRightMargin(0.15);
  //gPad->SetLogz();
  hpEff->Divide(hpAcc);
  hpEff->GetZaxis()->SetRangeUser(zmin,zmax);
  hpEff->GetXaxis()->SetTitle("|y^{#mu#mu#mu}(B_{c})|");
  hpEff->GetYaxis()->SetTitle("p_{T}^{#mu#mu#mu}(B_{c})");
  hpEff->SetTitle("Efficiency "+nameSuf);
  hpEff->Draw("COLZ0");
  l1->Draw("same"); l2->Draw("same"); l3->Draw("same"); l4->Draw("same"); l5->Draw("same"); l6->Draw("same");

  c2->SaveAs("figs/EfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+nameSuf+(TString)(ispp?"_pp":"_PbPb")+".pdf");
  c2->SaveAs("figs/EfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+nameSuf+(TString)(ispp?"_pp":"_PbPb")+".png");

}

void BuildEffMap(bool ispp = true, bool runAEtoys=true, bool secondStep=false, bool runMCclos=false, bool BDTuncorrFromM=false, bool integratePtBins=false){
  
  if(ispp) runMCclos=false;

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  //**************************************************************  
  //Grab the variations of the pT bias of MC, from first step r1 and r2
  vector<TH1F*> bias,bias_2ndStep;
  TFile *BiasFile = TFile::Open("../twoSteps/pTBiases.root","READ");
  if(secondStep || runAEtoys){
    for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++)
      bias.push_back((TH1F*)BiasFile->Get("pTbias_"+(TString)(ispp?"pp":"PbPb")+"_var"+(TString)to_string(v)));
    if(secondStep && runAEtoys){
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++)
	bias_2ndStep.push_back((TH1F*)BiasFile->Get("pTbias_"+(TString)(ispp?"pp":"PbPb")+"_var"+(TString)to_string(v)+"_2ndStep"));
    }
  }

  vector<TH1F*> biasMCclos = vector<TH1F*>();
  vector<TH1F*> biasMCclos_2ndStep = vector<TH1F*>();
  if(secondStep || runAEtoys){
    for(int t=0;t<_nMCclos;t++){
      biasMCclos.push_back(runMCclos?((TH1F*)BiasFile->Get("pTbias_PbPb_MCclosure_toy"+(TString)to_string(t))):NULL);
      if(runAEtoys)
        biasMCclos_2ndStep.push_back((secondStep && runMCclos)?((TH1F*)BiasFile->Get("pTbias_PbPb_MCclosure_toy"+(TString)to_string(t)+"_2ndStep")):NULL);
    }
  }
  
  //**************************************************************  
  //Create Tree and branches
  TFile *fileMC = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/Oniatree_MC_Bc_trimuons_21112019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/BcToJpsiMuNu_BCVEGPY_PYTHIA8_2018PbPb5TeV_HINPbPbAutumn18DR-00196_08092020_4200k_ONIATREE.root");
  TTree* T = (TTree*)fileMC->Get("hionia/myTree");
  int nentries = T->GetEntries();
  std::cout<<"nevents MC = "<<nentries<<"\n";
  //T->Print();

  //**************************************************************  
  //Get input branches

  UInt_t eventNb;
  UInt_t runNb;
  UInt_t LS;
  Short_t nPV;
  int Centrality;
  ULong64_t HLTriggers; 
  Short_t Reco_3mu_size;
  Short_t Reco_3mu_charge[300];
  Short_t Reco_3mu_QQ1_idx[300];
  Short_t Reco_3mu_QQ2_idx[300];
  Short_t Reco_3mu_QQss_idx[300];
  Short_t Reco_3mu_muW_idx[300];
  Short_t Reco_3mu_muW2_idx[300];
  Short_t Reco_3mu_mumi_idx[300];
  Short_t Reco_3mu_mupl_idx[300];
  float Reco_3mu_ctau[300];
  float Reco_3mu_ctauErr[300];
  float Reco_3mu_ctau3D[300]; 
  float Reco_3mu_ctauErr3D[300];
  float Reco_3mu_cosAlpha[300]; 
  float Reco_3mu_cosAlpha3D[300];
  float Reco_3mu_VtxProb[300]; 
  float Reco_QQ_VtxProb[100]; 
  float Reco_QQ_dca[100]; 
  Short_t Reco_QQ_sign[100]; 
  Short_t Reco_QQ_mumi_idx[100];
  Short_t Reco_QQ_mupl_idx[100];
  Short_t Reco_mu_charge[300]; 
  float Reco_mu_dxy[100]; 
  float Reco_mu_dz[100]; 
  float Reco_mu_dxyErr[100];
  float Reco_mu_dzErr[100]; 
  float Reco_mu_normChi2_global[100];
  float Reco_mu_normChi2_inner[100]; 
  int Reco_mu_SelType[100]; 
  int Reco_mu_nPixWMea[100];
  int Reco_mu_nTrkWMea[100];
  ULong64_t Reco_mu_trig[100];
  bool Reco_mu_highPurity[100];
  bool Reco_mu_InLooseAcc[100];
  bool Reco_mu_InTightAcc[100];
  float Reco_3mu_muW_dz_muonlessVtx[300]; 
  float Reco_3mu_mumi_dz_muonlessVtx[300];
  float Reco_3mu_mupl_dz_muonlessVtx[300];
  float Reco_3mu_muW_dxy_muonlessVtx[300];
  float Reco_3mu_mumi_dxy_muonlessVtx[300];
  float Reco_3mu_mupl_dxy_muonlessVtx[300];
  Short_t Reco_QQ_size;
  Short_t Reco_mu_size;
  Short_t Reco_QQ_whichGen[100];
  Short_t Reco_3mu_whichGen[300];
  Short_t Gen_QQ_size;
  Short_t Gen_Bc_size;
  Short_t Gen_Bc_QQ_idx[100];
  Short_t Gen_Bc_muW_idx[100];
  Short_t Gen_QQ_mumi_idx[100];
  Short_t Gen_QQ_mupl_idx[100];
  Short_t Gen_3mu_whichRec[100];

  TBranch *b_eventNb = T->GetBranch("eventNb");
  b_eventNb->SetAddress(&eventNb);

  if(ispp){
    TBranch *b_nPV = T->GetBranch("nPV");
    b_nPV->SetAddress(&nPV);
  } else {
    TBranch *b_Centrality = T->GetBranch("Centrality");
    b_Centrality->SetAddress(&Centrality);
  }

  TBranch *b_HLTriggers = T->GetBranch("HLTriggers");
  b_HLTriggers->SetAddress(&HLTriggers);

  TBranch *b_Reco_3mu_size = T->GetBranch("Reco_3mu_size");
  b_Reco_3mu_size->SetAddress(&Reco_3mu_size);

  TBranch *b_Reco_3mu_charge = T->GetBranch("Reco_3mu_charge");
  b_Reco_3mu_charge->SetAddress(&Reco_3mu_charge);

  TBranch *b_Reco_3mu_QQ1_idx = T->GetBranch("Reco_3mu_QQ1_idx");
  b_Reco_3mu_QQ1_idx->SetAddress(&Reco_3mu_QQ1_idx);

  TBranch *b_Reco_3mu_QQ2_idx = T->GetBranch("Reco_3mu_QQ2_idx");
  b_Reco_3mu_QQ2_idx->SetAddress(&Reco_3mu_QQ2_idx);

  TBranch *b_Reco_3mu_QQss_idx = T->GetBranch("Reco_3mu_QQss_idx");
  b_Reco_3mu_QQss_idx->SetAddress(&Reco_3mu_QQss_idx);

  TBranch *b_Reco_3mu_muW2_idx = T->GetBranch("Reco_3mu_muW2_idx");
  b_Reco_3mu_muW2_idx->SetAddress(&Reco_3mu_muW2_idx);
    
  TBranch *b_Reco_3mu_muW_idx = T->GetBranch("Reco_3mu_muW_idx");
  b_Reco_3mu_muW_idx->SetAddress(&Reco_3mu_muW_idx);

  TBranch *b_Reco_3mu_mumi_idx = T->GetBranch("Reco_3mu_mumi_idx"); 
  b_Reco_3mu_mumi_idx->SetAddress(&Reco_3mu_mumi_idx);

  TBranch *b_Reco_3mu_mupl_idx = T->GetBranch("Reco_3mu_mupl_idx");
  b_Reco_3mu_mupl_idx->SetAddress(&Reco_3mu_mupl_idx);

  TBranch *b_Reco_3mu_ctau = T->GetBranch("Reco_3mu_ctau");
  b_Reco_3mu_ctau->SetAddress(&Reco_3mu_ctau);

  TBranch *b_Reco_3mu_ctauErr = T->GetBranch("Reco_3mu_ctauErr");
  b_Reco_3mu_ctauErr->SetAddress(&Reco_3mu_ctauErr);

  TBranch *b_Reco_3mu_ctau3D = T->GetBranch("Reco_3mu_ctau3D");
  b_Reco_3mu_ctau3D->SetAddress(&Reco_3mu_ctau3D);

  TBranch *b_Reco_3mu_ctauErr3D = T->GetBranch("Reco_3mu_ctauErr3D");
  b_Reco_3mu_ctauErr3D->SetAddress(&Reco_3mu_ctauErr3D);

  TBranch *b_Reco_3mu_cosAlpha = T->GetBranch("Reco_3mu_cosAlpha");
  b_Reco_3mu_cosAlpha->SetAddress(&Reco_3mu_cosAlpha);

  TBranch *b_Reco_3mu_cosAlpha3D = T->GetBranch("Reco_3mu_cosAlpha3D");
  b_Reco_3mu_cosAlpha3D->SetAddress(&Reco_3mu_cosAlpha3D);

  TBranch *b_Reco_3mu_VtxProb = T->GetBranch("Reco_3mu_VtxProb");
  b_Reco_3mu_VtxProb->SetAddress(&Reco_3mu_VtxProb);

  TBranch *b_Reco_QQ_VtxProb = T->GetBranch("Reco_QQ_VtxProb");
  b_Reco_QQ_VtxProb->SetAddress(&Reco_QQ_VtxProb);

  TBranch *b_Reco_QQ_dca = T->GetBranch("Reco_QQ_dca");
  b_Reco_QQ_dca->SetAddress(&Reco_QQ_dca);

  TBranch *b_Reco_QQ_sign = T->GetBranch("Reco_QQ_sign");
  b_Reco_QQ_sign->SetAddress(&Reco_QQ_sign);

  TBranch *b_Reco_QQ_mumi_idx = T->GetBranch("Reco_QQ_mumi_idx");
  b_Reco_QQ_mumi_idx->SetAddress(&Reco_QQ_mumi_idx);

  TBranch *b_Reco_QQ_mupl_idx = T->GetBranch("Reco_QQ_mupl_idx");
  b_Reco_QQ_mupl_idx->SetAddress(&Reco_QQ_mupl_idx);

  TBranch *b_Reco_mu_charge = T->GetBranch("Reco_mu_charge");
  b_Reco_mu_charge->SetAddress(&Reco_mu_charge);

  TBranch *b_Reco_mu_dxy = T->GetBranch("Reco_mu_dxy");
  b_Reco_mu_dxy->SetAddress(&Reco_mu_dxy);

  TBranch *b_Reco_mu_dz = T->GetBranch("Reco_mu_dz");
  b_Reco_mu_dz->SetAddress(&Reco_mu_dz);

  TBranch *b_Reco_mu_dxyErr = T->GetBranch("Reco_mu_dxyErr");
  b_Reco_mu_dxyErr->SetAddress(&Reco_mu_dxyErr);

  TBranch *b_Reco_mu_dzErr = T->GetBranch("Reco_mu_dzErr");
  b_Reco_mu_dzErr->SetAddress(&Reco_mu_dzErr);

  TBranch *b_Reco_mu_normChi2_global = T->GetBranch("Reco_mu_normChi2_global");
  b_Reco_mu_normChi2_global->SetAddress(&Reco_mu_normChi2_global);

  TBranch *b_Reco_mu_normChi2_inner = T->GetBranch("Reco_mu_normChi2_inner");
  b_Reco_mu_normChi2_inner->SetAddress(&Reco_mu_normChi2_inner);

  TBranch *b_Reco_mu_SelType = T->GetBranch("Reco_mu_SelectionType");
  b_Reco_mu_SelType->SetAddress(&Reco_mu_SelType);

  TBranch *b_Reco_mu_nPixWMea = T->GetBranch("Reco_mu_nPixWMea");
  b_Reco_mu_nPixWMea->SetAddress(&Reco_mu_nPixWMea);

  TBranch *b_Reco_mu_nTrkWMea = T->GetBranch("Reco_mu_nTrkWMea");
  b_Reco_mu_nTrkWMea->SetAddress(&Reco_mu_nTrkWMea);

  TBranch *b_Reco_mu_trig = T->GetBranch("Reco_mu_trig");
  b_Reco_mu_trig->SetAddress(&Reco_mu_trig);

  TBranch *b_Reco_mu_highPurity = T->GetBranch("Reco_mu_highPurity");
  b_Reco_mu_highPurity->SetAddress(&Reco_mu_highPurity);

  TBranch *b_Reco_mu_InLooseAcc = T->GetBranch("Reco_mu_InLooseAcc");
  b_Reco_mu_InLooseAcc->SetAddress(&Reco_mu_InLooseAcc);

  TBranch *b_Reco_mu_InTightAcc = T->GetBranch("Reco_mu_InTightAcc");
  b_Reco_mu_InTightAcc->SetAddress(&Reco_mu_InTightAcc);

  if(ispp){
    TBranch *b_Reco_3mu_muW_dz_muonlessVtx = T->GetBranch("Reco_3mu_muW_dz_muonlessVtx");
    b_Reco_3mu_muW_dz_muonlessVtx->SetAddress(&Reco_3mu_muW_dz_muonlessVtx);

    TBranch *b_Reco_3mu_mumi_dz_muonlessVtx = T->GetBranch("Reco_3mu_mumi_dz_muonlessVtx");
    b_Reco_3mu_mumi_dz_muonlessVtx->SetAddress(&Reco_3mu_mumi_dz_muonlessVtx);

    TBranch *b_Reco_3mu_mupl_dz_muonlessVtx = T->GetBranch("Reco_3mu_mupl_dz_muonlessVtx");
    b_Reco_3mu_mupl_dz_muonlessVtx->SetAddress(&Reco_3mu_mupl_dz_muonlessVtx);

    TBranch *b_Reco_3mu_muW_dxy_muonlessVtx = T->GetBranch("Reco_3mu_muW_dxy_muonlessVtx");
    b_Reco_3mu_muW_dxy_muonlessVtx->SetAddress(&Reco_3mu_muW_dxy_muonlessVtx);

    TBranch *b_Reco_3mu_mumi_dxy_muonlessVtx = T->GetBranch("Reco_3mu_mumi_dxy_muonlessVtx");
    b_Reco_3mu_mumi_dxy_muonlessVtx->SetAddress(&Reco_3mu_mumi_dxy_muonlessVtx);

    TBranch *b_Reco_3mu_mupl_dxy_muonlessVtx = T->GetBranch("Reco_3mu_mupl_dxy_muonlessVtx");
    b_Reco_3mu_mupl_dxy_muonlessVtx->SetAddress(&Reco_3mu_mupl_dxy_muonlessVtx);
  }

  TBranch *b_Reco_3mu_4mom;
  TClonesArray *Reco_3mu_4mom = new TClonesArray();
  T->SetBranchAddress("Reco_3mu_4mom", &Reco_3mu_4mom, &b_Reco_3mu_4mom);

  TBranch *b_Reco_mu_4mom;
  TClonesArray *Reco_mu_4mom = new TClonesArray();
  T->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);

  TBranch *b_Reco_QQ_4mom;
  TClonesArray *Reco_QQ_4mom = new TClonesArray();
  T->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);

  TBranch *b_Reco_QQ_size = T->GetBranch("Reco_QQ_size");
  b_Reco_QQ_size->SetAddress(&Reco_QQ_size);

  TBranch *b_Reco_mu_size = T->GetBranch("Reco_mu_size");
  b_Reco_mu_size->SetAddress(&Reco_mu_size);

  TBranch *b_Gen_Bc_size = T->GetBranch("Gen_Bc_size");
  b_Gen_Bc_size->SetAddress(&Gen_Bc_size);

  TBranch *b_Gen_QQ_size = T->GetBranch("Gen_QQ_size");
  b_Gen_QQ_size->SetAddress(&Gen_QQ_size);

  TBranch *b_Gen_Bc_QQ_idx = T->GetBranch("Gen_Bc_QQ_idx");
  b_Gen_Bc_QQ_idx->SetAddress(&Gen_Bc_QQ_idx);

  TBranch *b_Gen_Bc_muW_idx = T->GetBranch("Gen_Bc_muW_idx");
  b_Gen_Bc_muW_idx->SetAddress(&Gen_Bc_muW_idx);

  TBranch *b_Gen_QQ_mumi_idx = T->GetBranch("Gen_QQ_mumi_idx");
  b_Gen_QQ_mumi_idx->SetAddress(&Gen_QQ_mumi_idx);

  TBranch *b_Gen_QQ_mupl_idx = T->GetBranch("Gen_QQ_mupl_idx");
  b_Gen_QQ_mupl_idx->SetAddress(&Gen_QQ_mupl_idx);

  TBranch *b_Gen_3mu_whichRec = T->GetBranch("Gen_3mu_whichRec");
  b_Gen_3mu_whichRec->SetAddress(&Gen_3mu_whichRec);

  TBranch *b_Gen_QQ_4mom;
  TClonesArray *Gen_QQ_4mom = new TClonesArray();
  T->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);

  TBranch *b_Gen_3mu_4mom;
  TClonesArray *Gen_3mu_4mom = new TClonesArray();
  T->SetBranchAddress("Gen_3mu_4mom", &Gen_3mu_4mom, &b_Gen_3mu_4mom);

  TBranch *b_Gen_mu_4mom;
  TClonesArray *Gen_mu_4mom = new TClonesArray();
  T->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  TBranch *b_Reco_QQ_whichGen = T->GetBranch("Reco_QQ_whichGen");
  b_Reco_QQ_whichGen->SetAddress(&Reco_QQ_whichGen);

  TBranch *b_Reco_3mu_whichGen = T->GetBranch("Reco_3mu_whichGen");
  b_Reco_3mu_whichGen->SetAddress(&Reco_3mu_whichGen);

  //counters, histos
  vector<float> ntot(_NanaBins+1, 0);
  float nacc = 0, nsel = 0, nsel2 = 0, nBDT23 = 0, nBDT3 = 0;
  TH2Poly *hp = _hp();
  TH2Poly *hp_coarser = _hp_coarser();
  TH2Poly *hp_all = (TH2Poly*) hp->Clone("hp_all");
  vector<TH2Poly*> hp_acc = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hp_sel = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hp_sel_noSF = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hp_sel_selectiveSFapplication = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hp_sel_simpleAverage = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hpcoarse_sel = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hpcoarse_inBDT23 = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hpcoarse_inBDT3 = vector<TH2Poly*>(_NcentBins+1);
  TH1F *BDT23effVsPt = new TH1F("BDT23effVsPt","BDT23effVsPt",22,6,35);
  TH1F *BDT3effVsPt = new TH1F("BDT3effVsPt","BDT3effVsPt",22,6,35);
  TH1F *selectedVsPt = new TH1F("selectedVsPt","selectedVsPt",22,6,35);
  TH1F *SFs = new TH1F("SFs","SFs",100,ispp?0.86:0.86,ispp?1.2:1.4);
  TH1F *SFs_selectiveSFapplication = new TH1F("SFs_selectiveSFapplication","SFs_selectiveSFapplication",100,ispp?0.86:0.86,ispp?1.2:1.4);
  TH1F *SFs_simpleAverage = new TH1F("SFs_simpleAverage","SFs_simpleAverage",100,ispp?0.86:0.86,ispp?1.2:1.4);
  vector<vector<float> > passing_oneBinned(_NanaBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<vector<float> > passing_oneBinned_centDep(_NcentBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<vector<float> > passing_oneBinned_noSF(_NanaBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<vector<float> > eff_oneBinned(_NanaBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<vector<float> > eff_oneBinned_centDep(_NcentBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<vector<float> > eff_oneBinned_noSF(_NanaBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<float> accepted_oneBinned(_NanaBins+1,0);
  vector<float> accepted_oneBinned_centDep(_NcentBins+1,0);
  vector<vector<vector<float> > > passing_oneB(_NanaBins+1,vector<vector<float> >(4,vector<float>(5,0))); //4 is: muid (or glb), trk (in pp), trg (or muidtrg), tot
  vector<vector<vector<float> > > passing_oneB_centDep(_NcentBins+1,vector<vector<float> >(4,vector<float>(5,0))); //4 is: muid (or glb), trk (in pp), trg (or muidtrg), tot
  //MC with biased pT distributions
  vector<vector<float> > accepted_oneB_biased(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0));
  vector<vector<float> > passing_oneB_biased(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0));
  vector<vector<float> > eff_oneB_biased(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 1));
  vector<vector<float> > accepted_oneB_centDep_biased(_NcentBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0));
  vector<vector<float> > passing_oneB_centDep_biased(_NcentBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0));
  vector<vector<float> > eff_oneB_centDep_biased(_NcentBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 1));
  vector<vector<float> > accepted_MCclos(_NanaBins+1, vector<float>(_nMCclos, 0));
  vector<vector<float> > passing_MCclos(_NanaBins+1, vector<float>(_nMCclos, 0));
  vector<vector<float> > eff_MCclos(_NanaBins+1, vector<float>(_nMCclos, 1));

  for(int centb=0;centb<=_NcentBins;centb++){
    hpcoarse_sel[centb] = (TH2Poly*) hp_coarser->Clone("hpcoarse_sel"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hpcoarse_inBDT23[centb] = (TH2Poly*) hp_coarser->Clone("hpcoarse_inBDT23"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hpcoarse_inBDT3[centb] = (TH2Poly*) hp_coarser->Clone("hpcoarse_inBDT3"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hp_acc[centb] = (TH2Poly*) hp->Clone("hp_acc"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hp_sel[centb] = (TH2Poly*) hp->Clone("hp_sel"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hp_sel_noSF[centb] = (TH2Poly*) hp->Clone("hp_sel_noSF"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hp_sel_selectiveSFapplication[centb] = (TH2Poly*) hp->Clone("hp_sel_selectiveSFapplication"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hp_sel_simpleAverage[centb] = (TH2Poly*) hp->Clone("hp_sel_simpleAverage"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
  }

  //Jpsi mass histo for JpsiChoiceWeight
  vector<TH1F*> JpsiM(_nChan(ispp)+1);
  vector<TH1F*> JpsiM_tight(_nChan(ispp)+1);
  TFile *f_JpsiM = TFile::Open("../BDT/JpsiMassDistr.root","READ");
  for(int k=0;k<=_nChan(ispp);k++){
    JpsiM[k] = (TH1F*)f_JpsiM->Get("JpsiMass_data_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)))+(TString)(ispp?"_pp":"_PbPb"));
    JpsiM_tight[k] = (TH1F*)f_JpsiM->Get("JpsiMassCentralEta_data_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)))+(TString)(ispp?"_pp":"_PbPb"));
  }

  //Map to find BDT values
  std::map<std::pair<UInt_t,int>,float> bdtval;

  //**************************************************************
  //Pre-selected tree to get BDT efficiency map
  auto preselFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  TTree* T_presel = (TTree*)preselFile->Get("signal_MC");

  UInt_t evNb; TBranch *b_evNb = T_presel->GetBranch("eventNb");
  b_evNb->SetAddress(&evNb);
  float Bc_Pt; TBranch *b_Bc_Pt = T_presel->GetBranch("Bc_Pt");
  b_Bc_Pt->SetAddress(&Bc_Pt);
  float Bc_Y; TBranch *b_Bc_Y = T_presel->GetBranch("Bc_Y");
  b_Bc_Y->SetAddress(&Bc_Y);
  float weight; TBranch *b_weight = T_presel->GetBranch(secondStep?"weight2":"weight");
  b_weight->SetAddress(&weight);
  float BDT; TBranch *b_BDT = T_presel->GetBranch(secondStep?"BDT2":"BDT");
  b_BDT->SetAddress(&BDT);
  float QQ_dca; TBranch *b_QQ_dca = T_presel->GetBranch("QQ_dca");
  b_QQ_dca->SetAddress(&QQ_dca);
  int Centr; 
  if(!ispp){
    TBranch *b_Centr = T_presel->GetBranch("Centrality");
    b_Centr->SetAddress(&Centr);
  }

  for(int j=0; j<T_presel->GetEntries(); j++){
    T_presel->GetEntry(j);
    nsel2 += weight;
    selectedVsPt->Fill(Bc_Pt, weight);

    bdtval[make_pair(evNb,100000*QQ_dca)] = BDT;

    int kinb = 0;
    int centBin = 0;
    if(!integratePtBins) kinb = (Bc_Pt<_BcPtmax[1])?1:2;    
    if(!ispp) centBin = (Centr<_Centmax[1])?1:2;

    for(int centb=0;centb<=_NcentBins;centb++){
      if(ispp && centb>0) break;
      if(centb>0 && centb!=centBin) continue;
      
      hpcoarse_sel[centb]->Fill(fabs(Bc_Y),Bc_Pt, weight);
      if(BDT>_BDTcuts(ispp,kinb,centBin,secondStep,BDTuncorrFromM)[1]){
	if(centb==0){
	  nBDT23 += weight;
	  BDT23effVsPt->Fill(Bc_Pt, weight);
	}
	hpcoarse_inBDT23[centb]->Fill(fabs(Bc_Y),Bc_Pt, weight);}
      if(BDT>_BDTcuts(ispp,kinb,centBin,secondStep,BDTuncorrFromM)[2]){
	if(centb==0){
	  nBDT3 += weight;
	  BDT3effVsPt->Fill(Bc_Pt, weight);
	}
	hpcoarse_inBDT3[centb]->Fill(fabs(Bc_Y),Bc_Pt, weight);}
    }
  }

  BDT23effVsPt->Divide(selectedVsPt);
  BDT3effVsPt->Divide(selectedVsPt);
  for(int centb=0;centb<=_NcentBins;centb++){
    if(ispp && centb>0) break;
    hpcoarse_inBDT23[centb]->Divide(hpcoarse_sel[centb]);
    hpcoarse_inBDT3[centb]->Divide(hpcoarse_sel[centb]);
  }
  float eff_BDT23 = (float)nBDT23/(float)nsel2;
  float eff_BDT3 = (float)nBDT3/(float)nsel2;

  
  float npass=0,npass2=0;
  //**************************************************************
  //loop on events for preselection efficiency
  for(int j=0; j<nentries; j++){//nentries
    if(j%100000==0){ cout<<"Scanned "<<100.*(double)j/nentries<<"% of entries"<<endl; }
    
    Reco_3mu_4mom->Clear();
    Reco_mu_4mom->Clear();
    Reco_QQ_4mom->Clear();
    Gen_3mu_4mom->Clear();
    Gen_mu_4mom->Clear();
    Gen_QQ_4mom->Clear();
    
    T->GetEntry(j);

    int centBin=0;
    if(!ispp){
      if((float)Centrality < 2*_Centmin[0] || (float)Centrality >= 2*_Centmax[0]) continue;
      for(int centb=1;centb<=_NcentBins;centb++){
	if((float)Centrality >= 2*_Centmin[centb] && (float)Centrality < 2*_Centmax[centb])
	  centBin = centb;
      }
      if (centBin==0) continue; //keep 0-90% centrality
    }

    for(int igen=0;igen<Gen_Bc_size;igen++){
      TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom->At(igen);

      float weight = _scaleMCsig[ispp];
      if(!ispp) {float weightNcoll = (float)findNcoll(Centrality); //for PbPb MC
	weight *= weightNcoll;}
      weight *= (secondStep?( getBias( bias[_nomMethVar] , gen3mu->Pt()) ):1);
      float wadd = (secondStep && runAEtoys)?( getBias(bias_2ndStep[_nomMethVar],gen3mu->Pt()) ):1;
      weight *= wadd;

      //for MC closure test                                                                                                                                                                                                                 
      vector<float> wt = vector<float>(_nMCclos, ispp?1:(_scaleMCsig[ispp]*((float)findNcoll(Centrality))));
      if(runMCclos){
	for(int t=0;t<_nMCclos;t++){
	  if(!secondStep && !runAEtoys) wt[t] *= MCclosurePTw(gen3mu->Pt(),t);
	  wt[t] *= secondStep?( getBias(biasMCclos[t],gen3mu->Pt()) ):1;
	  wt[t] *= (secondStep && runAEtoys)?( getBias(biasMCclos_2ndStep[t],gen3mu->Pt()) ):1;
	}
      }

      ntot[0] += weight;
      hp_all->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight);
      
      int genQQidx = Gen_Bc_QQ_idx[igen];
      TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom->At(genQQidx);
      TLorentzVector *genBc_muW = (TLorentzVector*) Gen_mu_4mom->At(Gen_Bc_muW_idx[igen]);
      TLorentzVector *genBc_mumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[genQQidx]);
      TLorentzVector *genBc_mupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[genQQidx]);
      
      for(int b=1;b<=_NanaBins;b++){
	if(inFidCuts(b,gen3mu->Pt(),gen3mu->Rapidity())){
	  ntot[b] += weight;
	}
      }
      
      if(!InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,_withTM)) continue;
      for(int b=1;b<=_NanaBins;b++){
	if(inFidCuts(b,gen3mu->Pt(),gen3mu->Rapidity())){
	  accepted_oneBinned[b] += weight;
	  accepted_oneBinned[0] += weight;
	  if(!ispp){
	    accepted_oneBinned_centDep[centBin] += weight;
	  }

	  //MC closure
	  if(runMCclos){
	    for(int t=0;t<_nMCclos;t++){
	      accepted_MCclos[0][t] += wt[t];
	      accepted_MCclos[b][t] += wt[t];
	    }
	  }

	  //biased MC
	  if(runAEtoys && !runMCclos){
	    if(!secondStep) cout<<"BEWARE! The pT bias might not be defined for all variations for the first step!"<<endl;
	    for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	      float wtmp = getBias( (secondStep?bias_2ndStep:bias)[v] , gen3mu->Pt()) * weight / wadd;
	      accepted_oneB_biased[0][v] += wtmp; 
	      accepted_oneB_biased[b][v] += wtmp;
	      if(!ispp){
		accepted_oneB_centDep_biased[centBin][v] += wtmp;
	      }
	    }
	  }

	}//end accepted in good ana bin
      }
      
      hp_acc[0]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight);
      if(!ispp)
	hp_acc[centBin]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight);
      
      int irec = Gen_3mu_whichRec[igen];
      if(irec<0) continue;

      int QQidx_2[2] = { Reco_3mu_QQ1_idx[irec], Reco_3mu_QQ2_idx[irec] };
      int muWidx_2[2] = { Reco_3mu_muW_idx[irec], Reco_3mu_muW2_idx[irec] };
      int mumiidx_2[2] = { Reco_3mu_mumi_idx[irec], Reco_QQ_mumi_idx[QQidx_2[1]] };
      int muplidx_2[2] = { Reco_3mu_mupl_idx[irec], Reco_QQ_mupl_idx[QQidx_2[1]] };
      
      for(int k=0; k<2; k++){
	float wei = weight;
	vector<float> weit = wt;

	int k2 = (k==0)?1:0;
	int QQidx = QQidx_2[k];
	int QQ2idx = QQidx_2[k2];
	int muWidx = muWidx_2[k];
	int mumiidx = mumiidx_2[k];
	int muplidx = muplidx_2[k];
	
	if(QQidx>-1){
	  TLorentzVector *recBc = (TLorentzVector*) Reco_3mu_4mom->At(irec);
	  TLorentzVector *recQQ = (TLorentzVector*) Reco_QQ_4mom->At(QQidx);
	  if(muWidx==-1 || mumiidx==-1 || muplidx==-1) {cout<<"!!!!!!! one muon has a -1 index !! Skip this candidate"<<endl; continue;}
	  TLorentzVector *recBc_muW = (TLorentzVector*) Reco_mu_4mom->At(muWidx);
	  TLorentzVector *recBc_mumi = (TLorentzVector*) Reco_mu_4mom->At(mumiidx);
	  TLorentzVector *recBc_mupl = (TLorentzVector*) Reco_mu_4mom->At(muplidx);
	  
	  float BcCandE = sqrt(pow(recQQ->P(),2) + pow(m_Jpsi,2) ) + recBc_muW->E() ;
	  float BcCandM = sqrt( pow( BcCandE ,2) - pow(recBc->P(),2) );
	  float QQM = recQQ->M();
	  
	  float muW_eta = recBc_muW->Eta();
	  float mumi_eta = recBc_mumi->Eta();
	  float mupl_eta = recBc_mupl->Eta();
	  float maxEta = max(fabs(muW_eta),max(fabs(mumi_eta),fabs(mupl_eta)));

	  float muW_pt = recBc_muW->Pt();
	  float mumi_pt = recBc_mumi->Pt();
	  float mupl_pt = recBc_mupl->Pt();
	  
	  bool goodTree = fabs(Reco_3mu_charge[irec])==1 && Reco_QQ_sign[QQidx]==0 && inLooseMassRange(QQM) // in Jpsi mass region
	    && (BcCandM < _mBcMax) && (BcCandM > _mBcMin) // in Bc mass region
	    && Reco_3mu_whichGen[irec]>-1
	    && Reco_QQ_whichGen[QQidx]>-1;
	  if(inJpsiMassSB(QQM, maxEta<1.5)) {
	    wei *= -1;
	    for(int t=0;t<_nMCclos;t++){
	      weit[t] *= -1;
	    }	    
	  }
	  else if(!(inJpsiMassRange(QQM, maxEta<1.5))){ wei = 0;}

	  if(!goodTree || wei == 0) continue;
	  
	  //**************************************************************
	  float Bc_ctauSignif = Reco_3mu_ctau[irec] / Reco_3mu_ctauErr[irec] ;
	  float Bc_ctauSignif3D = Reco_3mu_ctau3D[irec] / Reco_3mu_ctauErr3D[irec] ;
	  float Bc_alpha = TMath::ACos(Reco_3mu_cosAlpha[irec]);
	  float Bc_alpha3D = TMath::ACos(Reco_3mu_cosAlpha3D[irec]);
	  float Bc_VtxProb = Reco_3mu_VtxProb[irec];
	  float QQ_VtxProb = Reco_QQ_VtxProb[QQidx];
	  float QQ_dca = Reco_QQ_dca[QQidx];
	  float QQ2_VtxProb = (QQ2idx<Reco_QQ_size && QQ2idx>-1)?(Reco_QQ_VtxProb[QQ2idx]):0;
	  float QQ2_dca = (QQ2idx<Reco_QQ_size && QQ2idx>-1)?(Reco_QQ_dca[QQ2idx]):100;

	  bool muW_inLooseAcc = Reco_mu_InLooseAcc[muWidx];
	  bool muW_inTightAcc = Reco_mu_InTightAcc[muWidx];
	  bool muW_isGlb = (Reco_mu_SelType[muWidx]&2)>0;
	  bool muW_trig = (Reco_mu_trig[muWidx]&(ispp?8:4096))>0 ; //DoubleMu0 trigger = 2^3 for pp //HL_THIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 = 2^12 for PbPb
          bool muW_L2 = muW_trig && (ispp || (Reco_mu_trig[muWidx]&(65536))>0) ;
          bool muW_L3 = muW_trig && (ispp || (Reco_mu_trig[muWidx]&(131072))>0) ;
	  bool muW_isSoft = isSoft( false, muW_isGlb , (Reco_mu_SelType[muWidx]&8)>0 , ((Reco_mu_SelType[muWidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[muWidx] , Reco_mu_dxy[muWidx] , Reco_mu_dz[muWidx] , Reco_mu_nPixWMea[muWidx] , Reco_mu_nTrkWMea[muWidx] );

	  bool mumi_inLooseAcc = Reco_mu_InLooseAcc[mumiidx];
	  bool mumi_inTightAcc = Reco_mu_InTightAcc[mumiidx];
	  bool mumi_isGlb = (Reco_mu_SelType[mumiidx]&2)>0;
	  bool mumi_trig = (Reco_mu_trig[mumiidx]&(ispp?8:4096))>0;
          bool mumi_L2 = mumi_trig && (ispp || (Reco_mu_trig[mumiidx]&(65536))>0) ;
          bool mumi_L3 = mumi_trig && (ispp || (Reco_mu_trig[mumiidx]&(131072))>0) ;
	  bool mumi_isSoft = isSoft( false, mumi_isGlb , (Reco_mu_SelType[mumiidx]&8)>0  , ((Reco_mu_SelType[mumiidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[mumiidx] , Reco_mu_dxy[mumiidx] , Reco_mu_dz[mumiidx] , Reco_mu_nPixWMea[mumiidx] , Reco_mu_nTrkWMea[mumiidx] );

	  bool mupl_inLooseAcc = Reco_mu_InLooseAcc[muplidx];
	  bool mupl_inTightAcc = Reco_mu_InTightAcc[muplidx];
	  bool mupl_isGlb = (Reco_mu_SelType[muplidx]&2)>0;
	  bool mupl_trig = (Reco_mu_trig[muplidx]&(ispp?8:4096))>0;
          bool mupl_L2 = mupl_trig && (ispp || (Reco_mu_trig[muplidx]&(65536))>0) ;
          bool mupl_L3 = mupl_trig && (ispp || (Reco_mu_trig[muplidx]&(131072))>0) ;
	  bool mupl_isSoft = isSoft( false, mupl_isGlb , (Reco_mu_SelType[muplidx]&8)>0  , ((Reco_mu_SelType[muplidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[muplidx] , Reco_mu_dxy[muplidx] , Reco_mu_dz[muplidx] , Reco_mu_nPixWMea[muplidx] , Reco_mu_nTrkWMea[muplidx] );

	  float muW_dxy = ispp?(Reco_3mu_muW_dxy_muonlessVtx[irec]):(Reco_mu_dxy[muWidx]);
	  float mumi_dxy = ispp?(Reco_3mu_mumi_dxy_muonlessVtx[irec]):(Reco_mu_dxy[mumiidx]);
	  float mupl_dxy = ispp?(Reco_3mu_mupl_dxy_muonlessVtx[irec]):(Reco_mu_dxy[muplidx]);
	  float muW_dz = ispp?(Reco_3mu_muW_dz_muonlessVtx[irec]):(Reco_mu_dz[muWidx]);
	  float mumi_dz = ispp?(Reco_3mu_mumi_dz_muonlessVtx[irec]):(Reco_mu_dz[mumiidx]);
	  float mupl_dz = ispp?(Reco_3mu_mupl_dz_muonlessVtx[irec]):(Reco_mu_dz[muplidx]);

	  // cout<<(Bc_ctauSignif>1.5)<<" "<<(Bc_alpha<0.8)<<" "<<(Bc_alpha3D<0.8)<<" "<<(Bc_VtxProb>0.01)<<" "<<(QQ_dca<0.3)<<" "<<(muW_isSoft)<<" "<<( mumi_isSoft)<<" "<<(mupl_isSoft)<<" "<<( (muW_isGlb && muW_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) || (muW_isGlb && muW_inLooseAcc && mumi_isGlb && mumi_inLooseAcc) || (mumi_isGlb && mumi_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) )<<" "<<((muW_trig && mupl_trig && muW_inTightAcc && mupl_inTightAcc ) ||(muW_trig && mumi_trig && muW_inTightAcc && mumi_inTightAcc ) || (mumi_trig && mupl_trig && mumi_inTightAcc && mupl_inTightAcc ))<<" "<<(fabs(muW_dz)<(ispp?0.6:0.8) && fabs(mumi_dz)<(ispp?0.6:0.8) && fabs(mupl_dz)<(ispp?0.6:0.8))<<" "<<((HLTriggers&((ispp)?8:4096))>0)<<endl;
	  if(
	     //**************************************************************
	     // Pre-selections
	     Bc_ctauSignif>_ctauSignif_cut
	     && Bc_ctauSignif3D>_ctauSignif3D_cut //maybe not a good idea in pp due to pile-up along z. But actually not a big deal if we select nonprompt objects
	     && Bc_alpha<_alpha_cut(ispp)
	     && Bc_alpha3D<_alpha3D_cut(ispp) //combined with alpha cut, kills 0.3% of signal 
	     && Bc_VtxProb>_vtxProb_cut //kills 1.2% more signal compared to Bc_VtxProb>0.005
	     && QQ_VtxProb>_QQvtxProb_cut //drop the QQ_VtxProb, too correlated with Bc_VtxProb?
	     && QQ_dca<_QQdca_cut //keep this that kills 1.8% of signal and 10% of WRONSIGN/BCMASS
	     && muW_isSoft
	     && mumi_isSoft
	     && mupl_isSoft
	     && (_withTM?( (muW_isGlb && muW_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) || //only one muon can be tracker and out of LooseAcceptance
			  (muW_isGlb && muW_inLooseAcc && mumi_isGlb && mumi_inLooseAcc) || 
			  (mumi_isGlb && mumi_inLooseAcc && mupl_isGlb && mupl_inLooseAcc)
			  ):(
			     mumi_isGlb && mupl_isGlb && muW_isGlb
			     && mumi_inLooseAcc && mupl_inLooseAcc && muW_inLooseAcc           
			     ))
	     && ( ( muW_trig && mupl_trig && muW_inTightAcc && mupl_inTightAcc ) || //two muons among three must trigger //BEWARE ! Not sure if TightAcceptance should be put there
		  ( muW_trig && mumi_trig && muW_inTightAcc && mumi_inTightAcc ) ||
		  ( mumi_trig && mupl_trig && mumi_inTightAcc && mupl_inTightAcc ) //only this last option can be true for dimuon+trk
		  )
	     && fabs(muW_dz)<0.6 && fabs(mumi_dz)<0.6 && fabs(mupl_dz)<0.6
	     && (!ispp || (HLTriggers&(ispp?8:4096))>0) //the event must fire the trigger as well
	     ){

	    float QQ2_M = (Reco_mu_charge[muWidx]>0)?((*recBc_mumi+*recBc_muW).M()):((*recBc_mupl+*recBc_muW).M()); //QQ2 is the second OS pair
	    float QQ3_M = (Reco_mu_charge[muWidx]>0)?((*recBc_mupl+*recBc_muW).M()):((*recBc_mumi+*recBc_muW).M()); //QQ3 is the SS pair
	    float QQM_correction = m_Jpsi / QQM;

	    //corrected mass selection
	    float BcP = sqrt( pow(QQM_correction * recQQ->Px()+recBc_muW->Px() ,2) + pow(QQM_correction * recQQ->Py()+recBc_muW->Py() ,2) + pow(QQM_correction * recQQ->Pz()+recBc_muW->Pz() ,2) );
            float PperpTrimu = TMath::Sin(Bc_alpha3D) * BcP;
            float Bc_CorrM = sqrt(BcCandM*BcCandM + PperpTrimu*PperpTrimu) + PperpTrimu;
            if(Bc_CorrM>_BcCorrM_cut(ispp)) continue;

	    //selection on sum of deltaR's
	    float dRjpsi = (float) recBc_mumi->DeltaR(*recBc_mupl);
	    float dRmuWmi = (float) recBc_muW->DeltaR(*recBc_mumi);
	    float dRmuWpl = (float) recBc_muW->DeltaR(*recBc_mupl);
	    vector<float> dRcorr = dRcorrectedForQQM(dRjpsi,dRmuWmi,dRmuWpl,QQM_correction,mumi_pt,mupl_pt);
	    float dR_sum = dRcorr[0] + dRcorr[1] + dRcorr[2];
	    if(dR_sum>_dRsum_cut) continue;

	    float dR_sum_QQ2 = 8;
	    if(inJpsiMassRange(QQ2_M, maxEta<1.5) || inJpsiMassSB(QQ2_M, maxEta<1.5) ){
	      float QQ2M_correction = m_Jpsi / QQ2_M;
	      TLorentzVector *recBc_muW2 = (TLorentzVector*) Reco_mu_4mom->At(muWidx_2[k2]);
	      TLorentzVector *recBc_mumi2 = (TLorentzVector*) Reco_mu_4mom->At(mumiidx_2[k2]);
	      TLorentzVector *recBc_mupl2 = (TLorentzVector*) Reco_mu_4mom->At(muplidx_2[k2]);
	      float dRjpsi2 = (float) recBc_mumi2->DeltaR(*recBc_mupl2);
	      float dRmuWmi2 = (float) recBc_muW2->DeltaR(*recBc_mumi2);
	      float dRmuWpl2 = (float) recBc_muW2->DeltaR(*recBc_mupl2);
	      vector<float> dRcorr2 = dRcorrectedForQQM(dRjpsi2,dRmuWmi2,dRmuWpl2,QQ2M_correction,recBc_mumi2->Pt(),recBc_mupl2->Pt());
	      dR_sum_QQ2 = dRcorr2[0] + dRcorr2[1] + dRcorr2[2];
	    }

	    if(muW_trig && !(muW_L2 || muW_L3)) cout<<"!!!! PROBLEM with muW trigger consistency"<<endl;
	    if(mupl_trig && !(mupl_L2 || mupl_L3)) cout<<"!!!! PROBLEM with mupl trigger consistency"<<endl;
	    if(mumi_trig && !(mumi_L2 || mumi_L3)) cout<<"!!!! PROBLEM with mumi trigger consistency"<<endl;
		
	    //**** Deal with the Jpsi dimuon choice
	    if(
	       (inJpsiMassRange(QQ2_M, maxEta<1.5) || inJpsiMassSB(QQ2_M, maxEta<1.5))
	       && Reco_QQ_VtxProb[QQ2idx]>_QQvtxProb_cut && Reco_QQ_dca[QQ2idx]<_QQdca_cut && Reco_QQ_dca[QQ2idx]>0 && dR_sum_QQ2<_dRsum_cut
	       // (inJpsiMassRange(QQM, maxEta<1.5) && inJpsiMassSB(QQ2_M, maxEta<1.5)
	       // 	&& Reco_QQ_VtxProb[QQ2idx]>_QQvtxProb_cut && Reco_QQ_dca[QQ2idx]<_QQdca_cut && Reco_QQ_dca[QQ2idx]>0)
	       // || (inJpsiMassSB(QQM, maxEta<1.5) && inJpsiMassRange(QQ2_M, maxEta<1.5)
	       // 	   && Reco_QQ_VtxProb[QQ2idx]>_QQvtxProb_cut && Reco_QQ_dca[QQ2idx]<_QQdca_cut && Reco_QQ_dca[QQ2idx]>0)
	       ){

	      float bdtv = bdtval.find(make_pair(eventNb,100000*QQ_dca))->second;
	      int kbin = 0;
	      for(int k=0;k<_nChan(ispp);k++){
		if(bdtv>_BDTcuts(ispp,0,0,false)[k] && bdtv<_BDTcuts(ispp,0,0,false)[k+1]) kbin = k+1;
	      }

	      float binc_QQ1 = ( (maxEta<1.5)?JpsiM_tight:JpsiM )[kbin]->GetBinContent(( (maxEta<1.5)?JpsiM_tight:JpsiM )[kbin]->FindBin(QQM));
	      float binc_QQ2 = ( (maxEta<1.5)?JpsiM_tight:JpsiM )[kbin]->GetBinContent(( (maxEta<1.5)?JpsiM_tight:JpsiM )[kbin]->FindBin(QQ2_M));

	      if (binc_QQ1==0) {
		wei = 0;
		if(runMCclos){
		  for(int t=0;t<_nMCclos;t++){
		    weit[t] = 0;
		  }	    
		}
	      }
	      else {
		wei *= binc_QQ1 / (binc_QQ1+binc_QQ2);
		if(runMCclos){
		  for(int t=0;t<_nMCclos;t++){
		    weit[t] *= binc_QQ1 / (binc_QQ1+binc_QQ2);
		  }
		}
	      }
	    }

	    //**** Apply the scale factors to correct MC efficiencies
	    //H = tight criteria (tight acceptance, triggering muon...), L = loose criteria (loose acceptance, no trigger required...), T = fire L3 trigger (in PbPb case) -> implies H
	    //2 muons have to respect H, 1 has to respect L
	    //In PbPb, at least 1 triggering muon must pass L3
	    //muon 1 is mumi, muon 2 is mupl, muon3 is muW
	    double effcorr=1; //effcorr is the correction to be applied
	    double effcorr_simpleAverage=1, effcorr_selectiveSFapplication=1; 
	    vector<int> erIdx = {2,1,-2,-1,0}; //statlo, stathi, systlo, systhi, nominal
	    
	    for(int sftype=0; sftype<3; sftype++){ //nominal result is run in sftype==2
	      if(ispp && sftype==1) continue;
	      for(int id=0;id<5;id++){
		if(id!=4 && (integratePtBins || BDTuncorrFromM)) continue;
		int er = erIdx[id];
		int IDer = (sftype==0)?er:0; //id (PbPb) or glb (pp)
		int TKer = (sftype==1)?er:0; //tracker (PbPb)
		int TGer = (sftype==2)?er:0; //trigger (PbPb) or muidtrig (pp)

		double SF1_T,SF2_T,SF3_T; //Definition of the scale-factor for the 3 muons, L or H or T criteria
		double SF1_H,SF2_H,SF3_H; //[trk syst or stat][muid/glb][trg/muidtrg]
		double SF1_L,SF2_L,SF3_L; //[trk/glb syst or stat][muid]
		double L1,L2,L3,H1,H2,H3,T1,T2,T3;//Definition of the MC efficiencies

		if (!ispp){//Case PbPb
		  SF1_L = tnp_weight_muid_looseacceptance_pbpb(mumi_pt,mumi_eta,IDer) * tnp_weight_trk_looseacceptance_pbpb(mumi_eta,TKer);//Product of muid and trk, loose acceptance
		  SF2_L = tnp_weight_muid_looseacceptance_pbpb(mupl_pt,mupl_eta,IDer) * tnp_weight_trk_looseacceptance_pbpb(mupl_eta,TKer);
		  SF3_L = tnp_weight_muid_looseacceptance_pbpb(muW_pt,muW_eta,IDer) * tnp_weight_trk_looseacceptance_pbpb(muW_eta,TKer);
		  L1 = tnp_weight_muid_looseacceptance_pbpb(mumi_pt,mumi_eta,3) * tnp_weight_trk_looseacceptance_pbpb(mumi_eta,3); //When varying the SF with the index, this efficiency should not vary
		  L2 = tnp_weight_muid_looseacceptance_pbpb(mupl_pt,mupl_eta,3) * tnp_weight_trk_looseacceptance_pbpb(mupl_eta,3);
		  L3 = tnp_weight_muid_looseacceptance_pbpb(muW_pt,muW_eta,3) * tnp_weight_trk_looseacceptance_pbpb(muW_eta,3);

		  SF1_H = tnp_weight_muid_pbpb(mumi_pt,mumi_eta,IDer) * tnp_weight_trk_pbpb(mumi_eta,TKer) * tnp_weight_trg_pbpb(mumi_pt,mumi_eta,0,TGer);//Product of muid, trk, and L2 trg, tight acceptance
		  SF2_H = tnp_weight_muid_pbpb(mupl_pt,mupl_eta,IDer) * tnp_weight_trk_pbpb(mupl_eta,TKer) * tnp_weight_trg_pbpb(mupl_pt,mupl_eta,0,TGer);
		  SF3_H = tnp_weight_muid_pbpb(muW_pt,muW_eta,IDer) * tnp_weight_trk_pbpb(muW_eta,TKer) * tnp_weight_trg_pbpb(muW_pt,muW_eta,0,TGer);
		  H1 = tnp_weight_muid_pbpb(mumi_pt,mumi_eta,3) * tnp_weight_trk_pbpb(mumi_eta,3) * tnp_weight_trg_pbpb(mumi_pt,mumi_eta,0,3);
		  H2 = tnp_weight_muid_pbpb(mupl_pt,mupl_eta,3) * tnp_weight_trk_pbpb(mupl_eta,3) * tnp_weight_trg_pbpb(mupl_pt,mupl_eta,0,3);
		  H3 = tnp_weight_muid_pbpb(muW_pt,muW_eta,3) * tnp_weight_trk_pbpb(muW_eta,3) * tnp_weight_trg_pbpb(muW_pt,muW_eta,0,3);

		  SF1_T = tnp_weight_muid_pbpb(mumi_pt,mumi_eta,IDer) * tnp_weight_trk_pbpb(mumi_eta,TKer) * tnp_weight_trg_pbpb(mumi_pt,mumi_eta,1,TGer);//Product of muid, trk, and L3 trg, tight acceptance
		  SF2_T = tnp_weight_muid_pbpb(mupl_pt,mupl_eta,IDer) * tnp_weight_trk_pbpb(mupl_eta,TKer) * tnp_weight_trg_pbpb(mupl_pt,mupl_eta,1,TGer);
		  SF3_T = tnp_weight_muid_pbpb(muW_pt,muW_eta,IDer) * tnp_weight_trk_pbpb(muW_eta,TKer) * tnp_weight_trg_pbpb(muW_pt,muW_eta,1,TGer);
		  T1 = tnp_weight_muid_pbpb(mumi_pt,mumi_eta,3) * tnp_weight_trk_pbpb(mumi_eta,3) * tnp_weight_trg_pbpb(mumi_pt,mumi_eta,1,3);
		  T2 = tnp_weight_muid_pbpb(mupl_pt,mupl_eta,3) * tnp_weight_trk_pbpb(mupl_eta,3) * tnp_weight_trg_pbpb(mupl_pt,mupl_eta,1,3);
		  T3 = tnp_weight_muid_pbpb(muW_pt,muW_eta,3) * tnp_weight_trk_pbpb(muW_eta,3) * tnp_weight_trg_pbpb(muW_pt,muW_eta,1,3);
		  //Both the muons that respect H have to be triggering, one has to pass the  L3 filter, one has to pass the L2 filter

		  // This part might be redundant
		  if(!mumi_inTightAcc) {H1 = 0; T1 = 0; SF1_H = 1; SF1_T = 1;} //SF=1 is only to avoid NaN's
		  if(!mupl_inTightAcc) {H2 = 0; T2 = 0; SF2_H = 1; SF2_T = 1;}
		  if(!muW_inTightAcc) {H3 = 0; T3 = 0; SF3_H = 1; SF3_T = 1;}

		  double effMC = asymTrimuEff(L1, L2, L3, H1, H2, H3, T1, T2, T3);
		  double lh[] = {SF1_L*L1, SF2_L*L2, SF3_L*L3, SF1_H*H1, SF2_H*H2, SF3_H*H3,  SF1_T*T1, SF2_T*T2, SF3_T*T3};
		  for(int x=0;x<9;x++) //prevent eff>1 due to SF variation
		    if(lh[x]>1) lh[x]=1;
		  double effdata = asymTrimuEff(lh[0], lh[1], lh[2], lh[3], lh[4], lh[5], lh[6], lh[7], lh[8]);
		  effcorr = effdata/effMC;

		  if(er==0 && sftype==2){//only for nominal sf
		    effcorr_simpleAverage = (SF1_L*SF2_H*SF3_T*(H2>0 && T3>0) + SF3_L*SF1_H*SF2_T*(H1>0 && T2>0) + SF2_L*SF3_H*SF1_T*(H3>0 && T1>0)
					     + SF1_L*SF3_H*SF2_T*(H3>0 && T2>0) + SF3_L*SF2_H*SF1_T*(H2>0 && T1>0) + SF2_L*SF1_H*SF3_T*(H1>0 && T3>0))
		      / ((H2>0 && T3>0) + (H1>0 && T2>0) + (H3>0 && T1>0) + (H3>0 && T2>0) + (H2>0 && T1>0) + (H1>0 && T3>0));

		    // We do NOT want to require mumi_trig for the trigger efficiency on mumi to be taken into account !! (at least in the nominal method, but maybe we can do it in the simple method)
		    if(!mumi_inTightAcc || !mumi_trig) {H1 = 0; T1 = 0; SF1_H = 1; SF1_T = 1;} //SF=1 is only to avoid NaN's
		    if(!mumi_inTightAcc || !mumi_L3) {T1 = 0; SF1_T = 1;}
		    if(!mupl_inTightAcc || !mupl_trig) {H2 = 0; T2 = 0; SF2_H = 1; SF2_T = 1;}
		    if(!mupl_inTightAcc || !mupl_L3) {T1 = 0; SF1_T = 1;}
		    if(!muW_inTightAcc || !muW_trig) {H3 = 0; T3 = 0; SF3_H = 1; SF3_T = 1;}
		    if(!muW_inTightAcc || !muW_L3) {T3 = 0; SF3_T = 1;}
	      
		    double effMC_selectiveSFapplication = asymTrimuEff(L1, L2, L3, H1, H2, H3, T1, T2, T3);
		    double effdata_selectiveSFapplication = asymTrimuEff(SF1_L*L1, SF2_L*L2, SF3_L*L3, SF1_H*H1, SF2_H*H2, SF3_H*H3, SF1_T*T1, SF2_T*T2, SF3_T*T3);
		    effcorr_selectiveSFapplication = effdata_selectiveSFapplication/effMC_selectiveSFapplication;
		  }
		}//end PbPb

		else {//Case pp
		  SF1_H = tnp_weight_glb_tightacceptance_pp(mumi_pt,mumi_eta,IDer)*tnp_weight_muidtrg_tightacceptance_pp(mumi_pt,mumi_eta,TGer);//H is glb and muidtrg, tight acceptance
		  SF2_H = tnp_weight_glb_tightacceptance_pp(mupl_pt,mupl_eta,IDer)*tnp_weight_muidtrg_tightacceptance_pp(mupl_pt,mupl_eta,TGer);
		  SF3_H = tnp_weight_glb_tightacceptance_pp(muW_pt,muW_eta,IDer)*tnp_weight_muidtrg_tightacceptance_pp(muW_pt,muW_eta,TGer);
		  H1 = tnp_weight_glb_tightacceptance_pp(mumi_pt,mumi_eta,3)*tnp_weight_muidtrg_tightacceptance_pp(mumi_pt,mumi_eta,3);
		  H2 = tnp_weight_glb_tightacceptance_pp(mupl_pt,mupl_eta,3)*tnp_weight_muidtrg_tightacceptance_pp(mupl_pt,mupl_eta,3);
		  H3 = tnp_weight_glb_tightacceptance_pp(muW_pt,muW_eta,3)*tnp_weight_muidtrg_tightacceptance_pp(muW_pt,muW_eta,3);

		  SF1_L = tnp_weight_muid_looseacceptance_pp(mumi_pt,mumi_eta,TGer)*tnp_weight_glb_looseacceptance_pp(mumi_pt,mumi_eta,IDer);// L is glb and muid, loose acceptance
		  SF2_L = tnp_weight_muid_looseacceptance_pp(mupl_pt,mupl_eta,TGer)*tnp_weight_glb_looseacceptance_pp(mupl_pt,mupl_eta,IDer);
		  SF3_L = tnp_weight_muid_looseacceptance_pp(muW_pt,muW_eta,IDer)*tnp_weight_glb_looseacceptance_pp(muW_pt,muW_eta,IDer);
		  L1 = tnp_weight_muid_looseacceptance_pp(mumi_pt,mumi_eta,3)*tnp_weight_glb_looseacceptance_pp(mumi_pt,mumi_eta,3);
		  L2 = tnp_weight_muid_looseacceptance_pp(mupl_pt,mupl_eta,3)*tnp_weight_glb_looseacceptance_pp(mupl_pt,mupl_eta,3);
		  L3 = tnp_weight_muid_looseacceptance_pp(muW_pt,muW_eta,3)*tnp_weight_glb_looseacceptance_pp(muW_pt,muW_eta,3);

		  // This part might be redundant
		  if(!mumi_inTightAcc) {H1 = 0; SF1_H = 1;} //SF=1 is only to avoid NaN's
		  if(!mupl_inTightAcc) {H2 = 0; SF2_H = 1;}
		  if(!muW_inTightAcc) {H3 = 0; SF3_H = 1;}

		  //This happens somtimes when tnp_weight_glb_tightacceptance_pp > tnp_weight_glb_looseacceptance_pp
		  if(H1>L1) H1 *= tnp_weight_glb_looseacceptance_pp(mumi_pt,mumi_eta,3) / tnp_weight_glb_tightacceptance_pp(mumi_pt,mumi_eta,3);
		  if(H2>L2) H2 *= tnp_weight_glb_looseacceptance_pp(mupl_pt,mupl_eta,3) / tnp_weight_glb_tightacceptance_pp(mupl_pt,mupl_eta,3);
		  if(H3>L3) H3 *= tnp_weight_glb_looseacceptance_pp(muW_pt,muW_eta,3) / tnp_weight_glb_tightacceptance_pp(muW_pt,muW_eta,3);

		  double effMC = asymTrimuEff(L1, L2, L3, H1, H2, H3);
		  // cout<<"IDer, TGer = "<<IDer<<" "<<TGer<<endl;
		  // cout<<"l1, l2, l3, h1, h2, h3 = "<<SF1_L*L1<<" "<< SF2_L*L2<<" "<< SF3_L*L3<<" "<< SF1_H*H1<<" "<< SF2_H*H2<<" "<< SF3_H*H3<<endl;
		  double lh[] = {SF1_L*L1, SF2_L*L2, SF3_L*L3, SF1_H*H1, SF2_H*H2, SF3_H*H3};
		  for(int x=0;x<6;x++) //prevent eff>1 due to SF variation
		    if(lh[x]>1) lh[x]=1;
		  double effdata = asymTrimuEff(lh[0], lh[1], lh[2], lh[3], lh[4], lh[5]);
		  effcorr = effdata/effMC;
	      
		  if(er==0 && sftype==2){//only for nominal sf
		    effcorr_simpleAverage = (SF1_L*SF2_H*SF3_H*(H2>0 && H3>0) + SF3_L*SF1_H*SF2_H*(H1>0 && H2>0) + SF2_L*SF3_H*SF1_H*(H3>0 && H1>0))
		      / ((H2>0 && H3>0) + (H1>0 && H2>0) + (H3>0 && H1>0));

		    // We do NOT want to require mumi_trig for the trigger efficiency on mumi to be taken into account !! (at least in the nominal method, but maybe we can do it in the simple method)
		    if(!mumi_inTightAcc || !mumi_trig) {H1 = 0; SF1_H = 1;} //SF=1 is only to avoid NaN's
		    if(!mupl_inTightAcc || !mupl_trig) {H2 = 0; SF2_H = 1;}
		    if(!muW_inTightAcc || !muW_trig) {H3 = 0; SF3_H = 1;}

		    double effMC_selectiveSFapplication = asymTrimuEff(L1, L2, L3, H1, H2, H3);
		    double effdata_selectiveSFapplication = asymTrimuEff(SF1_L*L1, SF2_L*L2, SF3_L*L3, SF1_H*H1, SF2_H*H2, SF3_H*H3);
		    effcorr_selectiveSFapplication = effdata_selectiveSFapplication/effMC_selectiveSFapplication;
		  }	      
		}//end pp

		//**** Actually fill the 2D histogram
		for(int b=1;b<=_NanaBins;b++){
		  if(inFidCuts(b,gen3mu->Pt(),gen3mu->Rapidity())){
		    passing_oneB[b][sftype][id] += wei * effcorr;
		    passing_oneB[0][sftype][id] += wei * effcorr;
		    if(!ispp){
		      passing_oneB_centDep[centBin][sftype][id] += wei * effcorr;
		    }
		    if(er==0 && sftype==2){//nominal concerning SF variation
		      //MC closure
		      if(runMCclos){
			for(int t=0;t<_nMCclos;t++){
			  passing_MCclos[0][t] += weit[t] * effcorr;
			  passing_MCclos[b][t] += weit[t] * effcorr;
			}
		      }

		      if(runAEtoys && !runMCclos){
			//pT biased MC
			for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
			  float wtmp = getBias( (secondStep?bias_2ndStep:bias)[v] , gen3mu->Pt()) * wei * effcorr / wadd;
			  passing_oneB_biased[0][v] += wtmp;
			  passing_oneB_biased[b][v] += wtmp;
			  if(!ispp){
			    passing_oneB_centDep_biased[centBin][v] += wtmp;
			  }
			}
		      }
		    }
		  }
		}

	      }//end loop on Error index
	    }//end loop on sf error type

	    //last ran variation is the nominal (er=0), so effcorr has the nominal value
	    hp_sel_noSF[0]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei);
	    hp_sel[0]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei * effcorr);
	    hp_sel_simpleAverage[0]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei * effcorr_simpleAverage);
	    hp_sel_selectiveSFapplication[0]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei * effcorr_selectiveSFapplication);
	    if(!ispp){
	      hp_sel_noSF[centBin]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei);
	      hp_sel[centBin]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei * effcorr);
	      hp_sel_simpleAverage[centBin]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei * effcorr_simpleAverage);
	      hp_sel_selectiveSFapplication[centBin]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), wei * effcorr_selectiveSFapplication);
	    }

	    for(int b=1;b<=_NanaBins;b++){
	      if(inFidCuts(b,gen3mu->Pt(),gen3mu->Rapidity())){
		passing_oneBinned_noSF[b][0] += wei;
		passing_oneBinned_noSF[0][0] += wei;
	      }
	    }

	    SFs->Fill(effcorr, wei);
	    SFs_simpleAverage->Fill(effcorr_simpleAverage, wei);
	    SFs_selectiveSFapplication->Fill(effcorr_selectiveSFapplication, wei);

	  } //end if passes full selection
	      
	} //end if(QQ_isValid)
      } //end loop on 2/3 possible Jpsi dimuon choice
    } //end loop on Bc candidates
  } //end loop on entries



  //Calculate errors and systematics
  //pT bins
  for(int b=0;b<=_NanaBins;b++){
    //nominal
    passing_oneBinned[b][0] = passing_oneB[b][0][4]; //nominal has the same value for [muid][4] or [trg][4] or [trk][4]

    //stat errors
    float idStatErr = (fabs(passing_oneB[b][0][0]-passing_oneB[b][0][4]) + fabs(passing_oneB[b][0][1]-passing_oneB[b][0][4]))/2; //average of deviations of statLo [0] and statHi [1] from nominal [4]
    float trkStatErr;
    if(!ispp) trkStatErr = (fabs(passing_oneB[b][1][0]-passing_oneB[b][1][4]) + fabs(passing_oneB[b][1][1]-passing_oneB[b][1][4]))/2; //average of deviations of statLo [0] and statHi [1] from nominal [4]
    else trkStatErr = 0;
    float trgStatErr = (fabs(passing_oneB[b][2][0]-passing_oneB[b][2][4]) + fabs(passing_oneB[b][2][1]-passing_oneB[b][2][4]))/2; //average of deviations of statLo [0] and statHi [1] from nominal [4]
    passing_oneBinned[b][1] = sqrt(pow(idStatErr,2)+pow(trkStatErr,2)+pow(trgStatErr,2));

    //syst errors
    float idSystErr = (fabs(passing_oneB[b][0][2]-passing_oneB[b][0][4]) + fabs(passing_oneB[b][0][3]-passing_oneB[b][0][4]))/2; //average of deviations of systLo [2] and systHi [3] from nominal [4]
    float trkSystErr;
    if(!ispp) trkSystErr = (fabs(passing_oneB[b][1][2]-passing_oneB[b][1][4]) + fabs(passing_oneB[b][1][3]-passing_oneB[b][1][4]))/2; //average of deviations of systLo [2] and systHi [3] from nominal [4]
    else trkSystErr = 0;
    float trgSystErr = (fabs(passing_oneB[b][2][2]-passing_oneB[b][2][4]) + fabs(passing_oneB[b][2][3]-passing_oneB[b][2][4]))/2; //average of deviations of systLo [2] and systHi [3] from nominal [4]
    passing_oneBinned[b][2] = sqrt(pow(idSystErr,2)+pow(trkSystErr,2)+pow(trgSystErr,2));

    //total error
    passing_oneBinned[b][3] = sqrt(pow(passing_oneBinned[b][1],2)+pow(passing_oneBinned[b][2],2));

    //**************************************************************
    //One-binned efficiency
    for(int e=0;e<4;e++)
      eff_oneBinned[b][e] = passing_oneBinned[b][e]/accepted_oneBinned[b]; 
    eff_oneBinned_noSF[b][0] = passing_oneBinned_noSF[b][0]/accepted_oneBinned[b];

    //MC with biased pT
    if(runAEtoys && !runMCclos){
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	eff_oneB_biased[b][v] = passing_oneB_biased[b][v]/accepted_oneB_biased[b][v];
      }
    }

    //MC closure
    if(runMCclos){
      for(int t=0;t<_nMCclos;t++){
	eff_MCclos[b][t] = passing_MCclos[b][t]/accepted_MCclos[b][t];
	cout<<"MC closure: t,b,passing, accepted, eff = "<<t<<" "<<b<<" "<<passing_MCclos[b][t]<<" "<<accepted_MCclos[b][t]<<" "<<eff_MCclos[b][t] <<endl;
      }
    }

    if(b>0){
      // cout<<"b nom stat syst tot = "<<b<<" "<<passing_oneBinned[b][0]<<" "<<passing_oneBinned[b][1]<<" "<<passing_oneBinned[b][2]<<" "<<passing_oneBinned[b][3]<<endl;
      // cout<<"id stat syst = "<<idStatErr<<" "<<idSystErr<<endl;
      // cout<<"trk stat syst = "<<trkStatErr<<" "<<trkSystErr<<endl;
      // cout<<"trg stat syst = "<<trgStatErr<<" "<<trgSystErr<<endl;
      cout<<"ngen, accepted_oneBinned, passing_oneBinned, efficiency = " <<ntot[b]<<" "<<accepted_oneBinned[b]<<" "<<passing_oneBinned[b][0]<<" "<<eff_oneBinned[b][0]<<" "<<endl;
      cout<<"err passing_oneBinned, relerr efficiency = " <<passing_oneBinned[b][3]<<" "<<eff_oneBinned[b][3]/eff_oneBinned[b][0]<<" "<<endl;
    }
    cout<<"bin "<<b<<" efficiency with or without scale factors, and ratio = "<<eff_oneBinned[b][0]<<" "<<eff_oneBinned_noSF[b][0]<<" "<<eff_oneBinned[b][0]/eff_oneBinned_noSF[b][0]<<endl;
  }

  //Centrality bins
  if(!ispp){
    for(int b=1;b<=_NcentBins;b++){
      //nominal
      passing_oneBinned_centDep[b][0] = passing_oneB_centDep[b][0][4]; //nominal has the same value for [muid][4] or [trg][4] or [trk][4]

      //stat errors
      float idStatErr = (fabs(passing_oneB_centDep[b][0][0]-passing_oneB_centDep[b][0][4]) + fabs(passing_oneB_centDep[b][0][1]-passing_oneB_centDep[b][0][4]))/2; //average of deviations of statLo [0] and statHi [1] from nominal [4]
      float trkStatErr;
      if(!ispp) trkStatErr = (fabs(passing_oneB_centDep[b][1][0]-passing_oneB_centDep[b][1][4]) + fabs(passing_oneB_centDep[b][1][1]-passing_oneB_centDep[b][1][4]))/2; //average of deviations of statLo [0] and statHi [1] from nominal [4]
      else trkStatErr = 0;
      float trgStatErr = (fabs(passing_oneB_centDep[b][2][0]-passing_oneB_centDep[b][2][4]) + fabs(passing_oneB_centDep[b][2][1]-passing_oneB_centDep[b][2][4]))/2; //average of deviations of statLo [0] and statHi [1] from nominal [4]
      passing_oneBinned_centDep[b][1] = sqrt(pow(idStatErr,2)+pow(trkStatErr,2)+pow(trgStatErr,2));

      //syst errors
      float idSystErr = (fabs(passing_oneB_centDep[b][0][2]-passing_oneB_centDep[b][0][4]) + fabs(passing_oneB_centDep[b][0][3]-passing_oneB_centDep[b][0][4]))/2; //average of deviations of systLo [2] and systHi [3] from nominal [4]
      float trkSystErr;
      if(!ispp) trkSystErr = (fabs(passing_oneB_centDep[b][1][2]-passing_oneB_centDep[b][1][4]) + fabs(passing_oneB_centDep[b][1][3]-passing_oneB_centDep[b][1][4]))/2; //average of deviations of systLo [2] and systHi [3] from nominal [4]
      else trkSystErr = 0;
      float trgSystErr = (fabs(passing_oneB_centDep[b][2][2]-passing_oneB_centDep[b][2][4]) + fabs(passing_oneB_centDep[b][2][3]-passing_oneB_centDep[b][2][4]))/2; //average of deviations of systLo [2] and systHi [3] from nominal [4]
      passing_oneBinned_centDep[b][2] = sqrt(pow(idSystErr,2)+pow(trkSystErr,2)+pow(trgSystErr,2));

      //total error
      passing_oneBinned_centDep[b][3] = sqrt(pow(passing_oneBinned_centDep[b][1],2)+pow(passing_oneBinned_centDep[b][2],2));

      //**************************************************************
      //One-binned efficiency
      for(int e=0;e<4;e++)
	eff_oneBinned_centDep[b][e] = passing_oneBinned_centDep[b][e] / accepted_oneBinned_centDep[b];
      //MC with biased pT
      if(runAEtoys && !runMCclos){
	for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	  eff_oneB_centDep_biased[b][v] = passing_oneB_centDep_biased[b][v]/accepted_oneB_centDep_biased[b][v];
	}
      }

      if(b>0){
	// cout<<"b nom stat syst tot = "<<b<<" "<<passing_oneBinned_centDep[b][0]<<" "<<passing_oneBinned_centDep[b][1]<<" "<<passing_oneBinned_centDep[b][2]<<" "<<passing_oneBinned_centDep[b][3]<<endl;
	// cout<<"id stat syst = "<<idStatErr<<" "<<idSystErr<<endl;
	// cout<<"trk stat syst = "<<trkStatErr<<" "<<trkSystErr<<endl;
	// cout<<"trg stat syst = "<<trgStatErr<<" "<<trgSystErr<<endl;
	cout<<"ngen, accepted_oneBinned_centDep, passing_oneBinned_centDep, efficiency_centDep = " <<ntot[b]<<" "<<accepted_oneBinned_centDep[b]<<" "<<passing_oneBinned_centDep[b][0]<<" "<<eff_oneBinned_centDep[b][0]<<" "<<endl;
	cout<<"err passing_oneBinned_centDep, relerr efficiency_centDep = " <<passing_oneBinned_centDep[b][3]<<" "<<eff_oneBinned_centDep[b][3]/eff_oneBinned_centDep[b][0]<<" "<<endl;
      }
    }

    //fill the integrated bin of centralityDep, for consistency
    for(int e=0;e<4;e++)
      eff_oneBinned_centDep[0][e] = eff_oneBinned[0][e];
    if(runAEtoys && !runMCclos){
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++)
	eff_oneB_centDep_biased[0][v] = eff_oneB_biased[0][v];
    }

  }

  cout<<"ngen, accepted_oneBinned[0], passing_oneBinned[0], efficiency[0] = "<<ntot[0]<<" "<<accepted_oneBinned[0]<<" "<<passing_oneBinned[0][0]<<" "<<eff_oneBinned[0][0]<<" "<<endl;
  cout<<"err passing_oneBinned[0], relerr efficiency[0] = " <<passing_oneBinned[0][3]<<" "<<eff_oneBinned[0][3]/eff_oneBinned[0][0]<<" "<<endl;

  //**************************************************************
  //Lines for fiducial cuts 
  TLine *line1 = new TLine(_BcYmin[0],_BcPtmax[1],_BcYmin[1],_BcPtmax[1]);
  TLine *line2 = new TLine(_BcYmin[1],_BcPtmin[1],_BcYmin[1],_BcPtmax[1]);
  TLine *line3 = new TLine(_BcYmin[1],_BcPtmin[1],_BcYmax[1],_BcPtmin[1]);
  TLine *line4 = new TLine(_BcYmax[1],_BcPtmin[0],_BcYmax[1],_BcPtmax[0]);
  TLine *line5 = new TLine(_BcYmin[1],_BcPtmax[1],_BcYmax[1],_BcPtmax[1]);
  TLine *line6 = new TLine(_BcYmin[0],_BcPtmax[0],_BcYmax[0],_BcPtmax[0]);
  line1->SetLineWidth(4);  line1->SetLineColor(kBlack);
  line2->SetLineWidth(4);  line2->SetLineColor(kBlack);
  line3->SetLineWidth(4);  line3->SetLineColor(kBlack);
  line4->SetLineWidth(4);  line4->SetLineColor(kBlack);
  line5->SetLineWidth(3); line5->SetLineStyle(2); line5->SetLineColor(kBlack);
  line6->SetLineWidth(4);  line6->SetLineColor(kBlack);

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(50);

  //**************************************************************
  //Draw TH2Poly for efficiencies

  vector<TH2Poly*> hp_efficiency = vector<TH2Poly*>(_NcentBins);
  vector<TH2Poly*> hp_efficiency_noSF = vector<TH2Poly*>(_NcentBins);
  vector<TH2Poly*> hp_efficiency_selectiveSFapplication = vector<TH2Poly*>(_NcentBins);
  vector<TH2Poly*> hp_efficiency_simpleAverage = vector<TH2Poly*>(_NcentBins);
  if(!runMCclos){
    for(int centb=0;centb<=_NcentBins;centb++){
      if(ispp && centb>0) break;

      hp_efficiency[centb] = (TH2Poly*)hp_sel[centb]->Clone("hp_efficiency"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
      if(!integratePtBins && !BDTuncorrFromM) {
	drawEffMap(hp_efficiency[centb], hp_sel[centb], hp_acc[centb], line1, line2, line3, line4, line5, line6, ispp, ""+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb)))+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""), 0,0.45);
      }

      hp_efficiency_noSF[centb] = (TH2Poly*)hp_sel_noSF[centb]->Clone("hp_efficiency_noSF"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
      if(!integratePtBins && !BDTuncorrFromM){
	hp_efficiency_noSF[centb]->Divide(hp_efficiency[centb]);
	drawEffMap(hp_efficiency_noSF[centb], hp_sel_noSF[centb], hp_acc[centb], line1, line2, line3, line4, line5, line6, ispp, "_noSF_DivideNominal"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb)))+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""), ispp?0.84:0.83,ispp?1.16:1.17);
      }

      hp_efficiency_selectiveSFapplication[centb] = (TH2Poly*)hp_sel_selectiveSFapplication[centb]->Clone("hp_efficiency_selectiveSFapplication"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
      if(!integratePtBins && !BDTuncorrFromM){ 
	hp_efficiency_selectiveSFapplication[centb]->Divide(hp_efficiency[centb]);
	drawEffMap(hp_efficiency_selectiveSFapplication[centb], hp_sel_selectiveSFapplication[centb], hp_acc[centb], line1, line2, line3, line4, line5, line6, ispp, "_selectiveSFapplication_DivideNominal"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb)))+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""), ispp?0.92:0.89,ispp?1.08:1.11);
      }

      hp_efficiency_simpleAverage[centb] = (TH2Poly*)hp_sel_simpleAverage[centb]->Clone("hp_efficiency_simpleAverage"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
      if(!integratePtBins && !BDTuncorrFromM) {
	hp_efficiency_simpleAverage[centb]->Divide(hp_efficiency[centb]);
	drawEffMap(hp_efficiency_simpleAverage[centb], hp_sel_simpleAverage[centb], hp_acc[centb], line1, line2, line3, line4, line5, line6, ispp, "_simpleAverage_DivideNominal"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb)))+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""), ispp?0.9:0.84,ispp?1.1:1.16);
      }
    }
  }

  //**************************************************************
  //Draw scale factors
  if(!integratePtBins && !BDTuncorrFromM){
    TCanvas *c6 = new TCanvas("c6","c6",2500,2500);
    SFs->SetTitle("Scale factors for MC signal trimuons");
    SFs->SetLineWidth(2);
    SFs->GetXaxis()->SetTitle("Scale factors on efficiencies");
    SFs->GetYaxis()->SetTitle("Expected signal trimuons");
    SFs->GetYaxis()->SetRangeUser(0,1.2*SFs->GetMaximum());
    SFs->Draw("hist");
    cout<<"Mean of nominal SF's = "<<SFs->GetMean()<<endl;

    SFs_selectiveSFapplication->SetLineWidth(2);
    SFs_selectiveSFapplication->SetLineColor(kRed);
    SFs_selectiveSFapplication->Draw("histsame");
    cout<<"Mean of selectiveApplication SF's = "<<SFs_selectiveSFapplication->GetMean()<<endl;

    SFs_simpleAverage->SetLineWidth(2);
    SFs_simpleAverage->SetLineColor(kGreen+2);
    SFs_simpleAverage->Draw("histsame");
    cout<<"Mean of simpleAverage SF's = "<<SFs_simpleAverage->GetMean()<<endl;

    c6->SetLeftMargin(0.12);
    TLegend *leg = new TLegend(0.42,0.7,0.9,0.9);
    leg->SetTextSize(0.035);
    leg->SetHeader("MC trimuon scale factors");
    leg->SetBorderSize(0);
    leg->AddEntry(SFs,"nominal");
    leg->AddEntry(SFs_selectiveSFapplication,"SF only on passing muons");
    leg->AddEntry(SFs_simpleAverage,"Simple average of single-mu SF");
    leg->Draw("same");

    if(!runMCclos) c6->SaveAs("figs/ScaleFactorsComparison"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):"")+".pdf");
  }
  
  //**************************************************************
  //Grab acceptance map
  TFile acc_file("../acceptance/acceptanceMap.root","READ");
  vector<TH2Poly*> hp_acceptance = vector<TH2Poly*>(_NcentBins+1);
  vector<TH2Poly*> hp_acceff = vector<TH2Poly*>(_NcentBins+1);

  for(int centb=0;centb<=_NcentBins;centb++){
    if(ispp && centb>0) break;

    hp_acceptance[centb]= (TH2Poly*)acc_file.Get("hp_acceptance_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?(runAEtoys?"_3rdStep":"_2ndStep"):""));
    hp_acceff[centb] = (TH2Poly*)hp_acceptance[centb]->Clone("hp_acceff"+(TString)((centb==0)?"":("_centBin"+(TString)to_string(centb))));
    hp_acceff[centb]->SetDirectory(0);
    hp_acceptance[centb]->SetDirectory(0);
    if(!runMCclos) hp_acceff[centb]->Multiply(hp_efficiency[centb]);

    if(!integratePtBins && !BDTuncorrFromM && !runMCclos){
      TCanvas *c3 = new TCanvas("c3","c3",2500,2500);
      c3->Divide(2,2);
  
      c3->cd(1);
      hp_acceptance[centb]->Draw("COLZ");
      line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same"); line6->Draw("same");
      gPad->SetLogz();
      gPad->SetRightMargin(0.15);
      hp_acceptance[centb]->GetZaxis()->SetRangeUser(5e-5,1);

      TPaletteAxis *palette = (TPaletteAxis*)hp_acceptance[centb]->GetListOfFunctions()->FindObject("palette");
      // the following lines moe the paletter. Choose the values you need for the position.
      palette->SetX1NDC(0.86);
      palette->SetX2NDC(0.91);
      palette->SetY1NDC(0.1);
      palette->SetY2NDC(0.9);
      gPad->Modified();
      gPad->Update();
  
      c3->cd(2);
      hp_efficiency[centb]->SetTitle("Efficiency");
      hp_efficiency[centb]->Draw("COLZ0");
      gPad->SetRightMargin(0.15);
      line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same"); line6->Draw("same");
      c3->cd(3);
      hp_sel[centb]->Draw("COLZ0");
      line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same"); line6->Draw("same");

      c3->cd(4);
      hp_acceff[centb]->Draw("COLZ0");
      line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same"); line6->Draw("same");
      hp_acceff[centb]->SetTitle("Acceptance #times Efficiency");
      cout<<"smallest in fid cuts acc eff = "<<hp_acceff[centb]->GetBinContent(hp_acceff[centb]->FindBin(1.32,6.05))<<endl;
      gPad->SetLogz();
      gPad->SetRightMargin(0.15);
      hp_acceff[centb]->GetZaxis()->SetRangeUser(3e-5,0.45);

      TPaletteAxis *palette2 = (TPaletteAxis*)hp_acceff[centb]->GetListOfFunctions()->FindObject("palette");
      // the following lines moe the paletter. Choose the values you need for the position.
      palette2->SetX1NDC(0.86);
      palette2->SetX2NDC(0.91);
      palette2->SetY1NDC(0.1);
      palette2->SetY2NDC(0.9);
      gPad->Modified();
      gPad->Update();

      if(centb==0){
	c3->SaveAs("figs/AcceptanceEfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):"")+".pdf");
	c3->SaveAs("figs/AcceptanceEfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):"")+".png");
      }

      //Save AccEff Pad only
      c3->cd(4)->SetRightMargin(0.18);
      c3->cd(4)->SetTopMargin(0.04);
      hp_acceff[centb]->SetTitle("");
      hp_acceff[centb]->GetZaxis()->SetTitle("acceptance #times efficiency");
      hp_acceff[centb]->GetZaxis()->SetTitleOffset(1.5);

      TLatex CMStag;
      CMStag.SetNDC();
      CMStag.SetTextFont(42);
      CMStag.SetTextSize(0.041);
      CMStag.SetTextAlign(11);
      CMStag.DrawLatex(0.13,0.16,"#splitline{#font[61]{CMS}}{#font[52]{Preliminary}}");

      palette2->SetX1NDC(0.83);
      palette2->SetX2NDC(0.88);
      palette2->SetY1NDC(0.1);
      palette2->SetY2NDC(0.96);
      gPad->Modified();
      gPad->Update();

      c3->SaveAs("figs/AcceptanceEfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):"")+"_ForAccEffOnly.pdf");

      //**************************************************************
      //Draw TH2Poly for BDT efficiency
      TCanvas *c4 = new TCanvas("c4","c4",1500,1500);
      c4->Divide(2,2);

      c4->cd(1);
      hpcoarse_inBDT23[centb]->GetZaxis()->SetRangeUser(eff_BDT23-0.2, eff_BDT23+0.2);
      hpcoarse_inBDT23[centb]->SetTitle("Efficiency for BDT bin 2-3;|Y|;p_{T} [GeV]");
      hpcoarse_inBDT23[centb]->Draw("COLZ");

      c4->cd(2);
      hpcoarse_inBDT3[centb]->GetZaxis()->SetRangeUser(eff_BDT3-0.2, eff_BDT3+0.2);
      hpcoarse_inBDT3[centb]->SetTitle("Efficiency for BDT bin 3;|Y|;p_{T} [GeV]");
      hpcoarse_inBDT3[centb]->Draw("COLZ");

      c4->cd(3);
      BDT23effVsPt->GetYaxis()->SetRangeUser(eff_BDT23-0.2, eff_BDT23+0.2);
      BDT23effVsPt->SetTitle("Efficiency for BDT bin 2-3;p_{T} [GeV]");
      BDT23effVsPt->Draw("E");

      c4->cd(4);
      BDT3effVsPt->GetYaxis()->SetRangeUser(eff_BDT3-0.2, eff_BDT3+0.2);
      BDT3effVsPt->SetTitle("Efficiency for BDT bin 3;p_{T} [GeV]");
      BDT3effVsPt->Draw("E");

      if(centb==0)
	c4->SaveAs("figs/BDTefficiencyMap"+(TString)(_withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):"")+".pdf");
    }
  }

  //**************************************************************
  //Store maps
  TFile* fout = TFile::Open("AcceptanceEfficiencyMap.root","UPDATE");
  for(int centb=0;centb<=_NcentBins;centb++){
    if(ispp && centb>0) break;

    if(!runMCclos){
      hp_efficiency[centb]->Write("hp_efficiency_"+(TString)((centb==0)?"":("centBin"+(TString)to_string(centb)+"_"))+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
      hp_efficiency_noSF[centb]->Write("hp_efficiency_noSF_"+(TString)((centb==0)?"":("centBin"+(TString)to_string(centb)+"_"))+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
      hp_efficiency_simpleAverage[centb]->Write("hp_efficiency_simpleAverage_"+(TString)((centb==0)?"":("centBin"+(TString)to_string(centb)+"_"))+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
      hp_efficiency_selectiveSFapplication[centb]->Write("hp_efficiency_selectiveSFapplication_"+(TString)((centb==0)?"":("centBin"+(TString)to_string(centb)+"_"))+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));

      hp_acceff[centb]->Write("hp_acceff_"+(TString)((centb==0)?"":("centBin"+(TString)to_string(centb)+"_"))+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
      hpcoarse_inBDT23[centb]->Write("hpcoarse_inBDT23_"+(TString)((centb==0)?"":("centBin"+(TString)to_string(centb)+"_"))+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(integratePtBins?"_integratePtBins":"")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
      hpcoarse_inBDT3[centb]->Write("hpcoarse_inBDT3_"+(TString)((centb==0)?"":("centBin"+(TString)to_string(centb)+"_"))+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(integratePtBins?"_integratePtBins":"")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
    }

    fout->WriteObject(&eff_oneBinned,"efficiency_oneBinned"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
    if(!ispp) fout->WriteObject(&eff_oneBinned_centDep,"efficiency_oneBinned_centDep"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));
    if(runAEtoys) {
      fout->WriteObject(&eff_oneB_biased,"efficiency_oneBinned_biased"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?"_2ndStep":""));
      if(!ispp) fout->WriteObject(&eff_oneB_centDep_biased,"efficiency_oneBinned_centDep_biased"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?"_2ndStep":""));
    }
  }

  if(runMCclos){
    fout->WriteObject(&eff_MCclos,"efficiency_MCclosure"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));    
    fout->WriteObject(&passing_MCclos,"Nobs_MCclosure"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep":"_2ndStep")):""));    
  }
  fout->Close();
}
