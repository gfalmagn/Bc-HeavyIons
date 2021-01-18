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
#include "TPaletteAxis.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include <TVirtualFitter.h>
#include "../PbPb18/Utilities/EVENTUTILS.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"
#include "../acceptance/SgMuonAcceptanceCuts.h"
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

void drawEffMap(TH2Poly* hpEff, TH2Poly* hpSel, TH2Poly* hpAcc, TLine* l1, TLine* l2, TLine* l3, TLine* l4, TLine* l5, bool ispp, TString nameSuf, bool makeEffDiv, float zmin, float zmax){

  TCanvas *c2 = new TCanvas("c2","c2",3000,1500);
  c2->Divide(2,1);

  c2->cd(1);
  //  gPad->SetLogz();
  //hpSel->GetZaxis()->SetRangeUser(1,200);
  hpSel->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  hpSel->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  hpSel->SetTitle("Selected B_{c}'s");
  hpSel->Draw("COLZ0");
  l1->Draw("same"); l2->Draw("same"); l3->Draw("same"); l4->Draw("same"); l5->Draw("same");

  c2->cd(2)->SetRightMargin(0.15);
  //gPad->SetLogz();
  if(makeEffDiv) hpEff->Divide(hpAcc);
  else hpEff->Multiply(hpAcc);
  hpEff->GetZaxis()->SetRangeUser(zmin,zmax);
  hpEff->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  hpEff->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  hpEff->SetTitle("Efficiency "+nameSuf);
  hpEff->Draw("COLZ0");
  l1->Draw("same"); l2->Draw("same"); l3->Draw("same"); l4->Draw("same"); l5->Draw("same");
  //  cout<<hpEff->GetBinContent(hpEff->FindBin(1.3001,6.001))<<" "<<hpEff->GetBinContent(hpEff->FindBin(1.301,8.01))<<" "<<hpEff->GetBinContent(hpEff->FindBin(1.301,10.01))<<endl;

  c2->SaveAs("figs/EfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+nameSuf+(TString)(ispp?"_pp":"_PbPb")+".pdf");
  c2->SaveAs("figs/EfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+nameSuf+(TString)(ispp?"_pp":"_PbPb")+".png");

}

void BuildEffMap(bool ispp = true, bool BDTuncorrFromM=false, bool integratePtBins=false, bool runAEtoys=true){

  //**************************************************************  
  //Grab the variations of the pT bias of MC, from first step r1 and r2
  vector<TH1F*> bias;
  if(runAEtoys==false){
    TFile *BiasFile = TFile::Open("../twoSteps/pTBiases.root","READ");
    for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
      bias.push_back((TH1F*)BiasFile->Get("pTbias_"+(TString)(ispp?"pp":"PbPb")+"_var"+(TString)to_string(v)));
      //cout<<"variation value_at_11 = "<<v<<" "<<bias[v]->Eval(11)<<endl;
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

  //**************************************************************
  //Some variables and XS corrections
  std::map<bool, float> scaleMCsig = {{ true,  304800 * 2.54e0 * 0.668 / 3000000}, // Lumi_pp[nb-1] (from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM) * (XS_Bc_pp * BF((J/psi -> mu mu) mu nu))[nb] * (XS(5.02 TeV) / XS(7 TeV)) / nevents(uncut MC sample)
				      { false, 1.0 * 1.71641 * 2.54e0 * 0.668 * (7461 / 67.6) / 1960000} }; //Lumi_PbPb[nb-1] (from https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/948.html) * (XS_Bc_pp * BF((J/psi -> mu mu) mu nu))[nb] * (XS(5.02 TeV) / XS(7 TeV)) * ( XS^geom_PbPb / XS_Nucleon-Nucleon ) / nevents(uncut MC sample)
                        //weighted by Ncoll(centrality of given event) later, with an average Ncoll_MB = 392. The value of XS^geom was set to A^2 * XS_NN / Ncoll_MB, where XS_NN is taken from Glauber MC d'Enterria PRC 97.054910
                        //Assuming R_AA(Bc)=1

  //counters, histos
  float ntot = 0, nacc = 0, nsel = 0, nsel2 = 0, nBDT23 = 0, nBDT3 = 0;
  TH2Poly *hp = _hp();
  TH2Poly *hp_coarser = _hp_coarser();
  TH2Poly *hp_all = (TH2Poly*) hp->Clone("hp_all");
  TH2Poly *hp_acc = (TH2Poly*) hp->Clone("hp_acc");
  TH2Poly *hp_sel = (TH2Poly*) hp->Clone("hp_sel");
  TH2Poly *hp_sel_noSF = (TH2Poly*) hp->Clone("hp_sel_noSF");
  TH2Poly *hp_sel_selectiveSFapplication = (TH2Poly*) hp->Clone("hp_sel_selectiveSFapplication");
  TH2Poly *hp_sel_simpleAverage = (TH2Poly*) hp->Clone("hp_sel_simpleAverage");
  TH2Poly *hpcoarse_sel = (TH2Poly*) hp_coarser->Clone("hpcoarse_sel");
  TH2Poly *hpcoarse_inBDT23 = (TH2Poly*) hp_coarser->Clone("hpcoarse_inBDT23");
  TH2Poly *hpcoarse_inBDT3 = (TH2Poly*) hp_coarser->Clone("hpcoarse_inBDT3");
  TH1F *BDT23effVsPt = new TH1F("BDT23effVsPt","BDT23effVsPt",22,6,50);
  TH1F *BDT3effVsPt = new TH1F("BDT3effVsPt","BDT3effVsPt",22,6,50);
  TH1F *selectedVsPt = new TH1F("selectedVsPt","selectedVsPt",22,6,50);
  TH1F *SFs = new TH1F("SFs","SFs",100,ispp?0.86:0.86,ispp?1.2:1.4);
  TH1F *SFs_selectiveSFapplication = new TH1F("SFs_selectiveSFapplication","SFs_selectiveSFapplication",100,ispp?0.86:0.86,ispp?1.2:1.4);
  TH1F *SFs_simpleAverage = new TH1F("SFs_simpleAverage","SFs_simpleAverage",100,ispp?0.86:0.86,ispp?1.2:1.4);
  vector<vector<float> > passing_oneBinned(_NanaBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<vector<float> > eff_oneBinned(_NanaBins+1,vector<float>(4,0)); //nominal, statErr, systErr, totErr
  vector<float> accepted_oneBinned(_NanaBins+1,0);
  vector<vector<vector<float> > > passing_oneB(_NanaBins+1,vector<vector<float> >(4,vector<float>(5,0))); //4 is: muid (or glb), trk (in pp), trg (or muidtrg), tot
  //MC with biased pT distributions
  vector<vector<float> > eff_oneB_biased(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 1));
  vector<vector<float> > accepted_oneB_biased(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0));
  vector<vector<float> > passing_oneB_biased(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0));


  //Jpsi mass histo for JpsiChoiceWeight
  vector<TH1F*> JpsiM(_nChan(ispp)+1);
  vector<TH1F*> JpsiM_tight(_nChan(ispp)+1);
  TFile *f_JpsiM = TFile::Open("../BDT/JpsiMassDistr.root","READ");
  for(int k=0;k<=_nChan(ispp);k++){
    JpsiM[k] = (TH1F*)f_JpsiM->Get("JpsiMass_data_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)))+(TString)(ispp?"_pp":"_PbPb"));
    JpsiM_tight[k] = (TH1F*)f_JpsiM->Get("JpsiMassCentralEta_data_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)))+(TString)(ispp?"_pp":"_PbPb"));
  }
  

  //**************************************************************
  //loop on events
  for(int j=0; j<nentries; j++){//nentries
    if(j%100000==0){ cout<<"Scanned "<<100.*(double)j/nentries<<"% of entries"<<endl; }
    
    Reco_3mu_4mom->Clear();
    Reco_mu_4mom->Clear();
    Reco_QQ_4mom->Clear();
    Gen_3mu_4mom->Clear();
    Gen_mu_4mom->Clear();
    Gen_QQ_4mom->Clear();
    
    T->GetEntry(j);
    
    for(int igen=0;igen<Gen_Bc_size;igen++){

      float weight = scaleMCsig[ispp];
      if(!ispp) {float weightNcoll = (float)findNcoll(Centrality); //for PbPb MC
	weight *= weightNcoll;}
      
      ntot += weight;      
      TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom->At(igen);
      hp_all->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight);
      
      int genQQidx = Gen_Bc_QQ_idx[igen];
      TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom->At(genQQidx);
      TLorentzVector *genBc_muW = (TLorentzVector*) Gen_mu_4mom->At(Gen_Bc_muW_idx[igen]);
      TLorentzVector *genBc_mumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[genQQidx]);
      TLorentzVector *genBc_mupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[genQQidx]);
      
      if(!InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,_withTM)) continue;
      for(int b=1;b<=_NanaBins;b++){
	if(fabs(gen3mu->Rapidity())>_BcYmin[b] && fabs(gen3mu->Rapidity())<_BcYmax[b] && gen3mu->Pt()>_BcPtmin[b] && gen3mu->Pt()<_BcPtmax[b] ){
	  accepted_oneBinned[b] += weight;
	  accepted_oneBinned[0] += weight;

	  //biased MC
	  if(runAEtoys){
	    for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	      accepted_oneB_biased[0][v] += getBias( bias[v] , gen3mu->Pt()) * weight;
	      accepted_oneB_biased[b][v] += getBias( bias[v] , gen3mu->Pt()) * weight;
	    }
	  }

	}//end accepted in good ana bin
      }

      hp_acc->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight);
      
      int irec = Gen_3mu_whichRec[igen];

      if(irec<0) continue;
      int QQidx_2[2] = { Reco_3mu_QQ1_idx[irec], Reco_3mu_QQ2_idx[irec] };
      int muWidx_2[2] = { Reco_3mu_muW_idx[irec], Reco_3mu_muW2_idx[irec] };
      int mumiidx_2[2] = { Reco_3mu_mumi_idx[irec], Reco_QQ_mumi_idx[QQidx_2[1]] };
      int muplidx_2[2] = { Reco_3mu_mupl_idx[irec], Reco_QQ_mupl_idx[QQidx_2[1]] };
      
      for(int k=0; k<2; k++){
	int QQidx = QQidx_2[k];
	int QQ2idx = QQidx_2[(k==0)?1:0];
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
	  if(inJpsiMassSB(QQM, maxEta<1.5)) {weight *= -1;}
	  else if(!(inJpsiMassRange(QQM, maxEta<1.5))){ weight = 0;}
	    
	  if(!goodTree || weight == 0) continue;
	  
	  //**************************************************************
	  //Write the wanted variables into the chosen (goodTree=true) output tree
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

	  //cout<<(Bc_ctauSignif>1.5)<<" "<<(Bc_alpha<0.8)<<" "<<(Bc_alpha3D<0.8)<<" "<<(Bc_VtxProb>0.01)<<" "<<(QQ_dca<0.3)<<" "<<(muW_isSoft)<<" "<<( mumi_isSoft)<<" "<<(mupl_isSoft)<<" "<<( (muW_isGlb && muW_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) || (muW_isGlb && muW_inLooseAcc && mumi_isGlb && mumi_inLooseAcc) || (mumi_isGlb && mumi_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) )<<" "<<((muW_trig && mupl_trig && muW_inTightAcc && mupl_inTightAcc ) ||(muW_trig && mumi_trig && muW_inTightAcc && mumi_inTightAcc ) || (mumi_trig && mupl_trig && mumi_inTightAcc && mupl_inTightAcc ))<<" "<<(fabs(muW_dz)<(ispp?0.6:0.8) && fabs(mumi_dz)<(ispp?0.6:0.8) && fabs(mupl_dz)<(ispp?0.6:0.8))<<" "<<((HLTriggers&((ispp || i>=4)?8:4096))>0)<<endl;
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

            float PperpTrimu = TMath::Sin(Bc_alpha3D) * recBc->P();
            float Bc_CorrM = sqrt(BcCandM*BcCandM + PperpTrimu*PperpTrimu) + PperpTrimu;
            if(Bc_CorrM>_BcCorrM_cut(ispp)) continue;
	    if(!ispp && Centrality>180) continue; //keep 0-90% centrality                                                                                                                                                              

	    if(muW_trig && !(muW_L2 || muW_L3)) cout<<"!!!! PROBLEM with muW trigger consistency"<<endl;
	    if(mupl_trig && !(mupl_L2 || mupl_L3)) cout<<"!!!! PROBLEM with mupl trigger consistency"<<endl;
	    if(mumi_trig && !(mumi_L2 || mumi_L3)) cout<<"!!!! PROBLEM with mumi trigger consistency"<<endl;

	    float QQ2_M = (Reco_mu_charge[muWidx]>0)?((*recBc_mumi+*recBc_muW).M()):((*recBc_mupl+*recBc_muW).M()); //QQ2 is the second OS pair
	    float QQ3_M = (Reco_mu_charge[muWidx]>0)?((*recBc_mupl+*recBc_muW).M()):((*recBc_mumi+*recBc_muW).M()); //QQ3 is the SS pair
		
	    //**** Deal with the Jpsi dimuon choice
	    if((inJpsiMassRange(QQM, maxEta<1.5) && inJpsiMassSB(QQ2_M, maxEta<1.5)
		&& Reco_QQ_VtxProb[QQ2idx]>_QQvtxProb_cut && Reco_QQ_dca[QQ2idx]<_QQdca_cut && Reco_QQ_dca[QQ2idx]>0)
	       || (inJpsiMassSB(QQM, maxEta<1.5) && inJpsiMassRange(QQ2_M, maxEta<1.5)
		   && Reco_QQ_VtxProb[QQ2idx]>_QQvtxProb_cut && Reco_QQ_dca[QQ2idx]<_QQdca_cut && Reco_QQ_dca[QQ2idx]>0)
	       ){
	      
	      float binc_QQ1 = ( (maxEta<1.5)?JpsiM_tight:JpsiM )[0]->GetBinContent(( (maxEta<1.5)?JpsiM_tight:JpsiM )[0]->FindBin(QQM));
	      float binc_QQ2 = ( (maxEta<1.5)?JpsiM_tight:JpsiM )[0]->GetBinContent(( (maxEta<1.5)?JpsiM_tight:JpsiM )[0]->FindBin(QQ2_M));

	      if (binc_QQ1==0) weight = 0;
	      else weight *= binc_QQ1 / (binc_QQ1+binc_QQ2);
	    }


	    //**** Apply the scale factors to correct MC efficiencies
	    //H = tight criteria (tight acceptance, triggering muon...), L = loose criteria (loose acceptance, no trigger required...), T = fire L3 trigger (in PbPb case) -> implies H
	    //2 muons have to respect H, 1 has to respect L
	    //In PbPb, at least 1 triggering muon must pass L3
	    //muon 1 is mumi, muon 2 is mupl, muon3 is muW
	    double effcorr; //effcorr is the correction to be applied
	    double effcorr_simpleAverage, effcorr_selectiveSFapplication; 
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
		  if(fabs(gen3mu->Rapidity())>_BcYmin[b] && fabs(gen3mu->Rapidity())<_BcYmax[b] && gen3mu->Pt()>_BcPtmin[b] && gen3mu->Pt()<_BcPtmax[b] ){
		    passing_oneB[b][sftype][id] += weight * effcorr;
		    passing_oneB[0][sftype][id] += weight * effcorr;

		    if(er==0 && sftype==2 && runAEtoys){//nominal concerning SF variation
		      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
			passing_oneB_biased[0][v] += getBias( bias[v] , gen3mu->Pt()) * weight * effcorr;
			passing_oneB_biased[b][v] += getBias( bias[v] , gen3mu->Pt()) * weight * effcorr;
		      }
		    }

		  }
		}

	      }//end loop on Error index
	    }//end loop on sf error type

	    //last ran variation is the nominal (er=0), so effcorr has the nominal value
	    hp_sel_noSF->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight);
	    hp_sel->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight * effcorr);
	    hp_sel_simpleAverage->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight * effcorr_simpleAverage);
	    hp_sel_selectiveSFapplication->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight * effcorr_selectiveSFapplication);

	    SFs->Fill(effcorr, weight);
	    SFs_simpleAverage->Fill(effcorr_simpleAverage, weight);
	    SFs_selectiveSFapplication->Fill(effcorr_selectiveSFapplication, weight);

	  } //end if passes full selection
	      
	} //end if(QQ_isValid)
      } //end loop on 2/3 possible Jpsi dimuon choice
    } //end loop on Bc candidates
  } //end loop on entries

  //Calculate errors and systematics
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
    //MC with biased pT
    if(runAEtoys){
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	eff_oneB_biased[b][v] = passing_oneB_biased[b][v]/accepted_oneB_biased[b][v];
      }
    }

    if(b>0){
      // cout<<"b nom stat syst tot = "<<b<<" "<<passing_oneBinned[b][0]<<" "<<passing_oneBinned[b][1]<<" "<<passing_oneBinned[b][2]<<" "<<passing_oneBinned[b][3]<<endl;
      // cout<<"id stat syst = "<<idStatErr<<" "<<idSystErr<<endl;
      // cout<<"trk stat syst = "<<trkStatErr<<" "<<trkSystErr<<endl;
      // cout<<"trg stat syst = "<<trgStatErr<<" "<<trgSystErr<<endl;
      cout<<"accepted_oneBinned, passing_oneBinned, efficiency = " <<ntot<<" "<<accepted_oneBinned[b]<<" "<<passing_oneBinned[b][0]<<" "<<eff_oneBinned[b][0]<<" "<<endl;
      cout<<"err passing_oneBinned, relerr efficiency = " <<passing_oneBinned[b][3]<<" "<<eff_oneBinned[b][3]/eff_oneBinned[b][0]<<" "<<endl;
    }
  }

  cout<<"ntot, accepted_oneBinned[0], passing_oneBinned[0], efficiency[0] = " <<ntot<<" "<<accepted_oneBinned[0]<<" "<<passing_oneBinned[0][0]<<" "<<eff_oneBinned[0][0]<<" "<<endl;
  cout<<"err passing_oneBinned[0], relerr efficiency[0] = " <<passing_oneBinned[0][3]<<" "<<eff_oneBinned[0][3]/eff_oneBinned[0][0]<<" "<<endl;


  //**************************************************************
  //Pre-selected tree to get BDT efficiency map
  auto preselFile = TFile::Open("../BDT/BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  TTree* T_presel = (TTree*)preselFile->Get("signal_MC");

  float Bc_Pt; TBranch *b_Bc_Pt = T_presel->GetBranch("Bc_Pt");
  b_Bc_Pt->SetAddress(&Bc_Pt);
  float Bc_Y; TBranch *b_Bc_Y = T_presel->GetBranch("Bc_Y");
  b_Bc_Y->SetAddress(&Bc_Y);
  float weight; TBranch *b_weight = T_presel->GetBranch("weight");
  b_weight->SetAddress(&weight);
  float BDT; TBranch *b_BDT = T_presel->GetBranch("BDT");
  b_BDT->SetAddress(&BDT);

  for(int j=0; j<T_presel->GetEntries(); j++){//T_presel->GetEntries()
    T_presel->GetEntry(j);
    hpcoarse_sel->Fill(fabs(Bc_Y),Bc_Pt, weight);
    selectedVsPt->Fill(Bc_Pt, weight);
    nsel2 += weight;
    int kinb = 0;
    if(!integratePtBins) kinb = (Bc_Pt<_BcPtmax[1])?1:2;
    if(BDT>_BDTcuts(ispp,kinb,BDTuncorrFromM)[1]){ //should be adapted, when changing to centrality binning
      nBDT23 += weight;
      BDT23effVsPt->Fill(Bc_Pt, weight);
      hpcoarse_inBDT23->Fill(fabs(Bc_Y),Bc_Pt, weight);}
    if(BDT>_BDTcuts(ispp,kinb,BDTuncorrFromM)[2]){
      nBDT3 += weight;
      BDT3effVsPt->Fill(Bc_Pt, weight);
      hpcoarse_inBDT3->Fill(fabs(Bc_Y),Bc_Pt, weight);}
  }

  BDT23effVsPt->Divide(selectedVsPt);
  BDT3effVsPt->Divide(selectedVsPt);
  hpcoarse_inBDT23->Divide(hpcoarse_sel);
  hpcoarse_inBDT3->Divide(hpcoarse_sel);
  float eff_BDT23 = (float)nBDT23/(float)nsel2;
  float eff_BDT3 = (float)nBDT3/(float)nsel2;

  //**************************************************************
  //Lines for fiducial cuts 
  TLine *line1 = new TLine(_BcYmin[0],_BcPtmax[1],_BcYmin[1],_BcPtmax[1]);
  TLine *line2 = new TLine(_BcYmin[1],_BcPtmin[1],_BcYmin[1],_BcPtmax[1]);
  TLine *line3 = new TLine(_BcYmin[1],_BcPtmin[1],_BcYmax[1],_BcPtmin[1]);
  TLine *line4 = new TLine(_BcYmax[1],_BcPtmin[0],_BcYmax[1],_BcPtmax[0]);
  TLine *line5 = new TLine(_BcYmin[1],_BcPtmax[1],_BcYmax[1],_BcPtmax[1]);
  line1->SetLineWidth(4);  line1->SetLineColor(kBlack);
  line2->SetLineWidth(4);  line2->SetLineColor(kBlack);
  line3->SetLineWidth(4);  line3->SetLineColor(kBlack);
  line4->SetLineWidth(4);  line4->SetLineColor(kBlack);
  line5->SetLineWidth(3); line5->SetLineStyle(2); line5->SetLineColor(kBlack);

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(50);

  //**************************************************************
  //Draw TH2Poly for efficiencies

  TH2Poly* hp_efficiency = (TH2Poly*)hp_sel->Clone("hp_efficiency");
  if(!integratePtBins && !BDTuncorrFromM) 
    drawEffMap(hp_efficiency, hp_sel, hp_acc, line1, line2, line3, line4, line5, ispp, "", true,0,0.45);

  TH2Poly* hp_efficiency_noSF = (TH2Poly*)hp_efficiency->Clone("hp_efficiency_noSF");
  hp_efficiency_noSF->Divide(hp_sel_noSF);
 if(!integratePtBins && !BDTuncorrFromM) 
   drawEffMap(hp_efficiency_noSF, hp_sel_noSF, hp_acc, line1, line2, line3, line4, line5, ispp, "_noSF_DivideNominal", false,ispp?0.84:0.83,ispp?1.16:1.17);

  TH2Poly* hp_efficiency_selectiveSFapplication = (TH2Poly*)hp_efficiency->Clone("hp_efficiency_selectiveSFapplication");
  hp_efficiency_selectiveSFapplication->Divide(hp_sel_selectiveSFapplication);
  if(!integratePtBins && !BDTuncorrFromM) 
    drawEffMap(hp_efficiency_selectiveSFapplication, hp_sel_selectiveSFapplication, hp_acc, line1, line2, line3, line4, line5, ispp, "_selectiveSFapplication_DivideNominal", false,ispp?0.92:0.89,ispp?1.08:1.11);

  TH2Poly* hp_efficiency_simpleAverage = (TH2Poly*)hp_efficiency->Clone("hp_efficiency_simpleAverage");
  hp_efficiency_simpleAverage->Divide(hp_sel_simpleAverage);
  if(!integratePtBins && !BDTuncorrFromM) 
    drawEffMap(hp_efficiency_simpleAverage, hp_sel_simpleAverage, hp_acc, line1, line2, line3, line4, line5, ispp, "_simpleAverage_DivideNominal", false,ispp?0.9:0.84,ispp?1.1:1.16);

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

    c6->SaveAs("figs/ScaleFactorsComparison"+(TString)(ispp?"_pp":"_PbPb")+".pdf");
  }
  
  //**************************************************************
  //Grab acceptance map
  TFile acc_file("../acceptance/acceptanceMap.root","READ");
  TH2Poly* hp_acceptance = (TH2Poly*)acc_file.Get("hp_acceptance");
  hp_acceptance->SetDirectory(0);
  TH2Poly* hp_acceff = (TH2Poly*)hp_acceptance->Clone("hp_acceff");
  hp_acceff->SetDirectory(0);
  hp_acceff->Multiply(hp_efficiency);

  if(!integratePtBins && !BDTuncorrFromM){
    TCanvas *c3 = new TCanvas("c3","c3",2500,2500);
    c3->Divide(2,2);
  
    c3->cd(1);
    hp_acceptance->Draw("COLZ");
    line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same");
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    hp_acceptance->GetZaxis()->SetRangeUser(5e-5,1);

    TPaletteAxis *palette = (TPaletteAxis*)hp_acceptance->GetListOfFunctions()->FindObject("palette");
    // the following lines moe the paletter. Choose the values you need for the position.
    palette->SetX1NDC(0.86);
    palette->SetX2NDC(0.91);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();
  
    c3->cd(2);
    hp_efficiency->Draw("COLZ0");
    gPad->SetRightMargin(0.15);
    line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same");
    c3->cd(3);
    hp_sel->Draw("COLZ0");
    line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same");

    c3->cd(4);
    hp_acceff->Draw("COLZ0");
    line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same"); line5->Draw("same");
    hp_acceff->SetTitle("Acceptance #times  Efficiency");
    cout<<"smallest in fid cuts acc eff = "<<hp_acceff->GetBinContent(hp_acceff->FindBin(1.32,6.05))<<endl;
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    hp_acceff->GetZaxis()->SetRangeUser(3e-5,0.45);

    TPaletteAxis *palette2 = (TPaletteAxis*)hp_acceff->GetListOfFunctions()->FindObject("palette");
    // the following lines moe the paletter. Choose the values you need for the position.
    palette2->SetX1NDC(0.86);
    palette2->SetX2NDC(0.91);
    palette2->SetY1NDC(0.1);
    palette2->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();

    c3->SaveAs("figs/AcceptanceEfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".pdf");
    c3->SaveAs("figs/AcceptanceEfficiencyMap_tunedBins"+(TString)(_withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".png");

    //**************************************************************
    //Draw TH2Poly for BDT efficiency
    TCanvas *c4 = new TCanvas("c4","c4",1500,1500);
    c4->Divide(2,2);

    c4->cd(1);
    hpcoarse_inBDT23->GetZaxis()->SetRangeUser(eff_BDT23-0.2, eff_BDT23+0.2);
    hpcoarse_inBDT23->SetTitle("Efficiency for BDT bin 2-3;|Y|;p_{T} [GeV]");
    hpcoarse_inBDT23->Draw("COLZ");

    c4->cd(2);
    hpcoarse_inBDT3->GetZaxis()->SetRangeUser(eff_BDT3-0.2, eff_BDT3+0.2);
    hpcoarse_inBDT3->SetTitle("Efficiency for BDT bin 3;|Y|;p_{T} [GeV]");
    hpcoarse_inBDT3->Draw("COLZ");

    c4->cd(3);
    BDT23effVsPt->GetYaxis()->SetRangeUser(eff_BDT23-0.2, eff_BDT23+0.2);
    BDT23effVsPt->SetTitle("Efficiency for BDT bin 2-3;p_{T} [GeV]");
    BDT23effVsPt->Draw("E");

    c4->cd(4);
    BDT3effVsPt->GetYaxis()->SetRangeUser(eff_BDT3-0.2, eff_BDT3+0.2);
    BDT3effVsPt->SetTitle("Efficiency for BDT bin 3;p_{T} [GeV]");
    BDT3effVsPt->Draw("E");

    c4->SaveAs("figs/BDTefficiencyMap"+(TString)(_withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".pdf");
  }

  //**************************************************************
  //Store maps
  TFile* fout = TFile::Open("AcceptanceEfficiencyMap.root","UPDATE");
  hp_efficiency->Write("hp_efficiency_"+(TString)(ispp?"pp":"PbPb"));
  hp_efficiency_noSF->Write("hp_efficiency_noSF_"+(TString)(ispp?"pp":"PbPb"));
  hp_efficiency_simpleAverage->Write("hp_efficiency_simpleAverage_"+(TString)(ispp?"pp":"PbPb"));
  hp_efficiency_selectiveSFapplication->Write("hp_efficiency_selectiveSFapplication_"+(TString)(ispp?"pp":"PbPb"));

  hp_acceff->Write("hp_acceff_"+(TString)(ispp?"pp":"PbPb"));
  hpcoarse_inBDT23->Write("hpcoarse_inBDT23_"+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(integratePtBins?"_integratePtBins":""));
  hpcoarse_inBDT3->Write("hpcoarse_inBDT3_"+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(integratePtBins?"_integratePtBins":""));
  fout->WriteObject(&eff_oneBinned,"efficiency_oneBinned"+(TString)(ispp?"_pp":"_PbPb"));
  if(runAEtoys) fout->WriteObject(&eff_oneB_biased,"efficiency_oneBinned_biased"+(TString)(ispp?"_pp":"_PbPb"));
  fout->Close();
}
