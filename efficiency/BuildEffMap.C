#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include <TVirtualFitter.h>
#include "../PbPb18/Utilities/EVENTUTILS.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../acceptance/SgMuonAcceptanceCuts.h"
#include "../helpers/AccEff2DBinning.h"

void BuildEffMap(bool ispp = true, bool BDTuncorrFromM=false, bool integratePtBins=false){
  bool withTM = false;

  //**************************************************************  
  //Create Tree and branches
  TFile *fileMC = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/Oniatree_MC_Bc_trimuons_21112019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/BcToJpsiMuNu_BCVEGPY_PYTHIA8_2018PbPb5TeV_09012020_1948k_ONIATREE.root");
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
  float ntot = 0, nfid = 0, nacc = 0, nsel = 0, nsel2 = 0, nBDT23 = 0, nBDT3 = 0;
  TH2Poly *hp = _hp();
  TH2Poly *hp_coarser = _hp_coarser();
  TH2Poly *hp_all = (TH2Poly*) hp->Clone("hp_all");
  TH2Poly *hp_acc = (TH2Poly*) hp->Clone("hp_acc");
  TH2Poly *hp_sel = (TH2Poly*) hp->Clone("hp_sel");
  TH2Poly *hpcoarse_sel = (TH2Poly*) hp_coarser->Clone("hpcoarse_sel");
  TH2Poly *hpcoarse_inBDT23 = (TH2Poly*) hp_coarser->Clone("hpcoarse_inBDT23");
  TH2Poly *hpcoarse_inBDT3 = (TH2Poly*) hp_coarser->Clone("hpcoarse_inBDT3");
  TH1F *BDT23effVsPt = new TH1F("BDT23effVsPt","BDT23effVsPt",22,6,50);
  TH1F *BDT3effVsPt = new TH1F("BDT3effVsPt","BDT3effVsPt",22,6,50);
  TH1F *selectedVsPt = new TH1F("selectecVsPt","selectedVsPt",22,6,50);
  vector<float> passing_oneBinned(_NanaBins+1,0);
  vector<float> accepted_oneBinned(_NanaBins+1,0);
  vector<float> eff_oneBinned(_NanaBins+1,1);

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
      
      for(int b=1;b<=_NanaBins;b++)
	if(fabs(gen3mu->Rapidity())>_BcYmin[b] && fabs(gen3mu->Rapidity())<_BcYmax[b] && gen3mu->Pt()>_BcPtmin[b] && gen3mu->Pt()<_BcPtmax[b] ) nfid += weight;

      if(!InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM)) continue;
      for(int b=1;b<=_NanaBins;b++){
	if(fabs(gen3mu->Rapidity())>_BcYmin[b] && fabs(gen3mu->Rapidity())<_BcYmax[b] && gen3mu->Pt()>_BcPtmin[b] && gen3mu->Pt()<_BcPtmax[b] ){
	  accepted_oneBinned[b] += weight;
	  accepted_oneBinned[0] += weight;
	}
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
	  
	  bool goodTree = fabs(Reco_3mu_charge[irec])==1 && Reco_QQ_sign[QQidx]==0 && inLooseMassRange(QQM) // in Jpsi mass region
	    && (BcCandM < m_Bc + 1.0) && (BcCandM > 3.5) // in Bc mass region
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
	  bool muW_isSoft = isSoft( false, muW_isGlb , (Reco_mu_SelType[muWidx]&8)>0 , ((Reco_mu_SelType[muWidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[muWidx] , Reco_mu_dxy[muWidx] , Reco_mu_dz[muWidx] , Reco_mu_nPixWMea[muWidx] , Reco_mu_nTrkWMea[muWidx] );

	  bool mumi_inLooseAcc = Reco_mu_InLooseAcc[mumiidx];
	  bool mumi_inTightAcc = Reco_mu_InTightAcc[mumiidx];
	  bool mumi_isGlb = (Reco_mu_SelType[mumiidx]&2)>0;
	  bool mumi_trig = (Reco_mu_trig[mumiidx]&(ispp?8:4096))>0;
	  bool mumi_isSoft = isSoft( false, mumi_isGlb , (Reco_mu_SelType[mumiidx]&8)>0  , ((Reco_mu_SelType[mumiidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[mumiidx] , Reco_mu_dxy[mumiidx] , Reco_mu_dz[mumiidx] , Reco_mu_nPixWMea[mumiidx] , Reco_mu_nTrkWMea[mumiidx] );

	  bool mupl_inLooseAcc = Reco_mu_InLooseAcc[muplidx];
	  bool mupl_inTightAcc = Reco_mu_InTightAcc[muplidx];
	  bool mupl_isGlb = (Reco_mu_SelType[muplidx]&2)>0;
	  bool mupl_trig = (Reco_mu_trig[muplidx]&(ispp?8:4096))>0;
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
	     && Bc_alpha<_alpha_cut
	     && Bc_alpha3D<_alpha3D_cut //combined with alpha cut, kills 0.3% of signal 
	     && Bc_VtxProb>_vtxProb_cut //kills 1.2% more signal compared to Bc_VtxProb>0.005
	     && QQ_VtxProb>_QQvtxProb_cut //drop the QQ_VtxProb, too correlated with Bc_VtxProb?
	     && QQ_dca<_QQdca_cut //keep this that kills 1.8% of signal and 10% of WRONSIGN/BCMASS
	     && muW_isSoft
	     && mumi_isSoft
	     && mupl_isSoft
	     && (withTM?( (muW_isGlb && muW_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) || //only one muon can be tracker and out of LooseAcceptance
			  (muW_isGlb && muW_inLooseAcc && mumi_isGlb && mumi_inLooseAcc) || 
			  (mumi_isGlb && mumi_inLooseAcc && mupl_isGlb && mupl_inLooseAcc)
			  ):(
			     mumi_isGlb && mupl_isGlb && muW_isGlb
			     ))
	     && ( ( muW_trig && mupl_trig && muW_inTightAcc && mupl_inTightAcc ) || //two muons among three must trigger //BEWARE ! Not sure if TightAcceptance should be put there
		  ( muW_trig && mumi_trig && muW_inTightAcc && mumi_inTightAcc ) ||
		  ( mumi_trig && mupl_trig && mumi_inTightAcc && mupl_inTightAcc ) //only this last option can be true for dimuon+trk
		  )
	     && fabs(muW_dz)<(ispp?0.6:0.8) && fabs(mumi_dz)<(ispp?0.6:0.8) && fabs(mupl_dz)<(ispp?0.6:0.8)
	     && (!ispp || (HLTriggers&(ispp?8:4096))>0) //the event must fire the trigger as well
	     ){

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

	    for(int b=1;b<=_NanaBins;b++){
	      if(fabs(gen3mu->Rapidity())>_BcYmin[b] && fabs(gen3mu->Rapidity())<_BcYmax[b] && gen3mu->Pt()>_BcPtmin[b] && gen3mu->Pt()<_BcPtmax[b] ){
		passing_oneBinned[b] += weight;
		passing_oneBinned[0] += weight;
	      }
	    }

	    hp_sel->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt(), weight);

	  } //end if passes full selection
	      
	} //end if(QQ_isValid)
      } //end loop on 2/3 possible Jpsi dimuon choice
    } //end loop on Bc candidates
  } //end loop on entries

  cout<<"ntot, accepted_oneBinned[0], passing_oneBinned[0], efficiency = " <<ntot<<" "<<accepted_oneBinned[0]<<" "<<passing_oneBinned[0]<<" "<<passing_oneBinned[0]/accepted_oneBinned[0]<<" "<<endl;


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
  //One-binned efficiency
  for(int b=1;b<=_NanaBins;b++)
    eff_oneBinned[b] = passing_oneBinned[b]/accepted_oneBinned[b];

  //**************************************************************
  //Lines for fiducial cuts 
  TLine *line1 = new TLine(_BcYmin[1],_BcPtmax[2],_BcYmin[2],_BcPtmax[2]);
  TLine *line2 = new TLine(_BcYmin[2],_BcPtmin[1],_BcYmin[2],_BcPtmax[2]);
  TLine *line3 = new TLine(_BcYmin[2],_BcPtmin[1],_BcYmax[1],_BcPtmin[1]);
  TLine *line4 = new TLine(_BcYmax[1],_BcPtmin[1],_BcYmax[1],_BcPtmax[1]);
  line1->SetLineWidth(4);  line1->SetLineColor(kBlack);
  line2->SetLineWidth(4);  line2->SetLineColor(kBlack);
  line3->SetLineWidth(4);  line3->SetLineColor(kBlack);
  line4->SetLineWidth(4);  line4->SetLineColor(kBlack);

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(50);

  //**************************************************************
  //Draw TH2Poly 

  TCanvas *c2 = new TCanvas("c2","c2",3000,1500);
  c2->Divide(2,1);

  TH2Poly* hp_efficiency = (TH2Poly*)hp_sel->Clone("hp_efficiency");
  c2->cd(1);
  //  gPad->SetLogz();
  //hp_sel->GetZaxis()->SetRangeUser(1,200);
  hp_sel->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  hp_sel->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  hp_sel->SetTitle("Selected B_{c}'s");
  hp_sel->Draw("COLZ0");
  line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same");

  c2->cd(2);
  //  gPad->SetLogz();
  hp_efficiency->Divide(hp_acc);
  //  hp_efficiency->GetZaxis()->SetRangeUser(0,1);//1e-4,1);
  hp_efficiency->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  hp_efficiency->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  hp_efficiency->SetTitle("Efficiency");
  hp_efficiency->Draw("COLZ0");
  line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same");

  c2->SaveAs("figs/EfficiencyMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".pdf");
  c2->SaveAs("figs/EfficiencyMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".png");

  //**************************************************************
  //Grab acceptance map
  TFile acc_file("../acceptance/acceptanceMap.root","READ");
  TH2Poly* hp_acceptance = (TH2Poly*)acc_file.Get("hp_acceptance");
  hp_acceptance->SetDirectory(0);
  TH2Poly* hp_acceff = (TH2Poly*)hp_acceptance->Clone("hp_acceff");
  hp_acceff->SetDirectory(0);
  hp_acceff->Multiply(hp_efficiency);

  TCanvas *c3 = new TCanvas("c3","c3",2500,2500);
  c3->Divide(2,2);
  
  c3->cd(1);
  hp_acceptance->Draw("COLZ");
  line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same");
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
  line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same");
  c3->cd(3);
  hp_sel->Draw("COLZ0");
  line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same");

  c3->cd(4);
  hp_acceff->Draw("COLZ0");
  line1->Draw("same"); line2->Draw("same"); line3->Draw("same"); line4->Draw("same");
  hp_acceff->SetTitle("Acceptance #times  Efficiency");
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

  c3->SaveAs("figs/AcceptanceEfficiencyMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".pdf");
  c3->SaveAs("figs/AcceptanceEfficiencyMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".png");

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

  c4->SaveAs("figs/BDTefficiencyMap"+(TString)(withTM?"_withTrackerMu":"")+(TString)(ispp?"_pp":"_PbPb")+".pdf");

  //**************************************************************
  //Store maps
  TFile* fout = TFile::Open("AcceptanceEfficiencyMap.root","UPDATE");
  hp_acceff->Write("hp_acceff_"+(TString)(ispp?"pp":"PbPb"));
  hpcoarse_inBDT23->Write("hpcoarse_inBDT23_"+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(integratePtBins?"_integratePtBins":""));
  hpcoarse_inBDT3->Write("hpcoarse_inBDT3_"+(TString)(ispp?"pp":"PbPb")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(integratePtBins?"_integratePtBins":""));
  fout->WriteObject(&eff_oneBinned,"efficiency_oneBinned"+(TString)(ispp?"_pp":"_PbPb"));
  fout->Close();
}
