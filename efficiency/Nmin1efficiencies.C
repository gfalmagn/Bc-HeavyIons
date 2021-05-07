#include <iostream>
#include <fstream>
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
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"
#include "../helpers/SgMuonAcceptanceCuts.h"

void Nmin1efficiencies(bool ispp = true){

  //**************************************************************  
  //Create Tree and branches
  vector<TTree*> Trees; vector<int> nentries;

  TFile *fileData = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/DoubleMu_Run2017G-09Aug2019_UL2017_AOD_Run_306546_306826_OniaTree_TripleMuBc_25012021.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/data/HIDoubleMuon_HIDoubleMuonPsiPeri_Run2018A_AOD_OniaTree_Run_326381_327564_BcTrimuon_25012021.root","READ");
  Trees.push_back( (TTree*)fileData->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[0]->GetEntries() );
  std::cout<<"nevents data = "<<nentries[0]<<"\n";

  TFile *fileMC = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/Oniatree_MC_Bc_trimuons_21112019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/BcToJpsiMuNu_BCVEGPY_PYTHIA8_2018PbPb5TeV_HINPbPbAutumn18DR-00196_08092020_4200k_ONIATREE.root","READ");
  Trees.push_back( (TTree*)fileMC->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[1]->GetEntries() );
  std::cout<<"nevents MC = "<<nentries[1]<<"\n";

  TFile *fileMCb = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/NonPromptJpsi/MConiatree/crab_BJPsiMM_TuneCUETP8M1_5p02TeV_pythia8_05082019_wLambdabFor10_ptHatMinCombined_ONIATREE.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/NonPromptJpsi/BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_HINPbPbAutumn18DR_trimuons_oniatree_14022021.root","READ");
  Trees.push_back( (TTree*)fileMCb->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[2]->GetEntries() );
  std::cout<<"nevents B->J/psi MC = "<<nentries[2]<<"\n";

  TFile *fileMCpromptJ = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/PromptJpsi/Oniatree_MC_trimuons_PromptJpsi_ptHatMinCombined_05082019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/PromptJpsi/Jpsi_pThat-2_TuneCP5_HydjetDrumMB_HINPbPbAutumn18DR_trimuons_oniatree_25012021.root","READ");
  Trees.push_back( (TTree*)fileMCpromptJ->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[3]->GetEntries() );
  std::cout<<"nevents Prompt J/psi MC = "<<nentries[3]<<"\n";

  TFile *fileFlipJpsi = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/flipJpsi/DoubleMu_Run2017G-09Aug2019_UL2017_AOD_Run_306546_306826_OniaTree_flippedJpsi_25012021.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/flippedJpsi/HIDoubleMuon_HIDoubleMuonPsiPeri_Run2018A_AOD_OniaTree_Run_326381_327564_flippedJpsi_BcTrimuon_25012021.root","READ");
  Trees.push_back( (TTree*)fileFlipJpsi->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[4]->GetEntries() );
  std::cout<<"nevents flipped-Jpsi data = "<<nentries[4]<<"\n";

  const int nInputT = Trees.size();

  //**************************************************************  
  //Get input branches

  UInt_t INeventNb[nInputT]; TBranch *b_INeventNb[nInputT];
  UInt_t INrunNb[nInputT]; TBranch *b_INrunNb[nInputT];
  UInt_t INLS[nInputT]; TBranch *b_INLS[nInputT];
  Short_t INnPV[nInputT]; TBranch *b_INnPV[nInputT];
  int INCentrality[nInputT]; TBranch *b_INCentrality[nInputT];
  ULong64_t HLTriggers[nInputT]; TBranch *b_HLTriggers[nInputT];
  Short_t Reco_3mu_size[nInputT]; TBranch *b_Reco_3mu_size[nInputT];
  Short_t Reco_3mu_charge[nInputT][300]; TBranch *b_Reco_3mu_charge[nInputT];
  Short_t Reco_3mu_QQ1_idx[nInputT][300]; TBranch *b_Reco_3mu_QQ1_idx[nInputT];
  Short_t Reco_3mu_QQ2_idx[nInputT][300]; TBranch *b_Reco_3mu_QQ2_idx[nInputT];
  Short_t Reco_3mu_QQss_idx[nInputT][300]; TBranch *b_Reco_3mu_QQss_idx[nInputT];
  Short_t Reco_3mu_muW_idx[nInputT][300]; TBranch *b_Reco_3mu_muW_idx[nInputT];
  Short_t Reco_3mu_muW2_idx[nInputT][300]; TBranch *b_Reco_3mu_muW2_idx[nInputT];
  Short_t Reco_3mu_mumi_idx[nInputT][300]; TBranch *b_Reco_3mu_mumi_idx[nInputT];
  Short_t Reco_3mu_mupl_idx[nInputT][300]; TBranch *b_Reco_3mu_mupl_idx[nInputT];
  float Reco_3mu_ctau[nInputT][300]; TBranch *b_Reco_3mu_ctau[nInputT];
  float Reco_3mu_ctauErr[nInputT][300]; TBranch *b_Reco_3mu_ctauErr[nInputT];
  float Reco_3mu_ctau3D[nInputT][300]; TBranch *b_Reco_3mu_ctau3D[nInputT];
  float Reco_3mu_ctauErr3D[nInputT][300]; TBranch *b_Reco_3mu_ctauErr3D[nInputT];
  float Reco_3mu_cosAlpha[nInputT][300]; TBranch *b_Reco_3mu_cosAlpha[nInputT];
  float Reco_3mu_cosAlpha3D[nInputT][300]; TBranch *b_Reco_3mu_cosAlpha3D[nInputT];
  float Reco_3mu_VtxProb[nInputT][300]; TBranch *b_Reco_3mu_VtxProb[nInputT];
  //  float Reco_3mu_CorrM[nInputT][300]; TBranch *b_Reco_3mu_CorrM[nInputT];
  float Reco_QQ_VtxProb[nInputT][100]; TBranch *b_Reco_QQ_VtxProb[nInputT];
  float Reco_QQ_dca[nInputT][100]; TBranch *b_Reco_QQ_dca[nInputT];
  Short_t Reco_QQ_sign[nInputT][100]; TBranch *b_Reco_QQ_sign[nInputT];
  Short_t Reco_QQ_mumi_idx[nInputT][100]; TBranch *b_Reco_QQ_mumi_idx[nInputT];
  Short_t Reco_QQ_mupl_idx[nInputT][100]; TBranch *b_Reco_QQ_mupl_idx[nInputT];
  Short_t Reco_mu_charge[nInputT][300]; TBranch *b_Reco_mu_charge[nInputT];
  float Reco_mu_dxy[nInputT][100]; TBranch *b_Reco_mu_dxy[nInputT];
  float Reco_mu_dz[nInputT][100]; TBranch *b_Reco_mu_dz[nInputT];
  float Reco_mu_dxyErr[nInputT][100]; TBranch *b_Reco_mu_dxyErr[nInputT];
  float Reco_mu_dzErr[nInputT][100]; TBranch *b_Reco_mu_dzErr[nInputT];
  float Reco_mu_normChi2_global[nInputT][100]; TBranch *b_Reco_mu_normChi2_global[nInputT];
  float Reco_mu_normChi2_inner[nInputT][100]; TBranch *b_Reco_mu_normChi2_inner[nInputT];
  int Reco_mu_SelType[nInputT][100]; TBranch *b_Reco_mu_SelType[nInputT];  
  int Reco_mu_nPixWMea[nInputT][100]; TBranch *b_Reco_mu_nPixWMea[nInputT];  
  int Reco_mu_nTrkWMea[nInputT][100]; TBranch *b_Reco_mu_nTrkWMea[nInputT];  
  ULong64_t Reco_mu_trig[nInputT][100]; TBranch *b_Reco_mu_trig[nInputT];
  bool Reco_mu_highPurity[nInputT][100]; TBranch *b_Reco_mu_highPurity[nInputT];
  bool Reco_mu_InLooseAcc[nInputT][100]; TBranch *b_Reco_mu_InLooseAcc[nInputT];
  bool Reco_mu_InTightAcc[nInputT][100]; TBranch *b_Reco_mu_InTightAcc[nInputT];
  bool Reco_trk_InLooseAcc[nInputT][100]; TBranch *b_Reco_trk_InLooseAcc[nInputT];
  bool Reco_trk_InTightAcc[nInputT][100]; TBranch *b_Reco_trk_InTightAcc[nInputT];
  float Reco_trk_dxy[nInputT][100]; TBranch *b_Reco_trk_dxy[nInputT];
  float Reco_trk_dz[nInputT][100]; TBranch *b_Reco_trk_dz[nInputT];
  float Reco_trk_dxyError[nInputT][100]; TBranch *b_Reco_trk_dxyError[nInputT];
  float Reco_trk_dzError[nInputT][100]; TBranch *b_Reco_trk_dzError[nInputT];
  float Reco_3mu_muW_dz_muonlessVtx[nInputT][300]; TBranch *b_Reco_3mu_muW_dz_muonlessVtx[nInputT]; //only in pp
  float Reco_3mu_mumi_dz_muonlessVtx[nInputT][300]; TBranch *b_Reco_3mu_mumi_dz_muonlessVtx[nInputT]; //only in pp
  float Reco_3mu_mupl_dz_muonlessVtx[nInputT][300]; TBranch *b_Reco_3mu_mupl_dz_muonlessVtx[nInputT]; //only in pp
  float Reco_3mu_muW_dxy_muonlessVtx[nInputT][300]; TBranch *b_Reco_3mu_muW_dxy_muonlessVtx[nInputT]; 
  float Reco_3mu_mumi_dxy_muonlessVtx[nInputT][300]; TBranch *b_Reco_3mu_mumi_dxy_muonlessVtx[nInputT];
  float Reco_3mu_mupl_dxy_muonlessVtx[nInputT][300]; TBranch *b_Reco_3mu_mupl_dxy_muonlessVtx[nInputT];
  TClonesArray *Reco_3mu_4mom[nInputT]; TBranch *b_Reco_3mu_4mom[nInputT];
  TClonesArray *Reco_mu_4mom[nInputT]; TBranch *b_Reco_mu_4mom[nInputT];
  TClonesArray *Reco_trk_4mom[nInputT]; TBranch *b_Reco_trk_4mom[nInputT];
  TClonesArray *Reco_QQ_4mom[nInputT]; TBranch *b_Reco_QQ_4mom[nInputT];
  TClonesArray *Reco_QQ_mumi_4mom[nInputT]; TBranch *b_Reco_QQ_mumi_4mom[nInputT];
  TClonesArray *Reco_QQ_mupl_4mom[nInputT]; TBranch *b_Reco_QQ_mupl_4mom[nInputT];
  TClonesArray *Gen_QQ_4mom[nInputT]; TBranch *b_Gen_QQ_4mom[nInputT];
  TClonesArray *Gen_3mu_4mom[nInputT]; TBranch *b_Gen_3mu_4mom[nInputT];
  TClonesArray *Gen_Bc_4mom[nInputT]; TBranch *b_Gen_Bc_4mom[nInputT];
  float Gen_weight[nInputT]; TBranch *b_Gen_weight[nInputT];
  float pthatweight[nInputT]; TBranch *b_pthatweight[nInputT];
  Short_t Gen_QQ_size[nInputT]; TBranch *b_Gen_QQ_size[nInputT];
  Short_t Reco_QQ_size[nInputT]; TBranch *b_Reco_QQ_size[nInputT];
  Short_t Reco_mu_size[nInputT]; TBranch *b_Reco_mu_size[nInputT];
  int Gen_QQ_momId[nInputT][100]; TBranch *b_Gen_QQ_momId[nInputT];
  Short_t Reco_QQ_whichGen[nInputT][100]; TBranch *b_Reco_QQ_whichGen[nInputT];
  bool Reco_3mu_muW_isGenJpsiBro[nInputT][300]; TBranch *b_Reco_3mu_muW_isGenJpsiBro[nInputT];
  int Reco_3mu_muW_trueId[nInputT][300]; TBranch *b_Reco_3mu_muW_trueId[nInputT];
  Short_t Reco_3mu_whichGen[nInputT][300]; TBranch *b_Reco_3mu_whichGen[nInputT];
  Short_t Reco_QQ_flipJpsi[nInputT][300]; TBranch *b_Reco_QQ_flipJpsi[nInputT];

  for(int i=0;i<nInputT;i++){
    b_INeventNb[i] = Trees[i]->GetBranch("eventNb");
    b_INeventNb[i]->SetAddress(&INeventNb[i]);

    // if(i==0 || i==4 || i==5){//data only
    //   b_INrunNb[i] = Trees[i]->GetBranch("runNb");
    //   b_INrunNb[i]->SetAddress(&INrunNb[i]);

    //   b_INLS[i] = Trees[i]->GetBranch("LS");
    //   b_INLS[i]->SetAddress(&INLS[i]);
    // }

    // if(ispp){
    //   b_INnPV[i] = Trees[i]->GetBranch("nPV");
    //   b_INnPV[i]->SetAddress(&INnPV[i]);
    // } else 
    if(!ispp) {
      b_INCentrality[i] = Trees[i]->GetBranch("Centrality");
      b_INCentrality[i]->SetAddress(&INCentrality[i]);
    }

    b_HLTriggers[i] = Trees[i]->GetBranch("HLTriggers");
    b_HLTriggers[i]->SetAddress(&HLTriggers[i]);

    b_Reco_3mu_size[i] = Trees[i]->GetBranch("Reco_3mu_size");
    b_Reco_3mu_size[i]->SetAddress(&Reco_3mu_size[i]);

    b_Reco_3mu_charge[i] = Trees[i]->GetBranch("Reco_3mu_charge");
    b_Reco_3mu_charge[i]->SetAddress(&Reco_3mu_charge[i]);

    b_Reco_3mu_QQ1_idx[i] = Trees[i]->GetBranch("Reco_3mu_QQ1_idx");
    b_Reco_3mu_QQ1_idx[i]->SetAddress(&Reco_3mu_QQ1_idx[i]);

    b_Reco_3mu_QQ2_idx[i] = Trees[i]->GetBranch("Reco_3mu_QQ2_idx");
    b_Reco_3mu_QQ2_idx[i]->SetAddress(&Reco_3mu_QQ2_idx[i]);

    b_Reco_3mu_QQss_idx[i] = Trees[i]->GetBranch("Reco_3mu_QQss_idx");
    b_Reco_3mu_QQss_idx[i]->SetAddress(&Reco_3mu_QQss_idx[i]);

    if(i<4){
      b_Reco_3mu_muW2_idx[i] = Trees[i]->GetBranch("Reco_3mu_muW2_idx");
      b_Reco_3mu_muW2_idx[i]->SetAddress(&Reco_3mu_muW2_idx[i]);
    }
    
    b_Reco_3mu_muW_idx[i] = Trees[i]->GetBranch("Reco_3mu_muW_idx");
    b_Reco_3mu_muW_idx[i]->SetAddress(&Reco_3mu_muW_idx[i]);

    b_Reco_3mu_mumi_idx[i] = Trees[i]->GetBranch("Reco_3mu_mumi_idx"); 
    b_Reco_3mu_mumi_idx[i]->SetAddress(&Reco_3mu_mumi_idx[i]);

    b_Reco_3mu_mupl_idx[i] = Trees[i]->GetBranch("Reco_3mu_mupl_idx");
    b_Reco_3mu_mupl_idx[i]->SetAddress(&Reco_3mu_mupl_idx[i]);

    b_Reco_3mu_ctau[i] = Trees[i]->GetBranch("Reco_3mu_ctau");
    b_Reco_3mu_ctau[i]->SetAddress(&Reco_3mu_ctau[i]);

    b_Reco_3mu_ctauErr[i] = Trees[i]->GetBranch("Reco_3mu_ctauErr");
    b_Reco_3mu_ctauErr[i]->SetAddress(&Reco_3mu_ctauErr[i]);

    b_Reco_3mu_ctau3D[i] = Trees[i]->GetBranch("Reco_3mu_ctau3D");
    b_Reco_3mu_ctau3D[i]->SetAddress(&Reco_3mu_ctau3D[i]);

    b_Reco_3mu_ctauErr3D[i] = Trees[i]->GetBranch("Reco_3mu_ctauErr3D");
    b_Reco_3mu_ctauErr3D[i]->SetAddress(&Reco_3mu_ctauErr3D[i]);

    b_Reco_3mu_cosAlpha[i] = Trees[i]->GetBranch("Reco_3mu_cosAlpha");
    b_Reco_3mu_cosAlpha[i]->SetAddress(&Reco_3mu_cosAlpha[i]);

    b_Reco_3mu_cosAlpha3D[i] = Trees[i]->GetBranch("Reco_3mu_cosAlpha3D");
    b_Reco_3mu_cosAlpha3D[i]->SetAddress(&Reco_3mu_cosAlpha3D[i]);

    b_Reco_3mu_VtxProb[i] = Trees[i]->GetBranch("Reco_3mu_VtxProb");
    b_Reco_3mu_VtxProb[i]->SetAddress(&Reco_3mu_VtxProb[i]);

    b_Reco_QQ_VtxProb[i] = Trees[i]->GetBranch("Reco_QQ_VtxProb");
    b_Reco_QQ_VtxProb[i]->SetAddress(&Reco_QQ_VtxProb[i]);

    b_Reco_QQ_dca[i] = Trees[i]->GetBranch("Reco_QQ_dca");
    b_Reco_QQ_dca[i]->SetAddress(&Reco_QQ_dca[i]);

    b_Reco_QQ_sign[i] = Trees[i]->GetBranch("Reco_QQ_sign");
    b_Reco_QQ_sign[i]->SetAddress(&Reco_QQ_sign[i]);

    b_Reco_QQ_mumi_idx[i] = Trees[i]->GetBranch("Reco_QQ_mumi_idx");
    b_Reco_QQ_mumi_idx[i]->SetAddress(&Reco_QQ_mumi_idx[i]);

    b_Reco_QQ_mupl_idx[i] = Trees[i]->GetBranch("Reco_QQ_mupl_idx");
    b_Reco_QQ_mupl_idx[i]->SetAddress(&Reco_QQ_mupl_idx[i]);

    b_Reco_mu_charge[i] = Trees[i]->GetBranch("Reco_mu_charge");
    b_Reco_mu_charge[i]->SetAddress(&Reco_mu_charge[i]);

    b_Reco_mu_dxy[i] = Trees[i]->GetBranch("Reco_mu_dxy");
    b_Reco_mu_dxy[i]->SetAddress(&Reco_mu_dxy[i]);

    b_Reco_mu_dz[i] = Trees[i]->GetBranch("Reco_mu_dz");
    b_Reco_mu_dz[i]->SetAddress(&Reco_mu_dz[i]);

    b_Reco_mu_dxyErr[i] = Trees[i]->GetBranch("Reco_mu_dxyErr");
    b_Reco_mu_dxyErr[i]->SetAddress(&Reco_mu_dxyErr[i]);

    b_Reco_mu_dzErr[i] = Trees[i]->GetBranch("Reco_mu_dzErr");
    b_Reco_mu_dzErr[i]->SetAddress(&Reco_mu_dzErr[i]);

    b_Reco_mu_SelType[i] = Trees[i]->GetBranch("Reco_mu_SelectionType");
    b_Reco_mu_SelType[i]->SetAddress(&Reco_mu_SelType[i]);

    b_Reco_mu_nPixWMea[i] = Trees[i]->GetBranch("Reco_mu_nPixWMea");
    b_Reco_mu_nPixWMea[i]->SetAddress(&Reco_mu_nPixWMea[i]);

    b_Reco_mu_nTrkWMea[i] = Trees[i]->GetBranch("Reco_mu_nTrkWMea");
    b_Reco_mu_nTrkWMea[i]->SetAddress(&Reco_mu_nTrkWMea[i]);

    b_Reco_mu_trig[i] = Trees[i]->GetBranch("Reco_mu_trig");
    b_Reco_mu_trig[i]->SetAddress(&Reco_mu_trig[i]);

    b_Reco_mu_highPurity[i] = Trees[i]->GetBranch("Reco_mu_highPurity");
    b_Reco_mu_highPurity[i]->SetAddress(&Reco_mu_highPurity[i]);

    b_Reco_mu_InLooseAcc[i] = Trees[i]->GetBranch("Reco_mu_InLooseAcc");
    b_Reco_mu_InLooseAcc[i]->SetAddress(&Reco_mu_InLooseAcc[i]);

    b_Reco_mu_InTightAcc[i] = Trees[i]->GetBranch("Reco_mu_InTightAcc");
    b_Reco_mu_InTightAcc[i]->SetAddress(&Reco_mu_InTightAcc[i]);

    if(ispp || i==4){ //for pp and flipJpsi
      b_Reco_3mu_muW_dz_muonlessVtx[i] = Trees[i]->GetBranch("Reco_3mu_muW_dz_muonlessVtx");
      b_Reco_3mu_muW_dz_muonlessVtx[i]->SetAddress(&Reco_3mu_muW_dz_muonlessVtx[i]);

      b_Reco_3mu_mumi_dz_muonlessVtx[i] = Trees[i]->GetBranch("Reco_3mu_mumi_dz_muonlessVtx");
      b_Reco_3mu_mumi_dz_muonlessVtx[i]->SetAddress(&Reco_3mu_mumi_dz_muonlessVtx[i]);

      b_Reco_3mu_mupl_dz_muonlessVtx[i] = Trees[i]->GetBranch("Reco_3mu_mupl_dz_muonlessVtx");
      b_Reco_3mu_mupl_dz_muonlessVtx[i]->SetAddress(&Reco_3mu_mupl_dz_muonlessVtx[i]);

      b_Reco_3mu_muW_dxy_muonlessVtx[i] = Trees[i]->GetBranch("Reco_3mu_muW_dxy_muonlessVtx");
      b_Reco_3mu_muW_dxy_muonlessVtx[i]->SetAddress(&Reco_3mu_muW_dxy_muonlessVtx[i]);

      b_Reco_3mu_mumi_dxy_muonlessVtx[i] = Trees[i]->GetBranch("Reco_3mu_mumi_dxy_muonlessVtx");
      b_Reco_3mu_mumi_dxy_muonlessVtx[i]->SetAddress(&Reco_3mu_mumi_dxy_muonlessVtx[i]);

      b_Reco_3mu_mupl_dxy_muonlessVtx[i] = Trees[i]->GetBranch("Reco_3mu_mupl_dxy_muonlessVtx");
      b_Reco_3mu_mupl_dxy_muonlessVtx[i]->SetAddress(&Reco_3mu_mupl_dxy_muonlessVtx[i]);
    }

    Reco_3mu_4mom[i] = new TClonesArray();
    Trees[i]->SetBranchAddress("Reco_3mu_4mom", &Reco_3mu_4mom[i], &b_Reco_3mu_4mom[i]);

    Reco_mu_4mom[i] = new TClonesArray();
    Trees[i]->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom[i], &b_Reco_mu_4mom[i]);

    Reco_QQ_4mom[i] = new TClonesArray();
    Trees[i]->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom[i], &b_Reco_QQ_4mom[i]);

    b_Reco_QQ_size[i] = Trees[i]->GetBranch("Reco_QQ_size");
    b_Reco_QQ_size[i]->SetAddress(&Reco_QQ_size[i]);

    b_Reco_mu_size[i] = Trees[i]->GetBranch("Reco_mu_size");
    b_Reco_mu_size[i]->SetAddress(&Reco_mu_size[i]);

    if(i>0 && i<4){
      Gen_QQ_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom[i], &b_Gen_QQ_4mom[i]);

      b_Reco_QQ_whichGen[i] = Trees[i]->GetBranch("Reco_QQ_whichGen");
      b_Reco_QQ_whichGen[i]->SetAddress(&Reco_QQ_whichGen[i]);

      if(i==2 || i==3){
	if(ispp){
	  b_pthatweight[i] = Trees[i]->GetBranch("pthatweight");
	  b_pthatweight[i]->SetAddress(&pthatweight[i]);}
	else{
	  b_Gen_weight[i] = Trees[i]->GetBranch("Gen_weight");
	  b_Gen_weight[i]->SetAddress(&Gen_weight[i]);}	  
      }

      b_Reco_3mu_muW_isGenJpsiBro[i] = Trees[i]->GetBranch("Reco_3mu_muW_isGenJpsiBro");
      b_Reco_3mu_muW_isGenJpsiBro[i]->SetAddress(&Reco_3mu_muW_isGenJpsiBro[i]);

      b_Reco_3mu_whichGen[i] = Trees[i]->GetBranch("Reco_3mu_whichGen");
      b_Reco_3mu_whichGen[i]->SetAddress(&Reco_3mu_whichGen[i]);
    }

    if(i==4){
      b_Reco_QQ_flipJpsi[i] = Trees[i]->GetBranch("Reco_QQ_flipJpsi");
      b_Reco_QQ_flipJpsi[i]->SetAddress(&Reco_QQ_flipJpsi[i]);

      Reco_QQ_mumi_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom[i], &b_Reco_QQ_mumi_4mom[i]);

      Reco_QQ_mupl_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom[i], &b_Reco_QQ_mupl_4mom[i]);
    }
  }

  //**************************************************************
  //Some variables and XS corrections
  TFile* fXS = TFile::Open("/home/llr/cms/falmagne/Bc/MCnormalization/NonPromptJpsiXS_scalefactor.root");
  TF1* SFpp = (TF1*)fXS->Get("SFpp");
  std::pair<TF1*,float> scalepp (SFpp, L_pp * 0.06); // SF(data/MC)(Non prompt Jpsi XS) * Lumi_pp[pb-1] (preliminary, from https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/794.html) * BF(Jpsi-> mu mu)
  TF1* SFPbPb = (TF1*)fXS->Get("SFPbPb");
  std::pair<TF1*,float> scalePbPb (SFPbPb, Leq_PbPb* 0.06);// SF(data/MC)(Non prompt Jpsi XS PbPb) * A^2 * Lumi_PbPb[pb-1] / BF(Jpsi-> mu mu)
  std::map<bool, std::pair<TF1*,float> > scaleMCb = {{ true, scalepp }, 
						     { false, scalePbPb } }; 

  TFile* fXSprompt = TFile::Open("/home/llr/cms/falmagne/Bc/MCnormalization/PromptJpsiXS_scalefactor.root");
  TF1* SFppPrompt = (TF1*)fXSprompt->Get("SFpp");
  std::pair<TF1*,float> scaleppPrompt (SFppPrompt, L_pp * 0.06); // SF(data/MC)(Non prompt Jpsi XS) * Lumi_pp[pb-1] (preliminary, from https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/794.html) * BF(Jpsi-> mu mu)
  TF1* SFPbPbPrompt = (TF1*)fXSprompt->Get("SFPbPb");
  std::pair<TF1*,float> scalePbPbPrompt (SFPbPbPrompt, Leq_PbPb* 0.06);// SF(data/MC)(Non prompt Jpsi XS PbPb) * A^2 * Lumi_PbPb[pb-1] / BF(Jpsi-> mu mu)
  std::map<bool, std::pair<TF1*,float> > scaleMCprompt = {{ true, scaleppPrompt }, 
							  { false, scalePbPbPrompt } }; 
  //For Jpsi choice weight
  TFile* fJpsiM = TFile::Open("../BDT/JpsiMassDistr.root","READ");
  TH1F* JpsiM = (TH1F*)fJpsiM->Get("JpsiMass_data_allBDTbins_"+(TString)(ispp?"pp":"PbPb"));
  TH1F* JpsiM_tight = (TH1F*)fJpsiM->Get("JpsiMassCentralEta_data_allBDTbins_"+(TString)(ispp?"pp":"PbPb"));

  //**************************************************************
  //Samples and variables to test
  int samNum = 7;//hardcoded
  TString samName[] = {"signal","fakeJpsi","dataHighMass","bMC ","bMCcorr","PromptMC","flipJpsi"}; //name of studied sample
  vector<int> samToDo[] = {{1,2}, {0}, {3,4}, {5},{6}}; //indices of samples that must be checked for a given input tree

  int effNum = 16;//hardcoded
  TString effToTest[] = {"empty","3global","TrimuVProb","DimuVProb","TrimuVProbDimuVProb","alpha3D","alpha2D","alpha2Dalpha3D","tauSignif3D","tauSignif2D","tauSignif3DtauSignif2D","dca","CorrM","dRsum","dz",""}; //the first empty one is necessary for subtetlies with the Jpsi choice (and associated weight)
  for(int e=2;e<15;e++) effToTest[15] += effToTest[e]; //full efficiency
  vector<vector<float> > npass(samNum, vector<float>(effNum,0) );
  vector<vector<float> > ntot(samNum, vector<float>(effNum,0) );

  //**************************************************************
  //For each input tree, redirect data/MC events to the right tree, and write into branches
  bool goodTree = false;

  for(int iIpt=0; iIpt<nInputT; iIpt++){//nInputT
    for(int j=0; j<nentries[iIpt]; j++){//nentries[iIpt]

      if(j%200000==0){ cout<<"Scanned "<<100.*(double)j/nentries[iIpt]<<"% of tree #"<<iIpt<<endl; }

      b_Reco_3mu_size[iIpt]->GetEntry(j);
      if(Reco_3mu_size[iIpt]==0) continue;

      Reco_3mu_4mom[iIpt]->Clear();
      Reco_mu_4mom[iIpt]->Clear();
      Reco_QQ_4mom[iIpt]->Clear();
      if(iIpt==4){
	Reco_QQ_mumi_4mom[iIpt]->Clear();
	Reco_QQ_mupl_4mom[iIpt]->Clear();}

      Trees[iIpt]->GetEntry(j);

      for(int BcNb=0;BcNb<Reco_3mu_size[iIpt];BcNb++){
	if(fabs(Reco_3mu_charge[iIpt][BcNb])!=1) continue; //no wrongsign here

	Short_t QQidx_3[3] = { Reco_3mu_QQ1_idx[iIpt][BcNb], Reco_3mu_QQ2_idx[iIpt][BcNb], Reco_3mu_QQss_idx[iIpt][BcNb] };
	Short_t muW3idx = (Reco_3mu_mumi_idx[iIpt][BcNb]!=Reco_3mu_muW2_idx[iIpt][BcNb])?(Reco_3mu_mumi_idx[iIpt][BcNb]):(Reco_3mu_mupl_idx[iIpt][BcNb]); //for the third dimuon choice, the muWidx=old_muWmi if old_muWmi!=muW2idx, and muWidx=old_muWpl otherwise
	Short_t muWidx_3[3] = { Reco_3mu_muW_idx[iIpt][BcNb], Reco_3mu_muW2_idx[iIpt][BcNb], muW3idx };
	Short_t mumiidx_3[3] = { Reco_3mu_mumi_idx[iIpt][BcNb], Reco_QQ_mumi_idx[iIpt][QQidx_3[1]], Reco_QQ_mumi_idx[iIpt][QQidx_3[2]] };
	Short_t muplidx_3[3] = { Reco_3mu_mupl_idx[iIpt][BcNb], Reco_QQ_mupl_idx[iIpt][QQidx_3[1]], Reco_QQ_mupl_idx[iIpt][QQidx_3[2]] };
	if(iIpt==4){
	  if(mumiidx_3[0] != mumiidx_3[1] && mumiidx_3[0] != muplidx_3[1]) muWidx_3[1] = mumiidx_3[0]; //Reco_3mu_muW2_idx is badly defined for flipJpsi
	  if(muplidx_3[0] != mumiidx_3[1] && muplidx_3[0] != muplidx_3[1]) muWidx_3[1] = muplidx_3[0];
	}

	int nCandidatePairs = 2;

	//**************************************************************
	//Loop on the 2 or 3 possible Jpsi dimuon choices 
	int thisTrimuGaveABc = -1;
	for(int k=0; k<nCandidatePairs; k++){
          int k2 = (k==0)?1:0;
	  Short_t QQidx = QQidx_3[k];
	  Short_t QQ2idx = QQidx_3[k2];
	  Short_t muWidx = muWidx_3[k];
	  Short_t mumiidx = mumiidx_3[k];
	  Short_t muplidx = muplidx_3[k];

	  if(QQidx>-1){
	    TLorentzVector *recBc = (TLorentzVector*) Reco_3mu_4mom[iIpt]->At(BcNb);
	    TLorentzVector *recQQ = (TLorentzVector*) Reco_QQ_4mom[iIpt]->At(QQidx);
	    if(recQQ==NULL) continue;
	    if(muWidx==-1 || mumiidx==-1 || muplidx==-1) {cout<<"!!!!!!! one muon has a -1 index !! Skip this candidate"<<endl; continue;}
	    if((muWidx>=Reco_mu_size[iIpt] || mumiidx>=Reco_mu_size[iIpt] || muplidx>=Reco_mu_size[iIpt])) {cout<<"!!!!!!! Muon index > muons vectro size ! Skip this candidate"<<endl; continue;}
	    TLorentzVector *recBc_mumi = (TLorentzVector*) ((iIpt<4)?(Reco_mu_4mom[iIpt]->At(mumiidx)):(Reco_QQ_mumi_4mom[iIpt]->At(QQidx)));
	    TLorentzVector *recBc_mupl = (TLorentzVector*) ((iIpt<4)?(Reco_mu_4mom[iIpt]->At(muplidx)):(Reco_QQ_mupl_4mom[iIpt]->At(QQidx)));
	    TLorentzVector *recBc_muW = (TLorentzVector*) ((iIpt==4 && k==1)?(
									      (Reco_mu_charge[iIpt][muWidx]>0)?(Reco_QQ_mupl_4mom[iIpt]->At(QQidx_3[0])):(Reco_QQ_mumi_4mom[iIpt]->At(QQidx_3[0]))
									      ):(
										 Reco_mu_4mom[iIpt]->At(muWidx)));

	    float BcCandE = sqrt(pow(recQQ->P(),2) + pow(m_Jpsi,2) ) + recBc_muW->E() ;
	    float BcCandM = sqrt( pow( BcCandE ,2) - pow(recBc->P(),2) );
	    float QQCandM = recQQ->M();

	    float muW_eta = recBc_muW->Eta();
	    float mumi_eta = recBc_mumi->Eta();
	    float mupl_eta = recBc_mupl->Eta();
	    float maxEta = max(fabs(muW_eta),max(fabs(mumi_eta),fabs(mupl_eta)));
	    float weight = 1;

	    int nstart = samToDo[iIpt][0]; int nstop = samToDo[iIpt][samToDo[iIpt].size()-1];
	    for(int i=nstart; i<=nstop; i++){
	      //**************************************************************
	      //What sample should I check?
	      bool goodSam = false;
	      switch(i) {
	      case 0: //MC SIGNAL
		goodSam = Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mBcMax) && (BcCandM > _mBcMin) // in Bc mass region
		  && Reco_3mu_whichGen[iIpt][BcNb]>-1
		  && Reco_QQ_whichGen[iIpt][QQidx]>-1;
		weight = _scaleMCsig[ispp];
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) {weight *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ weight = 0;}
		break;

	      case 1: //Jpsi mass sidebands
		goodSam = Reco_QQ_sign[iIpt][QQidx]==0 && inJpsiMassSB(QQCandM, maxEta<1.5) 
		  && (BcCandM < _mBcMax) && (BcCandM > _mBcMin);
		weight = 1;
		break;

	      case 2: //HighMassRegion in data
		goodSam = Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM)//(inJpsiMassRange(QQCandM, maxEta<1.5) || inJpsiMassSB(QQCandM, maxEta<1.5))
		  && (BcCandM > m_Bc) && (BcCandM < _mMax);
		if(inJpsiMassRange(QQCandM, maxEta<1.5)){
		  weight = 1;}
		else if(inJpsiMassSB(QQCandM, maxEta<1.5)){
		  weight = -1;} // simple background subtraction (considering linear background) to keep only true J/psi's
		else{
		  weight = 0; goodSam = false;}
		break;

	      case 3: //MC NonPrompt J/psi
		goodSam = Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mBcMax) && (BcCandM > _mBcMin) // in Bc mass region
		  && Reco_3mu_whichGen[iIpt][BcNb]==-1
		  && Reco_QQ_whichGen[iIpt][QQidx]>-1;
		if(goodSam){
		  TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom[iIpt]->At(Reco_QQ_whichGen[iIpt][QQidx]);
		  weight = ((ispp)?pthatweight[iIpt]:Gen_weight[iIpt]) * scaleMCb[ispp].second * (scaleMCb[ispp].first)->Eval(genQQ->Pt());
		}
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) { weight *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ weight = 0; goodSam = false;}
		break;

	      case 4: //MC NonPrompt B->J/psi X only (correlated)
		goodSam = Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mBcMax) && (BcCandM > _mBcMin) // in Bc mass region
		  && Reco_3mu_whichGen[iIpt][BcNb]==-1
		  && Reco_QQ_whichGen[iIpt][QQidx]>-1
		  && Reco_3mu_muW_isGenJpsiBro[iIpt][BcNb];
		if(goodSam){
		  TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom[iIpt]->At(Reco_QQ_whichGen[iIpt][QQidx]);
		  weight = ((ispp)?pthatweight[iIpt]:Gen_weight[iIpt]) * scaleMCb[ispp].second * (scaleMCb[ispp].first)->Eval(genQQ->Pt());
		}
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) { weight *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ weight = 0; goodSam = false;}
		break;

	      case 5: //MC Prompt J/psi
		goodSam = Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mBcMax) && (BcCandM > _mBcMin) // in Bc mass region
		  && Reco_3mu_whichGen[iIpt][BcNb]==-1
		  && Reco_QQ_whichGen[iIpt][QQidx]>-1;
		if(goodSam){
		  TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom[iIpt]->At(Reco_QQ_whichGen[iIpt][QQidx]);
		  weight = ((ispp)?pthatweight[iIpt]:Gen_weight[iIpt]) * scaleMCprompt[ispp].second * (scaleMCprompt[ispp].first)->Eval(genQQ->Pt());
		}
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) { weight *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ weight = 0; goodSam = false;}
		break;

	      case 6: //Jpsi flipping
		goodSam = Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM)
		  && (BcCandM < _mBcMax) && (BcCandM > _mBcMin);
		if(inJpsiMassRange(QQCandM, maxEta<1.5)){
		  weight = 1/13.;} //7 versions of Jpsi flipping are added
		else if(inJpsiMassSB(QQCandM, maxEta<1.5)){
		  weight = -1/13.;} // simple background subtraction (considering linear background) to keep only true J/psi's
		else{
		  weight = 0; goodSam = false;}
		break;

	      }

	      if(!goodSam) continue;
	      if (thisTrimuGaveABc==i) continue; //two dimuon choices cannot give two candidates in the same output tree //Redundant check !

	      if(!ispp && (iIpt>0 && iIpt<4) ) {float weightNcoll = (float)findNcoll(INCentrality[iIpt]); //for PbPb MC
		weight *= weightNcoll;}
	      
	      //**************************************************************
	      //Write the wanted variables into the chosen (goodSam=true) output tree
	      float Bc_ctauSignif = Reco_3mu_ctau[iIpt][BcNb] / Reco_3mu_ctauErr[iIpt][BcNb] ;
	      float Bc_ctauSignif3D = Reco_3mu_ctau3D[iIpt][BcNb] / Reco_3mu_ctauErr3D[iIpt][BcNb] ;
	      float Bc_alpha = TMath::ACos(Reco_3mu_cosAlpha[iIpt][BcNb]);
	      float Bc_alpha3D = TMath::ACos(Reco_3mu_cosAlpha3D[iIpt][BcNb]);
	      float Bc_VtxProb = Reco_3mu_VtxProb[iIpt][BcNb];
	      float QQ_VtxProb = Reco_QQ_VtxProb[iIpt][QQidx];
	      float QQ_dca = Reco_QQ_dca[iIpt][QQidx];
	      float QQ2_VtxProb = (QQ2idx<Reco_QQ_size[iIpt] && QQ2idx>-1)?(Reco_QQ_VtxProb[iIpt][QQ2idx]):0;
	      float QQ2_dca = (QQ2idx<Reco_QQ_size[iIpt] && QQ2idx>-1)?(Reco_QQ_dca[iIpt][QQ2idx]):100;
	      	      
	      bool muW_inLooseAcc = Reco_mu_InLooseAcc[iIpt][muWidx];
	      bool muW_inTightAcc = Reco_mu_InTightAcc[iIpt][muWidx];
	      bool muW_isGlb = (Reco_mu_SelType[iIpt][muWidx]&2)>0;
	      bool muW_trig = (Reco_mu_trig[iIpt][muWidx]&(ispp?8:4096))>0; //DoubleMu0 trigger = 2^3 for pp //HL_THIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 = 2^12 for PbPb
	      bool muW_isSoft = isSoft(false, muW_isGlb , (Reco_mu_SelType[iIpt][muWidx]&8)>0, ((Reco_mu_SelType[iIpt][muWidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[iIpt][muWidx] , Reco_mu_dxy[iIpt][muWidx] , Reco_mu_dz[iIpt][muWidx] , Reco_mu_nPixWMea[iIpt][muWidx] , Reco_mu_nTrkWMea[iIpt][muWidx] );
	    
	      bool mumi_inLooseAcc = Reco_mu_InLooseAcc[iIpt][mumiidx];
	      bool mumi_inTightAcc = Reco_mu_InTightAcc[iIpt][mumiidx];
	      bool mumi_isGlb = (Reco_mu_SelType[iIpt][mumiidx]&2)>0;
	      bool mumi_trig = (Reco_mu_trig[iIpt][mumiidx]&(ispp?8:4096))>0;
	      bool mumi_isSoft = isSoft( false, mumi_isGlb , (Reco_mu_SelType[iIpt][mumiidx]&8)>0 , ((Reco_mu_SelType[iIpt][mumiidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[iIpt][mumiidx] , Reco_mu_dxy[iIpt][mumiidx] , Reco_mu_dz[iIpt][mumiidx] , Reco_mu_nPixWMea[iIpt][mumiidx] , Reco_mu_nTrkWMea[iIpt][mumiidx] );
	    
	      bool mupl_inLooseAcc = Reco_mu_InLooseAcc[iIpt][muplidx];
	      bool mupl_inTightAcc = Reco_mu_InTightAcc[iIpt][muplidx];
	      bool mupl_isGlb = (Reco_mu_SelType[iIpt][muplidx]&2)>0;
	      bool mupl_trig = (Reco_mu_trig[iIpt][muplidx]&(ispp?8:4096))>0;
	      bool mupl_isSoft = isSoft( false, mupl_isGlb , (Reco_mu_SelType[iIpt][muWidx]&8)>0 , ((Reco_mu_SelType[iIpt][muplidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[iIpt][muplidx] , Reco_mu_dxy[iIpt][muplidx] , Reco_mu_dz[iIpt][muplidx] , Reco_mu_nPixWMea[iIpt][muplidx] , Reco_mu_nTrkWMea[iIpt][muplidx] );
	    
	      float muW_dxy = (ispp || iIpt==4)?(Reco_3mu_muW_dxy_muonlessVtx[iIpt][BcNb]):(Reco_mu_dxy[iIpt][muWidx]);
	      float mumi_dxy = (ispp || iIpt==4)?(Reco_3mu_mumi_dxy_muonlessVtx[iIpt][BcNb]):(Reco_mu_dxy[iIpt][mumiidx]);
	      float mupl_dxy = (ispp || iIpt==4)?(Reco_3mu_mupl_dxy_muonlessVtx[iIpt][BcNb]):(Reco_mu_dxy[iIpt][muplidx]);
	      float muW_dz = (ispp || iIpt==4)?(Reco_3mu_muW_dz_muonlessVtx[iIpt][BcNb]):(Reco_mu_dz[iIpt][muWidx]);
	      float mumi_dz = (ispp || iIpt==4)?(Reco_3mu_mumi_dz_muonlessVtx[iIpt][BcNb]):(Reco_mu_dz[iIpt][mumiidx]);
	      float mupl_dz = (ispp || iIpt==4)?(Reco_3mu_mupl_dz_muonlessVtx[iIpt][BcNb]):(Reco_mu_dz[iIpt][muplidx]);
	      float muW_Pt = recBc_muW->Pt();
	      float mumi_Pt = recBc_mumi->Pt();
	      float mupl_Pt = recBc_mupl->Pt();
	    
	      float QQM_correction = m_Jpsi / QQCandM;
	      float QQ2M_correction = m_Jpsi / QQCandM;

	      //corrected mass
	      float Bc_Pt = sqrt( pow(QQM_correction * recQQ->Px()+recBc_muW->Px() ,2) + pow(QQM_correction * recQQ->Py()+recBc_muW->Py() ,2) );
	      float Bc_P = sqrt( pow(Bc_Pt ,2) + pow(QQM_correction * recQQ->Pz()+recBc_muW->Pz() ,2) );
	      float PperpTrimu = TMath::Sin(Bc_alpha3D) * Bc_P;
	      float Bc_CorrM = sqrt(BcCandM*BcCandM + PperpTrimu*PperpTrimu) + PperpTrimu;

	      //sum of deltaR's
	      float dRjpsi = (float) recBc_mumi->DeltaR(*recBc_mupl);
	      float dRmuWmi = (float) recBc_muW->DeltaR(*recBc_mumi);
	      float dRmuWpl = (float) recBc_muW->DeltaR(*recBc_mupl);
	      vector<float> dRcorr = dRcorrectedForQQM(dRjpsi,dRmuWmi,dRmuWpl,QQM_correction,mumi_Pt,mupl_Pt);
	      float dR_jpsi = dRcorr[0];
	      float dR_muWmi = dRcorr[1];		
	      float dR_muWpl = dRcorr[2]; 
	      float dR_sum = dR_muWmi + dR_muWpl + dR_jpsi;

	      float dR_sum_QQ2;
	      TLorentzVector *recBc_muW2 = (TLorentzVector*) ((iIpt==4)?(
									 (Reco_mu_charge[iIpt][muWidx_3[k2]]>0)?(Reco_QQ_mupl_4mom[iIpt]->At(QQidx)):(Reco_QQ_mumi_4mom[iIpt]->At(QQidx))
									 ):(
									    Reco_mu_4mom[iIpt]->At(min(muWidx_3[k2],(Short_t)(Reco_mu_size[iIpt]-1))))
							      );
	      TLorentzVector *recBc_mumi2 = (TLorentzVector*) ((iIpt<4)?(Reco_mu_4mom[iIpt]->At(mumiidx_3[k2])
									 ):(
									    (Reco_mu_charge[iIpt][muWidx_3[k2]]<0)?(Reco_mu_4mom[iIpt]->At(muWidx_3[k2])):(Reco_QQ_mumi_4mom[iIpt]->At(QQidx)))
							       );
	      TLorentzVector *recBc_mupl2 = (TLorentzVector*) ((iIpt<4)?(Reco_mu_4mom[iIpt]->At(muplidx_3[k2])
									 ):(
									    (Reco_mu_charge[iIpt][muWidx_3[k2]]>0)?(Reco_mu_4mom[iIpt]->At(muWidx_3[k2])):(Reco_QQ_mupl_4mom[iIpt]->At(QQidx)))
							       );
	      if(!(iIpt==4 && (recBc_mumi2==NULL || recBc_mupl2==NULL))) {
		float dRjpsi2 = (float) recBc_mumi2->DeltaR(*recBc_mupl2);
		float dRmuWmi2 = (float) recBc_muW2->DeltaR(*recBc_mumi2);
		float dRmuWpl2 = (float) recBc_muW2->DeltaR(*recBc_mupl2);
		vector<float> dRcorr2 = dRcorrectedForQQM(dRjpsi2,dRmuWmi2,dRmuWpl2,QQ2M_correction,recBc_mumi2->Pt(),recBc_mupl2->Pt());
		dR_sum_QQ2 = dRcorr2[0] + dRcorr2[1] + dRcorr2[2];
	      }
	      else{
		dR_sum_QQ2 = dR_sum;
	      }

	      //mass shift for the TRUEJPSI (high mass control region) bkg
	      if(i==2){
		float BcM_correction = 1-1.2/BcCandM;
		BcCandM *= BcM_correction;
		Bc_CorrM = sqrt(BcCandM*BcCandM + PperpTrimu*PperpTrimu) + PperpTrimu;
		float dRmuWplmi = TMath::ACos(1- BcM_correction*BcM_correction* (1-TMath::Cos(dR_muWpl)) ) + 
		  TMath::ACos(1- BcM_correction*BcM_correction* (1-TMath::Cos(dR_muWmi)) );
		dR_sum = dRmuWplmi + dR_jpsi;
	      }

	      //cout<<(Bc_ctauSignif>1.5)<<" "<<(Bc_alpha<0.8)<<" "<<(Bc_alpha3D<0.8)<<" "<<(Bc_VtxProb>0.01)<<" "<<(QQ_dca<0.3)<<" "<<(muW_isSoft)<<" "<<( mumi_isSoft)<<" "<<(mupl_isSoft)<<" "<<( (muW_isGlb && muW_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) || (muW_isGlb && muW_inLooseAcc && mumi_isGlb && mumi_inLooseAcc) || (mumi_isGlb && mumi_inLooseAcc && mupl_isGlb && mupl_inLooseAcc) )<<" "<<(mumi_isGlb && mupl_isGlb && muW_isGlb && mumi_inLooseAcc && mupl_inLooseAcc && muW_inLooseAcc)<<" "<<muW_isGlb<<" "<<mupl_isGlb<<" "<<mumi_isGlb<<" "<<((muW_trig && mupl_trig && muW_inTightAcc && mupl_inTightAcc ) ||(muW_trig && mumi_trig && muW_inTightAcc && mumi_inTightAcc ) || (mumi_trig && mupl_trig && mumi_inTightAcc && mupl_inTightAcc ))<<" "<<(fabs(muW_dz)<(ispp?0.6:0.8) && fabs(mumi_dz)<(ispp?0.6:0.8) && fabs(mupl_dz)<(ispp?0.6:0.8))<<" "<<((HLTriggers[iIpt]&((ispp || i>=4)?8:4096))>0)<<endl;
	    
	      for(int ieff=0;ieff<effNum;ieff++){
		TString probeEff = effToTest[ieff];

		if(
		   //**************************************************************
		   // Pre-selections
		   muW_isSoft
		   && mumi_isSoft
		   && mupl_isSoft
		   && ( (muW_isGlb && muW_inLooseAcc && mupl_isGlb && mupl_inLooseAcc && looseAcc(mumi_Pt,mumi_eta,true) ) || //only one muon can be tracker and out of LooseAcceptance     
			(muW_isGlb && muW_inLooseAcc && mumi_isGlb && mumi_inLooseAcc && looseAcc(mupl_Pt,mupl_eta,true) ) ||
			(mumi_isGlb && mumi_inLooseAcc && mupl_isGlb && mupl_inLooseAcc && looseAcc(muW_Pt,muW_eta,true) )
			)
		   && ( ( muW_trig && mupl_trig && muW_inTightAcc && mupl_inTightAcc ) || //two muons among three must trigger //BEWARE ! Not sure if TightAcceptance should be put there
			( muW_trig && mumi_trig && muW_inTightAcc && mumi_inTightAcc ) ||
			( mumi_trig && mupl_trig && mumi_inTightAcc && mupl_inTightAcc ) //only this last option can be true for dimuon+trk
			)
		   && (!ispp || (HLTriggers[iIpt]&((ispp || i>=4)?8:4096))>0) //the event must fire the trigger as well //BEWARE ! this should be re-established for PbPb when samples are re-run
		   //!!!!!!!!!!!!!!!!!!!!!! Here, should reject the tracks that are global muons, for dimuon+track
		   && (ispp || INCentrality[iIpt]<180) //keep 0-90% centrality

		   && (_withTM || probeEff.Contains("3global") || (mumi_isGlb && mupl_isGlb && muW_isGlb
								   && mumi_inLooseAcc && mupl_inLooseAcc && muW_inLooseAcc))
		   && (probeEff.Contains("tauSignif2D") || Bc_ctauSignif>_ctauSignif_cut)
		   && (probeEff.Contains("tauSignif3D") || Bc_ctauSignif3D>_ctauSignif3D_cut) //maybe not a good idea in pp due to pile-up along z. But actually not a big deal if we select nonprompt objects
		   && (probeEff.Contains("alpha2D") || Bc_alpha<_alpha_cut(ispp))
		   && (probeEff.Contains("alpha3D") || Bc_alpha3D<_alpha3D_cut(ispp)) //combined with alpha cut, kills 1.4% of signal in pp (1.1% in PbPb) 
		   && (probeEff.Contains("TrimuVProb") || Bc_VtxProb>_vtxProb_cut) //kills 1.2% more signal compared to Bc_VtxProb>0.005
		   && (probeEff.Contains("DimuVProb") || QQ_VtxProb>_QQvtxProb_cut) //drop the QQ_VtxProb, too correlated with Bc_VtxProb?
		   && (probeEff.Contains("dca") || QQ_dca<_QQdca_cut) //keep this that kills 1.8% of signal and 10% of WRONSIGN/BCMASS
		   && (probeEff.Contains("dz") || (fabs(muW_dz)<0.6 && fabs(mumi_dz)<0.6 && fabs(mupl_dz)<0.6))
		   && (probeEff.Contains("CorrM") || Bc_CorrM<_BcCorrM_cut(ispp))
		   && (probeEff.Contains("dRsum") || dR_sum<_dRsum_cut)
		   ){
		  float QQ_M = QQCandM;
		  float QQ2_M = (Reco_mu_charge[iIpt][muWidx]>0)?((*recBc_mumi+*recBc_muW).M()):((*recBc_mupl+*recBc_muW).M()); //QQ2 is the second OS pair
		    
		  if(ieff==0){ //Jpsi choice and weight only at first pass, with full preselection
		    // if(iIpt==0){
		    //   //**** Deal with the Jpsi dimuon choice
		    //   if(i==2 && inJpsiMassRange(QQ_M, maxEta<1.5) && inJpsiMassRange(QQ2_M, maxEta<1.5) 
		    // 	 && (probeEff.Contains("dca") || QQ2_dca<QQ_dca) && (probeEff.Contains("DimuVProb") || QQ2_VtxProb>_QQvtxProb_cut))
		    // 	break; //In data, if both OS dimuons are in Jpsi peak region, choose the one with best dca (so no bias on the Jpsi mass)
		    //   if(i==1 && inJpsiMassSB(QQ_M, maxEta<1.5) && inJpsiMassSB(QQ2_M, maxEta<1.5)
		    // 	 && (probeEff.Contains("dca") || QQ2_dca<QQ_dca) && (probeEff.Contains("DimuVProb") || QQ2_VtxProb>_QQvtxProb_cut))
		    // 	break; //In data sidebands, if both OS dimuons are in Jpsi sidebands, choose the one with best dca (so no bias on the Jpsi mass)
		    // }

		    //**** Add the Jpsi choice weight
		    if((inJpsiMassSB(QQ2_M, maxEta<1.5) || inJpsiMassRange(QQ2_M, maxEta<1.5))
		       && (probeEff.Contains("dca") || QQ2_dca<QQ_dca) && (probeEff.Contains("DimuVProb") || QQ2_VtxProb>_QQvtxProb_cut) && (probeEff.Contains("dRsum") || dR_sum_QQ2>_dRsum_cut)
		       ){
		      //*** loose SB           
		      float binc_QQ1 = JpsiM->GetBinContent(JpsiM->FindBin(QQ_M));
		      float binc_QQ2 = JpsiM->GetBinContent(JpsiM->FindBin(QQ2_M));

		      //*** tight SB
		      if(maxEta<1.5){
			binc_QQ1 = JpsiM_tight->GetBinContent(JpsiM_tight->FindBin(QQ_M));
			binc_QQ2 = JpsiM_tight->GetBinContent(JpsiM_tight->FindBin(QQ2_M));
		      }

		      if (binc_QQ1==0) weight = 0;
		      else weight *= binc_QQ1 / (binc_QQ1+binc_QQ2);
		      //cout<<"Jpsi choice weight = "<<binc_QQ1 / (binc_QQ1+binc_QQ2)<<endl;
		    }
		  }

		  //cout<<"adding this weight to ntot: "<<weight<<endl;
		  ntot[i][ieff] += weight;

		  if(
		     (_withTM || (mumi_isGlb && mupl_isGlb && muW_isGlb
				  && mumi_inLooseAcc && mupl_inLooseAcc && muW_inLooseAcc))
		     && Bc_ctauSignif>_ctauSignif_cut
		     && Bc_ctauSignif3D>_ctauSignif3D_cut //maybe not a good idea in pp due to pile-up along z. But actually not a big deal if we select nonprompt objects
		     && Bc_alpha<_alpha_cut(ispp)
		     && Bc_alpha3D<_alpha3D_cut(ispp) //combined with alpha cut, kills 1.4% of signal in pp (1.1% in PbPb) 
		     && Bc_VtxProb>_vtxProb_cut //kills 1.2% more signal compared to Bc_VtxProb>0.005
		     && QQ_VtxProb>_QQvtxProb_cut //drop the QQ_VtxProb, too correlated with Bc_VtxProb?
		     && QQ_dca<_QQdca_cut //keep this that kills 1.8% of signal and 10% of WRONSIGN/BCMASS
		     && fabs(muW_dz)<0.6 && fabs(mumi_dz)<0.6 && fabs(mupl_dz)<0.6
		     && Bc_CorrM<_BcCorrM_cut(ispp)
		     && dR_sum<_dRsum_cut
		     ){
		    npass[i][ieff] += weight;
		  }
		
		} //end if passes full selection
	      }//end efficiency to test
	    } //end loop on checked samples 
	  } //end if(QQ_isValid)
	} //end loop on 2/3 possible Jpsi dimuon choice
      } //end loop on Bc candidates
    } //end loop on entries
    // cout<<"Jpsi mass signal region nall,npass,efficiency = "<<nall<<" "<<npass<<" "<<npass/nall<<endl;
    // cout<<"Jpsi mass sidebands nall,npass,efficiency = "<<nall<<" "<<nside<<" "<<nside/nall<<endl;
  } //end loop on input trees

  ofstream outfile;
  outfile.open("Nmin1efficiencies_"+(TString)(ispp?"pp":"PbPb")+".txt");
  for(int isam=0;isam<samNum;isam++){
    outfile<<endl<<"sample "<<samName[isam]<<endl;
    for(int ieff=0;ieff<effNum;ieff++){
      //outfile<<npass[isam][ieff]<<" "<<ntot[isam][ieff]<<endl;
      outfile<<"inefficiency of cut "<<effToTest[ieff]<<" = "<< (1-npass[isam][ieff]/ntot[isam][ieff]) <<endl;
    }
  }
  outfile.close();

}
