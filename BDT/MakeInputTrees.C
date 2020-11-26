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
#include "../acceptance/SgMuonAcceptanceCuts.h"

void MakeInputTrees(bool ispp = true){

  //**************************************************************  
  //Create Tree and branches
  vector<TTree*> Trees; vector<int> nentries;
  
  TFile *fileData = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/DoubleMu_Run2017G_AOD_Run_306546_306826_OniaTree_TripleMuBc_07082019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/data/HIDoubleMuon_HIDoubleMuonPsiPeri_Run2018A_AOD_OniaTree_Run_326381_327564_BcTrimuon_16122019.root","READ");
  Trees.push_back( (TTree*)fileData->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[0]->GetEntries() );
  std::cout<<"nevents data = "<<nentries[0]<<"\n";
  //  Trees[0]->Print();

  TFile *fileMC = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/Oniatree_MC_Bc_trimuons_21112019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/BcToJpsiMuNu_BCVEGPY_PYTHIA8_2018PbPb5TeV_HINPbPbAutumn18DR-00196_08092020_4200k_ONIATREE.root","READ");
  Trees.push_back( (TTree*)fileMC->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[1]->GetEntries() );
  std::cout<<"nevents MC = "<<nentries[1]<<"\n";
  //Trees[1]->Print();

  TFile *fileMCb = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/NonPromptJpsi/MConiatree/crab_BJPsiMM_TuneCUETP8M1_5p02TeV_pythia8_05082019_wLambdabFor10_ptHatMinCombined_ONIATREE.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/NonPromptJpsi/BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_trimuons_oniatree_09012020.root","READ");
  Trees.push_back( (TTree*)fileMCb->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[2]->GetEntries() );
  std::cout<<"nevents B->J/psi MC = "<<nentries[2]<<"\n";
  //Trees[2]->Print();

  TFile *fileMCpromptJ = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/PromptJpsi/Oniatree_MC_trimuons_PromptJpsi_ptHatMinCombined_05082019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/PromptJpsi/Jpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_trimuons_oniatree_09012020.root","READ");
  Trees.push_back( (TTree*)fileMCpromptJ->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[3]->GetEntries() );
  std::cout<<"nevents Prompt J/psi MC = "<<nentries[3]<<"\n";
  //Trees[3]->Print();

  TFile *fileDimuTrk = TFile::Open("/data_CMS/cms/falmagne/tuples/pp17/Bc/DimuonTrk/forTrimuon/DoubleMu_Run2017G_AOD_Run_306546_306826_OniaTree_DimuonTrkBc_TrimuonWay_21112019.root","READ");
  Trees.push_back( (TTree*)fileDimuTrk->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[4]->GetEntries() );
  std::cout<<"nevents dimuon+track data = "<<nentries[4]<<"\n";
  //Trees[4]->Print();

  TFile *fileFlipJpsi = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/flipJpsi/DoubleMu_Run2017G_AOD_Run_306546_306826_OniaTree_TripleMuBc_flippedJpsi_13112019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/flippedJpsi/HIDoubleMuon_HIDoubleMuonPsiPeri_Run2018A_AOD_OniaTree_Run_326381_327564_flippedJpsi_BcTrimuon_27122019.root","READ");
  Trees.push_back( (TTree*)fileFlipJpsi->Get("hionia/myTree") );
  nentries.push_back( (int)Trees[5]->GetEntries() );
  std::cout<<"nevents flipped-Jpsi data = "<<nentries[5]<<"\n";
  //  Trees[5]->Print();

  if(ispp){
    TFile *fileFlipJpsibMC = TFile::Open("/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/NonPromptJpsi/flipJpsi/BJPsiMM_TuneCUETP8M1_5p02TeV_pythia8_09092020_ptHatMin2_ONIATREE_flipJpsi.root","READ");
    Trees.push_back( (TTree*)fileFlipJpsibMC->Get("hionia/myTree") );
    nentries.push_back( (int)Trees[6]->GetEntries() );
    std::cout<<"nevents flipped-Jpsi b->Jpsi MC = "<<nentries[6]<<"\n";
    //  Trees[6]->Print();
  }

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

    if(i==0 || i==4 || i==5){//data only
      b_INrunNb[i] = Trees[i]->GetBranch("runNb");
      b_INrunNb[i]->SetAddress(&INrunNb[i]);

      b_INLS[i] = Trees[i]->GetBranch("LS");
      b_INLS[i]->SetAddress(&INLS[i]);
    }

    if(ispp){
      b_INnPV[i] = Trees[i]->GetBranch("nPV");
      b_INnPV[i]->SetAddress(&INnPV[i]);
    } else if(i!=4) {
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

    if(i<4){
      b_Reco_3mu_QQ2_idx[i] = Trees[i]->GetBranch("Reco_3mu_QQ2_idx");
      b_Reco_3mu_QQ2_idx[i]->SetAddress(&Reco_3mu_QQ2_idx[i]);

      b_Reco_3mu_QQss_idx[i] = Trees[i]->GetBranch("Reco_3mu_QQss_idx");
      b_Reco_3mu_QQss_idx[i]->SetAddress(&Reco_3mu_QQss_idx[i]);

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

    b_Reco_mu_normChi2_global[i] = Trees[i]->GetBranch("Reco_mu_normChi2_global");
    b_Reco_mu_normChi2_global[i]->SetAddress(&Reco_mu_normChi2_global[i]);

    b_Reco_mu_normChi2_inner[i] = Trees[i]->GetBranch("Reco_mu_normChi2_inner");
    b_Reco_mu_normChi2_inner[i]->SetAddress(&Reco_mu_normChi2_inner[i]);

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

    if(ispp || i==5){
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

    if(i==1){//MC signal only
      Gen_3mu_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Gen_3mu_4mom", &Gen_3mu_4mom[i], &b_Gen_3mu_4mom[i]);
      
      Gen_Bc_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Gen_Bc_4mom", &Gen_Bc_4mom[i], &b_Gen_Bc_4mom[i]);
    }

    if(i==4){
      Reco_trk_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Reco_trk_4mom", &Reco_trk_4mom[i], &b_Reco_trk_4mom[i]);

      b_Reco_trk_InLooseAcc[i] = Trees[i]->GetBranch("Reco_trk_InLooseAcc");
      b_Reco_trk_InLooseAcc[i]->SetAddress(&Reco_trk_InLooseAcc[i]);

      b_Reco_trk_InTightAcc[i] = Trees[i]->GetBranch("Reco_trk_InTightAcc");
      b_Reco_trk_InTightAcc[i]->SetAddress(&Reco_trk_InTightAcc[i]);

      b_Reco_trk_dxy[i] = Trees[i]->GetBranch("Reco_trk_dxy");
      b_Reco_trk_dxy[i]->SetAddress(&Reco_trk_dxy[i]);

      b_Reco_trk_dz[i] = Trees[i]->GetBranch("Reco_trk_dz");
      b_Reco_trk_dz[i]->SetAddress(&Reco_trk_dz[i]);

      b_Reco_trk_dxyError[i] = Trees[i]->GetBranch("Reco_trk_dxyError");
      b_Reco_trk_dxyError[i]->SetAddress(&Reco_trk_dxyError[i]);

      b_Reco_trk_dzError[i] = Trees[i]->GetBranch("Reco_trk_dzError");
      b_Reco_trk_dzError[i]->SetAddress(&Reco_trk_dzError[i]);
    }

    Reco_mu_4mom[i] = new TClonesArray();
    Trees[i]->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom[i], &b_Reco_mu_4mom[i]);

    Reco_QQ_4mom[i] = new TClonesArray();
    Trees[i]->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom[i], &b_Reco_QQ_4mom[i]);

    b_Reco_QQ_size[i] = Trees[i]->GetBranch("Reco_QQ_size");
    b_Reco_QQ_size[i]->SetAddress(&Reco_QQ_size[i]);

    b_Reco_mu_size[i] = Trees[i]->GetBranch("Reco_mu_size");
    b_Reco_mu_size[i]->SetAddress(&Reco_mu_size[i]);

    if((i>0 && i<4) || i>=6){
      Gen_QQ_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom[i], &b_Gen_QQ_4mom[i]);

      b_Gen_QQ_momId[i] = Trees[i]->GetBranch("Gen_QQ_momId");
      b_Gen_QQ_momId[i]->SetAddress(&Gen_QQ_momId[i]);

      b_Reco_QQ_whichGen[i] = Trees[i]->GetBranch("Reco_QQ_whichGen");
      b_Reco_QQ_whichGen[i]->SetAddress(&Reco_QQ_whichGen[i]);

      if(i==6){
	b_Gen_QQ_size[i] = Trees[i]->GetBranch("Gen_QQ_size");
	b_Gen_QQ_size[i]->SetAddress(&Gen_QQ_size[i]);
      }

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

      b_Reco_3mu_muW_trueId[i] = Trees[i]->GetBranch("Reco_3mu_muW_trueId");
      b_Reco_3mu_muW_trueId[i]->SetAddress(&Reco_3mu_muW_trueId[i]);

      b_Reco_3mu_whichGen[i] = Trees[i]->GetBranch("Reco_3mu_whichGen");
      b_Reco_3mu_whichGen[i]->SetAddress(&Reco_3mu_whichGen[i]);
    }

    if(i==5 || i==6){
      b_Reco_QQ_flipJpsi[i] = Trees[i]->GetBranch("Reco_QQ_flipJpsi");
      b_Reco_QQ_flipJpsi[i]->SetAddress(&Reco_QQ_flipJpsi[i]);

      Reco_QQ_mumi_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom[i], &b_Reco_QQ_mumi_4mom[i]);

      Reco_QQ_mupl_4mom[i] = new TClonesArray();
      Trees[i]->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom[i], &b_Reco_QQ_mupl_4mom[i]);
    }
  }

  //**************************************************************
  //Create the output file and trees
  TFile out_file( "BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+(TString)(_withTM?"_withTM":"")+".root","RECREATE");
  std::vector<TTree*> out_trees;
  out_trees.push_back(new TTree("bkgWRONGSIGN","tree with wrongsign background"));
  out_trees.push_back(new TTree("bkgBCMASS","tree with background from Jpsi mass sidebands"));
  out_trees.push_back(new TTree("bkgTRUEJPSI","tree with background from true Jpsi but wrong trimuon mass"));
  out_trees.push_back(new TTree("sigRegion","tree with data signal region"));
  out_trees.push_back(new TTree("signal_MC","tree with signal from MC"));
  out_trees.push_back(new TTree("bToJpsi_MC","MC tree with Jpsi from b-decays"));
  out_trees.push_back(new TTree("PromptJpsi_MC","MC tree with Prompt Jpsi"));
  out_trees.push_back(new TTree("dimuonTrk","tree with dimuon-track reconstructed the trimuon-way (muW is a track)"));
  out_trees.push_back(new TTree("flipJpsi","tree with trimuon formed from a Jpsi whose PV-SV and momentum have been flipped"));
  if(ispp) out_trees.push_back(new TTree("flipJpsibMC","tree with trimuon formed from a Jpsi whose PV-SV and momentum have been flipped, from the b->Jpsi MC"));

  //**************************************************************
  //Create the ouput branches   
  int ntrees = out_trees.size();
  UInt_t eventNb[ntrees];
  UInt_t runNb[ntrees];
  UInt_t LS[ntrees];
  int nPV[ntrees];
  int Centrality[ntrees];
  float Bc_M[ntrees];
  int Bc_charge[ntrees];
  float Bc_CorrM[ntrees];
  //  float Bc_M_muWismu[ntrees];
  float Bc_M_withQQM[ntrees];
  float Bc_M_muWisK[ntrees];
  //float Bc_M_muWisPi[ntrees];
  float Bc_P[ntrees];
  float Bc_Pt[ntrees];
  float genBc_Pt[ntrees];
  float gen3mu_Pt[ntrees];
  float Bc_Y[ntrees];
  float Bc_phi[ntrees];
  float Bc_ctau[ntrees];
  float Bc_ctauSignif[ntrees];
  float Bc_ctau3D[ntrees];
  float Bc_ctauSignif3D[ntrees];
  float Bc_alpha[ntrees];
  float Bc_alpha3D[ntrees];
  float Bc_VtxProb[ntrees];

  float Bc_M_shiftedM[ntrees];
  float Bc_CorrM_shiftedM[ntrees];
  float dR_sum_shiftedM[ntrees];
  float dR_jpsiOverMuW_shiftedM[ntrees];
  float dR_jpsiMuW_shiftedM[ntrees];

  float QQ_M[ntrees];
  float QQ_P[ntrees];
  float QQ_Pt[ntrees];
  float QQ_Y[ntrees];
  float QQ_phi[ntrees];
  float QQ2_M[ntrees];
  float QQ3_M[ntrees];
  float QQ_VtxProb[ntrees];
  float QQ_dca[ntrees];
  float QQ2_VtxProb[ntrees];
  float QQ2_dca[ntrees];
  int QQ_charge[ntrees];
  float muW_P[ntrees];
  float muW_Pt[ntrees];
  float muW_eta[ntrees];
  float muW_phi[ntrees];
  float mumi_P[ntrees];
  float mumi_Pt[ntrees];
  float mumi_eta[ntrees];
  float mumi_phi[ntrees];
  float mupl_P[ntrees];
  float mupl_Pt[ntrees];
  float mupl_eta[ntrees];
  float mupl_phi[ntrees];
  bool muW_isJpsiBro[ntrees];
  int muW_trueId[ntrees];
  bool muW_inLooseAcc[ntrees];
  bool mumi_inLooseAcc[ntrees];
  bool mupl_inLooseAcc[ntrees];
  bool muW_inTightAcc[ntrees];
  bool mumi_inTightAcc[ntrees];
  bool mupl_inTightAcc[ntrees];
  bool muW_isGlb[ntrees];
  bool mumi_isGlb[ntrees];
  bool mupl_isGlb[ntrees];
  bool muW_isSoft[ntrees];
  bool mumi_isSoft[ntrees];
  bool mupl_isSoft[ntrees];
  bool muW_trig[ntrees];
  bool mumi_trig[ntrees];
  bool mupl_trig[ntrees];
  float muW_normChi2_inner[ntrees];
  float muW_normChi2_glb[ntrees];

  float dR_muWmi[ntrees];
  float dR_muWpl[ntrees];
  float dR_jpsi[ntrees];
  float dR_sum[ntrees];
  float dR_jpsiMuW[ntrees];
  float dR_jpsiOverMuW[ntrees];
  float MuonDxySignif_sum[ntrees];
  float muW_dxySignif[ntrees];
  float mumi_dxySignif[ntrees];
  float mupl_dxySignif[ntrees];
  float MuonDzSignif_sum[ntrees];
  float muW_dzSignif[ntrees];
  float mumi_dzSignif[ntrees];
  float mupl_dzSignif[ntrees];
  float QQmuW_ptImbal[ntrees];
  int bkgType[ntrees];
  int QQ_momId[ntrees];
  float w_simple[ntrees];
  float w_unblind[ntrees];
  //  float weightJpsiChoice[ntrees];
  float weightNcoll[ntrees];
  int flipJpsi[ntrees];

  for(int i=0;i<ntrees;i++){
    out_trees[i]->Branch("eventNb",&eventNb[i],"eventNb/i");
    out_trees[i]->Branch("runNb",&runNb[i],"runNb/i");
    out_trees[i]->Branch("LS",&LS[i],"LS/i");
    if(ispp){
      out_trees[i]->Branch("nPV",&nPV[i],"nPV/I");
    }else{
      out_trees[i]->Branch("Centrality",&Centrality[i],"Centrality/I");
    }
    out_trees[i]->Branch("Bc_M",&Bc_M[i],"Bc_M/F");
    out_trees[i]->Branch("Bc_charge",&Bc_charge[i],"Bc_charge/I");
    out_trees[i]->Branch("Bc_CorrM",&Bc_CorrM[i],"Bc_CorrM/F");
    out_trees[i]->Branch("Bc_M_withQQM",&Bc_M_withQQM[i],"Bc_M_withQQM/F");
    out_trees[i]->Branch("Bc_M_muWisK",&Bc_M_muWisK[i],"Bc_M_muWisK/F");
    //    out_trees[i]->Branch("Bc_M_muWisPi",&Bc_M_muWisPi[i],"Bc_M_muWisPi/F");
    //out_trees[i]->Branch("Bc_M_muWismu",&Bc_M_muWismu[i],"Bc_M_muWismu/F");
    out_trees[i]->Branch("Bc_P",&Bc_P[i],"Bc_P/F");
    out_trees[i]->Branch("Bc_Pt",&Bc_Pt[i],"Bc_Pt/F");
    if(i==4){
      out_trees[i]->Branch("genBc_Pt",&genBc_Pt[i],"genBc_Pt/F");
      out_trees[i]->Branch("gen3mu_Pt",&gen3mu_Pt[i],"gen3mu_Pt/F");
    }
    out_trees[i]->Branch("Bc_Y",&Bc_Y[i],"Bc_Y/F");
    out_trees[i]->Branch("Bc_phi",&Bc_phi[i],"Bc_phi/F");
    out_trees[i]->Branch("Bc_ctau",&Bc_ctau[i],"Bc_ctau/F");
    out_trees[i]->Branch("Bc_ctauSignif",&Bc_ctauSignif[i],"Bc_ctauSignif/F");
    out_trees[i]->Branch("Bc_ctau3D",&Bc_ctau3D[i],"Bc_ctau3D/F");
    out_trees[i]->Branch("Bc_ctauSignif3D",&Bc_ctauSignif3D[i],"Bc_ctauSignif3D/F");
    out_trees[i]->Branch("Bc_alpha",&Bc_alpha[i],"Bc_alpha/F");
    out_trees[i]->Branch("Bc_alpha3D",&Bc_alpha3D[i],"Bc_alpha3D/F");
    out_trees[i]->Branch("Bc_VtxProb",&Bc_VtxProb[i],"Bc_VtxProb/F");

    out_trees[i]->Branch("QQ_M",&QQ_M[i],"QQ_M/F");
    out_trees[i]->Branch("QQ_P",&QQ_P[i],"QQ_P/F");
    out_trees[i]->Branch("QQ_Pt",&QQ_Pt[i],"QQ_Pt/F");
    out_trees[i]->Branch("QQ_Y",&QQ_Y[i],"QQ_Y/F");
    out_trees[i]->Branch("QQ_phi",&QQ_phi[i],"QQ_phi/F");
    out_trees[i]->Branch("QQ_VtxProb",&QQ_VtxProb[i],"QQ_VtxProb/F");
    out_trees[i]->Branch("QQ_dca",&QQ_dca[i],"QQ_dca/F");
    out_trees[i]->Branch("QQ_charge",&QQ_charge[i],"QQ_charge/I");
    out_trees[i]->Branch("QQ_momId",&QQ_momId[i],"QQ_momId/I");
    out_trees[i]->Branch("QQmuW_ptImbal",&QQmuW_ptImbal[i],"QQmuW_ptImbal/F");
    if(i!=7){
      out_trees[i]->Branch("QQ2_M",&QQ2_M[i],"QQ2_M/F");
      out_trees[i]->Branch("QQ3_M",&QQ3_M[i],"QQ3_M/F");
      out_trees[i]->Branch("QQ2_VtxProb",&QQ2_VtxProb[i],"QQ2_VtxProb/F");
      out_trees[i]->Branch("QQ2_dca",&QQ2_dca[i],"QQ2_dca/F");
    }

    out_trees[i]->Branch("muW_P",&muW_P[i],"muW_P/F");
    out_trees[i]->Branch("muW_Pt",&muW_Pt[i],"muW_Pt/F");
    out_trees[i]->Branch("muW_eta",&muW_eta[i],"muW_eta/F");
    out_trees[i]->Branch("muW_phi",&muW_phi[i],"muW_phi/F");
    out_trees[i]->Branch("mumi_P",&mumi_P[i],"mumi_P/F");
    out_trees[i]->Branch("mumi_Pt",&mumi_Pt[i],"mumi_Pt/F");
    out_trees[i]->Branch("mumi_eta",&mumi_eta[i],"mumi_eta/F");
    out_trees[i]->Branch("mumi_phi",&mumi_phi[i],"mumi_phi/F");
    out_trees[i]->Branch("mupl_P",&mupl_P[i],"mupl_P/F");
    out_trees[i]->Branch("mupl_Pt",&mupl_Pt[i],"mupl_Pt/F");
    out_trees[i]->Branch("mupl_eta",&mupl_eta[i],"mupl_eta/F");
    out_trees[i]->Branch("mupl_phi",&mupl_phi[i],"mupl_phi/F");
    out_trees[i]->Branch("muW_isJpsiBro",&muW_isJpsiBro[i],"muW_isJpsiBro/O");
    out_trees[i]->Branch("muW_trueId",&muW_trueId[i],"muW_trueId/I");
    out_trees[i]->Branch("muW_normChi2_inner",&muW_normChi2_inner[i],"muW_normChi2_inner/F");
    out_trees[i]->Branch("muW_normChi2_glb",&muW_normChi2_glb[i],"muW_normChi2_glb/F");
    out_trees[i]->Branch("muW_inLooseAcc",&muW_inLooseAcc[i],"muW_inLooseAcc/O");
    out_trees[i]->Branch("mupl_inLooseAcc",&mupl_inLooseAcc[i],"mupl_inLooseAcc/O");
    out_trees[i]->Branch("mumi_inLooseAcc",&mumi_inLooseAcc[i],"mumi_inLooseAcc/O");
    out_trees[i]->Branch("muW_inTightAcc",&muW_inTightAcc[i],"muW_inTightAcc/O");
    out_trees[i]->Branch("mupl_inTightAcc",&mupl_inTightAcc[i],"mupl_inTightAcc/O");
    out_trees[i]->Branch("mumi_inTightAcc",&mumi_inTightAcc[i],"mumi_inTightAcc/O");
    out_trees[i]->Branch("muW_isGlb",&muW_isGlb[i],"muW_isGlb/O");
    out_trees[i]->Branch("mupl_isGlb",&mupl_isGlb[i],"mupl_isGlb/O");
    out_trees[i]->Branch("mumi_isGlb",&mumi_isGlb[i],"mumi_isGlb/O");
    out_trees[i]->Branch("muW_isSoft",&muW_isSoft[i],"muW_isSoft/O");
    out_trees[i]->Branch("mupl_isSoft",&mupl_isSoft[i],"mupl_isSoft/O");
    out_trees[i]->Branch("mumi_isSoft",&mumi_isSoft[i],"mumi_isSoft/O");
    out_trees[i]->Branch("muW_trig",&muW_trig[i],"muW_trig/O");
    out_trees[i]->Branch("mupl_trig",&mupl_trig[i],"mupl_trig/O");
    out_trees[i]->Branch("mumi_trig",&mumi_trig[i],"mumi_trig/O");

    out_trees[i]->Branch("dR_muWmi",&dR_muWmi[i],"dR_muWmi/F");
    out_trees[i]->Branch("dR_muWpl",&dR_muWpl[i],"dR_muWpl/F");
    out_trees[i]->Branch("dR_jpsi",&dR_jpsi[i],"dR_jpsi/F");
    out_trees[i]->Branch("dR_jpsiMuW",&dR_jpsiMuW[i],"dR_jpsiMuW/F");
    out_trees[i]->Branch("dR_sum",&dR_sum[i],"dR_sum/F");
    out_trees[i]->Branch("dR_jpsiOverMuW",&dR_jpsiOverMuW[i],"dR_jpsiOverMuW/F");
    out_trees[i]->Branch("MuonDxySignif_sum",&MuonDxySignif_sum[i],"MuonDxySignif_sum/F");
    out_trees[i]->Branch("muW_dxySignif",&muW_dxySignif[i],"muW_dxySignif/F");
    out_trees[i]->Branch("mumi_dxySignif",&mumi_dxySignif[i],"mumi_dxySignif/F");
    out_trees[i]->Branch("mupl_dxySignif",&mupl_dxySignif[i],"mupl_dxySignif/F");
    out_trees[i]->Branch("MuonDzSignif_sum",&MuonDzSignif_sum[i],"MuonDzSignif_sum/F");
    out_trees[i]->Branch("muW_dzSignif",&muW_dzSignif[i],"muW_dzSignif/F");
    out_trees[i]->Branch("mumi_dzSignif",&mumi_dzSignif[i],"mumi_dzSignif/F");
    out_trees[i]->Branch("mupl_dzSignif",&mupl_dzSignif[i],"mupl_dzSignif/F");

    //shifted variables for TRUEJPSI in BDT training
    out_trees[i]->Branch("Bc_M_shiftedM",&Bc_M_shiftedM[i],"Bc_M_shiftedM/F");
    out_trees[i]->Branch("Bc_CorrM_shiftedM",&Bc_CorrM_shiftedM[i],"Bc_CorrM_shiftedM/F");
    out_trees[i]->Branch("dR_sum_shiftedM",&dR_sum_shiftedM[i],"dR_sum_shiftedM/F");
    out_trees[i]->Branch("dR_jpsiOverMuW_shiftedM",&dR_jpsiOverMuW_shiftedM[i],"dR_jpsiOverMuW_shiftedM/F");
    out_trees[i]->Branch("dR_jpsiMuW_shiftedM",&dR_jpsiMuW_shiftedM[i],"dR_jpsiMuW_shiftedM/F");
    
    out_trees[i]->Branch("bkgType",&bkgType[i],"bkgType/I");
    out_trees[i]->Branch("w_simple",&w_simple[i],"w_simple/F");
    if(!ispp && i==3) out_trees[i]->Branch("w_unblind",&w_unblind[i],"w_unblind/F");
    //    out_trees[i]->Branch("weightJpsiChoice",&weightJpsiChoice[i],"weightJpsiChoice/F");
    if(!ispp && (i>3 && i<7) )
      out_trees[i]->Branch("weightNcoll",&weightNcoll[i],"weightNcoll/F");
    if(i>=8)
      out_trees[i]->Branch("flipJpsi",&flipJpsi[i],"flipJpsi/I");
  }

  //**************************************************************
  //Some variables and XS corrections
  std::map<bool, float> scaleMCsig = {{ true,  304800 * 2.54e0 * 0.668 / 3000000}, // Lumi_pp[nb-1] (from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM) * (XS_Bc_pp * BF((J/psi -> mu mu) mu nu))[nb] * (XS(5.02 TeV) / XS(7 TeV)) / nevents(uncut MC sample)
				      { false, 1.0 * 1.6054 * 2.54e0 * 0.668 * (7656 / 67.6) / 4200000} }; //Lumi_PbPb[nb-1] (from https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/948.html) * (XS_Bc_pp * BF((J/psi -> mu mu) mu nu))[nb] * (XS(5.02 TeV) / XS(7 TeV)) * ( XS^geom_PbPb / XS_Nucleon-Nucleon ) / nevents(uncut MC sample)
                        //weighted by Ncoll(centrality of given event) later, with an average Ncoll_MB = 382. The value of XS^geom was set to A^2 * XS_NN / Ncoll_MB, where XS_NN is taken from Glauber MC d'Enterria PRC 97.054910
                        //Assuming R_AA(Bc)=1

  TFile* fXS = TFile::Open("/home/llr/cms/falmagne/Bc/MCnormalization/NonPromptJpsiXS_scalefactor.root");
  TF1* SFpp = (TF1*)fXS->Get("SFpp");
  std::pair<TF1*,float> scalepp (SFpp, 304800 * 0.06/1000); // SF(data/MC)(Non prompt Jpsi XS) * Lumi_pp[nb-1] (preliminary, from https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/794.html) * BF(Jpsi-> mu mu) / MC in pb->nb 
  TF1* SFPbPb = (TF1*)fXS->Get("SFPbPb");
  std::pair<TF1*,float> scalePbPb (SFPbPb, 208*208 * 1.6054 * 0.06/1000);// SF(data/MC)(Non prompt Jpsi XS PbPb) * A^2 * Lumi_PbPb[nb-1] * MC in pb->nb / BF(Jpsi-> mu mu)
  std::map<bool, std::pair<TF1*,float> > scaleMCb = {{ true, scalepp }, 
						     { false, scalePbPb } }; 

  TFile* fXSprompt = TFile::Open("/home/llr/cms/falmagne/Bc/MCnormalization/PromptJpsiXS_scalefactor.root");
  TF1* SFppPrompt = (TF1*)fXSprompt->Get("SFpp");
  std::pair<TF1*,float> scaleppPrompt (SFppPrompt, 304800 * 0.06/1000); // SF(data/MC)(Non prompt Jpsi XS) * Lumi_pp[nb-1] (preliminary, from https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/794.html) * BF(Jpsi-> mu mu) / MC in pb->nb 
  TF1* SFPbPbPrompt = (TF1*)fXSprompt->Get("SFPbPb");
  std::pair<TF1*,float> scalePbPbPrompt (SFPbPbPrompt, 208*208 * 1.6054 * 0.06/1000);// SF(data/MC)(Non prompt Jpsi XS PbPb) * A^2 * Lumi_PbPb[nb-1] * BF(Jpsi-> mu mu) / MC in pb->nb
  std::map<bool, std::pair<TF1*,float> > scaleMCprompt = {{ true, scaleppPrompt }, 
							  { false, scalePbPbPrompt } }; 
  float nall=0,npass=0,nside=0;

  //**************************************************************
  //For each input tree, redirect data/MC events to the right tree, and write into branches
  bool goodTree = false;

  for(int iIpt=0; iIpt<nInputT; iIpt++){//nInputT
    if(!ispp && iIpt==4) {
      out_trees[7]->Fill(); //dummy 1 event
      continue;} //no dimuon+track sample yet in PbPb

    for(int j=0; j<nentries[iIpt]; j++){//nentries[iIpt]

      if(j%200000==0){ cout<<"Scanned "<<100.*(double)j/nentries[iIpt]<<"% of tree #"<<iIpt<<endl; }

      b_Reco_3mu_size[iIpt]->GetEntry(j);
      if(Reco_3mu_size[iIpt]==0) continue;

      //Certain trees are filled only from certain files
      int ntreestart = 0; int ntreestop = 3;
      if(iIpt==1) {ntreestart = 4; ntreestop = 5;} 
      else if(iIpt>=2) {ntreestart = 3+iIpt; ntreestop = ntreestart;}

      for(int i=ntreestart; i<=ntreestop; i++){ //only trees with data events for iIpt=0, MC events for iIpt=1, b->Jpsi MC events for iIpt=2, prompt Jpsi MC events for iIpt=3, dimuon+track events for iIpt=4

	Reco_3mu_4mom[iIpt]->Clear();
	Reco_mu_4mom[iIpt]->Clear();
	if(iIpt==4) Reco_trk_4mom[iIpt]->Clear();
	Reco_QQ_4mom[iIpt]->Clear();
	if(iIpt>=5){
	  Reco_QQ_mumi_4mom[iIpt]->Clear();
	  Reco_QQ_mupl_4mom[iIpt]->Clear();}
	if(iIpt==1){
	  Gen_Bc_4mom[iIpt]->Clear();
	  Gen_3mu_4mom[iIpt]->Clear();
	}

	Trees[iIpt]->GetEntry(j);

	for(int BcNb=0;BcNb<Reco_3mu_size[iIpt];BcNb++){
	  Short_t genBcIdx = 0;
	  if(iIpt==1) genBcIdx = Reco_3mu_whichGen[iIpt][BcNb];

	  Short_t QQidx_3[3] = { Reco_3mu_QQ1_idx[iIpt][BcNb], Reco_3mu_QQ2_idx[iIpt][BcNb], Reco_3mu_QQss_idx[iIpt][BcNb] };
	  if(fabs(Reco_3mu_charge[iIpt][BcNb])==3 && iIpt!=4){ //This fix for wrongsign can be removed after rerunning the oniatree
	    bool foundit = false;
	    if(Reco_3mu_mumi_idx[iIpt][BcNb]!=Reco_QQ_mumi_idx[iIpt][Reco_3mu_QQ2_idx[iIpt][BcNb]] && Reco_3mu_mumi_idx[iIpt][BcNb]!=Reco_QQ_mupl_idx[iIpt][Reco_3mu_QQ2_idx[iIpt][BcNb]]) {
	      Reco_3mu_muW2_idx[iIpt][BcNb] = Reco_3mu_mumi_idx[iIpt][BcNb]; foundit=true;
	    }
	    if(Reco_3mu_mupl_idx[iIpt][BcNb]!=Reco_QQ_mumi_idx[iIpt][Reco_3mu_QQ2_idx[iIpt][BcNb]] && Reco_3mu_mupl_idx[iIpt][BcNb]!=Reco_QQ_mupl_idx[iIpt][Reco_3mu_QQ2_idx[iIpt][BcNb]]) {
	      Reco_3mu_muW2_idx[iIpt][BcNb] = Reco_3mu_mupl_idx[iIpt][BcNb]; foundit=true;
	    }
	    if(!foundit) cout<<"Did not find muW for wrongsign Bc"<<endl;
	  }
	  Short_t muW3idx = (Reco_3mu_mumi_idx[iIpt][BcNb]!=Reco_3mu_muW2_idx[iIpt][BcNb])?(Reco_3mu_mumi_idx[iIpt][BcNb]):(Reco_3mu_mupl_idx[iIpt][BcNb]); //for the third dimuon choice, the muWidx=old_muWmi if old_muWmi!=muW2idx, and muWidx=old_muWpl otherwise
	  if (iIpt!=4 && (muW3idx==Reco_3mu_muW_idx[iIpt][BcNb] || muW3idx==Reco_3mu_muW2_idx[iIpt][BcNb])) cout<<"!!!!!!!!! Wrong assignement of muW3idx !"<<endl;
	  Short_t muWidx_3[3] = { Reco_3mu_muW_idx[iIpt][BcNb], Reco_3mu_muW2_idx[iIpt][BcNb], muW3idx };
	  Short_t mumiidx_3[3] = { Reco_3mu_mumi_idx[iIpt][BcNb], Reco_QQ_mumi_idx[iIpt][QQidx_3[1]], Reco_QQ_mumi_idx[iIpt][QQidx_3[2]] };
	  Short_t muplidx_3[3] = { Reco_3mu_mupl_idx[iIpt][BcNb], Reco_QQ_mupl_idx[iIpt][QQidx_3[1]], Reco_QQ_mupl_idx[iIpt][QQidx_3[2]] };

	  int nCandidatePairs = (fabs(Reco_3mu_charge[iIpt][BcNb])==1)?2:3; //use candidates for all three Jpsi-dimuon choices if Bc_charge is wrong, and only the two OS pairs if Bc_charge is right
	  if(iIpt>=4) nCandidatePairs = 1;

	  //**************************************************************
	  //Loop on the 2 or 3 possible Jpsi dimuon choices 
	  int thisTrimuGaveABc = -1;
	  for(int k=0; k<nCandidatePairs; k++){
	    Short_t QQidx = QQidx_3[k];
	    Short_t QQ2idx = QQidx_3[(k==0 && iIpt!=4)?1:0];
	    Short_t muWidx = muWidx_3[k];
	    Short_t mumiidx = mumiidx_3[k];
	    Short_t muplidx = muplidx_3[k];

	    //cout<<"Dimuon choice #"<<k<<endl;
	    //can add pre-selection here
	    if(QQidx>-1){
	      //cout<<"Good dimuon"<<endl;

	      TLorentzVector *recBc = (TLorentzVector*) Reco_3mu_4mom[iIpt]->At(BcNb);
	      TLorentzVector *recQQ = (TLorentzVector*) Reco_QQ_4mom[iIpt]->At(QQidx);
	      if(muWidx==-1 || mumiidx==-1 || muplidx==-1) {cout<<"!!!!!!! one muon has a -1 index !! Skip this candidate"<<endl; continue;}
	      TLorentzVector *recBc_muW = (TLorentzVector*) ((iIpt!=4)?Reco_mu_4mom:Reco_trk_4mom)[iIpt]->At(muWidx); //muW is a track in case of dimuon+track
	      TLorentzVector *recBc_mumi = (TLorentzVector*) ((iIpt<5)?(Reco_mu_4mom[iIpt]->At(mumiidx)):(Reco_QQ_mumi_4mom[iIpt]->At(QQidx)));
	      TLorentzVector *recBc_mupl = (TLorentzVector*) ((iIpt<5)?(Reco_mu_4mom[iIpt]->At(muplidx)):(Reco_QQ_mupl_4mom[iIpt]->At(QQidx)));

	      float BcCandE = sqrt(pow(recQQ->P(),2) + pow(m_Jpsi,2) ) + recBc_muW->E() ;
	      float BcCandM = sqrt( pow( BcCandE ,2) - pow(recBc->P(),2) );
	      float QQCandM = recQQ->M();

	      muW_eta[i] = recBc_muW->Eta();
	      mumi_eta[i] = recBc_mumi->Eta();
	      mupl_eta[i] = recBc_mupl->Eta();
	      float maxEta = max(fabs(muW_eta[i]),max(fabs(mumi_eta[i]),fabs(mupl_eta[i])));

	      //**************************************************************
	      //What tree should we write the Bc candidate into?
	      goodTree = false;
	      switch(i) {
	      case 0: //WRONGSIGN
		goodTree = ( fabs(Reco_3mu_charge[iIpt][BcNb])!=1 || fabs(Reco_QQ_sign[iIpt][QQidx])!=0 ) && inJpsiMassRange(QQCandM, maxEta<1.5)
		  && (BcCandM < _mMax) && (BcCandM > 3.3);
		w_simple[i] = 3; //There are 2*6=12 combinations of 3 muons charges that give a total |charge|=1 (2 chosen OS pairs and 6 combinations), and 3*2=6 combinations giving |charge|=3. This means wrongsign deserves a weight=12/6=2 to be comparable to combinatorics background in data. //If we choose only one Bc candidate per trimuon, there are 6 combinations giving |charge|=1 and 2 combinations giving |charge|=3, so the ratio is 6/2=3
		break;

	      case 1: //BCMASS (Jpsi mass sidebands)
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inJpsiMassSB(QQCandM, maxEta<1.5) 
		  && (BcCandM < _mMax) && (BcCandM > 3.3);
		w_simple[i] = 1;
		break;

	      case 2: //TRUEJPSI
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM)//(inJpsiMassRange(QQCandM, maxEta<1.5) || inJpsiMassSB(QQCandM, maxEta<1.5))
		  && (BcCandM < _mMax) && (BcCandM > m_Bc);
		if(inJpsiMassRange(QQCandM, maxEta<1.5)){
		  w_simple[i] = 1;}
		else if(inJpsiMassSB(QQCandM, maxEta<1.5)){
		  w_simple[i] = -1;} // simple background subtraction (considering linear background) to keep only true J/psi's
		else{
		  w_simple[i] = 0;}
		break;

	      case 3: //SIGNAL region
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inJpsiMassRange(QQCandM, maxEta<1.5)
		  && (BcCandM < _mMax) && (BcCandM > 3.3);
		w_simple[i] = 1;
		break;

	      case 4: //MC SIGNAL
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mMax) && (BcCandM > 3.3) // in Bc mass region
		  && Reco_3mu_whichGen[iIpt][BcNb]>-1
		  && Reco_QQ_whichGen[iIpt][QQidx]>-1;
		w_simple[i] = scaleMCsig[ispp];
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) {w_simple[i] *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ w_simple[i] = 0;}
		break;

	      case 5: //MC B->J/psi X
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mMax) && (BcCandM > 3.3) // in Bc mass region
		  && Reco_3mu_whichGen[iIpt][BcNb]==-1
		  && Reco_QQ_whichGen[iIpt][QQidx]>-1;
		if(goodTree){
		  if(iIpt==1){ w_simple[i] = scaleMCsig[ispp];} //include Bc candidates from signal MC, unmatched to a gen Bc
		  else if(iIpt==2){
		    TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom[iIpt]->At(Reco_QQ_whichGen[iIpt][QQidx]);
		    w_simple[i] = ((ispp)?pthatweight[iIpt]:Gen_weight[iIpt]) * scaleMCb[ispp].second * (scaleMCb[ispp].first)->Eval(genQQ->Pt());
		  }
		}
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) { w_simple[i] *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ w_simple[i] = 0;}
		break;

	      case 6: //MC Prompt J/psi
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mMax) && (BcCandM > 3.3) // in Bc mass region
		  && Reco_3mu_whichGen[iIpt][BcNb]==-1
		  && Reco_QQ_whichGen[iIpt][QQidx]>-1;
		if(goodTree){
		  TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom[iIpt]->At(Reco_QQ_whichGen[iIpt][QQidx]);
		  w_simple[i] = ((ispp)?pthatweight[iIpt]:Gen_weight[iIpt]) * scaleMCprompt[ispp].second * (scaleMCprompt[ispp].first)->Eval(genQQ->Pt());
		}
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) { w_simple[i] *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ w_simple[i] = 0;}
		break;

	      case 7: //dimuon+track
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM)//(inJpsiMassRange(QQCandM, maxEta<1.5) || inJpsiMassSB(QQCandM, maxEta<1.5))
		  && (BcCandM < _mMax) && (BcCandM > 3.3);
		if(inJpsiMassRange(QQCandM, maxEta<1.5)){
		  w_simple[i] = 1;}
		else if(inJpsiMassSB(QQCandM, maxEta<1.5)){
		  w_simple[i] = -1;} // simple background subtraction (considering linear background) to keep only true J/psi's
		else{
		  w_simple[i] = 0;}
		break;

	      case 8: //Jpsi flipping
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM)//(inJpsiMassRange(QQCandM, maxEta<1.5) || inJpsiMassSB(QQCandM, maxEta<1.5))
		  && (BcCandM < _mMax) && (BcCandM > 3.3);
		if(inJpsiMassRange(QQCandM, maxEta<1.5)){
		  w_simple[i] = 1/7.;} //7 versions of Jpsi flipping are added
		else if(inJpsiMassSB(QQCandM, maxEta<1.5)){
		  w_simple[i] = -1/7.;} // simple background subtraction (considering linear background) to keep only true J/psi's
		else{
		  w_simple[i] = 0;}
		break;

	      case 9: //MC B->J/psi X , flipJpsi (only in pp)
		goodTree = fabs(Reco_3mu_charge[iIpt][BcNb])==1 && Reco_QQ_sign[iIpt][QQidx]==0 && inLooseMassRange(QQCandM) // in Jpsi mass region
		  && (BcCandM < _mMax) && (BcCandM > 3.3) // in Bc mass region
		  && Gen_QQ_size[iIpt]>0;
		if(goodTree){
		  TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom[iIpt]->At((Reco_QQ_whichGen[iIpt][QQidx]>-1)?Reco_QQ_whichGen[iIpt][QQidx]:0);
		  w_simple[i] = 1.17 * scaleMCb[ispp].second * (scaleMCb[ispp].first)->Eval(genQQ->Pt()); //1.17 is to compensate for the pthat>10 private sample missing in this flipJpsi b MC
		}
		if(inJpsiMassSB(QQCandM, maxEta<1.5)) { w_simple[i] *= -1;}
		else if(!(inJpsiMassRange(QQCandM, maxEta<1.5))){ w_simple[i] = 0;}
		break;

	      }

	      if(!goodTree || w_simple[i] == 0) continue;
	      if (thisTrimuGaveABc==i) continue; //two dimuon choices cannot give two candidates in the same output tree //Redundant check !

	      if(!ispp && (i>3 && i<7) ) {weightNcoll[i] = (float)findNcoll(INCentrality[iIpt]); //for PbPb MC
		w_simple[i] *= weightNcoll[i];}

	      //**************************************************************
	      //Write the wanted variables into the chosen (goodTree=true) output tree
	      Bc_ctauSignif[i] = Reco_3mu_ctau[iIpt][BcNb] / Reco_3mu_ctauErr[iIpt][BcNb] ;
	      Bc_ctauSignif3D[i] = Reco_3mu_ctau3D[iIpt][BcNb] / Reco_3mu_ctauErr3D[iIpt][BcNb] ;
	      Bc_alpha[i] = TMath::ACos(Reco_3mu_cosAlpha[iIpt][BcNb]);
	      Bc_alpha3D[i] = TMath::ACos(Reco_3mu_cosAlpha3D[iIpt][BcNb]);
	      Bc_VtxProb[i] = Reco_3mu_VtxProb[iIpt][BcNb];
	      QQ_VtxProb[i] = Reco_QQ_VtxProb[iIpt][QQidx];
	      QQ_dca[i] = Reco_QQ_dca[iIpt][QQidx];
	      QQ2_VtxProb[i] = (QQ2idx<Reco_QQ_size[iIpt] && QQ2idx>-1)?(Reco_QQ_VtxProb[iIpt][QQ2idx]):0;
	      QQ2_dca[i] = (QQ2idx<Reco_QQ_size[iIpt] && QQ2idx>-1)?(Reco_QQ_dca[iIpt][QQ2idx]):100;

	      muW_inLooseAcc[i] = ((iIpt!=4)?Reco_mu_InLooseAcc:Reco_trk_InLooseAcc)[iIpt][muWidx];
	      muW_inTightAcc[i] = ((iIpt!=4)?Reco_mu_InTightAcc:Reco_trk_InTightAcc)[iIpt][muWidx];
	      muW_isGlb[i] = (iIpt==4)?true:( (Reco_mu_SelType[iIpt][muWidx]&2)>0 );
	      muW_trig[i] = (iIpt==4)?false:( (Reco_mu_trig[iIpt][muWidx]&(ispp?8:4096))>0 ); //DoubleMu0 trigger = 2^3 for pp //HL_THIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 = 2^12 for PbPb
	      muW_isSoft[i] = (iIpt==4)?(isSoft(false, 1,1,1,1,Reco_trk_dxy[iIpt][muWidx], Reco_trk_dz[iIpt][muWidx], 20,20)):(
				         isSoft(false, muW_isGlb[i] , (Reco_mu_SelType[iIpt][muWidx]&8)>0, ((Reco_mu_SelType[iIpt][muWidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[iIpt][muWidx] , Reco_mu_dxy[iIpt][muWidx] , Reco_mu_dz[iIpt][muWidx] , Reco_mu_nPixWMea[iIpt][muWidx] , Reco_mu_nTrkWMea[iIpt][muWidx] ));

	      mumi_inLooseAcc[i] = Reco_mu_InLooseAcc[iIpt][mumiidx];
	      mumi_inTightAcc[i] = Reco_mu_InTightAcc[iIpt][mumiidx];
	      mumi_isGlb[i] = (Reco_mu_SelType[iIpt][mumiidx]&2)>0;
	      mumi_trig[i] = (Reco_mu_trig[iIpt][mumiidx]&(ispp?8:4096))>0;
	      mumi_isSoft[i] = isSoft( false, mumi_isGlb[i] , (Reco_mu_SelType[iIpt][mumiidx]&8)>0 , ((Reco_mu_SelType[iIpt][mumiidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[iIpt][mumiidx] , Reco_mu_dxy[iIpt][mumiidx] , Reco_mu_dz[iIpt][mumiidx] , Reco_mu_nPixWMea[iIpt][mumiidx] , Reco_mu_nTrkWMea[iIpt][mumiidx] );

	      mupl_inLooseAcc[i] = Reco_mu_InLooseAcc[iIpt][muplidx];
	      mupl_inTightAcc[i] = Reco_mu_InTightAcc[iIpt][muplidx];
	      mupl_isGlb[i] = (Reco_mu_SelType[iIpt][muplidx]&2)>0;
	      mupl_trig[i] = (Reco_mu_trig[iIpt][muplidx]&(ispp?8:4096))>0;
	      mupl_isSoft[i] = isSoft( false, mupl_isGlb[i] , (Reco_mu_SelType[iIpt][muWidx]&8)>0 , ((Reco_mu_SelType[iIpt][muplidx]&((int)pow(2,12)))>0) , Reco_mu_highPurity[iIpt][muplidx] , Reco_mu_dxy[iIpt][muplidx] , Reco_mu_dz[iIpt][muplidx] , Reco_mu_nPixWMea[iIpt][muplidx] , Reco_mu_nTrkWMea[iIpt][muplidx] );

	      float muW_dxy = (iIpt==4)?(Reco_trk_dxy[iIpt][muWidx]):( (ispp || iIpt>=5)?(Reco_3mu_muW_dxy_muonlessVtx[iIpt][BcNb]):(Reco_mu_dxy[iIpt][muWidx]) );
	      float mumi_dxy = (ispp || iIpt>=5)?(Reco_3mu_mumi_dxy_muonlessVtx[iIpt][BcNb]):(Reco_mu_dxy[iIpt][mumiidx]);
	      float mupl_dxy = (ispp || iIpt>=5)?(Reco_3mu_mupl_dxy_muonlessVtx[iIpt][BcNb]):(Reco_mu_dxy[iIpt][muplidx]);
	      float muW_dz = (iIpt==4)?(Reco_trk_dz[iIpt][muWidx]):( (ispp || iIpt>=5)?(Reco_3mu_muW_dz_muonlessVtx[iIpt][BcNb]):(Reco_mu_dz[iIpt][muWidx]) );
	      float mumi_dz = (ispp || iIpt>=5)?(Reco_3mu_mumi_dz_muonlessVtx[iIpt][BcNb]):(Reco_mu_dz[iIpt][mumiidx]);
	      float mupl_dz = (ispp || iIpt>=5)?(Reco_3mu_mupl_dz_muonlessVtx[iIpt][BcNb]):(Reco_mu_dz[iIpt][muplidx]);
	      muW_Pt[i] = recBc_muW->Pt();
	      mumi_Pt[i] = recBc_mumi->Pt();
	      mupl_Pt[i] = recBc_mupl->Pt();
	      
	      //cout<<(Bc_ctauSignif[i]>1.5)<<" "<<(Bc_alpha[i]<0.8)<<" "<<(Bc_alpha3D[i]<0.8)<<" "<<(Bc_VtxProb[i]>0.01)<<" "<<(QQ_dca[i]<0.3)<<" "<<(muW_isSoft[i])<<" "<<( mumi_isSoft[i])<<" "<<(mupl_isSoft[i])<<" "<<( (muW_isGlb[i] && muW_inLooseAcc[i] && mupl_isGlb[i] && mupl_inLooseAcc[i]) || (muW_isGlb[i] && muW_inLooseAcc[i] && mumi_isGlb[i] && mumi_inLooseAcc[i]) || (mumi_isGlb[i] && mumi_inLooseAcc[i] && mupl_isGlb[i] && mupl_inLooseAcc[i]) )<<" "<<(mumi_isGlb[i] && mupl_isGlb[i] && muW_isGlb[i] && mumi_inLooseAcc[i] && mupl_inLooseAcc[i] && muW_inLooseAcc[i])<<" "<<muW_isGlb[i]<<" "<<mupl_isGlb[i]<<" "<<mumi_isGlb[i]<<" "<<((muW_trig[i] && mupl_trig[i] && muW_inTightAcc[i] && mupl_inTightAcc[i] ) ||(muW_trig[i] && mumi_trig[i] && muW_inTightAcc[i] && mumi_inTightAcc[i] ) || (mumi_trig[i] && mupl_trig[i] && mumi_inTightAcc[i] && mupl_inTightAcc[i] ))<<" "<<(fabs(muW_dz)<(ispp?0.6:0.8) && fabs(mumi_dz)<(ispp?0.6:0.8) && fabs(mupl_dz)<(ispp?0.6:0.8))<<" "<<((HLTriggers[iIpt]&((ispp || i>=4)?8:4096))>0)<<endl;
	      if(
		 //**************************************************************
		 // Pre-selections
		 Bc_ctauSignif[i]>_ctauSignif_cut
		 && Bc_ctauSignif3D[i]>_ctauSignif3D_cut //maybe not a good idea in pp due to pile-up along z. But actually not a big deal if we select nonprompt objects
		 && Bc_alpha[i]<_alpha_cut(ispp)
		 && Bc_alpha3D[i]<_alpha3D_cut(ispp) //combined with alpha cut, kills 1.4% of signal in pp (1.1% in PbPb) 
		 && Bc_VtxProb[i]>((iIpt==4)?_vtxProb_cutLoose:_vtxProb_cut) //kills 1.2% more signal compared to Bc_VtxProb>0.005
		 && QQ_VtxProb[i]>((iIpt==4)?_QQvtxProb_cutTight:_QQvtxProb_cut) //drop the QQ_VtxProb, too correlated with Bc_VtxProb?
		 && QQ_dca[i]<_QQdca_cut //keep this that kills 1.8% of signal and 10% of WRONSIGN/BCMASS
		 && muW_isSoft[i]
		 && mumi_isSoft[i]
		 && mupl_isSoft[i]
		 && (iIpt!=4 || muW_inLooseAcc[i]) //force track inLooseAcceptance in case of dimuon+trk
		 && (_withTM?( (muW_isGlb[i] && muW_inLooseAcc[i] && mupl_isGlb[i] && mupl_inLooseAcc[i] && looseAcc(mumi_Pt[i],mumi_eta[i],true) ) || //only one muon can be tracker and out of LooseAcceptance     
			      (muW_isGlb[i] && muW_inLooseAcc[i] && mumi_isGlb[i] && mumi_inLooseAcc[i] && looseAcc(mupl_Pt[i],mupl_eta[i],true) ) ||
			      (mumi_isGlb[i] && mumi_inLooseAcc[i] && mupl_isGlb[i] && mupl_inLooseAcc[i] && looseAcc(muW_Pt[i],muW_eta[i],true) )
			      ):(
				 mumi_isGlb[i] && mupl_isGlb[i] && muW_isGlb[i]
				 && mumi_inLooseAcc[i] && mupl_inLooseAcc[i] && muW_inLooseAcc[i]
				 ))
		 && ( ( muW_trig[i] && mupl_trig[i] && muW_inTightAcc[i] && mupl_inTightAcc[i] ) || //two muons among three must trigger //BEWARE ! Not sure if TightAcceptance should be put there
		      ( muW_trig[i] && mumi_trig[i] && muW_inTightAcc[i] && mumi_inTightAcc[i] ) ||
		      ( mumi_trig[i] && mupl_trig[i] && mumi_inTightAcc[i] && mupl_inTightAcc[i] ) //only this last option can be true for dimuon+trk
		      )
		 && fabs(muW_dz)<0.6 && fabs(mumi_dz)<0.6 && fabs(mupl_dz)<0.6
		 && (!ispp || (HLTriggers[iIpt]&((ispp || i>=4)?8:4096))>0) //the event must fire the trigger as well //BEWARE ! this should be re-established for PbPb when samples are re-run
		 //!!!!!!!!!!!!!!!!!!!!!! Here, should reject the tracks that are global muons, for dimuon+track
		 ){
		
		QQ_M[i] = QQCandM;
		if(iIpt!=4){
		  QQ2_M[i] = (Reco_mu_charge[iIpt][muWidx]>0)?((*recBc_mumi+*recBc_muW).M()):((*recBc_mupl+*recBc_muW).M()); //QQ2 is the second OS pair
		  QQ3_M[i] = (Reco_mu_charge[iIpt][muWidx]>0)?((*recBc_mupl+*recBc_muW).M()):((*recBc_mumi+*recBc_muW).M()); //QQ3 is the SS pair
		}
		
		// //**** Measure efficiency of dimuon mass signal region //needs to remove previous Jpsi mass cuts 
		// if(k==0){
		//   nall += w_simple[i];
		//   if(inJpsiMassRange(QQ_M[i], maxEta<1.5) || inJpsiMassRange(QQ2_M[i], maxEta<1.5)) npass +=w_simple[i];
		//   if(!inJpsiMassRange(QQ_M[i], maxEta<1.5) && !inJpsiMassRange(QQ2_M[i], maxEta<1.5) && (inJpsiMassSB(QQ_M[i], maxEta<1.5) || inJpsiMassSB(QQ2_M[i], maxEta<1.5)) ) nside +=w_simple[i];
		// }
		//**** Deal with the Jpsi dimuon choice
		//		weightJpsiChoice[i] = 1;
		if((i==0 || i==2 || i==3) && inJpsiMassRange(QQ_M[i], maxEta<1.5) && inJpsiMassRange(QQ2_M[i], maxEta<1.5) 
		   && QQ2_dca[i]<QQ_dca[i] && QQ2_VtxProb[i]>_QQvtxProb_cut)
		  continue; //In data, if both OS dimuons are in Jpsi peak region, choose the one with best dca (so no bias on the Jpsi mass)
		if(i==1 && inJpsiMassSB(QQ_M[i], maxEta<1.5) && inJpsiMassSB(QQ2_M[i], maxEta<1.5)
		   && QQ2_dca[i]<QQ_dca[i] && QQ2_VtxProb[i]>_QQvtxProb_cut)
		   continue; //In data sidebands, if both OS dimuons are in Jpsi sidebands, choose the one with best dca (so no bias on the Jpsi mass)

		// if(//(i==1 || i==2) && 
		//    inJpsiMassSB(QQ_M[i], maxEta<1.5) && inJpsiMassRange(QQ2_M[i], maxEta<1.5) && Reco_QQ_dca[iIpt][QQ2idx]<0.3) weightJpsiChoice[i] = (ispp?0.213:0.662); //if QQ1 passes all BCMASS sample cuts, check if the other dimuon QQ2 passes the SIGNAL REGION cuts. If yes, apply a weight corresponding to the integrated proba of a data event not to contain a true Jpsi
		// if(i!=0 && //(i==0 || i==2 || i==3) && 
		//    inJpsiMassRange(QQ_M[i], maxEta<1.5) && inJpsiMassSB(QQ2_M[i], maxEta<1.5) && Reco_QQ_dca[iIpt][QQ2idx]<0.3) weightJpsiChoice[i] = (1-(ispp?0.213:0.662)); //if QQ1 passes all SIGNAL sample cuts, check if the other dimuon QQ2 passes the jpsi mass sidebands cuts. If yes, apply a weight corresponding to the integrated proba of a data event to contain a true Jpsi
		// w_simple[i] *= weightJpsiChoice[i];

		// if(!ispp && i==5 && iIpt==2 && recQQ->Pt()>6.5 && recQQ->Pt()<7.5)
		//   cout<< "weight, NcollWeight, Gen_weight, A^2*L_PbPb*BF, dataOverMC(QQ_Pt), weight/(NcollWeight*Gen_weight*lumi*dataOverMC) = "<<w_simple[i]<<" "<<weightNcoll[i]<<" "<<Gen_weight[iIpt]<<" "<< scaleMCb[ispp].second<<" "<<(scaleMCb[ispp].first)->Eval(recQQ->Pt())<<" "<<w_simple[i]/(weightNcoll[i]*Gen_weight[iIpt]*scaleMCb[ispp].second*(scaleMCb[ispp].first)->Eval(recQQ->Pt()))<<endl;
		// if(ispp && i==5 && iIpt==2 && recQQ->Pt()>6.5 && recQQ->Pt()<7.5)
		//   cout<< "weight, pthat_weight, A^2*L_PbPb*BF, dataOverMC(QQ_Pt), weight/(NcollWeight*Gen_weight*lumi*dataOverMC) = "<<w_simple[i]<<" "<<pthatweight[iIpt]<<" "<< scaleMCb[ispp].second<<" "<<(scaleMCb[ispp].first)->Eval(recQQ->Pt())<<" "<<w_simple[i]/(pthatweight[iIpt]*scaleMCb[ispp].second*(scaleMCb[ispp].first)->Eval(recQQ->Pt()))<<endl;
		
		//cout<<"Passes all selections!!!!!!!!!"<<endl;
		//**** Finish filling the output variables
		eventNb[i] = INeventNb[iIpt];
		if(iIpt==0 || iIpt==4 || iIpt==5){
		  runNb[i] = INrunNb[iIpt];
		  LS[i] = INLS[iIpt];}
		else{
		  runNb[i] = 0;
		  LS[i] = 0;
		}
		if(ispp){ nPV[i] = INnPV[iIpt];//(int)round(INnPV[iIpt]);
		}else{ Centrality[i] = INCentrality[iIpt];}

		bkgType[i] = i;
		if(iIpt==1 && i==5) bkgType[i] = 9; //new bkgType for unmatched trimuons from signal MC, merged in the b->Jpsi background //BEWARE, this was =7 before introduction of dimuon+track!
		if(iIpt==6) bkgType[i] = 10; //flipJpsi b MC

		float QQM_correction = m_Jpsi / recQQ->M();

		Bc_M[i] = BcCandM;
		//Bc_M[i] = sqrt( pow(  QQM_correction*recQQ->E() + recBc_muW->E()  ,2) - pow( Bc_P[i] ,2) ); //sqrt( (QQ_E + muW_E)^2 - Bc_P^2 )
		Bc_charge[i] = Reco_3mu_charge[iIpt][BcNb];
		//Bc_P[i] = recBc->P();
		//Bc_Pt[i] = recBc->Pt();
		Bc_Pt[i] = sqrt( pow(QQM_correction * recQQ->Px()+recBc_muW->Px() ,2) + pow(QQM_correction * recQQ->Py()+recBc_muW->Py() ,2) );
		Bc_P[i] = sqrt( pow(Bc_Pt[i] ,2) + pow(QQM_correction * recQQ->Pz()+recBc_muW->Pz() ,2) );
		Bc_Y[i] = recBc->Rapidity();//0.5* TMath::Log((BcCandE + recBc->Pz())/(BcCandE - recBc->Pz()));
		Bc_phi[i] = recBc->Phi();
		Bc_ctau[i] = Reco_3mu_ctau[iIpt][BcNb] * recBc->Pt()/Bc_Pt[i];
		Bc_ctau3D[i] = Reco_3mu_ctau3D[iIpt][BcNb] * recBc->P()/Bc_P[i];
		float PperpTrimu = TMath::Sin(Bc_alpha3D[i]) * Bc_P[i];
		Bc_CorrM[i] = sqrt(Bc_M[i]*Bc_M[i] + PperpTrimu*PperpTrimu) + PperpTrimu;

		Bc_M_withQQM[i] = recBc->M();
		Bc_M_muWisK[i] = sqrt( pow(   sqrt(pow(recBc_muW->P(),2) + pow(m_K,2) ) + sqrt(pow(recQQ->P(),2) + pow(m_Jpsi,2) )  ,2) - pow(Bc_P[i],2) );
		//Bc_M_muWismu[i] = sqrt( pow(   sqrt(pow(recBc_muW->P(),2) + pow(m_mu,2) ) + sqrt(pow(recQQ->P(),2) + pow(m_Jpsi,2) )  ,2) - pow(Bc_P[i],2) );
		//Bc_M_muWisPi[i] = sqrt( pow(   sqrt(pow(recBc_muW->P(),2) + pow(m_Pi,2) ) + sqrt(pow(recQQ->P(),2) + pow(m_Jpsi,2) )  ,2) - pow(Bc_P[i],2) );

		QQ_P[i] = QQM_correction * recQQ->P();
		QQ_Pt[i] = QQM_correction * recQQ->Pt();
		QQ_Y[i] = recQQ->Rapidity();//0.5* TMath::Log((sqrt(m_Jpsi*m_Jpsi + QQ_P[i]*QQ_P[i]) + QQM_correction * recQQ->Pz())/(sqrt(m_Jpsi*m_Jpsi + QQ_P[i]*QQ_P[i]) - QQM_correction * recQQ->Pz()));
		QQ_phi[i] = recQQ->Phi();
		QQ_charge[i] = Reco_QQ_sign[iIpt][QQidx];
		QQ_momId[i] = (iIpt==0)?0:Gen_QQ_momId[iIpt][Reco_QQ_whichGen[iIpt][QQidx]];

		dR_jpsi[i] = (float) recBc_mumi->DeltaR(*recBc_mupl);
		float dRjpsi_corr = (dR_jpsi[i]<=0)?1:( TMath::ACos(1- QQM_correction*QQM_correction* (1-TMath::Cos(dR_jpsi[i])) ) / dR_jpsi[i] );
		dR_jpsi[i] = dRjpsi_corr * dR_jpsi[i];
		if(dR_jpsi[i]<=0) dR_jpsi[i] = 1e-8; if(dR_jpsi[i]>4) dR_jpsi[i] = 4; //check that deltaR is in [0,4]

		float dRmuWmi = (float) recBc_muW->DeltaR(*recBc_mumi);
		float dRmuWpl = (float) recBc_muW->DeltaR(*recBc_mupl);
		float mumiCorr = 0.5*(dRjpsi_corr-1) * mupl_Pt[i]/(mupl_Pt[i]+mumi_Pt[i]);
		float muplCorr = mumiCorr*mumi_Pt[i]/mupl_Pt[i];
		dR_muWmi[i] = dRmuWmi*( 1+ mumiCorr*(1+ (pow(dR_jpsi[i]/dRjpsi_corr ,2) - pow( dRmuWpl,2))/pow(dRmuWmi ,2) ) );
		dR_muWpl[i] = dRmuWpl*( 1+ muplCorr*(1+ (pow(dR_jpsi[i]/dRjpsi_corr ,2) - pow( dRmuWmi,2))/pow(dRmuWpl ,2) ) );
		if(dR_muWmi[i]<=0) dR_muWmi[i] = 1e-8; if(dR_muWmi[i]>4) dR_muWmi[i] = 4; if(std::isnan(dR_muWmi[i])) dR_muWmi[i] = dRmuWmi; //check that deltaR is in [0,4]
		if(dR_muWpl[i]<=0) dR_muWpl[i] = 1e-8; if(dR_muWpl[i]>4) dR_muWpl[i] = 4; if(std::isnan(dR_muWpl[i])) dR_muWpl[i] = dRmuWpl; //check that deltaR is in [0,4]

		dR_jpsiMuW[i] = (float) recQQ->DeltaR(*recBc_muW);
		if(dR_jpsiMuW[i]<=0) dR_jpsiMuW[i] = 1e-8; if(dR_jpsiMuW[i]>4) dR_jpsiMuW[i] = 4; //check that deltaR is in [0,4]
		dR_sum[i] = dR_muWmi[i] + dR_muWpl[i] + dR_jpsi[i];
		dR_jpsiOverMuW[i] = dR_jpsi[i] / (dR_muWmi[i] + dR_muWpl[i]) ;
		if (dR_jpsiOverMuW[i]>1) dR_jpsiOverMuW[i] = 1.;

		//mass shift for the TRUEJPSI (high mass control region) bkg
		Bc_M_shiftedM[i] = Bc_M[i];
		Bc_CorrM_shiftedM[i] = Bc_CorrM[i];
		dR_jpsiOverMuW_shiftedM[i] = dR_jpsiOverMuW[i];
		dR_jpsiMuW_shiftedM[i] = dR_jpsiMuW[i];
		dR_sum_shiftedM[i] = dR_sum[i];
		if(i==2){
		  float BcM_correction = 1-1.2/Bc_M[i];
		  Bc_M_shiftedM[i] = Bc_M[i] * BcM_correction;
		  Bc_CorrM_shiftedM[i] = sqrt(Bc_M_shiftedM[i]*Bc_M_shiftedM[i] + PperpTrimu*PperpTrimu) + PperpTrimu;
		  float dRmuWplmi = TMath::ACos(1- BcM_correction*BcM_correction* (1-TMath::Cos(dR_muWpl[i])) ) + 
		                    TMath::ACos(1- BcM_correction*BcM_correction* (1-TMath::Cos(dR_muWmi[i])) );
		  dR_jpsiOverMuW_shiftedM[i] = dR_jpsi[i] / dRmuWplmi;
		  if (dR_jpsiOverMuW_shiftedM[i]>1) dR_jpsiOverMuW_shiftedM[i] = 1.; if(std::isnan(dR_jpsiOverMuW_shiftedM[i])) dR_jpsiOverMuW_shiftedM[i] = dR_jpsiOverMuW[i];
		  dR_jpsiMuW_shiftedM[i] = TMath::ACos(1- BcM_correction*BcM_correction* (1-TMath::Cos(dR_jpsiMuW[i])) );
		  if(dR_jpsiMuW_shiftedM[i]<0) dR_jpsiMuW_shiftedM[i] = 0; if(dR_jpsiMuW_shiftedM[i]>4) dR_jpsiMuW_shiftedM[i] = 4; if(std::isnan(dR_jpsiMuW_shiftedM[i])) dR_jpsiMuW_shiftedM[i] = dR_jpsiMuW[i]; //check that deltaR is in [0,4]
		  dR_sum_shiftedM[i] = dRmuWplmi + dR_jpsi[i];
		}

		muW_isJpsiBro[i] = (iIpt!=2 && iIpt!=3 && iIpt!=6)?false:Reco_3mu_muW_isGenJpsiBro[iIpt][BcNb];
		muW_trueId[i] = (iIpt!=2 && iIpt!=3 && iIpt!=6)?0:Reco_3mu_muW_trueId[iIpt][BcNb];
		muW_normChi2_inner[i] = (iIpt==4)?0:Reco_mu_normChi2_inner[iIpt][muWidx];
		muW_normChi2_glb[i] = (iIpt==4)?0:Reco_mu_normChi2_global[iIpt][muWidx];

		muW_P[i] = recBc_muW->P();
		muW_phi[i] = recBc_muW->Phi();
		mumi_P[i] = recBc_mumi->P();
		mumi_phi[i] = recBc_mumi->Phi();
		mupl_P[i] = recBc_mupl->P();
		mupl_phi[i] = recBc_mupl->Phi();

		muW_dxySignif[i] = fabs( muW_dxy / ((iIpt==4)?(Reco_trk_dxyError[iIpt][muWidx]):(Reco_mu_dxyErr[iIpt][muWidx]) ) );
		muW_dzSignif[i] = fabs( muW_dz / ((iIpt==4)?(Reco_trk_dzError[iIpt][muWidx]):(Reco_mu_dzErr[iIpt][muWidx]) ) );
		mumi_dxySignif[i] = fabs( mumi_dxy / Reco_mu_dxyErr[iIpt][mumiidx] );
		mumi_dzSignif[i] = fabs( mumi_dz / Reco_mu_dzErr[iIpt][mumiidx] );
		mupl_dxySignif[i] = fabs( mupl_dxy / Reco_mu_dxyErr[iIpt][muplidx] );
		mupl_dzSignif[i] = fabs( mupl_dz / Reco_mu_dzErr[iIpt][muplidx] );	      
		MuonDxySignif_sum[i] = muW_dxySignif[i] + mumi_dxySignif[i] + mupl_dxySignif[i];
		MuonDzSignif_sum[i] = muW_dzSignif[i] + mumi_dzSignif[i] + mupl_dzSignif[i];
	      
		QQmuW_ptImbal[i] = fabs((QQ_Pt[i]-muW_Pt[i])/(QQ_Pt[i]+muW_Pt[i]));

		if(iIpt>=5) flipJpsi[i] = Reco_QQ_flipJpsi[iIpt][QQidx];

		if(iIpt==1){
		  if(genBcIdx>-1){
		  TLorentzVector *genBc = (TLorentzVector*) Gen_Bc_4mom[iIpt]->At(genBcIdx);
		  TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom[iIpt]->At(genBcIdx);
		  genBc_Pt[i] = genBc->Pt();
		  gen3mu_Pt[i] = gen3mu->Pt();}
		  else{ 
		    genBc_Pt[i] = 0;
		    gen3mu_Pt[i] = 0;}
		}

		if(Bc_CorrM_shiftedM[i] > _BcCorrM_cut(ispp) || dR_sum_shiftedM[i]>7.5) continue; //cuts no signal in PbPb, and 0.1% in pp
		if(!ispp && Centrality[i]>180) continue; //keep 0-90% centrality
	    	     
		if(!(std::isnan(w_simple[i]))){
		  if(!ispp && i==3 && Bc_M[i]<_mBcMax){ //whether to blind 3/4 of the events of signal region
		    w_unblind[i] = w_simple[i];
		    if((j%4)==0) w_simple[i] *= 4;
		    else w_simple[i]=0;
		  }

		  thisTrimuGaveABc = i;
		  out_trees[i]->Fill();
		}
	      } //end if passes full selection

	    } //end if(QQ_isValid)
	  } //end loop on 2/3 possible Jpsi dimuon choice
	} //end loop on Bc candidates
      } //end loop on selection of output trees 
    } //end loop on entries
    // cout<<"Jpsi mass signal region nall,npass,efficiency = "<<nall<<" "<<npass<<" "<<npass/nall<<endl;
    // cout<<"Jpsi mass sidebands nall,npass,efficiency = "<<nall<<" "<<nside<<" "<<nside/nall<<endl;
  } //end loop on input trees

  //**************************************************************
  //Save trees and close file 
  for(int i=0;i<ntrees;i++){
    out_trees[i]->Print(); cout<<endl<<endl<<endl;
    out_trees[i]->AutoSave();
  }

  // out_file.Close();

}
