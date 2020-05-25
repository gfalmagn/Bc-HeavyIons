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

void addBDTvariable(bool ispp=true){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  TMVA::Tools::Instance();

  bool useVarCorrWMass = true;

  //need to start from untouched saved file
  //gSystem->Exec("cp BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+"_copystep1.root");  
  auto fullFile = TFile::Open("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","UPDATE");
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
  float QQ_M[ntrees];
  float QQ2_M[ntrees];
  float QQ3_M[ntrees];
  int bkgType[ntrees];
  int muW_isJpsiBro[ntrees];
  int muW_trueId[ntrees];
  float w_simple2[ntrees];
  UInt_t eventNb[ntrees];
  float BDT[ntrees];
  float BDTprob[ntrees];
  float BDTrarity[ntrees];

  //Which BDT is used, and what fitting strategy
  TString weightFile = "BDTfiner"+(TString)(ispp?"":"_PbPb")+"_withJpsiMC"+(TString)(useVarCorrWMass?"":"_dropVarCorrWMass"); //BDTfinerShallow
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi"};
  TString prettyName[] = {"WRONGSIGN","J/Psi sidebands","High mass control","signal region","MC signal expectation",
			  "MC NonPromptJpsi","MC PromptJpsi","dimuon+track (misID)","flipped J/Psi"};
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }

  //*******************************************
  //Get BDT variable for signal region
  //Declare Reader
  vector<TMVA::Reader*> reader1; //2 readers for the 2 half samples (train/test and test/train)
  vector<TMVA::Reader*> reader2;
  vector<TBranch*> b_BDT, b_BDTprob, b_BDTrarity;

  //*******************************************
  //For the signal region and the tmva output, fill the branches for BDT
  for(int iT=0; iT<(int)T.size(); iT++){
    std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    b_BDT.push_back( T[iT]->Branch("BDT",&BDT[iT],"BDT/F") );
    b_BDTprob.push_back( T[iT]->Branch("BDTprob",&BDTprob[iT],"BDTprob/F") );
    b_BDTrarity.push_back( T[iT]->Branch("BDTrarity",&BDTrarity[iT],"BDTrarity/F") );
		     
    //Prepare the analyzed tree for BDT reading 
    //first reader (train:eventNb even, test:eventNb odd)
    reader1.push_back(new TMVA::Reader( "Color:!Silent" ));

    //Declare the used variables -- must be the same as in training
    if(useVarCorrWMass) reader1[iT]->AddVariable("Bc_CorrM_shiftedM", &Bc_CorrM[iT]);
    reader1[iT]->AddSpectator("Bc_ctauSignif", &Bc_ctauSignif[iT]);
    reader1[iT]->AddVariable("Bc_log_ctauSignif3D := TMath::Log(Bc_ctauSignif3D+1e-8)", &Bc_log_ctauSignif3D[iT]);
    reader1[iT]->AddSpectator("Bc_alpha", &Bc_alpha[iT]);
    reader1[iT]->AddVariable("Bc_alpha3D", &Bc_alpha3D[iT]);
    reader1[iT]->AddVariable("Bc_log1min_VtxProb := TMath::Log(1-Bc_VtxProb+1e-6)", &Bc_log1min_VtxProb[iT]);
    reader1[iT]->AddSpectator("QQ_log1min_VtxProb := TMath::Log(1-QQ_VtxProb+1e-6)", &QQ_log1min_VtxProb[iT]);
    reader1[iT]->AddVariable("QQ_log_dca := TMath::Log(QQ_dca+1e-8)", &QQ_log_dca[iT]);
    if(useVarCorrWMass) {
      reader1[iT]->AddVariable("dR_sum_shiftedM", &dR_sum[iT]);
      reader1[iT]->AddVariable("dR_jpsiOverMuW_shiftedM", &dR_jpsiOverMuW[iT]);
      reader1[iT]->AddSpectator("dR_jpsiMuW_shiftedM", &dR_jpsiMuW[iT]);
    }
    reader1[iT]->AddVariable("log_MuonDxySignif_sum := TMath::Log(MuonDxySignif_sum+1e-8)", &log_MuonDxySignif_sum[iT]);
    reader1[iT]->AddSpectator("MuonDzSignif_sum", &MuonDzSignif_sum[iT]);
    reader1[iT]->AddVariable("muonsGlbInLooseAcc := (muW_isGlb && muW_inLooseAcc) + (mumi_isGlb && mumi_inLooseAcc) + (mupl_isGlb && mupl_inLooseAcc)",&muonsGlbInLooseAcc[iT]);
    reader1[iT]->AddVariable("muonsInTightAcc := muW_inTightAcc + mumi_inTightAcc + mupl_inTightAcc",&muonsInTightAcc[iT]);
    reader1[iT]->AddVariable("QQmuW_ptImbal", &QQmuW_ptImbal[iT]);
    reader1[iT]->AddSpectator("Bc_M", &Bc_M[iT]);
    reader1[iT]->AddSpectator("Bc_Pt", &Bc_Pt[iT]);
    reader1[iT]->AddSpectator("QQ_M", &QQ_M[iT]);
    // reader1[iT]->AddSpectator("QQ2_M", &QQ2_M[iT]);
    // reader1[iT]->AddSpectator("QQ3_M", &QQ3_M[iT]);
    reader1[iT]->AddSpectator("bkgType", &bkgType[iT]);
    reader1[iT]->AddSpectator("muW_isJpsiBro", &muW_isJpsiBro[iT]);
    reader1[iT]->AddSpectator("muW_trueId", &muW_trueId[iT]);
    reader1[iT]->AddSpectator("w_simple2", &w_simple2[iT]);

    //Book MVA methods
    reader1[iT]->BookMVA("BDT","dataset/weights/BcSigBkgClassification_"+weightFile+"_1stHalf.weights.xml");

    //second reader (train:eventNb odd, test:eventNb even)
    reader2.push_back(new TMVA::Reader( "Color:!Silent" ));

    //Declare the used variables -- must be the same as in training
    if(useVarCorrWMass) reader2[iT]->AddVariable("Bc_CorrM_shiftedM", &Bc_CorrM[iT]);
    reader2[iT]->AddSpectator("Bc_ctauSignif", &Bc_ctauSignif[iT]);
    reader2[iT]->AddVariable("Bc_log_ctauSignif3D := TMath::Log(Bc_ctauSignif3D+1e-8)", &Bc_log_ctauSignif3D[iT]);
    reader2[iT]->AddSpectator("Bc_alpha", &Bc_alpha[iT]);
    reader2[iT]->AddVariable("Bc_alpha3D", &Bc_alpha3D[iT]);
    reader2[iT]->AddVariable("Bc_log1min_VtxProb := TMath::Log(1-Bc_VtxProb+1e-6)", &Bc_log1min_VtxProb[iT]);
    reader2[iT]->AddSpectator("QQ_log1min_VtxProb := TMath::Log(1-QQ_VtxProb+1e-6)", &QQ_log1min_VtxProb[iT]);
    reader2[iT]->AddVariable("QQ_log_dca := TMath::Log(QQ_dca+1e-8)", &QQ_log_dca[iT]);
    if(useVarCorrWMass) {
      reader2[iT]->AddVariable("dR_sum_shiftedM", &dR_sum[iT]);
      reader2[iT]->AddVariable("dR_jpsiOverMuW_shiftedM", &dR_jpsiOverMuW[iT]);
      reader2[iT]->AddSpectator("dR_jpsiMuW_shiftedM", &dR_jpsiMuW[iT]);
    }
    reader2[iT]->AddVariable("log_MuonDxySignif_sum := TMath::Log(MuonDxySignif_sum+1e-8)", &log_MuonDxySignif_sum[iT]);
    reader2[iT]->AddSpectator("MuonDzSignif_sum", &MuonDzSignif_sum[iT]);
    reader2[iT]->AddVariable("muonsGlbInLooseAcc := (muW_isGlb && muW_inLooseAcc) + (mumi_isGlb && mumi_inLooseAcc) + (mupl_isGlb && mupl_inLooseAcc)",&muonsGlbInLooseAcc[iT]);
    reader2[iT]->AddVariable("muonsInTightAcc := muW_inTightAcc + mumi_inTightAcc + mupl_inTightAcc",&muonsInTightAcc[iT]);
    reader2[iT]->AddVariable("QQmuW_ptImbal", &QQmuW_ptImbal[iT]);
    reader2[iT]->AddSpectator("Bc_M", &Bc_M[iT]);
    reader2[iT]->AddSpectator("Bc_Pt", &Bc_Pt[iT]);
    reader2[iT]->AddSpectator("QQ_M", &QQ_M[iT]);
    // reader2[iT]->AddSpectator("QQ2_M", &QQ2_M[iT]);
    // reader2[iT]->AddSpectator("QQ3_M", &QQ3_M[iT]);
    reader2[iT]->AddSpectator("bkgType", &bkgType[iT]);
    reader2[iT]->AddSpectator("muW_isJpsiBro", &muW_isJpsiBro[iT]);
    reader2[iT]->AddSpectator("muW_trueId", &muW_trueId[iT]);
    reader2[iT]->AddSpectator("w_simple2", &w_simple2[iT]);

    //Book MVA methods //only difference between the two reader is here
    reader2[iT]->BookMVA("BDT","dataset/weights/BcSigBkgClassification_"+weightFile+"_2ndHalf.weights.xml");

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
    T[iT]->SetBranchAddress("dR_sum_shiftedM", &dR_sum[iT]);
    T[iT]->SetBranchAddress("dR_jpsiOverMuW_shiftedM", &dR_jpsiOverMuW[iT]);
    T[iT]->SetBranchAddress("dR_jpsiMuW_shiftedM", &dR_jpsiMuW[iT]);
    T[iT]->SetBranchAddress("MuonDxySignif_sum", &MuonDxySignif_sum[iT]);
    T[iT]->SetBranchAddress("MuonDzSignif_sum", &MuonDzSignif_sum[iT]);
    T[iT]->SetBranchAddress("Bc_M", &Bc_M[iT]);
    T[iT]->SetBranchAddress("Bc_Pt", &Bc_Pt[iT]);
    T[iT]->SetBranchAddress("QQ_M", &QQ_M[iT]);
    // if(iT<8){
    //   T[iT]->SetBranchAddress("QQ2_M", &QQ2_M[iT]);
    //   T[iT]->SetBranchAddress("QQ3_M", &QQ3_M[iT]);}
    T[iT]->SetBranchAddress("muW_isJpsiBro", &muW_isJpsiBro[iT]);
    T[iT]->SetBranchAddress("muW_trueId", &muW_trueId[iT]);
    T[iT]->SetBranchAddress("bkgType", &bkgType[iT]);
    T[iT]->SetBranchAddress("w_simple2", &w_simple2[iT]);
    T[iT]->SetBranchAddress("eventNb", &eventNb[iT]);

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){//T[iT]->GetEntries()

      T[iT]->GetEntry(j);

      muonsInTightAcc[iT] = muW_inTightAcc[iT] + mumi_inTightAcc[iT] + mupl_inTightAcc[iT];
      muonsGlbInLooseAcc[iT] = (muW_isGlb[iT] && muW_inLooseAcc[iT]) + (mumi_isGlb[iT] && mumi_inLooseAcc[iT]) + (mupl_isGlb[iT] && mupl_inLooseAcc[iT]);
      Bc_log1min_VtxProb[iT] = TMath::Log(1-Bc_VtxProb[iT]+1e-6);
      QQ_log1min_VtxProb[iT] = TMath::Log(1-QQ_VtxProb[iT]+1e-6);
      QQ_log_dca[iT] = TMath::Log(QQ_dca[iT]+1e-8);
      Bc_log_ctauSignif3D[iT] = TMath::Log(Bc_ctauSignif3D[iT]+1e-8);
      log_MuonDxySignif_sum[iT] = TMath::Log(MuonDxySignif_sum[iT]+1e-8);

      // Classifier response
      if((eventNb[iT]%2)==1){ //need events trained on (eventNb%2)==0
	BDT[iT] = reader1[iT]->EvaluateMVA("BDT");
	BDTprob[iT] = reader1[iT]->GetProba("BDT");
	BDTrarity[iT] = reader1[iT]->GetRarity("BDT");
      }
      else{ //need events trained on (eventNb%2)==1
	BDT[iT] = reader2[iT]->EvaluateMVA("BDT");
	BDTprob[iT] = reader2[iT]->GetProba("BDT");
	BDTrarity[iT] = reader2[iT]->GetRarity("BDT");
      }

      //      cout<<BDT[iT]<<endl;
      b_BDT[iT]->Fill();
      b_BDTprob[iT]->Fill();
      b_BDTrarity[iT]->Fill();
    }
    //END event loop

    //T[iT]->Print(); 
    T[iT]->Write("",TObject::kOverwrite); //overwrite, or two versions of the trees are saved 
  } 
  //END loop on trees
  
  delete fullFile;
  
}
