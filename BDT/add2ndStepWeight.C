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
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"

void add2ndStepWeight(bool ispp=true, bool withTM=false){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  //**************************************************************  
  //Grab the variations of the pT bias of MC, from first step r1 and r2
  TFile *BiasFile = TFile::Open("../twoSteps/pTBiases.root","READ");
  //var1 is the variation corresponding to pT^{n+m log(pT)}, fixed m
  TH1F* biasPTMC = (TH1F*)BiasFile->Get("pTbias_"+(TString)(ispp?"pp":"PbPb")+"_var"+(TString)to_string(_nomMethVar));
    
  auto fullFile = TFile::Open("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","UPDATE");
  int ntrees = 9;//ispp?10:9;

  //Initialization of variables
  // float Bc_M[ntrees];
  // float Bc_Pt[ntrees];
  float gen3mu_Pt[ntrees];
  float w_simple2[ntrees];
  float weight[ntrees];
  float weight2[ntrees];

  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi","flipJpsibMC"};
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }

  vector<TBranch*> b_weight2;

  //*******************************************
  //For the signal region and the tmva output, fill the branches for BDT
  for(int iT=0; iT<(int)ntrees; iT++){
    std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    b_weight2.push_back( T[iT]->Branch("weight2",&weight2[iT],"weight2/F") );
    T[iT]->SetBranchAddress("weight", &weight[iT]);
    if(iT==4) T[iT]->SetBranchAddress("gen3mu_Pt", &gen3mu_Pt[iT]);

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){//T[iT]->GetEntries()

      T[iT]->GetEntry(j);

      if(iT==4){
	float addw = getBias( biasPTMC , gen3mu_Pt[iT]);
	if(addw<=1e-7) addw = 1;
	weight2[iT] = weight[iT] * addw;
	b_weight2[iT]->Fill();
      }
    }
    //END event loop

    T[iT]->Write("",TObject::kOverwrite); //overwrite, or two versions of the trees are saved 
  } 
  //END loop on trees
  
  delete fullFile;
  
}
