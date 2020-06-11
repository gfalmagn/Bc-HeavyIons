#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/CrossValidation.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodBase.h"
#include "TMVA/Types.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

void Classifier(bool ispp = true, bool firstHalf=true, int kinBin=0){
  bool useVar_CorrWMass = true;

  //Declare Factory
  TMVA::Tools::Instance();
  //  ROOT::R::TRInterface &r = ROOT::R::TRInterface::Instance();

  auto inputFile = TFile::Open("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root");
  auto outputFile = TFile::Open("TMVAOutputCV_"+(TString)(ispp?"pp":"PbPb")+".root", "RECREATE");

  //Declare DataLoader(s)

  TMVA::DataLoader loader("dataset");
  //TMVA::DataLoader* transformed_loader1 = loader->VarTransform("VT(3)"); //keeps only variables different by 3 sigmas between sig and bkg

  if(useVar_CorrWMass) 
    loader.AddVariable("Bc_CorrM_shiftedM","B_{c} corrected mass","GeV",'F'); //"I" for discrete variable, 'F' is default //Different from Bc_CorrM only for bkgType==2 //_BcMcorr
  loader.AddSpectator("Bc_ctauSignif","Significance of non-zero B_{c} transverse lifetime","",'F');
  loader.AddVariable("Bc_log_ctauSignif3D := TMath::Log(Bc_ctauSignif3D+1e-8)","Significance of non-zero B_{c} 3D lifetime","",'F');
  loader.AddSpectator("Bc_alpha","B_{c} #alpha","",'F');
  if(_withTM) loader.AddVariable("Bc_alpha3D","B_{c} #alpha 3D","",'F');
  loader.AddVariable("Bc_log1min_VtxProb := TMath::Log(1-Bc_VtxProb+1e-6)","B_{c} vertex probability","",'F');
  loader.AddSpectator("QQ_log1min_VtxProb := TMath::Log(1-QQ_VtxProb+1e-6)","J/#psi vertex probability","",'F');
  if(_withTM) loader.AddVariable("QQ_log_dca := TMath::Log(QQ_dca+1e-8)","J/#psi dimuon distance of closest approach","mm",'F');
  if(useVar_CorrWMass){
    loader.AddVariable("dR_sum_shiftedM","Sum of #Delta R between the three pairs of muons","",'F');//Different from Bc_CorrM only for bkgType==2 //_BcMcorr
    loader.AddVariable("dR_jpsiOverMuW_shiftedM","Ratio of #Delta R of J/#psi over (#Delta R (#mu_{W}-#mu^{-}) + #Delta R (#mu_{W}-#mu^{+}))","",'F');
    loader.AddSpectator("dR_jpsiMuW_shiftedM","#Delta R(J/#psi,#mu_{W})","",'F');
  }
  loader.AddSpectator("MuonDzSignif_sum","Sum of the significances of longitudinal displacements of the three muons to the PV","",'F');
  if(_withTM){
    loader.AddVariable("log_MuonDxySignif_sum := TMath::Log(MuonDxySignif_sum+1e-8)","Sum of the significances of transverse displacements of the three muons to the PV","",'F');
    loader.AddVariable("muonsGlbInLooseAcc := (muW_isGlb && muW_inLooseAcc) + (mumi_isGlb && mumi_inLooseAcc) + (mupl_isGlb && mupl_inLooseAcc)","Number of muons global and in loose acceptance","",'I');
    loader.AddVariable("muonsInTightAcc := muW_inTightAcc + mumi_inTightAcc + mupl_inTightAcc","Number of muons in tight acceptance","",'I');
  }
  loader.AddVariable("QQmuW_ptImbal","J/#psi-#mu_{W} p_{T} imbalance","",'F'); //"I" for discrete variable, 'F' is default
  loader.AddSpectator("Bc_M","B_{c} mass","GeV",'F'); //JUST FOR CORRELATIONS, SHOULD NOT BE A VARIABLE FOR ACTUAL TRAINING
  loader.AddSpectator("Bc_Pt","B_{c} P_{T}","GeV",'F');
  loader.AddSpectator("QQ_M","J/#psi mass","GeV",'F');
  // loader.AddSpectator("QQ2_M","J/#psi 2 mass","GeV",'F');
  // loader.AddSpectator("QQ3_M","J/#psi 3 mass","GeV",'F');
  loader.AddSpectator("bkgType","background type","",'I');
  loader.AddSpectator("muW_isJpsiBro","muW and Jpsi from same parent decay","",'I');
  loader.AddSpectator("muW_trueId","muW gen pdgID","",'I');
  loader.AddSpectator("w_simple2","w_simple2","",'F');

  //Setup Dataset(s)

  bool addMCjpsi = true;

  TTree *tsig, *tbkgWS, *tbkgBCM, *tbkgTRUEJ, *tbkgBTOJPSI, *tbkgPROMPTJPSI, *tbkgFLIPJ;
  inputFile->GetObject("signal_MC", tsig);
  inputFile->GetObject("bkgWRONGSIGN", tbkgWS);
  inputFile->GetObject("bkgBCMASS", tbkgBCM);
  inputFile->GetObject("bkgTRUEJPSI", tbkgTRUEJ);
  inputFile->GetObject("bToJpsi_MC", tbkgBTOJPSI);
  inputFile->GetObject("PromptJpsi_MC", tbkgPROMPTJPSI);
  inputFile->GetObject("flipJpsi", tbkgFLIPJ);

  TString fidCut = "fabs(Bc_Y)>"+(TString)to_string(_BcYmin[kinBin])+" && fabs(Bc_Y)<"+(TString)to_string(_BcYmax[kinBin])+" && Bc_Pt>"+(TString)to_string(_BcPtmin[kinBin])+" && Bc_Pt<"+(TString)to_string(_BcPtmax[kinBin]); //fiducial cut
  TString cutSig = "Bc_M<6.25 && "+fidCut;
  TString cutBkg = "((Bc_M<6.25) || (bkgType==2)) && "+fidCut;

  loader.AddTree(tsig,"Signal",  1.,(TCut)(cutSig+" && (eventNb%2)=="+(TString)(firstHalf?"0":"1")), TMVA::Types::kTraining); //signal weight  = 1
  loader.AddTree(tsig,"Signal",  1.,(TCut)(cutSig+" && (eventNb%2)=="+(TString)(firstHalf?"1":"0")), TMVA::Types::kTesting);  //signal weight  = 1
  loader.AddTree(tbkgWS, "Background", 1.,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"0":"1")), TMVA::Types::kTraining);   //background weight = 1 
  loader.AddTree(tbkgWS, "Background", 1.,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"1":"0")), TMVA::Types::kTesting);   //background weight = 1 
  loader.AddTree(tbkgBCM, "Background", 1.,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"0":"1")), TMVA::Types::kTraining); //background weight = 1 
  loader.AddTree(tbkgBCM, "Background", 1.,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"1":"0")), TMVA::Types::kTesting);  //background weight = 1 
  // loader.AddTree(tbkgTRUEJ, "Background", 1.,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"0":"1")), TMVA::Types::kTraining);//this sample is somewhat redundant, but keep it for more complete description of bkg
  // loader.AddTree(tbkgTRUEJ, "Background", 1.,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"1":"0")), TMVA::Types::kTesting);
  if(ispp){
    loader.AddTree(tbkgFLIPJ, "Background", ispp?0.65:0,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"0":"1")), TMVA::Types::kTraining);
    loader.AddTree(tbkgFLIPJ, "Background", ispp?0.65:0,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"1":"0")), TMVA::Types::kTesting);
  }
  if(addMCjpsi){
    loader.AddTree(tbkgBTOJPSI, "Background", ispp?1.7:1.5,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"0":"1")), TMVA::Types::kTraining);
    loader.AddTree(tbkgBTOJPSI, "Background", ispp?1.7:1.5,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"1":"0")), TMVA::Types::kTesting); 
    loader.AddTree(tbkgPROMPTJPSI, "Background", ispp?1.:1.5,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"0":"1")), TMVA::Types::kTraining);
    loader.AddTree(tbkgPROMPTJPSI, "Background", ispp?1.:1.5,(TCut)(cutBkg+" && (eventNb%2)=="+(TString)(firstHalf?"1":"0")), TMVA::Types::kTesting); 
  }
  loader.SetWeightExpression( "w_simple2" );

  loader.PrepareTrainingAndTestTree("","", "SplitMode=Random:SplitSeed=101:NormMode=None:!V" ); //nTrain_Signal=0 is the default (split in half for training and testing) //SplitSeed=100 ensures to have the same splitting each time tmva is ran. Use SplitSeed=0 for a seed changing every time.
  

  // // Setup cross-validation 
  // TMVA::CrossValidation cv(&loader);
  // TMVA::CrossValidation *cv = new TMVA::CrossValidation(loader);

  // cv->BookMethod(TMVA::Types::kBDT, "BDTfinerShallow"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":""),
  // 		 "!V:NTrees=900:MinNodeSize=6%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30:CreateMVAPdfs");

  // // Run cross-validation and print results 
  // cv->SetNumFolds(3);
  // cv->Evaluate();
  // const TMVA::CrossValidationResult &results = cv->GetResults();
  // results.Print();
  // results.Draw();




  //  ****************** Simple Factory
  TMVA::Factory factory("BcSigBkgClassification", outputFile,
  			"!V:ROC:!Correlations:!Silent:Color:DrawProgressBar:AnalysisType=Classification" ); //Transformations=I;D;P;G,D:

  //Boosted Decision Trees
  if(_withTM){
    // factory.BookMethod(&loader,TMVA::Types::kBDT, "BDT"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"") +(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
    // 		     "!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
    factory.BookMethod(&loader,TMVA::Types::kBDT, 
		       "BDTfiner"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+"__withTM"+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
		       (ispp?"!V:NTrees=700:MinNodeSize=3.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs":
			"!V:NTrees=500:MinNodeSize=4.%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs"));
    // auto method = factory.BookMethod(&loader,TMVA::Types::kBDT, "BDTfinerShallow"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
    // 				   (ispp?"!V:NTrees=1300:MinNodeSize=3.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs":
    // 				    "!V:NTrees=900:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs") );
  }
   
  else{
    factory.BookMethod(&loader,TMVA::Types::kBDT, 
		       "BDT"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"") +(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
		       (TString)(ispp?"!V:NTrees=120:MinNodeSize=10%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=15:CreateMVAPdfs":"!V:NTrees=120:MinNodeSize=10%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs") );
    factory.BookMethod(&loader,TMVA::Types::kBDT, 
		       "BDTfiner"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
		       (TString)(ispp?"!V:NTrees=350:MinNodeSize=10%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=10:CreateMVAPdfs":"!V:NTrees=350:MinNodeSize=10%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs"));
    auto method = factory.BookMethod(&loader,TMVA::Types::kBDT, 
				     "BDTfinerShallow"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
				     (TString)(ispp?"!V:NTrees=600:MinNodeSize=7%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs":"!V:NTrees=600:MinNodeSize=7%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFrac\
tion=0.5:SeparationType=GiniIndex:nCuts=15:CreateMVAPdfs"));
				     }
  

  //  method->OptimizeTuningParameters();

  // factory.BookMethod(TMVA::Types::kRXGB, "RXGBdefault", "!V:NRounds=10:MaxDepth=6:Eta=0.3");
  // factory.BookMethod(TMVA::Types::kRXGB, "RXGB", "!V:NRounds=80:MaxDepth=2:Eta=1");

  //Multi-Layer Perceptron (Neural Network)
  // factory.BookMethod(&loader, TMVA::Types::kMLP, "MLP",
  // 		     "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=100:HiddenLayers=N+5:TestRate=5:!UseRegulator:CreateMVAPdfs" );




  //Train Methods
  factory.TrainAllMethods();

  //Test and Evaluate Methods
  factory.TestAllMethods();
  factory.EvaluateAllMethods();


  //Plot ROC Curve

  //We enable JavaScript visualisation for the plots

  //%jsroot on

  auto c1 = factory.GetROCCurve(&loader);
  c1->Draw();

  outputFile->Close();
  //if (!gROOT->IsBatch() && !firstHalf) TMVA::TMVAGui(  "TMVAOutputCV_"+(TString)(ispp?"pp":"PbPb")+".root" );

  // // Clean up
  //  delete factory;
  //  r.SetVerbose(1);
  //  delete dataloader;
 
}



void ClassifierSigBkg(bool ispp=true){
  Classifier(ispp,true,0);
  Classifier(ispp,false,0);
  Classifier(ispp,true,1);
  Classifier(ispp,false,1);
}
