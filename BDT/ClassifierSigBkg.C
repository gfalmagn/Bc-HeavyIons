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

void Classifier(bool ispp = true, bool firstHalf=true, int kinBin=0, bool secondStep=false){ //kinBin==-1 gathers the two half samples and adds mass variable, to plot correlation matrix
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
    loader.AddVariable("Bc_CorrM","Corrected mass","GeV",'F'); //"I" for discrete variable, 'F' is default //Different from Bc_CorrM only for bkgType==2 //_BcMcorr
  loader.AddSpectator("Bc_ctauSignif","Significance of 2D lifetime","",'F');
  loader.AddVariable("Bc_log_ctauSignif3D := TMath::Log(Bc_ctauSignif3D)","Significance of 3D lifetime","",'F'); 
  loader.AddVariable("Bc_alpha","#alpha 2D","",'F');
  if(_withTM) 
    loader.AddVariable("Bc_alpha3D","#alpha 3D","",'F'); 
  else loader.AddSpectator("Bc_alpha3D","#alpha 3D","",'F'); 
  loader.AddSpectator("log_MuonDxyzSignif_sum := TMath::Log(MuonDxyzSignif_sum+1e-2)","Sum of the significances of 3D displacements of the three muons to the PV","",'F'); 
  loader.AddSpectator("log_MuonDxyzSignif_min := TMath::Log(MuonDxyzSignif_min+1e-2)","Min of the significances of 3D displacements of the three muons to the PV","",'F'); 
  loader.AddVariable("log_muW_dxyzSignif := TMath::Log(muW_dxyzSignif+1e-2)","Significance of d_{PV,3D}(#mu_{W})","",'F'); 
  loader.AddVariable("Bc_log_VtxProb := TMath::Log(Bc_VtxProb)","B_{c} vertex probability","",'F');
  loader.AddSpectator("QQ_log_VtxProb := TMath::Log(QQ_VtxProb)","J/#psi vertex probability","",'F');
  if(_withTM) loader.AddVariable("QQ_log_dca := TMath::Log(QQ_dca+1e-2)","J/#psi dca","mm",'F');
  if(useVar_CorrWMass){
    loader.AddVariable("dR_sum","Sum(#Delta R(#mu_{i}#mu_{j}) )","",'F');//Different from Bc_CorrM only for bkgType==2 //_BcMcorr
    loader.AddVariable("dR_jpsiOverMuW","#DeltaR(J/#psi)/(#DeltaR(#mu_{W}#mu^{-})+#DeltaR(#mu_{W}#mu^{+}))","",'F');
    loader.AddSpectator("dR_jpsiMuW","#Delta R(J/#psi,#mu_{W})","",'F');
  }
  if(_withTM){
    loader.AddVariable("log_MuonDxySignif_sum := TMath::Log(MuonDxySignif_sum+1e-2)","Sum of the significances of transverse displacements of the three muons to the PV","",'F'); 
    loader.AddVariable("muonsGlbInLooseAcc := (muW_isGlb && muW_inLooseAcc) + (mumi_isGlb && mumi_inLooseAcc) + (mupl_isGlb && mupl_inLooseAcc)","Number of muons global and in loose acceptance","",'I');
    loader.AddVariable("muonsInTightAcc := muW_inTightAcc + mumi_inTightAcc + mupl_inTightAcc","Number of muons in tight acceptance","",'I');
  }
  loader.AddVariable("QQmuW_ptImbal","J/#psi-#mu_{W} p_{T} imbalance","",'F'); //"I" for discrete variable, 'F' is default


  if(kinBin==-1)
    loader.AddVariable("Bc_M","B_{c} mass","GeV",'F'); //JUST FOR CORRELATIONS, SHOULD NOT BE A VARIABLE FOR ACTUAL TRAINING
  else
    loader.AddSpectator("Bc_M","Trimu mass","GeV",'F');
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

  int kBin = (kinBin==-1)?0:kinBin;
  TString fidCut = "fabs(Bc_Y)>"+(TString)to_string(_BcYmin[kBin])+" && fabs(Bc_Y)<"+(TString)to_string(_BcYmax[kBin])+" && Bc_Pt>"+(TString)to_string(_BcPtmin[kBin])+" && Bc_Pt<"+(TString)to_string(_BcPtmax[kBin]); //fiducial cut
  TString cutSig = "Bc_M<6.2 && Bc_M>3.5 && "+fidCut;
  TString cutBkg = "((Bc_M<6.2 && Bc_M>3.5) || (bkgType==2)) && "+fidCut;
  TString cutHalf1 = " && (eventNb%"+(TString)to_string((kinBin==-1)?10:2)+")=="+(TString)(firstHalf?"0":"1"); //9/10th if kinBin==-1
  TString cutHalf2 = " && (eventNb%"+(TString)to_string((kinBin==-1)?10:2)+")!="+(TString)(firstHalf?"0":"1");

  loader.AddTree(tsig,"Signal",  1.,(TCut)(cutSig+cutHalf1), TMVA::Types::kTraining); //signal weight  = 1
  loader.AddTree(tsig,"Signal",  1.,(TCut)(cutSig+cutHalf2), TMVA::Types::kTesting);  //signal weight  = 1
  loader.AddTree(tbkgWS, "Background", 1.,(TCut)(cutBkg+cutHalf1), TMVA::Types::kTraining);   //background weight = 1 
  loader.AddTree(tbkgWS, "Background", 1.,(TCut)(cutBkg+cutHalf2), TMVA::Types::kTesting);   //background weight = 1 
  loader.AddTree(tbkgBCM, "Background", 1.,(TCut)(cutBkg+cutHalf1), TMVA::Types::kTraining); //background weight = 1 
  loader.AddTree(tbkgBCM, "Background", 1.,(TCut)(cutBkg+cutHalf2), TMVA::Types::kTesting);  //background weight = 1 
  // loader.AddTree(tbkgTRUEJ, "Background", 1.,(TCut)(cutBkg+cutHalf1), TMVA::Types::kTraining);//this sample is somewhat redundant, but keep it for more complete description of bkg
  // loader.AddTree(tbkgTRUEJ, "Background", 1.,(TCut)(cutBkg+cutHalf2), TMVA::Types::kTesting);
  loader.AddTree(tbkgFLIPJ, "Background", ispp?0.7:0.7,(TCut)(cutBkg+cutHalf1), TMVA::Types::kTraining);
  loader.AddTree(tbkgFLIPJ, "Background", ispp?0.7:0.7,(TCut)(cutBkg+cutHalf2), TMVA::Types::kTesting);
  if(addMCjpsi){
    loader.AddTree(tbkgBTOJPSI, "Background", ispp?1.75:2.6,(TCut)(cutBkg+(TString)" && muW_isJpsiBro"+cutHalf1), TMVA::Types::kTraining);
    loader.AddTree(tbkgBTOJPSI, "Background", ispp?1.75:2.6,(TCut)(cutBkg+(TString)" && muW_isJpsiBro"+cutHalf2), TMVA::Types::kTesting); 
    loader.AddTree(tbkgBTOJPSI, "Background", 0.6,(TCut)(cutBkg+(TString)" && !muW_isJpsiBro"+cutHalf1), TMVA::Types::kTraining);
    loader.AddTree(tbkgBTOJPSI, "Background", 0.6,(TCut)(cutBkg+(TString)" && !muW_isJpsiBro"+cutHalf2), TMVA::Types::kTesting); 
    loader.AddTree(tbkgPROMPTJPSI, "Background", 0.6,(TCut)(cutBkg+cutHalf1), TMVA::Types::kTraining);
    loader.AddTree(tbkgPROMPTJPSI, "Background", 0.6,(TCut)(cutBkg+cutHalf2), TMVA::Types::kTesting); 
  }
  
  if(!secondStep)
    loader.SetWeightExpression( "w_simple2" );
  else{
    loader.SetSignalWeightExpression( "weight2" );
    loader.SetBackgroundWeightExpression( "weight" );
  }

  loader.PrepareTrainingAndTestTree("","", "SplitMode=Random:SplitSeed=101:NormMode=None:!V" ); //nTrain_Signal=0 is the default (split in half for training and testing) //SplitSeed=100 ensures to have the same splitting each time tmva is ran. Use SplitSeed=0 for a seed changing every time. //SplitMode is actually dummy because the train and test samples are defined before
  

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
		       "BDTfiner"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+"_withTM"+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
		       (ispp?"!V:NTrees=700:MinNodeSize=3.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs":
			"!V:NTrees=500:MinNodeSize=4.%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs"));
    // auto method = factory.BookMethod(&loader,TMVA::Types::kBDT, "BDTfinerShallow"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
    // 				   (ispp?"!V:NTrees=1300:MinNodeSize=3.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs":
    // 				    "!V:NTrees=900:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40:CreateMVAPdfs") );
  }
   
  else{
    factory.BookMethod(&loader,TMVA::Types::kBDT, 
		       "BDT"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"") +(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
		       (TString)(ispp?"!V:NTrees=120:MinNodeSize=10%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=15:CreateMVAPdfs":"!V:NTrees=120:MinNodeSize=10%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs") );
    factory.BookMethod(&loader,TMVA::Types::kBDT, 
		       "BDTfiner"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
		       (TString)(ispp?"!V:NTrees=350:MinNodeSize=10%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=10:CreateMVAPdfs":"!V:NTrees=350:MinNodeSize=10%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs"));
    auto method = factory.BookMethod(&loader,TMVA::Types::kBDT, 
				     "BDTfinerShallow"+(TString)(ispp?"":"_PbPb") +(TString)(addMCjpsi?"_withJpsiMC":"")+(TString)(useVar_CorrWMass?"":"_dropVarCorrWMass")+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(kinBin)+(TString)(firstHalf?"_1stHalf":"_2ndHalf"),
				     (TString)(ispp?"!V:NTrees=600:MinNodeSize=7%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30:CreateMVAPdfs":"!V:NTrees=600:MinNodeSize=7%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=15:CreateMVAPdfs"));
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



void ClassifierSigBkg(bool ispp=true, bool secondStep=false, bool forCorrMatrix=false){
  if(forCorrMatrix){
    Classifier(ispp,true,-1); //first half is 9/10th in this case
  }
  else{
    Classifier(ispp,true,1,secondStep);
    cout<<endl<<endl;
    Classifier(ispp,false,1,secondStep);
    cout<<endl<<endl;
    Classifier(ispp,true,2,secondStep);
    cout<<endl<<endl;
    Classifier(ispp,false,2,secondStep);
  }
}
