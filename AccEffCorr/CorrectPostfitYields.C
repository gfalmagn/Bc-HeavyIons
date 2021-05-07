#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"

//first syst name corresponds to the suffix on InputForCombine...root, and the second is the additional suffix on fitDiagnostics...root
void CorrectPostfitYields(bool ispp = true, bool secondStep=false, TString metafitSyst="", TString systExt="", bool verbose=true){ 

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  bool RescalePredToData = true; //rescale the sig+bkg in each trimuon mass bin to the actual data yield in this bin
  bool limitAccEffCorrErrors = true;
  bool truncateLargeCorr = false; //impose a ceiling on the large corrections

  const int nproc = 5;
  TString procNameDef[] = {"BcSig","data_obs","FakeJpsi",(ispp?"JpsiMC":"NPJpsi"),"flipJpsi"};
  vector<vector<TString> > systName(nproc);
  if(ispp)
    systName = {{""},{""}, {"","_JpsiSBUp","_JpsiSBDown"}, 
		{"","_flipJSameSideUp","_flipJSameSideDown","_wPromptMCUp","_wPromptMCDown"}   ,   {"","_flipJSameSideUp","_flipJSameSideDown","_wPromptMCUp","_wPromptMCDown"} };
  else 
    systName = {{""},{""}, {"","_JpsiSBUp","_JpsiSBDown"},
		{"","_FlipJorMCUp","_FlipJorMCDown","_UncorrNPJUp","_UncorrNPJDown"}   ,   {"","_FlipJorMCUp","_FlipJorMCDown","_UncorrNPJUp","_UncorrNPJDown"} }; //PbPb

  vector<pair<int,int> > allProc;
  for(int i=0; i<nproc; i++){
    for(int sys=0; sys<(int)systName[i].size(); sys++){
      allProc.push_back(make_pair(i,sys));
    }
  }

  vector<vector<vector<TH1F*> > > h_AECorr(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_AECorr_BDT23(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_AECorr_BDT3(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_postfit(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_postfit_cent(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_prefit(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_AEweighted(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_AEweighted_BDT23(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_AEweighted_BDT3(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<float> > > Yields_pref(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1))); //prefit
  vector<vector<vector<float> > > Yields_postf(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1))); //postfit, except for data
  vector<vector<vector<float> > > Yields_postf_cent(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1))); //postfit, except for data
  vector<vector<float> > corrYield(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_MC(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_MCv2(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_BDT23(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_BDT3(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYieldErr(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYieldErr_BDT23(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYieldErr_BDT3(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<float> r1r2Corr(1,0);
  vector<float> r1r2Corr_cent(1,0);
  vector<float> rsig(_NanaBins+1,0), rsig_true(_NanaBins+1,0), rsig_errl(_NanaBins+1,0), rsig_errh(_NanaBins+1,0);
  vector<vector<float> > rsig_relerr(_NanaBins+1,vector<float>(3,0));
  vector<float> rsig_cent(_NanaBins+1,0), rsig_errl_cent(_NanaBins+1,0), rsig_errh_cent(_NanaBins+1,0);
  vector<vector<float> > rsig_relerr_cent(_NanaBins+1,vector<float>(3,0));
  vector<vector<TH1F*> > h_BcM_BkgSubData(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_BkgSubData_BDT23(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_BkgSubData_BDT3(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_Bkg(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_Data(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_DataMinBkg(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_SigMC(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<vector<TH1F*> > h_BcM_SigMCv2(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1));
  vector<TCanvas*> can(_NanaBins+1);
  vector<TLegend*> leg(_NanaBins+1);

  //********************************************************
  // Extract mass histograms before and after AccEff corrections, for all samples, including systematic variations (allProc)
  //********************************************************  
  auto histFile = TFile::Open("../templateFit/InputForCombine_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+metafitSyst+".root"); 
  cout<<"Extract (mass-dependent) acceptance x efficiency corrections from file "<<"../templateFit/InputForCombine_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+metafitSyst+".root"<<endl;

  for(int p=0; p<(int)allProc.size(); p++){
    int MainProc = allProc[p].first;
    int sys = allProc[p].second;
    for(int b=0;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){
	h_AECorr[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffCorr");
	h_AECorr_BDT23[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffCorr_BDTeff23");
	h_AECorr_BDT3[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffCorr_BDTeff3");
	h_BcM_AEweighted[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffWeighted"); //grab BcM, corrected by event-by-event acc-eff
	h_BcM_AEweighted_BDT23[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffWeighted_BDTeff23"); //grab BcM, corrected by event-by-event acc-eff
	h_BcM_AEweighted_BDT3[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffWeighted_BDTeff3"); //grab BcM, corrected by event-by-event acc-eff

	if(sys==0) {
	  h_BcM_prefit[MainProc][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM"); //prefit histo
	  h_BcM_prefit[MainProc][b][k]->SetDirectory(0);
	}
	if(p==1) {
	  h_BcM_postfit[MainProc][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM"); //grab data BcM histo
	  h_BcM_postfit[MainProc][b][k]->SetDirectory(0);
	  if(!ispp){
	    h_BcM_postfit_cent[MainProc][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+(TString)((b==0)?"":("Cent"+(TString)to_string(b)))+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM"); //grab data BcM histo
	    h_BcM_postfit_cent[MainProc][b][k]->SetDirectory(0);
	  }
	}

	h_AECorr[p][b][k]->SetDirectory(0);
	h_AECorr_BDT23[p][b][k]->SetDirectory(0);
	h_AECorr_BDT3[p][b][k]->SetDirectory(0);
	h_BcM_AEweighted[p][b][k]->SetDirectory(0);
	h_BcM_AEweighted_BDT23[p][b][k]->SetDirectory(0);
	h_BcM_AEweighted_BDT3[p][b][k]->SetDirectory(0);
      }
    }
  }
  
  histFile->Close();

  //********************************************************
  //Extract POSTFIT SHAPES and SHAPE NUISANCE PARAMETERS from combine output
  //******************************************************** 
  TString normFileName = "../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_2bins"+(TString)(secondStep?"_2ndStep":"")+metafitSyst+systExt+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
  TString normFileIntName = "../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_integrated"+(TString)(secondStep?"_2ndStep":"")+metafitSyst+systExt+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
  TString normFileCentName = "../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_centBins"+(TString)(secondStep?"_2ndStep":"")+metafitSyst+systExt+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
  cout<<"Extract normalisations from files "<<normFileName<<" "<<(ispp?"":normFileCentName)<<" "<<normFileIntName<<" "<<endl;
  auto normFile = TFile::Open(normFileName);
  auto normFileInt = TFile::Open(normFileIntName);
  auto normFileCent= TFile::Open(ispp?normFileName:normFileCentName);

  for(int proc=0; proc<nproc; proc++){
    if(proc==1) continue; //data is not a postfit shape
    for(int b=0;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){
	h_BcM_postfit[proc][b][k] = (TH1F*)h_AECorr[1][b][k]->Clone("Bc_M_postfit_"+procNameDef[proc]); //Do this to recover proper X axis //not essential, but prettier	
	h_BcM_postfit[proc][b][k]->SetDirectory(0);
	if(!ispp){
	  h_BcM_postfit_cent[proc][b][k] = (TH1F*)h_AECorr[1][b][k]->Clone("Bc_M_postfit_"+procNameDef[proc]); //Do this to recover proper X axis //not essential, but prettier	
	  h_BcM_postfit_cent[proc][b][k]->SetDirectory(0);
	}

	//fetch postfit shapes                                                                           
	//cout<<"fetching histogram shapes_fit_s/BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[proc]<<endl;
	TH1F* h_tmp = (TH1F*)((b==0)?normFileInt:normFile)->Get("shapes_fit_s/BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/"+procNameDef[proc]) ;
	for(int bin=0;bin<=h_BcM_postfit[proc][b][k]->GetNbinsX();bin++){
	  h_BcM_postfit[proc][b][k]->SetBinContent(bin, h_tmp->GetBinContent(bin));
	  h_BcM_postfit[proc][b][k]->SetBinError(bin, h_tmp->GetBinError(bin));
	}
	delete h_tmp;

	if(!ispp){
	  TH1F* h_tmp = (TH1F*)((b==0)?normFileInt:normFileCent)->Get("shapes_fit_s/BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Cent"+(TString)to_string(b)))+"/"+procNameDef[proc]) ;
	  for(int bin=0;bin<=h_BcM_postfit_cent[proc][b][k]->GetNbinsX();bin++){
	    h_BcM_postfit_cent[proc][b][k]->SetBinContent(bin, h_tmp->GetBinContent(bin));
	    h_BcM_postfit_cent[proc][b][k]->SetBinError(bin, h_tmp->GetBinError(bin));
	  }
	  delete h_tmp;
	}

      }
    }
  }

  //get values of nuisance parameters
  RooArgList fittedPars = ((RooFitResult*)normFile->Get("fit_s"))->floatParsFinal();
  RooArgList fittedParsInt = ((RooFitResult*)normFileInt->Get("fit_s"))->floatParsFinal();
  RooArgList fittedParsCent = ispp?RooArgList():( ((RooFitResult*)normFileCent->Get("fit_s"))->floatParsFinal() );
  float JpsiSB = ((RooRealVar*)fittedPars.find("JpsiSB"))->getValV(); //getError()
  float flipJ = ((RooRealVar*)fittedPars.find(ispp?"flipJSameSide":"FlipJorMC"))->getValV(); //getError()
  float JMC = ispp?( ((RooRealVar*)fittedPars.find(ispp?"wPromptMC":"UncorrNPJ"))->getValV() ):0; //getError()
  r1r2Corr[0] = ((RooFitResult*)normFile->Get("fit_s"))->correlation("r1","r2");
  if(!ispp) r1r2Corr_cent[0] = ((RooFitResult*)normFileCent->Get("fit_s"))->correlation("r1","r2");

  //Fill yields pre/post-fit (integrals of BcM distros)
  for(int proc=0; proc<nproc; proc++){
    for(int b=0;b<=_NanaBins;b++){
      Yields_pref[proc][b][0] = 0;
      Yields_postf[proc][b][0] = 0;
      for(int k=1;k<=_nChan(ispp);k++){
	Yields_pref[proc][b][k] = h_BcM_prefit[proc][b][k]->Integral(1,_nbinMSR(ispp)[k-1]) ;
	Yields_pref[proc][b][0] += Yields_pref[proc][b][k];
	Yields_postf[proc][b][k] = h_BcM_postfit[proc][b][k]->Integral(1,_nbinMSR(ispp)[k-1]) ;
	Yields_postf[proc][b][0] += Yields_postf[proc][b][k];
        if(!ispp){
	  Yields_postf_cent[proc][b][k] = h_BcM_postfit_cent[proc][b][k]->Integral(1,_nbinMSR(ispp)[k-1]) ;
	  Yields_postf_cent[proc][b][0] += Yields_postf_cent[proc][b][k];
	}

      }
    }
  }

  for(int b=0;b<=_NanaBins;b++){
    rsig[b] = ((RooRealVar*)((b==0)?fittedParsInt:fittedPars).find("r"+(TString)to_string(max(b,1))))->getValV();
    rsig_errl[b] = ((RooRealVar*)((b==0)?fittedParsInt:fittedPars).find("r"+(TString)to_string(max(b,1))))->getErrorLo();
    rsig_errh[b] = ((RooRealVar*)((b==0)?fittedParsInt:fittedPars).find("r"+(TString)to_string(max(b,1))))->getErrorHi();
    rsig_true[b] = Yields_postf[0][b][0]/Yields_pref[0][b][0];
    rsig_relerr[b][1] = -rsig_errl[b]/rsig[b];
    rsig_relerr[b][2] =  rsig_errh[b]/rsig[b];
    rsig_relerr[b][0] =  (rsig_relerr[b][1]+rsig_relerr[b][2])/2;

    if(!ispp){
      rsig_cent[b] = ((RooRealVar*)((b==0)?fittedParsInt:fittedParsCent).find("r"+(TString)to_string(max(b,1))))->getValV();
      rsig_errl_cent[b] = ((RooRealVar*)((b==0)?fittedParsInt:fittedParsCent).find("r"+(TString)to_string(max(b,1))))->getErrorLo();
      rsig_errh_cent[b] = ((RooRealVar*)((b==0)?fittedParsInt:fittedParsCent).find("r"+(TString)to_string(max(b,1))))->getErrorHi();
      rsig_relerr_cent[b][1] = -rsig_errl_cent[b]/rsig_cent[b];
      rsig_relerr_cent[b][2] =  rsig_errh_cent[b]/rsig_cent[b];
      rsig_relerr_cent[b][0] =  (rsig_relerr_cent[b][1]+rsig_relerr_cent[b][2])/2;
    }
  }

  cout<<"Shape morphing parameters (JpsiSB, flipJSameSide/FlipJorMC, wPromptMC/UncorrNPJ) are found to be "<<JpsiSB<<" "<<flipJ<<" "<<JMC<<endl;
  cout<<"Signal normalisation modifiers for the 2 pt bins are found to be "<<rsig[1]<<" "<<rsig[2]<<endl;
  if(!ispp) cout<<"Signal normalisation modifiers for the 2 centrality bins are found to be "<<rsig_cent[1]<<" "<<rsig_cent[2]<<endl;
  cout<<"           versus actual ratio of postfit and prefit signal normalisations = "<<rsig_true[1]<<" "<<rsig_true[2]<<endl;

  //only for writing it out
  RooArgSet *Yields = (RooArgSet*)normFile->Get("norm_fit_s");
  RooArgSet *YieldsInt = (RooArgSet*)normFileInt->Get("norm_fit_s");
  RooArgSet *YieldsCent = NULL;
  if(!ispp) YieldsCent = (RooArgSet*)normFileCent->Get("norm_fit_s");
  vector<float> nsig = vector<float>(_NanaBins+1,0);
  vector<float> nsig_cent = vector<float>(_NanaBins+1,0);
  for(int b=0;b<=_NanaBins;b++){
    for(int k=1;k<=_nChan(ispp);k++){
      nsig[b] += ((b==0)?YieldsInt:Yields)->getRealValue("BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+"/BcSig");
      if(!ispp)
	nsig_cent[b] += ((b==0)?YieldsInt:YieldsCent)->getRealValue("BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Cent"+(TString)to_string(b)))+"/BcSig");
    }
  }
  normFile->Close();

  if(metafitSyst=="" && systExt==""){
    //********************************************************
    //Subtract background from data, considering AccEff-corrected yields
    //******************************************************** 
    for(int b=0;b<=_NanaBins;b++){
      if(verbose) cout<<endl<<"b = "<<b<<endl;
      corrYield[b][0] = 0;
      corrYield_BDT23[b][0] = 0;
      corrYield_BDT3[b][0] = 0;
      corrYield_MC[b][0] = 0;
      corrYield_MCv2[b][0] = 0;

      for(int k=1;k<=_nChan(ispp);k++){
	if(verbose) cout<<endl<<"k = "<<k<<endl;
	corrYield[b][k]=0;
	corrYield_BDT23[b][k]=0;
	corrYield_BDT3[b][k]=0;
	corrYield_MC[b][k]=0;
	corrYield_MCv2[b][k]=0;

	h_BcM_BkgSubData[b][k] = (TH1F*)h_BcM_AEweighted[1][b][k]->Clone("BackgroundSubtractedData_AEcorr_BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"));
	h_BcM_BkgSubData_BDT23[b][k] = (TH1F*)h_BcM_AEweighted_BDT23[1][b][k]->Clone("BackgroundSubtractedData_AEcorr_BDTeff23_BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"));
	h_BcM_BkgSubData_BDT3[b][k] = (TH1F*)h_BcM_AEweighted_BDT3[1][b][k]->Clone("BackgroundSubtractedData_AEcorr_BDTeff3_BDT"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"));

	h_BcM_Data[b][k] = new TH1F("corrData_bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"), "Data-background, corrected by AccEff, bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b))), _nbinMSR(ispp)[k-1], _mBcMin, _mBcMax);    
	h_BcM_Bkg[b][k] = new TH1F("corrBkg_bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"), "Scaled background, corrected by AccEff, bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b))), _nbinMSR(ispp)[k-1], _mBcMin, _mBcMax);    
	h_BcM_DataMinBkg[b][k] = new TH1F("corrDataMinBkg_bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"), "Data-background, corrected by AccEff, bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b))), _nbinMSR(ispp)[k-1], _mBcMin, _mBcMax);    
	h_BcM_SigMC[b][k] = new TH1F("corrSigMC_bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"), "Scaled signal MC corrected by AccEff, bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b))), _nbinMSR(ispp)[k-1], _mBcMin, _mBcMax);    
	h_BcM_SigMCv2[b][k] = new TH1F("corrSigMCv2_bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b)))+metafitSyst+systExt+(TString)(ispp?"pp":"PbPb"), "Signal MC (v2) corrected by AccEff, bin"+(TString)(to_string(k))+(TString)((b==0)?"":("Kin"+(TString)to_string(b))), _nbinMSR(ispp)[k-1], _mBcMin, _mBcMax);    
      
	for(int bin=1;bin<=h_BcM_BkgSubData[b][k]->GetNbinsX();bin++){ //forget underflow here?

	  if(verbose) cout<<"bin #"<<bin<<endl;      
	  //Rescale the sig+bkg prediction to exactly fit the data in this bin
	  float predY = 0;
	  vector<float> yield;
	  for(int proc=0; proc<nproc; proc++){
	    yield.push_back( h_BcM_postfit[proc][b][k]->GetBinContent(bin) );
	    if(yield[proc]<0) yield[proc] = 0;
	    if(proc!=1) predY += yield[proc];
	  }
	  if(RescalePredToData){
	    for(int proc=0; proc<nproc; proc++){
	      if(proc==1) continue; //could maybe maintain proc==0 here, to check if it changed sthg on corrYield_MC
	      yield[proc] *= yield[1] / predY; //correct by data over total prediction
	    }
	  }

	  vector<float> binBkg(_nChan(ispp),0);
	  vector<float> binBkgErr(_nChan(ispp),0);

	  //lambda function to grab the AE correction for a given process (for this precise bin)
	  auto AEcorr = [h_AECorr,h_AECorr_BDT23,h_AECorr_BDT3,b,k,bin,limitAccEffCorrErrors,truncateLargeCorr](int process, int BDTeffMode)->float{
	    int pout = process;
	    vector<vector<vector<TH1F*> > > h_AECorr_tmp = (BDTeffMode==0)?h_AECorr:((BDTeffMode==1)?h_AECorr_BDT23:h_AECorr_BDT3);
	    if(limitAccEffCorrErrors){
	      if(pout>2 && pout<5
		 && fabs(h_AECorr_tmp[pout][b][k]->GetBinContent(bin) - h_AECorr_tmp[2][b][k]->GetBinContent(bin)) < h_AECorr_tmp[pout][b][k]->GetBinError(bin) //error larger than difference between nominal and systematic shape
		 && h_AECorr_tmp[pout][b][k]->GetBinError(bin) > h_AECorr_tmp[2][b][k]->GetBinError(bin) ) //larger error for systematic than nominal
		pout = 2; //keep the nominal AE for this bin and this systematic
	      if(pout>5 && pout<11
		 && fabs(h_AECorr_tmp[pout][b][k]->GetBinContent(bin) - h_AECorr_tmp[5][b][k]->GetBinContent(bin)) < h_AECorr_tmp[pout][b][k]->GetBinError(bin) //error larger than difference between nominal and systematic shape
		 && h_AECorr_tmp[pout][b][k]->GetBinError(bin) > h_AECorr_tmp[5][b][k]->GetBinError(bin) ) //larger error for systematic than nominal
		pout = 5; //keep the nominal AE for this bin and this systematic
	      if(pout>10 && pout<16
		 && fabs(h_AECorr_tmp[pout][b][k]->GetBinContent(bin) - h_AECorr_tmp[10][b][k]->GetBinContent(bin)) < h_AECorr_tmp[pout][b][k]->GetBinError(bin) //error larger than difference between nominal and systematic shape
		 && h_AECorr_tmp[pout][b][k]->GetBinError(bin) > h_AECorr_tmp[10][b][k]->GetBinError(bin) ) //larger error for systematic than nominal
		pout = 10;
	    }
	    float res = h_AECorr_tmp[pout][b][k]->GetBinContent(bin);
	    if(truncateLargeCorr){
	      float maxcorr = (500/(float)k)/((BDTeffMode==0)?1.:((BDTeffMode==1)?0.8:0.4)); //600 is a bit arbitrary//changed from 480
	      if(res>maxcorr) {
		if(BDTeffMode==0) cout<<"behead to 600/k in BDTeffMode "<<BDTeffMode<<" process "<<pout<<" this large correction: "<<res<<endl;
		res = maxcorr; }
	    }
	    return res;
	  };

	  for(int m=0;m<3;m++){
	    //1 shape morphing parameter for Jpsi sidebands
	    float p_JSB = 0;//min((float)1,fabs(JpsiSB));
	    binBkg[m] += yield[2] * ( (1-p_JSB) * AEcorr(2,m) + p_JSB * AEcorr( 3+((JpsiSB>0)?0:1) ,m) );
	    binBkgErr[m] = sqrt(pow(binBkgErr[m],2) + pow(h_BcM_AEweighted[2][b][k]->GetBinError(bin),2)  * yield[1] / predY );
	    //cout<<"bkg#2 corrected yield = "<<yield[2] * ( (1-p_JSB) * AEcorr(2,m) + p_JSB * AEcorr( 3+((JpsiSB>0)?0:1) ,m) )<<endl;
	    //2 shape morphing parameters for Jpsi MC and for flipJpsi
	    float p_flipJ = 0;//min((float)1,fabs(flipJ));
	    float p_JMC = 0;//min((float)1,fabs(JMC));
	    binBkg[m] += yield[3] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(5,m) + p_flipJ * (1-p_JMC/2) * AEcorr(6+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(8+((JMC>0)?0:1) ,m) );
	    binBkgErr[m] = sqrt(pow(binBkgErr[m],2) + pow(h_BcM_AEweighted[5][b][k]->GetBinError(bin),2)  * yield[1] / predY );
	    //cout<<"bkg#3 corrected yield = "<<yield[3] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(5,m) + p_flipJ * (1-p_JMC/2) * AEcorr(6+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(8+((JMC>0)?0:1),m) )<<endl;
	    binBkg[m] += yield[4] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(10,m) + p_flipJ * (1-p_JMC/2) * AEcorr(11+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(13+((JMC>0)?0:1),m) );
	    binBkgErr[m] = sqrt(pow(binBkgErr[m],2) + pow(h_BcM_AEweighted[10][b][k]->GetBinError(bin),2)  * yield[1] / predY );
	    //cout<<"bkg#4 corrected yield = "<<yield[4] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(10,m) + p_flipJ * (1-p_JMC/2) * AEcorr(11+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(13+((JMC>0)?0:1),m) )<<endl;
	  }
	
	  if(verbose) cout<<"corrected data vs corrected bkg: "<<h_BcM_BkgSubData[b][k]->GetBinContent(bin)<<" "<<binBkg[0]<<endl;
	  if(verbose) cout<<"data-bkg vs MC signal v1 and v2: "<<h_BcM_BkgSubData[b][k]->GetBinContent(bin) - binBkg[0]<<" "<<yield[0] * h_AECorr[0][b][k]->GetBinContent(bin)<<" "<<rsig[b]*h_BcM_AEweighted[0][b][k]->GetBinContent(bin)<<endl;
	  if(verbose) cout<<" detail of MC val = "<<yield[0] <<" "<< h_AECorr[0][b][k]->GetBinContent(bin)<<endl;

	  h_BcM_BkgSubData[b][k]->SetBinContent(bin, h_BcM_BkgSubData[b][k]->GetBinContent(bin) - binBkg[0]); //substract corrected background from corrected data
	  h_BcM_BkgSubData[b][k]->SetBinError(bin, sqrt(pow(h_BcM_BkgSubData[b][k]->GetBinContent(bin),2) + pow(binBkgErr[0],2)) ); //forget the shape morphing subtleties in the error
	  h_BcM_BkgSubData_BDT23[b][k]->SetBinContent(bin, h_BcM_BkgSubData_BDT23[b][k]->GetBinContent(bin) - binBkg[1]); //substract corrected background from corrected data
	  h_BcM_BkgSubData_BDT23[b][k]->SetBinError(bin, sqrt(pow(h_BcM_BkgSubData_BDT23[b][k]->GetBinContent(bin),2) + pow(binBkgErr[1],2)) ); //forget the shape morphing subtleties in the error
	  h_BcM_BkgSubData_BDT3[b][k]->SetBinContent(bin, h_BcM_BkgSubData_BDT3[b][k]->GetBinContent(bin) - binBkg[2]); //substract corrected background from corrected data
	  h_BcM_BkgSubData_BDT3[b][k]->SetBinError(bin, sqrt(pow(h_BcM_BkgSubData_BDT3[b][k]->GetBinContent(bin),2) + pow(binBkgErr[2],2)) ); //forget the shape morphing subtleties in the error

	  if(h_BcM_BkgSubData[b][k]->GetBinCenter(bin)<_mBcMax){ 
	    corrYield[b][k] += h_BcM_BkgSubData[b][k]->GetBinContent(bin);
	    corrYield_BDT23[b][k] += h_BcM_BkgSubData_BDT23[b][k]->GetBinContent(bin);
	    corrYield_BDT3[b][k] += h_BcM_BkgSubData_BDT3[b][k]->GetBinContent(bin);
	    corrYield_MC[b][k] += yield[0] * h_AECorr[0][b][k]->GetBinContent(bin);
	    corrYield_MCv2[b][k] += rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinContent(bin);

	    corrYieldErr[b][k] += pow(h_BcM_BkgSubData[b][k]->GetBinError(bin),2);
	    corrYieldErr_BDT23[b][k] += pow(h_BcM_BkgSubData_BDT23[b][k]->GetBinError(bin),2);
	    corrYieldErr_BDT3[b][k] += pow(h_BcM_BkgSubData_BDT3[b][k]->GetBinError(bin),2);

	    //Fill histos for display
	    h_BcM_Data[b][k]->SetBinContent(bin , h_BcM_AEweighted[1][b][k]->GetBinContent(bin) );
	    h_BcM_Data[b][k]->SetBinError(bin , h_BcM_AEweighted[1][b][k]->GetBinError(bin) );
	    h_BcM_Bkg[b][k]->SetBinContent(bin , binBkg[0] );
	    h_BcM_Bkg[b][k]->SetBinError(bin , binBkgErr[0] );
	    h_BcM_DataMinBkg[b][k]->SetBinContent(bin , h_BcM_BkgSubData[b][k]->GetBinContent(bin) );
	    h_BcM_DataMinBkg[b][k]->SetBinError(bin , sqrt( pow(h_BcM_AEweighted[1][b][k]->GetBinError(bin),2) + pow(binBkgErr[0],2) ) );
	    h_BcM_SigMC[b][k]->SetBinContent(bin , yield[0] * h_AECorr[0][b][k]->GetBinContent(bin) );
	    h_BcM_SigMC[b][k]->SetBinError(bin , rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinError(bin));
	    h_BcM_SigMCv2[b][k]->SetBinContent(bin , rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinContent(bin) );
	    h_BcM_SigMCv2[b][k]->SetBinError(bin , rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinError(bin));
	  }

	}//end loop on mass bins
      
	if(verbose) cout<<"corrected yield (in mass signal region) for analysis bin #"<<b<<" and BDT bin "<<k<<" is "<<corrYield[b][k]<<"\n    vs MC v1 and v2: "<<corrYield_MC[b][k]<<" "<<corrYield_MCv2[b][k]<<endl;

	corrYield[b][0] += corrYield[b][k];
	if(k>1) corrYield_BDT23[b][0] += corrYield_BDT23[b][k];
	if(k>2) corrYield_BDT3[b][0] += corrYield_BDT3[b][k];
	corrYield_MC[b][0] += corrYield_MC[b][k];
	corrYield_MCv2[b][0] += corrYield_MCv2[b][k];

	corrYieldErr[b][0] += corrYieldErr[b][k];
	corrYieldErr[b][k] = sqrt(corrYieldErr[b][k]);
	if(k>1) corrYieldErr_BDT23[b][0] += corrYieldErr_BDT23[b][k];
	corrYieldErr_BDT23[b][k] = sqrt(corrYieldErr_BDT23[b][k]);
	if(k>2) corrYieldErr_BDT3[b][0] += corrYieldErr_BDT3[b][k];
	corrYieldErr_BDT3[b][k] = sqrt(corrYieldErr_BDT3[b][k]);
      
	if(k>2){
	  corrYieldErr[b][0] = sqrt(corrYieldErr[b][0]);
	  corrYieldErr_BDT23[b][0] = sqrt(corrYieldErr_BDT23[b][0]);
	  corrYieldErr_BDT3[b][0] = sqrt(corrYieldErr_BDT3[b][0]);
	}

      }//end loop on BDT bins
      cout<<"corrected yield (in mass signal region) for analysis bin #"<<b<<" is "<<corrYield[b][0]<<"\n    vs MC v1 and v2:"<<corrYield_MC[b][0]<<" "<<corrYield_MCv2[b][0]<<endl;
      cout<<"    vs BDTeff23 and BDTeff3: "<<corrYield_BDT23[b][0]<<" "<<corrYield_BDT3[b][0]<<endl;

    }//end loop on analysis bins

    //********************************************************
    //DRAWING
    //******************************************************** 
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TLine *vert = new TLine(_mBcMax,0,_mBcMax,1000);
    vert->SetLineWidth(3);
    vert->SetLineStyle(1);
    TLine *hor = new TLine(_mBcMax,0,_mBcMax,1000);
    hor->SetLineWidth(ispp?10:12);
    hor->SetLineStyle(1);
    TLine *vert0 = new TLine(_mBcMax,0,_mBcMax,1000);
    vert0->SetLineWidth(2);
    vert0->SetLineStyle(7);
    TLatex BDTbinT, title2;

    for(int b=1;b<=_NanaBins;b++){

      float maxC1 = 1.15 * h_BcM_Data[b][1]->GetMaximum();
      float minC2 = 1.5 * min(0., max(h_BcM_DataMinBkg[b][1]->GetMinimum(), -1000/1.5));
      float maxC2 = 1.5 * max(h_BcM_DataMinBkg[b][1]->GetMaximum(),h_BcM_DataMinBkg[b][2]->GetMaximum());
    
      can[b] = new TCanvas("c"+(TString)to_string(b),"c1",3000,2000);
      // can[b]->SetRightMargin(0.02);
      // can[b]->SetLeftMargin(0.04);
      can[b]->SetBottomMargin(0.15);
      can[b]->Divide(3,2,0.,0.);
    
      for(int k=1;k<=_nChan(ispp);k++){

	if(!ispp){
	  maxC1 = ((k==1)?1.15:2.3) * h_BcM_Data[b][k]->GetMaximum();
	  minC2 = 1.15 * min(0., h_BcM_DataMinBkg[b][k]->GetMinimum());
	  maxC2 = ((k==1)?1.15:2.3) * h_BcM_DataMinBkg[b][k]->GetMaximum();
	  if(fabs(minC2)>maxC2) maxC2=fabs(minC2);
	}

	can[b]->cd(k);
	gPad->SetTickx(2);
	if(k==_nChan(ispp)) gPad->SetRightMargin(0.001); 
	if(k>1 && ispp) h_BcM_Data[b][k]->GetYaxis()->SetTickLength(0);
	if(k>0 && !ispp) gPad->SetLeftMargin(0.085);
	if(k<_nChan(ispp)) gPad->SetRightMargin(0.001);
	h_BcM_Data[b][k]->GetYaxis()->SetRangeUser(0, maxC1);
	//h_BcM_Data[b][k]->SetTitle("#splitline{Data VS background,}{corrected by AccEff, bin"+(TString)to_string(b)+"}");
	h_BcM_Data[b][k]->SetLineWidth(3);
	h_BcM_Data[b][k]->SetLineColor(kRed);
	h_BcM_Bkg[b][k]->SetLineWidth(3);
	h_BcM_Bkg[b][k]->SetLineColor(kBlack);
	h_BcM_Data[b][k]->Draw("E");
	h_BcM_Bkg[b][k]->Draw("histEsame");
	h_BcM_Data[b][k]->Draw("Esame");
	gPad->Update();
	if(k>1)
	  vert->DrawLine((ispp?1.002:0.97)*gPad->GetUxmin(),gPad->GetUymin() - 0.15*(gPad->GetUymax()-gPad->GetUymin()),(ispp?1.002:0.97)*gPad->GetUxmin(),gPad->GetUymax());

	can[b]->cd(k+_nChan(ispp));
	gPad->SetTickx(2);
	if(k==_nChan(ispp)) gPad->SetRightMargin(0.001);
	if(k>1 && ispp) h_BcM_DataMinBkg[b][k]->GetYaxis()->SetTickLength(0);
	if(k>0 && !ispp) gPad->SetLeftMargin(0.085);
	if(k<_nChan(ispp)) gPad->SetRightMargin(0.001);
	h_BcM_DataMinBkg[b][k]->GetYaxis()->SetRangeUser(minC2, maxC2);

	//      h_BcM_DataMinBkg[b][k]->SetTitle("#splitline{Data#minusbackground VS signal MC,}{corrected by AccEff, bin"+(TString)to_string(b)+"};M(B_{c}) [GeV]");
	if(k==_nChan(ispp)) h_BcM_DataMinBkg[b][k]->GetXaxis()->SetTitle("M(B_{c}) [GeV]");
	h_BcM_DataMinBkg[b][k]->SetLineWidth(3);
	h_BcM_DataMinBkg[b][k]->SetLineColor(kBlue);
	h_BcM_DataMinBkg[b][k]->GetXaxis()->SetTitleSize(0.04);
	h_BcM_SigMC[b][k]->SetLineWidth(3);
	h_BcM_SigMC[b][k]->SetLineColor(kGreen);
	h_BcM_SigMCv2[b][k]->SetLineWidth(3);
	h_BcM_SigMCv2[b][k]->SetLineColor(kGreen+2);
	h_BcM_DataMinBkg[b][k]->Draw("E");
	h_BcM_SigMC[b][k]->Draw("histsame");
	h_BcM_SigMCv2[b][k]->Draw("histsame");
	gPad->Update();

	if(k>1)
	  vert->DrawLine((ispp?1.002:0.97)*gPad->GetUxmin(),gPad->GetUymin() - 0.15*(gPad->GetUymax()-gPad->GetUymin()) ,(ispp?1.002:0.97)*gPad->GetUxmin(),gPad->GetUymax());
	vert0->DrawLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
	hor->DrawLine(((k>0)?0.9:1.)*gPad->GetUxmin(),gPad->GetUymax()-0.005,((k<_nChan(ispp))?1.1:1.)*gPad->GetUxmax(),gPad->GetUymax()-0.005);
	BDTbinT.SetNDC();
	BDTbinT.SetTextSize(0.05);
	BDTbinT.DrawLatex(0.45,0.04,"BDT bin"+(TString)to_string(k));

      }

      can[b]->cd();
      TPad *pad = new TPad("pad","",0,0,1,1); //transparent pad to overlay text
      pad->Draw();
      pad->cd();
      pad->SetFillStyle(4000);

      title2.SetNDC();
      title2.SetTextSize(0.038);
      title2.DrawLatex(0.13,0.85,"#splitline{Data VS background,}{corrected by AccEff}");

      title2.SetNDC();
      title2.SetTextSize(0.038);
      title2.DrawLatex(0.13,0.39,"#splitline{Data#minusbackground VS signal MC,}{corrected by AccEff}");
    
      can[b]->cd(3);
      gPad->Update();
      leg[b] = new TLegend(0.12,0.4,0.85,0.97);
      leg[b]->SetBorderSize(0);
      leg[b]->SetFillStyle(0);
      leg[b]->SetTextSize(0.075);
      leg[b]->SetHeader(Form("%s, %.0f<p_{T}<%.0f GeV",ispp?"pp":"PbPb",_BcPtmin[b],_BcPtmax[b]));
      leg[b]->AddEntry(h_BcM_Data[b][1], "data","le");
      leg[b]->AddEntry(h_BcM_Bkg[b][1], "total background");
      leg[b]->AddEntry(h_BcM_DataMinBkg[b][1], "Data - background","le");
      leg[b]->AddEntry(h_BcM_SigMC[b][1], "MC signal");
      leg[b]->AddEntry(h_BcM_SigMCv2[b][1], "MC signal (v2)");
      leg[b]->DrawClone("same");

      can[b]->SaveAs("figs/DataMinBkg_vsMC_wAccEffCorrections_eventByEvent_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+"_bin"+(TString)to_string(b)+metafitSyst+systExt+(TString)(limitAccEffCorrErrors?"_limitAccEffCorrErrors":"")+(TString)(truncateLargeCorr?"_truncateLargeCorr":"")+".pdf");
      can[b]->SaveAs("figs/DataMinBkg_vsMC_wAccEffCorrections_eventByEvent_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+"_bin"+(TString)to_string(b)+metafitSyst+systExt+(TString)(limitAccEffCorrErrors?"_limitAccEffCorrErrors":"")+(TString)(truncateLargeCorr?"_truncateLargeCorr":"")+".png");
    }
  }

  TFile * outf = new TFile("corrected_yields"+(TString)(secondStep?"_2ndStep":"")+".root","UPDATE");
  outf->WriteObject(&Yields_pref,"Yields_prefit"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);
  outf->WriteObject(&Yields_postf,"Yields_postfit"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYield,"corrYield"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYield_BDT23,"corrYield_BDTeff23"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYield_BDT3,"corrYield_BDTeff3"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYieldErr,"corrYieldErr"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYieldErr_BDT23,"corrYieldErr_BDTeff23"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYieldErr_BDT3,"corrYieldErr_BDTeff3"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYield_MC,"corrYield_MC"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&corrYield_MCv2,"corrYield_MCv2"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&nsig,"nsig"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&rsig,"rsig"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&rsig_relerr,"rsig_relerr"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);  
  outf->WriteObject(&r1r2Corr,"r1r2Correlation"+(TString)(ispp?"_pp":"_PbPb")+metafitSyst+systExt);

  if(!ispp){
    outf->WriteObject(&Yields_postf_cent,"Yields_postfit_centralityDep_PbPb"+metafitSyst+systExt);  
    outf->WriteObject(&nsig_cent,"nsig_centralityDep_PbPb"+metafitSyst+systExt);  
    outf->WriteObject(&rsig_cent,"rsig_centralityDep_PbPb"+metafitSyst+systExt);  
    outf->WriteObject(&rsig_relerr_cent,"rsig_relerr_centralityDep_PbPb"+metafitSyst+systExt);  
    outf->WriteObject(&r1r2Corr_cent,"r1r2Correlation_centralityDep_PbPb"+metafitSyst+systExt);
  }
  outf->Close();
}
