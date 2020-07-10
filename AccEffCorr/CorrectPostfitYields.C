#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

void CorrectPostfitYields(bool ispp = true, bool BDTuncorrFromM=false, int skipBDTbins=0){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  bool useFlipJpsi = ispp;
  bool flipJpsiSameSide = false; // whether to keep only events with flipJpsi angle on same |eta| side
  bool bToJpsiOnly = false;//ispp;
  bool RescalePredToData = true; //rescale the sig+bkg in each trimuon mass bin to the actual data yield in this bin
  bool limitAccEffCorrErrors = false;

  const int nproc = 5;
  TString procNameDef[] = {"BcSig","data_obs","FakeJpsi",(ispp?"JpsiMC":"NPJpsi"),(ispp?"flipJpsi":"PromptJpsi")};
  vector<vector<TString> > systName(nproc);
  if(ispp)
    systName = {{""},{""}, {"","_JpsiSBUp","_JpsiSBDown"}, 
		{"","_flipJSameSideUp","_flipJSameSideDown","_wPromptMCUp","_wPromptMCDown"}   ,   {"","_flipJSameSideUp","_flipJSameSideDown","_wPromptMCUp","_wPromptMCDown"} };
  else 
    systName = {{""},{""}, {"","_JpsiSBUp","_JpsiSBDown"},
		{"","_PromptOrFlipJUp","_PromptOrFlipJDown","_bJpsiFracUp","_bJpsiFracDown"}   ,   {"","_PromptOrFlipJUp","_PromptOrFlipJDown","_bJpsiFracUp","_bJpsiFracDown"} }; //PbPb

  vector<pair<int,int> > allProc;
  for(int i=0; i<nproc; i++){
    for(int sys=0; sys<systName[i].size(); sys++){
      allProc.push_back(make_pair(i,sys));
    }
  }

  vector<vector<vector<TH1F*> > > h_AECorr(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_AECorr_BDT23(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_AECorr_BDT3(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_postfit(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_prefit(nproc, vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_AEweighted(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_AEweighted_BDT23(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<TH1F*> > > h_BcM_AEweighted_BDT3(allProc.size(), vector<vector<TH1F*> >(_NanaBins+1, vector<TH1F*>(_nChan(ispp)+1)));  
  vector<vector<vector<float> > > Yields_pref(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1))); //prefit
  vector<vector<vector<float> > > Yields_postf(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1))); //postfit, except for data
  vector<vector<vector<float> > > YieldsCorr(nproc, vector<vector<float> >(_NanaBins+1, vector<float>(_nChan(ispp)+1))); //postfit, except for data
  vector<vector<float> > corrYield(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_MC(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_MCv2(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_BDT23(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<vector<float> > corrYield_BDT3(_NanaBins+1, vector<float>(_nChan(ispp)+1, 0));
  vector<float> r1r2Corr(1,0);
  vector<float> rsig(_NanaBins+1,0), rsig_true(_NanaBins+1,0), rsig_errl(_NanaBins+1,0), rsig_errh(_NanaBins+1,0), rsig_relerr(_NanaBins+1,0);
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
  auto histFile = TFile::Open("../templateFit/InputForCombine_"+(TString)(BDTuncorrFromM?"BDTuncorrFromM_":"")+(TString)(ispp?"pp":"PbPb")+".root"); 
  cout<<"Extract (mass-dependent) acceptance x efficiency corrections from file "<<"../templateFit/InputForCombine_"+(TString)(BDTuncorrFromM?"BDTuncorrFromM_":"")+(TString)(ispp?"pp":"PbPb")+".root"<<endl;

  for(int p=0; p<allProc.size(); p++){
    int MainProc = allProc[p].first;
    int sys = allProc[p].second;
    for(int b=1;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){
	h_AECorr[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffCorr");
	h_AECorr_BDT23[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffCorr_BDTeff23");
	h_AECorr_BDT3[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffCorr_BDTeff3");
	h_BcM_AEweighted[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffWeighted"); //grab BcM, corrected by event-by-event acc-eff
	h_BcM_AEweighted_BDT23[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffWeighted_BDTeff23"); //grab BcM, corrected by event-by-event acc-eff
	h_BcM_AEweighted_BDT3[p][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM_AccEffWeighted_BDTeff3"); //grab BcM, corrected by event-by-event acc-eff

	if(sys==0) {
	  h_BcM_prefit[MainProc][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM"); //prefit histo
	  h_BcM_prefit[MainProc][b][k]->SetDirectory(0);
	}
	if(p==1) {
	  h_BcM_postfit[MainProc][b][k] = (TH1F*)histFile->Get("BDT"+(TString)to_string(k)+"Kin"+(TString)to_string(b)+"/"+procNameDef[MainProc]+systName[MainProc][sys]+"/BcM"); //grab data BcM histo
	  h_BcM_postfit[MainProc][b][k]->SetDirectory(0);
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
  TString normFileName = "../templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(bToJpsiOnly?"bToJpsi":"NonPromptJpsi")+(TString)(useFlipJpsi?(flipJpsiSameSide?"_flipJpsiSameSide":"_flipJpsi"):(ispp?"":"_PromptJpsi"))+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)((skipBDTbins==2)?"_ignoreBin1-2":((skipBDTbins==1)?"_ignoreBin1":""))+(TString)(ispp?"_pp":"_PbPb")+"_2bins"+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
  cout<<"Extract normalisations from file "<<normFileName<<endl;
  auto normFile = TFile::Open(normFileName);

  for(int proc=0; proc<nproc; proc++){
    if(proc==1) continue; //data is not a postfit shape
    for(int b=1;b<=_NanaBins;b++){
      for(int k=1;k<=_nChan(ispp);k++){

	h_BcM_postfit[proc][b][k] = (TH1F*)h_AECorr[1][b][k]->Clone("Bc_M_postfit_"+procNameDef[proc]); //Do this to recover proper X axis //not essential, but prettier	
	h_BcM_postfit[proc][b][k]->SetDirectory(0);

	//fetch postfit shapes                                                                           
	//cout<<"fetching histogram shapes_fit_s/BDT"+(TString)(to_string(k))+"Kin"+(TString)to_string(b)+"/"+procNameDef[proc]<<endl;
	TH1F* h_tmp = (TH1F*)normFile->Get("shapes_fit_s/BDT"+(TString)(to_string(k))+"Kin"+(TString)to_string(b)+"/"+procNameDef[proc]) ;
	for(int bin=0;bin<=h_BcM_postfit[proc][b][k]->GetNbinsX();bin++){
	  h_BcM_postfit[proc][b][k]->SetBinContent(bin, h_tmp->GetBinContent(bin));
	  h_BcM_postfit[proc][b][k]->SetBinError(bin, h_tmp->GetBinError(bin));
	}
	delete h_tmp;

      }
    }
  }

  //get values of nuisance parameters
  RooArgList fittedPars = ((RooFitResult*)normFile->Get("fit_s"))->floatParsFinal();
  float JpsiSB = ((RooRealVar*)fittedPars.find("JpsiSB"))->getValV(); //getError()
  float flipJ = ((RooRealVar*)fittedPars.find(ispp?"flipJSameSide":"PromptOrFlipJ"))->getValV(); //getError()
  float JMC = ((RooRealVar*)fittedPars.find(ispp?"wPromptMC":"bJpsiFrac"))->getValV(); //getError()
  r1r2Corr[0] = ((RooFitResult*)normFile->Get("fit_s"))->correlation("r1","r2");

  //Fill yields pre/post-fit (integrals of BcM distros)
  for(int proc=0; proc<nproc; proc++){
    for(int b=1;b<=_NanaBins;b++){
      Yields_pref[proc][b][0] = 0;
      Yields_postf[proc][b][0] = 0;
      for(int k=1;k<=_nChan(ispp);k++){
	Yields_pref[proc][b][k] = h_BcM_prefit[proc][b][k]->Integral(1,_nbinMSR(ispp)) ;
	Yields_postf[proc][b][k] = h_BcM_postfit[proc][b][k]->Integral(1,_nbinMSR(ispp)) ;
	Yields_pref[proc][b][0] += Yields_pref[proc][b][k];
	Yields_postf[proc][b][0] += Yields_postf[proc][b][k];
      }
    }
    for(int k=0;k<=_nChan(ispp);k++){
      Yields_pref[proc][0][k] = 0;
      Yields_postf[proc][0][k] = 0;
      for(int b=1;b<=_NanaBins;b++){
	Yields_pref[proc][0][k] += Yields_pref[proc][b][k];
	Yields_postf[proc][0][k] += Yields_postf[proc][b][k];
      }
    }
  }

  for(int b=1;b<=_NanaBins;b++){
    rsig[b] = ((RooRealVar*)fittedPars.find("r"+(TString)to_string(b)))->getValV();
    rsig_errl[b] = ((RooRealVar*)fittedPars.find("r"+(TString)to_string(b)))->getError();
    rsig_errh[b] = ((RooRealVar*)fittedPars.find("r"+(TString)to_string(b)))->getError();
    rsig_true[b] = Yields_postf[0][b][0]/Yields_pref[0][b][0];
    rsig_relerr[b] = rsig_errh[b]/rsig[b];
  }
  
  cout<<"Shape morphing parameters (JpsiSB, flipJSameSide/PromptOrFlipJ, wPromptMC/bJpsiFrac) are found to be "<<JpsiSB<<" "<<flipJ<<" "<<JMC<<endl;
  cout<<"Signal normalisation modifiers for the 2 pt bins are found to be "<<rsig[1]<<" "<<rsig[2]<<endl;
  cout<<"           versus actual ratio of postfit and prefit signal normalisations = "<<rsig_true[1]<<" "<<rsig_true[2]<<endl;
  normFile->Close();

  //********************************************************
  //Subtract background from data, considering AccEff-corrected yields
  //******************************************************** 
  for(int b=1;b<=_NanaBins;b++){
    cout<<endl<<"b = "<<b<<endl;
    corrYield[b][0] = 0;
    corrYield_BDT23[b][0] = 0;
    corrYield_BDT3[b][0] = 0;
    corrYield_MC[b][0] = 0;
    corrYield_MCv2[b][0] = 0;

    for(int k=1;k<=_nChan(ispp);k++){
      cout<<endl<<"k = "<<k<<endl;
      corrYield[b][k]=0;
      corrYield_BDT23[b][k]=0;
      corrYield_BDT3[b][k]=0;
      corrYield_MC[b][k]=0;
      corrYield_MCv2[b][k]=0;

      h_BcM_BkgSubData[b][k] = (TH1F*)h_BcM_AEweighted[1][b][k]->Clone("BackgroundSubtractedData_AEcorr_BDT"+(TString)(to_string(k))+"Kin"+(TString)to_string(b));
      h_BcM_BkgSubData_BDT23[b][k] = (TH1F*)h_BcM_AEweighted_BDT23[1][b][k]->Clone("BackgroundSubtractedData_AEcorr_BDTeff23_BDT"+(TString)(to_string(k))+"Kin"+(TString)to_string(b));
      h_BcM_BkgSubData_BDT3[b][k] = (TH1F*)h_BcM_AEweighted_BDT3[1][b][k]->Clone("BackgroundSubtractedData_AEcorr_BDTeff3_BDT"+(TString)(to_string(k))+"Kin"+(TString)to_string(b));

      h_BcM_Data[b][k] = new TH1F("corrData_bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), "Data-background, corrected by AccEff, bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), _nbinMSR(ispp), 3.5, 6.2);    
      h_BcM_Bkg[b][k] = new TH1F("corrBkg_bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), "Scaled background, corrected by AccEff, bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), _nbinMSR(ispp), 3.5, 6.2);    
      h_BcM_DataMinBkg[b][k] = new TH1F("corrDataMinBkg_bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), "Data-background, corrected by AccEff, bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), _nbinMSR(ispp), 3.5, 6.2);    
      h_BcM_SigMC[b][k] = new TH1F("corrSigMC_bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), "Scaled signal MC corrected by AccEff, bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), _nbinMSR(ispp), 3.5, 6.2);    
      h_BcM_SigMCv2[b][k] = new TH1F("corrSigMCv2_bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), "Signal MC (v2) corrected by AccEff, bin"+(TString)(to_string(k))+"Kin"+(TString)to_string(b), _nbinMSR(ispp), 3.5, 6.2);    
      
      for(int bin=1;bin<=h_BcM_BkgSubData[b][k]->GetNbinsX();bin++){ //forget underflow here?

	cout<<"bin #"<<bin<<endl;      
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

	//lambda function to grab the AE correction for a given process (for this precise bin)
	auto AEcorr = [h_AECorr,h_AECorr_BDT23,h_AECorr_BDT3,b,k,bin,limitAccEffCorrErrors](int process, int BDTeffMode)->float{
	  int pout = process;
	  vector<vector<vector<TH1F*> > > h_AECorr_tmp = (BDTeffMode==0)?h_AECorr:((BDTeffMode==1)?h_AECorr_BDT23:h_AECorr_BDT3);
	  if(limitAccEffCorrErrors){
	    if(pout>5 && pout<11
	       && fabs(h_AECorr_tmp[pout][b][k]->GetBinContent(bin) - h_AECorr_tmp[5][b][k]->GetBinContent(bin)) < h_AECorr_tmp[pout][b][k]->GetBinError(bin) //error larger than difference between nominal and systematic shape
	       && h_AECorr_tmp[pout][b][k]->GetBinError(bin) > h_AECorr_tmp[5][b][k]->GetBinError(bin) ) //larger error for systematic than nominal
	      pout = 5; //keep the nominal AE for this bin and this systematic
	    if(pout>10 && pout<16
	       && fabs(h_AECorr_tmp[pout][b][k]->GetBinContent(bin) - h_AECorr_tmp[10][b][k]->GetBinContent(bin)) < h_AECorr_tmp[pout][b][k]->GetBinError(bin) //error larger than difference between nominal and systematic shape
	       && h_AECorr_tmp[pout][b][k]->GetBinError(bin) > h_AECorr_tmp[10][b][k]->GetBinError(bin) ) //larger error for systematic than nominal
	      pout = 10;
	  }
	  return h_AECorr_tmp[pout][b][k]->GetBinContent(bin);
	};

	for(int m=0;m<3;m++){
	  //1 shape morphing parameter for Jpsi sidebands
	  float p_JSB = min((float)1,fabs(JpsiSB));
	  binBkg[m] += yield[2] * ( (1-p_JSB) * AEcorr(2,m) + p_JSB * AEcorr( 3+((JpsiSB>0)?0:1) ,m) );
	  //cout<<"bkg#2 corrected yield = "<<yield[2] * ( (1-p_JSB) * AEcorr(2,m) + p_JSB * AEcorr( 3+((JpsiSB>0)?0:1) ,m) )<<endl;
	  //2 shape morphing parameters for Jpsi MC and for flipJpsi
	  float p_flipJ = min((float)1,fabs(flipJ));
	  float p_JMC = min((float)1,fabs(JMC));
	  binBkg[m] += yield[3] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(5,m) + p_flipJ * (1-p_JMC/2) * AEcorr(6+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(8+((JMC>0)?0:1) ,m) );
	  //cout<<"bkg#3 corrected yield = "<<yield[3] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(5,m) + p_flipJ * (1-p_JMC/2) * AEcorr(6+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(8+((JMC>0)?0:1),m) )<<endl;
	  binBkg[m] += yield[4] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(10,m) + p_flipJ * (1-p_JMC/2) * AEcorr(11+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(13+((JMC>0)?0:1),m) );
	  //cout<<"bkg#4 corrected yield = "<<yield[4] * ( (1-p_flipJ) * (1-p_JMC) * AEcorr(10,m) + p_flipJ * (1-p_JMC/2) * AEcorr(11+((flipJ>0)?0:1),m) + p_JMC * (1-p_flipJ/2) * AEcorr(13+((JMC>0)?0:1),m) )<<endl;
	}
	
	cout<<"corrected data vs corrected bkg: "<<h_BcM_BkgSubData[b][k]->GetBinContent(bin)<<" "<<binBkg[0]<<endl;
	cout<<"data-bkg vs MC signal v1 and v2: "<<h_BcM_BkgSubData[b][k]->GetBinContent(bin) - binBkg[0]<<" "<<yield[0] * h_AECorr[0][b][k]->GetBinContent(bin)<<" "<<rsig[b]*h_BcM_AEweighted[0][b][k]->GetBinContent(bin)<<endl;
	h_BcM_BkgSubData[b][k]->SetBinContent(bin, h_BcM_BkgSubData[b][k]->GetBinContent(bin) - binBkg[0]); //substract corrected background from corrected data
	h_BcM_BkgSubData_BDT23[b][k]->SetBinContent(bin, h_BcM_BkgSubData_BDT23[b][k]->GetBinContent(bin) - binBkg[1]); //substract corrected background from corrected data
	h_BcM_BkgSubData_BDT3[b][k]->SetBinContent(bin, h_BcM_BkgSubData_BDT3[b][k]->GetBinContent(bin) - binBkg[2]); //substract corrected background from corrected data
	if(h_BcM_BkgSubData[b][k]->GetBinCenter(bin)<6.2){ 
	  corrYield[b][k] += h_BcM_BkgSubData[b][k]->GetBinContent(bin);
	  corrYield_BDT23[b][k] += h_BcM_BkgSubData_BDT23[b][k]->GetBinContent(bin);
	  corrYield_BDT3[b][k] += h_BcM_BkgSubData_BDT3[b][k]->GetBinContent(bin);
	  corrYield_MC[b][k] += yield[0] * h_AECorr[0][b][k]->GetBinContent(bin);
	  corrYield_MCv2[b][k] += rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinContent(bin);

	  //Fill histos for display
	  h_BcM_Data[b][k]->SetBinContent(bin , h_BcM_AEweighted[1][b][k]->GetBinContent(bin) );
	  h_BcM_Data[b][k]->SetBinError(bin , h_BcM_AEweighted[1][b][k]->GetBinError(bin) );
	  h_BcM_Bkg[b][k]->SetBinContent(bin , binBkg[0] );
	  h_BcM_DataMinBkg[b][k]->SetBinContent(bin , h_BcM_BkgSubData[b][k]->GetBinContent(bin) );
	  h_BcM_DataMinBkg[b][k]->SetBinError(bin , h_BcM_AEweighted[1][b][k]->GetBinError(bin) );
	  h_BcM_SigMC[b][k]->SetBinContent(bin , yield[0] * h_AECorr[0][b][k]->GetBinContent(bin) );
	  h_BcM_SigMC[b][k]->SetBinError(bin , rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinError(bin));
	  h_BcM_SigMCv2[b][k]->SetBinContent(bin , rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinContent(bin) );
	  h_BcM_SigMCv2[b][k]->SetBinError(bin , rsig_true[b] * h_BcM_AEweighted[0][b][k]->GetBinError(bin));
	}

      }//end loop on mass bins
      
      cout<<"corrected yield (in mass signal region) for analysis bin #"<<b<<" and BDT bin "<<k<<" is "<<corrYield[b][k]<<"\n    vs MC v1 and v2: "<<corrYield_MC[b][k]<<" "<<corrYield_MCv2[b][k]<<endl;

      corrYield[b][0] += corrYield[b][k];
      if(k>1) corrYield_BDT23[b][0] += corrYield_BDT23[b][k];
      if(k>2) corrYield_BDT3[b][0] += corrYield_BDT3[b][k];
      corrYield_MC[b][0] += corrYield_MC[b][k];
      corrYield_MCv2[b][0] += corrYield_MCv2[b][k];
    }//end loop on BDT bins
    cout<<"corrected yield (in mass signal region) for analysis bin #"<<b<<" is "<<corrYield[b][0]<<"\n    vs MC v1 and v2:"<<corrYield_MC[b][0]<<" "<<corrYield_MCv2[b][0]<<endl;
    cout<<"    vs BDTeff23 and BDTeff3: "<<corrYield_BDT23[b][0]<<" "<<corrYield_BDT3[b][0]<<endl;

  }//end loop on analysis bins

  //********************************************************
  //DRAWING
  //******************************************************** 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLine *vert = new TLine(6.2,0,6.2,1000);
  vert->SetLineWidth(3);
  vert->SetLineStyle(1);
  TLine *hor = new TLine(6.2,0,6.2,1000);
  hor->SetLineWidth(ispp?10:12);
  hor->SetLineStyle(1);
  TLine *vert0 = new TLine(6.2,0,6.2,1000);
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
      h_BcM_Bkg[b][k]->Draw("histsame");
      gPad->Update();
      if(k>1)
	vert->DrawLine((ispp?1.002:0.97)*gPad->GetUxmin(),gPad->GetUymin() - 0.15*(gPad->GetUymax()-gPad->GetUymin()),(ispp?1.002:0.97)*gPad->GetUxmin(),gPad->GetUymax());

      can[b]->cd(k+_nChan(ispp));
      gPad->SetTickx(2);
      if(k==_nChan(ispp)) gPad->SetRightMargin(0.001);
      if(k>1 && ispp) h_BcM_DataMinBkg[b][k]->GetYaxis()->SetTickLength(0);
      if(k>0 && !ispp) gPad->SetLeftMargin(0.085);
      if(k<_nChan(ispp)) gPad->SetRightMargin(0.001);
      h_BcM_DataMinBkg[b][k]->GetYaxis()->SetRangeUser(minC2 , maxC2);
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

    can[b]->SaveAs("DataMinBkg_vsMC_wAccEffCorrections_eventByEvent_"+(TString)(ispp?"pp":"PbPb")+"_bin"+(TString)to_string(b)+".pdf");
    can[b]->SaveAs("DataMinBkg_vsMC_wAccEffCorrections_eventByEvent_"+(TString)(ispp?"pp":"PbPb")+"_bin"+(TString)to_string(b)+".png");
  }

  TFile * outf = new TFile("corrected_yields.root","UPDATE");
  outf->WriteObject(&Yields_pref,"Yields_prefit"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&Yields_postf,"Yields_postfit"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&corrYield,"corrYield"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&corrYield_BDT23,"corrYield_BDTeff23"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&corrYield_BDT3,"corrYield_BDTeff3"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&corrYield_MC,"corrYield_MC"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&corrYield_MCv2,"corrYield_MCv2"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&rsig_relerr,"rsig_relerr"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
  outf->WriteObject(&r1r2Corr,"r1r2Correlation"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"));  
}
