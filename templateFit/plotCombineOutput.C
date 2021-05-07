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
#include "TCanvas.h"
#include "TPad.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"

void plotCombineOutput(bool ispp=true, bool secondStep=false, bool bkgOnly=false, TString metafitSyst="", TString systExt=""){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  int ntrees = 9;

  //************************ meta-parameters (mostly outdated)
  bool bToJpsiOnly = !ispp;
  int nParam = 3; //effective number of parameters from the template fit, !!!hard-coded and approximate
  int varyBDTbin = 0;
  if(metafitSyst=="_BDTbinsUp") varyBDTbin = 1; 
  else if(metafitSyst=="_BDTbinsDown") varyBDTbin = -1;

  //************************ For significance vs BDT plot
  float BDTcut_l[_NanaBins+1][_NcentBins+1][_nChan(ispp)+2]; 
  for(int b=0;b<=_NanaBins;b++){
    for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){ //centrality bins only for PbPb and kinBin=0
      vector<float> BDTcut = _BDTcuts(ispp,b,cb,secondStep,(metafitSyst=="_BDTuncorrFromM"),varyBDTbin); //vector of BDT cut values put into array
      //  if(ignore1stBin) BDTcut.erase(0);
      BDTcut_l[b][cb][0] = -1;
      for(int k=0;k<=_nChan(ispp);k++)
	BDTcut_l[b][cb][k+1] = BDTcut[k];
    }
  }

  //************************ Mass histos definitions
  vector<vector<vector<vector<TH1F*> > > > h_BcM(ntrees, vector<vector<vector<TH1F*> > >(_NanaBins+1, vector<vector<TH1F*> >(_NcentBins+1, vector<TH1F*>(_nChan(ispp)+1))));
  vector<vector<vector<TH1F*> > > h_BcM_exp(_NanaBins+1, vector<vector<TH1F* > >(_NcentBins+1, vector<TH1F*>(_nChan(ispp)+1))); //expected
  vector<vector<vector<TGraph*> > > h_BcM_pull(_NanaBins+1, vector<vector<TGraph*> >(_NcentBins+1, vector<TGraph*>(_nChan(ispp)+1)));
  vector<vector<vector<vector<TH1F*> > > > h_BcM_prefit(ntrees, vector<vector<vector<TH1F*> > >(_NanaBins+1, vector<vector<TH1F*> >(_NcentBins+1, vector<TH1F*>(_nChan(ispp)+1))));
  vector<vector<vector<vector<TH1F*> > > > h_BcM_postfit(ntrees, vector<vector<vector<TH1F*> > >(_NanaBins+1, vector<vector<TH1F*> >(_NcentBins+1, vector<TH1F*>(_nChan(ispp)+1))));
  vector<vector<vector<vector<TH1F*> > > > h_BcM_DrawPostfit(ntrees, vector<vector<vector<TH1F*> > >(_NanaBins+1, vector<vector<TH1F*> >(_NcentBins+1, vector<TH1F*>(_nChan(ispp)+1))));
  for(int i=0; i<ntrees; i++){
    for(int b=0;b<=_NanaBins;b++){
      for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){ //centrality bins only for PbPb and kinBin=0
	for(int k=0;k<=_nChan(ispp);k++){
	  h_BcM_prefit[i][b][cb].push_back(new TH1F( "Bc_M_prefit_"+(TString)to_string(i)+"_"+(TString)to_string(b)+"_"+(TString)to_string(cb)+"_"+(TString)to_string(k), "Bc mass", _nbinM(ispp)[(k>0)?(k-1):0], _Mbinning(ispp, (k>0)?(k-1):0) ));
	  h_BcM_postfit[i][b][cb].push_back(new TH1F( "Bc_M_postfit_"+(TString)to_string(i)+"_"+(TString)to_string(b)+"_"+(TString)to_string(cb)+"_"+(TString)to_string(k), "Bc mass", _nbinM(ispp)[(k>0)?(k-1):0], _Mbinning(ispp, (k>0)?(k-1):0) ));
	  h_BcM_DrawPostfit[i][b][cb].push_back(new TH1F( "Bc_M_DrawPostfit_"+(TString)to_string(i)+"_"+(TString)to_string(b)+"_"+(TString)to_string(cb)+"_"+(TString)to_string(k), "Bc mass", _nbinM(ispp)[(k>0)?(k-1):0], _Mbinning(ispp, (k>0)?(k-1):0) ));
	}
      }
    }
  }

  //scale down control-region bins because they are wider
  vector<float> CRbinwRatio;
  for(int k=0;k<=2;k++) CRbinwRatio.push_back( ((_mBcMax-_mBcMin)/_nbinMSR(ispp)[k]) * (_nbinMCR(ispp)[k]/(_mMax-_mBcMax)) );

  //************************ Tree names, and which trees do we draw in the plots
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi"};
  TString prettyName[] = {"wrong-sign","fake J/#psi","high-mass control","signal region","MC signal",
  			  "MC B#rightarrow J/#psi X","MC PromptJpsi","dimuon+track (misID)","pivoted J/#psi"};
  TString procNameDef[] = {"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig",(ispp?"JpsiMC":"NPJpsi"),"PromptJpsi","dimuTrk","flipJpsi"};

  bool drawShape[] = {true,true,false,true,true,true,//nonprompt
		      false,//prompt//!ispp
		      false,//dimutrk
		      true//flipJpsi
  }; //mostly for drawing WrongSign on top
  bool usedForFit[] = {false,true,false,true,true,true,//nonprompt
		       false,//prompt//!ispp
		       false,//dimutrk
		       true//flipJpsi
  };
  int nPiled = -1; 
  for(int j=0;j<ntrees;j++) nPiled+= usedForFit[j];
  int idxPiled[] = {1,5,8,4}; //order of piled histos when plotting

  //************************ Efficiencies and significance
  vector<vector<TH1F*> > effSig(_NanaBins+1),rejBkg(_NanaBins+1),sigSignif(_NanaBins+1),sigSignifAsim(_NanaBins+1),sigPurity(_NanaBins+1);
  for(int b=0;b<=_NanaBins;b++){
    for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){ //centrality bins only for PbPb and kinBin=0
      effSig[b].push_back(new TH1F( "effSig_KinBin"+(TString)to_string(b)+(TString)to_string(cb), "signal efficiency", _nChan(ispp)+1, BDTcut_l[b][cb] ));
      rejBkg[b].push_back(new TH1F( "rejBkg_KinBin"+(TString)to_string(b)+(TString)to_string(cb), "background rejection", _nChan(ispp)+1, BDTcut_l[b][cb] ));
      sigSignif[b].push_back(new TH1F( "sigSignif_KinBin"+(TString)to_string(b)+(TString)to_string(cb), "signal significance S/sqrt(S+B)", _nChan(ispp)+1, BDTcut_l[b][cb] ));
      sigPurity[b].push_back(new TH1F( "sigPurity_KinBin"+(TString)to_string(b)+(TString)to_string(cb), "signal purity S/(S+B)", _nChan(ispp)+1, BDTcut_l[b][cb] ));
      sigSignifAsim[b].push_back(new TH1F( "sigSignifAsim_KinBin"+(TString)to_string(b)+(TString)to_string(cb), "signal significance (Asimov improved) S/sqrt(S+B)", _nChan(ispp)+1, BDTcut_l[b][cb] ));
    }
  }

  //********************************************************
  //Extract BcM HISTOGRAMS used as combine input + Record PREFIT NORMALISATIONS
  //********************************************************
  vector<vector<vector<vector<float> > > > yields_prefit(ntrees, vector<vector<vector<float> > >(_NanaBins+1, vector<vector<float> >(_NcentBins+1,vector<float>(_nChan(ispp)+1))));

  auto histFile = TFile::Open("InputForCombine_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+metafitSyst+".root");
  for(int i=0; i<ntrees; i++){
    if(!drawShape[i]) continue;
    if(i==7) continue; //no dimutrk for now
    for(int b=0;b<=_NanaBins;b++){
      for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){ //centrality bins only for PbPb and kinBin=0
	for(int k=0;k<=_nChan(ispp);k++){
	  int kk = (k==0)?1:k;
	  TString s_bin = (cb==0)?((b==0)?"":("Kin"+(TString)to_string(b))):("Cent"+(TString)to_string(cb));
	  h_BcM[i][b][cb][k] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+s_bin+"/"+procNameDef[i]+"/BcM"+(TString)(_preFitCorrAE?"_AccEffWeighted":""));

	  //**************** Record PREFIT NORMALISATIONS
	  yields_prefit[i][b][cb][k] = h_BcM[i][b][cb][k]->Integral(); //without compensating for larger control-region bin width

	  //compensate for larger bins in the high mass CR
	  for(int crbin=_nbinMSR(ispp)[kk-1]+1;crbin<=_nbinMSR(ispp)[kk-1]+_nbinMCR(ispp)[kk-1];crbin++){
	    h_BcM[i][b][cb][k]->SetBinContent(crbin, h_BcM[i][b][cb][k]->GetBinContent(crbin) * CRbinwRatio[kk-1]); 
	    h_BcM[i][b][cb][k]->SetBinError(crbin, h_BcM[i][b][cb][k]->GetBinError(crbin) * CRbinwRatio[kk-1]); 
	  }

	  if(k>1) { //manually sum histos in the BDT-integrated histo
	    AddTH1(h_BcM[i][b][cb][0], h_BcM[i][b][cb][k]);
	    yields_prefit[i][b][cb][0] += yields_prefit[i][b][cb][k];
	  }

	}
      }
    }
  }

  //********************************************************
  //DRAWING Bc mass for all trees and all BDT bins, PREFIT
  //********************************************************
  gStyle->SetOptStat(0);

  Color_t cols[] = {kCyan, kMagenta+1, kGreen, kRed-2, kBlue, kOrange-7, kOrange, kGreen+4, kGreen+1};

  cout<<"\n----------------------\n Drawing PREFIT mass plots"<<endl;
  for(int b=0;b<=_NanaBins;b++){
    for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
      cout<<endl<<"\n **** "<<(TString)((b==0)?"Kin-integrated":("Kinematic Bin #"+(TString)(to_string(b))))<<" "<<(TString)((cb==0)?"":(", Centrality Bin #"+(TString)(to_string(cb))))<<endl;
      TCanvas* c3 = new TCanvas("c3", "Bc mass, in BDT bins, prefit", _nChan(ispp)*1500, 1500);
      c3->Divide(_nChan(ispp),1);
      TCanvas* c4 = new TCanvas("c4", "Bc mass, no BDT cuts, prefit", 1500, 1500);

      for(int k=0;k<=_nChan(ispp);k++){
	if(k==0) {
	  c4->cd();
	  c4->SetRightMargin(0.04);}
	else {
	  c3->cd(k);
	  c3->cd(k)->SetRightMargin(0.04);}

	cout<<"\n **** For BDTval in ["<<_BDTcut_s(ispp,b,cb,secondStep,(metafitSyst=="_BDTuncorrFromM"),varyBDTbin)[max(k-1,0)]<<", "<<_BDTcut_s(ispp,b,cb,secondStep,(metafitSyst=="_BDTuncorrFromM"),varyBDTbin)[(k==0)?_nChan(ispp):k]<<"]"<<endl;
    
	//Formating all histos
	for(int i=0; i<ntrees; i++){
	  if(i==7) continue;
	  if(!drawShape[i]) continue;
	  h_BcM[i][b][cb][k]->SetFillStyle(usedForFit[i]?1001:0);
	  h_BcM[i][b][cb][k]->SetFillColor(cols[i]);
	  h_BcM[i][b][cb][k]->SetLineColor(((i==4)?kBlack:(cols[i]+2)));
	  h_BcM[i][b][cb][k]->SetLineWidth(3);
	}

	//Building the PILED (summed) histograms
	for(int i=0; i<nPiled; i++){
	  int iold = idxPiled[i];
	  h_BcM_prefit[i][b][cb][k] = (TH1F*)h_BcM[iold][b][cb][k]->Clone();
      
	  if(i>0){
	    h_BcM_prefit[i][b][cb][k]->Add(h_BcM_prefit[i-1][b][cb][k]);
	  }

	  h_BcM_prefit[i][b][cb][k]->SetFillStyle(3002);
	  h_BcM_prefit[i][b][cb][k]->SetFillColor(cols[iold]);
	  h_BcM_prefit[i][b][cb][k]->SetLineColor(((iold==4)?kBlack:(cols[iold]+2)));
	}

	//DRAW data signal region
	h_BcM[3][b][cb][k]->SetTitle("B_{c} candidates mass, "+(TString)((k==0)?"no BDT cuts":( (TString)((k==1)?"low":((k==2)?"medium":"high")) + " BDT"//#in "+(TString)((k==1)?"]-#infty":("["+_BDTcut_s(ispp,b,cb,secondStep,(metafitSyst=="_BDTuncorrFromM"),varyBDTbin)[k-1]))+","+(TString)((k==_nChan(ispp))?"+#infty[":(_BDTcut_s(ispp,b,cb,secondStep,(metafitSyst=="_BDTuncorrFromM"),varyBDTbin)[k]+"]"))
												) ));
	h_BcM[3][b][cb][k]->GetXaxis()->SetTitle("trimuon mass [GeV]");
	h_BcM[3][b][cb][k]->GetYaxis()->SetTitle("candidates/("+(TString)Form("%.2f",(_mBcMax-_mBcMin)/_nbinMSR(ispp)[(k>0)?(k-1):0]) +" GeV)");
	h_BcM[3][b][cb][k]->GetYaxis()->SetRangeUser(1e-4, 1.1*((h_BcM[3][b][cb][k]->GetMaximum() > h_BcM_prefit[nPiled-1][b][cb][k]->GetMaximum())?h_BcM[3][b][cb][k]->GetMaximum():h_BcM_prefit[nPiled-1][b][cb][k]->GetMaximum())  );
	h_BcM[3][b][cb][k]->GetYaxis()->SetTitleOffset(1.2);
	h_BcM[3][b][cb][k]->Draw("E");
      
	//DRAW PILED histos
	auto legend = new TLegend(0.63,0.59,0.96,0.9);
	for(int i=nPiled-1; i>=0; i--){
	  int iold = idxPiled[i];
	  //cout<<"--- Drawing piled backgrounds and expected signal: "<<prettyName[iold]<<endl;
	  h_BcM_prefit[i][b][cb][k]->Draw("histsame");      
	  legend->AddEntry(h_BcM_prefit[i][b][cb][k],  prettyName[iold], "lf");
	}

	//DRAW data signal region again
	//cout<<"--- Drawing again the data: "<<prettyName[3]<<endl;
	h_BcM[3][b][cb][k]->Draw("Esame");
	legend->AddEntry(h_BcM[3][b][cb][k], prettyName[3] ,"l");
    
	//DRAW remaining NON-PILED histos
	for(int i=0; i<ntrees; i++){
	  if(drawShape[i] && !usedForFit[i]){
	    //cout<<"--- Drawing other (non-piled) background: "<<prettyName[i]<<endl;
	    h_BcM[i][b][cb][k]->SetLineWidth(4);      
	    h_BcM[i][b][cb][k]->SetLineStyle(2);      
	    h_BcM[i][b][cb][k]->Draw("histsame");      
	    legend->AddEntry(h_BcM[i][b][cb][k],  prettyName[i], "l");
	  }
	}

	legend->SetTextSize(0.039);
	legend->Draw("same");
      
	//DRAW CMS Preliminary
	TLatex CMStag;
	CMStag.SetNDC();
	CMStag.SetTextFont(42);
	CMStag.SetTextSize(0.042);
	CMStag.DrawLatex(0.13,0.85,"#font[61]{CMS "+(TString)(ispp?Form("pp} (%.0f pb^{-1}, 5.02 TeV)",L_pp):Form("PbPb} (%.2f nb^{-1}, 5.02 TeV)",L_PbPb*1e3)));
	CMStag.DrawLatex(0.13,0.80,"#font[52]{Preliminary}");

	//Save canvas
	//((k==0)?c4:(c3->cd(k)))->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+(TString)((k==0)?"_noBDTcut":"_BDTbin"+(TString)(to_string(k)))+"_PREFIT"+metafitSyst+systExt+".pdf");
	//((k==0)?c4:(c3->cd(k)))->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+(TString)((k==0)?"_noBDTcut":"_BDTbin"+(TString)(to_string(k)))+"_PREFIT"+metafitSyst+systExt+".png");
	if(k==_nChan(ispp)) {
	  c3->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_allBDTbins_PREFIT"+metafitSyst+systExt+".png");
	  c3->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_allBDTbins_PREFIT"+metafitSyst+systExt+".pdf");}

      }
    }
  }

  //********************************************************
  //Extract results (POSTFIT NORMALISATIONS & SHAPES) from combine output
  //********************************************************

  //************************ Get files and yields
  TString normFileName = "./CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_2bins"+(TString)(secondStep?"_2ndStep":"")+metafitSyst+systExt+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
  TString normIntFileName = "./CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_integrated"+(TString)(secondStep?"_2ndStep":"")+metafitSyst+systExt+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
  auto normFile = TFile::Open(normFileName);
  auto normIntFile = TFile::Open(normIntFileName); TFile* normCentFile;
  cout<<"Extract normalisations from files: \n"<<normFileName<<"\n"<<normIntFileName<<endl;
  if(!ispp) {
    TString normCentFileName = "./CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(ispp?"pp":"PbPb")+"_centBins"+(TString)(secondStep?"_2ndStep":"")+metafitSyst+systExt+(TString)(_preFitCorrAE?"_wAccEff":"")+".root";
    normCentFile = TFile::Open(normCentFileName);
    cout<<normCentFileName<<endl;}
  RooArgSet *Yields = (RooArgSet*)normFile->Get("norm_fit_"+(TString)(bkgOnly?"b":"s"));
  RooArgSet *YieldsInt = (RooArgSet*)normIntFile->Get("norm_fit_"+(TString)(bkgOnly?"b":"s"));
  RooArgSet *YieldsCent = ispp?NULL:((RooArgSet*)normCentFile->Get("norm_fit_"+(TString)(bkgOnly?"b":"s")));

  vector<vector<vector<vector<float> > > > yields(ntrees, vector<vector<vector<float> > >(_NanaBins+1, vector<vector<float> >(_NcentBins+1,vector<float>(_nChan(ispp)+1))));
  for(int i=0; i<ntrees; i++){
    if(i==7) continue;
    if(!drawShape[i]) continue;
    for(int b=0;b<=_NanaBins;b++){
      for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
	for(int k=0;k<=_nChan(ispp);k++){

	  //Clone trimuon mass histos prefit into postfit. Processes not used in fitting will not change histos postfit
	  h_BcM_postfit[i][b][cb][k] = (TH1F*)h_BcM[i][b][cb][k]->Clone(); //Do this to recover proper X axis
	  h_BcM_postfit[i][b][cb][k]->SetDirectory(0);

	  if(i==3 || !usedForFit[i])
	    {yields[i][b][cb][k] = yields_prefit[i][b][cb][k];}
	  else{
	  
	    TString s_bin = (cb==0)?((b==0)?"":("Kin"+(TString)to_string(b))):("Cent"+(TString)to_string(cb));
	    if(k==0) yields[i][b][cb][k] = 0;
	    else{ 
	      //fetch postfit yields
	      yields[i][b][cb][k] = ((b==0 && cb==0)?YieldsInt:((cb==0)?Yields:YieldsCent))->getRealValue("BDT"+(TString)(to_string(k))+s_bin+"/"+procNameDef[i]);
	      yields[i][b][cb][0] += yields[i][b][cb][k]; 

	      //fetch postfit shapes
	      TH1F* h_tmp = (k==0)?NULL:( (TH1F*)((b==0 && cb==0)?normIntFile:((cb==0)?normFile:normCentFile))->Get("shapes_fit_"+(TString)(bkgOnly?"b":"s")+"/BDT"+(TString)(to_string(k))+s_bin+"/"+procNameDef[i]) );
	      //compensate for larger bins in the high mass CR
	      for(int crbin=_nbinMSR(ispp)[k-1]+1;crbin<=_nbinMSR(ispp)[k-1]+_nbinMCR(ispp)[k-1];crbin++){
		h_tmp->SetBinContent(crbin, h_tmp->GetBinContent(crbin) * CRbinwRatio[k-1]); 
		h_tmp->SetBinError(crbin, h_tmp->GetBinError(crbin) * CRbinwRatio[k-1]); 
	      }
	    
	      //cout<<"got shapes_fit_"+(TString)(bkgOnly?"b":"s")+"/BDT"+(TString)(to_string(k))+"Kin"+(TString)to_string(b)+"/"+procNameDef[i]<<endl;
	      for(int bin=0;bin<=h_BcM_postfit[i][b][cb][k]->GetNbinsX();bin++){
		h_BcM_postfit[i][b][cb][k]->SetBinContent(bin, (k==0)?0:(h_tmp->GetBinContent(bin)));
		h_BcM_postfit[i][b][cb][k]->SetBinError(bin, (k==0)?0:(h_tmp->GetBinError(bin)));
	      }
	      delete h_tmp;
	    }
	  }

	  //inclusive on BDT bins
	  if(k>0) AddTH1(h_BcM_postfit[i][b][cb][0], h_BcM_postfit[i][b][cb][k]);
	
	}//end loop on BDT bins
      } //end loop on centBin
    } //end loop on kinBin

    if(i!=3 && usedForFit[i]) {

      //extrapolate the yields for the bins not used in fit
      for(int b=0;b<=_NanaBins;b++){
	for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
	  if(systExt=="_noBDT1"){
	    
	    yields[i][b][cb][1] *= yields[i][b][cb][(systExt=="_noBDT12")?3:2]/yields_prefit[i][b][cb][(systExt=="_noBDT12")?3:2];
	    yields[i][b][cb][0] += yields[i][b][cb][1];}
	  if(systExt=="_noBDT12"){
	    yields[i][b][cb][2] *= yields[i][b][cb][3]/yields_prefit[i][b][cb][3];
	    yields[i][b][cb][0] += yields[i][b][cb][2];}	
	}
      }
    }

    for(int b=0;b<=_NanaBins;b++){
      for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
	for(int k=0;k<=_nChan(ispp);k++) 
	  cout<<"yield tree "<<i<<" kinBin "<<b<<" centBin "<<cb<<" BDTbin "<<k<<" = "<<yields[i][b][cb][k] <<"    ratio to prefit = "<<yields[i][b][cb][k]/yields_prefit[i][b][cb][k]<<endl; 
      }
    }
  }//end loop on trees

  normFile->Close();
  normIntFile->Close();
  if(!ispp) normCentFile->Close();
 
  //********************************************************
  //Measure SIGNAL SIGNIFICANCE and yields
  //********************************************************
  vector<vector<vector<float> > > nBkg = vector<vector<vector<float> > >(_NanaBins+1, vector<vector<float> >(_NcentBins+1, vector<float>(_nChan(ispp)+1))); //initialized to zeros
  vector<vector<vector<float> > > nSig = vector<vector<vector<float> > >(_NanaBins+1, vector<vector<float> >(_NcentBins+1, vector<float>(_nChan(ispp)+1)));
  vector<vector<vector<float> > > nData = vector<vector<vector<float> > >(_NanaBins+1, vector<vector<float> >(_NcentBins+1, vector<float>(_nChan(ispp)+1)));
  
  for(int b=0;b<=_NanaBins;b++){
    for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
      
      for(int k=0;k<=_nChan(ispp);k++){    

	//Fill sig and bkg
	nData[b][cb][k] = yields[3][b][cb][k];
	nSig[b][cb][k] = yields[4][b][cb][k];
	nBkg[b][cb][k] = 0;
	for(int i=0; i<ntrees; i++){
	  if(!usedForFit[i] || i==3 || i==4) continue; //only yields of histos used in combine fit
	  nBkg[b][cb][k] += yields[i][b][cb][k];
	}

	//Fill histos of efficiencies and significance
	effSig[b][cb]->SetBinContent(k+1, nSig[b][cb][k]/nSig[b][cb][0]);
	rejBkg[b][cb]->SetBinContent(k+1, 1- nBkg[b][cb][k]/nBkg[b][cb][0]);
	sigSignif[b][cb]->SetBinContent(k+1, (nSig[b][cb][k]+nBkg[b][cb][k]>0)?
				    (nSig[b][cb][k]/sqrt(nBkg[b][cb][k]+nSig[b][cb][k]) ):0 );
	sigPurity[b][cb]->SetBinContent(k+1, (nSig[b][cb][k]+nBkg[b][cb][k]>0)?
				    (nSig[b][cb][k]/(nBkg[b][cb][k]+nSig[b][cb][k]) ):0 );
	sigSignifAsim[b][cb]->SetBinContent(k+1, (nSig[b][cb][k]+nBkg[b][cb][k]>0 && ((2*nSig[b][cb][k]+nBkg[b][cb][k])*log(1+nSig[b][cb][k]/(nSig[b][cb][k]+nBkg[b][cb][k])) > nSig[b][cb][k]) )?
					(sqrt(2*( (2*nSig[b][cb][k]+nBkg[b][cb][k])*log(1+nSig[b][cb][k]/(nSig[b][cb][k]+nBkg[b][cb][k])) - nSig[b][cb][k]))   ):0);
      }
    }
  }
  
  //********************************************************
  //Make PULL graphs, and measure chi2
  //********************************************************
  float pullMax = 3.5;
  vector<vector<vector<int> > > ndf = vector<vector<vector<int> > >(_NanaBins+1, vector<vector<int> >(_NcentBins+1, vector<int>(_nChan(ispp)+1,0)));
  for(int b=0;b<=_NanaBins;b++){ 
    for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
        for(int k=0;k<=_nChan(ispp);k++) ndf[b][cb][k] = _nbinM(ispp)[(k>0)?(k-1):0]-nParam; }
  }
  vector<vector<vector<float> > > chi2 = vector<vector<vector<float> > >(_NanaBins+1, vector<vector<float> >(_NcentBins+1, vector<float>(_nChan(ispp)+1,0)));

  for(int b=0;b<=_NanaBins;b++){
    for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
      for(int k=0;k<=_nChan(ispp);k++){
	//adapt to mass binning!
	const int n = _nbinM(ispp)[(k==0)?0:(k-1)];

	double x[n];
	double y[n];

	bool seen = false;
	for(int i=0; i<ntrees; i++){
	  if(!usedForFit[i] || i==3) continue;
	  if(!seen) {h_BcM_exp[b][cb][k] = (TH1F*) h_BcM_postfit[i][b][cb][k]->Clone();
	    seen=true;}
	  else h_BcM_exp[b][cb][k]->Add(h_BcM_postfit[i][b][cb][k]); 
	}

	for(int bin=1; bin<=n; bin++){
	  x[bin-1] = h_BcM_exp[b][cb][k]->GetXaxis()->GetBinCenter(bin);

	  float binc = h_BcM[3][b][cb][k]->GetBinContent(bin);
	  float bine = h_BcM[3][b][cb][k]->GetBinError(bin);
	  float binc2 = h_BcM_exp[b][cb][k]->GetBinContent(bin);
	  y[bin-1] = (bine<=0)?0:( (binc-binc2)/bine );
	  if(y[bin-1]>=pullMax) y[bin-1] = pullMax*0.99999; //for plotting
	  if(y[bin-1]<=-pullMax) y[bin-1] = -pullMax*0.99999; //for plotting
	  chi2[b][cb][k] += pow(y[bin-1],2);
	  if(binc==0) ndf[b][cb][k] -= 1;
	}
	//      cout<<"chi2/ndf data vs expectation = "<<chi2[b][cb][k]/ndf[b][cb][k]<<endl;//" "<<chi2[b][cb][k]<<" "<<ndf[b][cb][k]<<endl;

	h_BcM_pull[b][cb][k] = new TGraphErrors(n,x,y);
	h_BcM_pull[b][cb][k]->SetMarkerColor(kBlack);    
	h_BcM_pull[b][cb][k]->SetFillColor(kBlack);    
	h_BcM_pull[b][cb][k]->SetMarkerSize(1.);    
      }
    }
  }

  //dummy histo to draw axis in the pull graph's pad
  float padratio = 5.5, padcorr=0.9;
  vector<TH1F*> dummy;
  for(int k=0;k<=_nChan(ispp);k++){
    dummy.push_back( new TH1F("dummy"+(TString)to_string(k), "", _nbinM(ispp)[(k>0)?(k-1):0], _Mbinning(ispp,(k>0)?(k-1):0)) );
    dummy[k]->GetXaxis()->SetTickLength(padcorr*padratio * h_BcM_exp[1][0][1]->GetXaxis()->GetTickLength());
    dummy[k]->GetXaxis()->SetLabelSize(padcorr*padratio * h_BcM_exp[1][0][1]->GetLabelSize("X"));
    dummy[k]->GetYaxis()->SetLabelSize(padcorr*padratio * h_BcM_exp[1][0][1]->GetLabelSize("Y"));
    dummy[k]->GetYaxis()->SetRangeUser(-pullMax,pullMax);  
    dummy[k]->GetYaxis()->SetNdivisions(205);//7
    //  dummy[k]->GetXaxis()->SetRangeUser(h_BcM_exp[1][1]->GetBinLowEdge(1), h_BcM_exp[1][1]->GetBinLowEdge(n+1));
    dummy[k]->GetXaxis()->SetTitle("trimuon mass [GeV]");    
    dummy[k]->GetYaxis()->SetTitle("pull");    
    dummy[k]->GetYaxis()->SetTitleSize(padcorr*padratio * h_BcM_exp[1][0][1]->GetXaxis()->GetTitleSize()); 
    dummy[k]->GetXaxis()->SetTitleSize(padcorr*padratio * h_BcM_exp[1][0][1]->GetXaxis()->GetTitleSize()); 
    dummy[k]->GetYaxis()->SetTitleOffset(0.19); 
    dummy[k]->GetXaxis()->SetLabelOffset(0.03); 
  }

  //********************************************************
  //DRAWING Bc mass for all trees and all BDT bins, POSTFIT
  //********************************************************
  cout<<"\n----------------------\n Drawing POSTFIT mass plots"<<endl;
  for(int b=0;b<=_NanaBins;b++){
    for(int cb=0;cb<=((ispp || b>0)?0:_NcentBins);cb++){
      
      cout<<endl<<"\n **** For Kinematic Bin "<<(TString)((b==0)?"integrated":("#"+(TString)to_string(b)))<< (TString)((cb>0)?(" centrality bin#"+(TString)to_string(cb)):"")<<endl;
      TCanvas* c6 = new TCanvas("c6", "Bc mass, in BDT bins, postfit", _nChan(ispp)*1500, 1500);
      c6->Divide(_nChan(ispp),2);
      TCanvas* c7 = new TCanvas("c7", "Bc mass, no BDT cuts, postfit", 1500, 1500);
      c7->Divide(1,2);

      for(int k=0;k<=_nChan(ispp);k++){

	if(k==0) {  
	  c7->cd(2)->SetPad(0.,0.,1.,1/padratio);//lower pad for pull
	  gPad->SetTickx(2);
	  gPad->SetGridy();
	  gPad->SetMargin(0.1,0.04, 0.38,0.);//left,right,bottom,top
	  dummy[k]->DrawClone();
	  c7->cd(1)->SetPad(0.,1/padratio,1.,1.); //upper pad for trimuon mass distro
	  gPad->SetMargin(0.1,0.04, 0.,0.1);
	  c7->cd(1);
	}
	else {
	  c6->cd(_nChan(ispp)+k)->SetPad((k-1)/(float)_nChan(ispp),0.,k/(float)_nChan(ispp),1/(padratio+1));//lower pad for pull
	  gPad->SetTickx(2);
	  gPad->SetGridy();
	  gPad->SetMargin(0.1,0.04, 0.38,0.);//left,right,bottom,top
	  dummy[k]->DrawClone();
	  c6->cd(k)->SetPad((k-1)/(float)_nChan(ispp),1/(padratio+1),k/(float)_nChan(ispp),1.);//upper pad for trimuon mass distro
	  gPad->SetMargin(0.1,0.04, 0.,0.1);
	  c6->cd(k);
	}

	cout<<"\n **** For BDTval in ["<<_BDTcut_s(ispp,b,cb,secondStep,(metafitSyst=="_BDTuncorrFromM"),varyBDTbin)[max(k-1,0)]<<", "<<_BDTcut_s(ispp,b,cb,secondStep,(metafitSyst=="_BDTuncorrFromM"),varyBDTbin)[(k==0)?_nChan(ispp):k]<<"]"<<endl;
	cout<<"background rejection = "<<rejBkg[b][cb]->GetBinContent(k+1)<<endl;
	if(!bkgOnly) cout<<"MC signal efficiency = "<<effSig[b][cb]->GetBinContent(k+1)<<endl;
	if(!bkgOnly) cout<<"Number of expected (MC postfit norm) signal events = "<<nSig[b][cb][k]<<endl;
	cout<<"Number of events in signal region = "<<nData[b][cb][k]<<endl;
	if(!bkgOnly) cout<<"Theoretical significance S/sqrt(S+B) = "<<sigSignif[b][cb]->GetBinContent(k+1)<<endl;
	if(!bkgOnly) cout<<"Theoretical significance (Asimov improved) S/sqrt(S+B) = "<<sigSignifAsim[b][cb]->GetBinContent(k+1)<<endl;

	auto legend = new TLegend(0.63,0.59,0.96,0.9);

	//DRAW data signal region
	h_BcM[3][b][cb][k]->Draw("E");

	//Building the PILED (summed) histograms
	for(int i=0; i<nPiled; i++){
	  int iold = idxPiled[i];
	  h_BcM_DrawPostfit[i][b][cb][k] = (TH1F*)h_BcM_postfit[iold][b][cb][k]->Clone();//"Bc_M_DrawPostfit_"+treeName[iold]);    
	
	  if(i>0){
	    h_BcM_DrawPostfit[i][b][cb][k]->Add(h_BcM_DrawPostfit[i-1][b][cb][k]);
	  }

	  h_BcM_DrawPostfit[i][b][cb][k]->SetFillStyle(3002);
	  h_BcM_DrawPostfit[i][b][cb][k]->SetFillColor(cols[iold]);
	  h_BcM_DrawPostfit[i][b][cb][k]->SetLineColor(((iold==4)?kBlack:(cols[iold]+2)));
	}

	//DRAW PILED histos
	for(int i=nPiled-1; i>=0; i--){
	  int iold = idxPiled[i];
	  //cout<<"--- Drawing piled backgrounds and expected signal: "<<prettyName[iold]<<endl;
	  h_BcM_DrawPostfit[i][b][cb][k]->Draw("histsame");      
	  legend->AddEntry(h_BcM_DrawPostfit[i][b][cb][k],  prettyName[iold], "lf");
	}

	//DRAW data signal region again
	h_BcM[3][b][cb][k]->GetYaxis()->SetRangeUser(1e-4, 1.1*((h_BcM[3][b][cb][k]->GetMaximum() > h_BcM_DrawPostfit[nPiled-1][b][cb][k]->GetMaximum())?h_BcM[3][b][cb][k]->GetMaximum():h_BcM_DrawPostfit[nPiled-1][b][cb][k]->GetMaximum())  );
	//cout<<"--- Drawing signal region: "<<prettyName[3]<<endl;
	h_BcM[3][b][cb][k]->Draw("Esame");
	legend->AddEntry(h_BcM[3][b][cb][k], prettyName[3] ,"l");

	//DRAW remaining NON-PILED histos
	for(int i=0; i<ntrees; i++){
	  if(drawShape[i] && !usedForFit[i]){
	    //cout<<"--- Drawing other (non-piled) background: "<<prettyName[i]<<endl;
	    h_BcM[i][b][cb][k]->SetLineWidth(4);      
	    h_BcM[i][b][cb][k]->SetLineStyle(2);      
	    h_BcM[i][b][cb][k]->Draw("histsame");      
	    legend->AddEntry(h_BcM[i][b][cb][k],  prettyName[i], "l");
	  }
	}

	legend->SetTextSize(0.039);
	legend->Draw("same");

	//DRAW values of efficiencies and significance
	char eff1[100]; sprintf(eff1,"f_{signal} = %.3f",effSig[b][cb]->GetBinContent(k+1));
	char eff2[100]; sprintf(eff2,"f_{background} = %.3f",1-rejBkg[b][cb]->GetBinContent(k+1));
	char eff3[100]; sprintf(eff3,"S/#sqrt{S+B} = %.1f",sigSignif[b][cb]->GetBinContent(k+1));
	char eff5[100]; sprintf(eff5,"purity = %.3f",sigPurity[b][cb]->GetBinContent(k+1));
	char eff4[100]; sprintf(eff4,"N_{signal}^{postfit} = %.0f",nSig[b][cb][k]);
	char s_chi2[100]; sprintf(s_chi2,"#chi^{2}/ndf = %.2f",chi2[b][cb][k]/ndf[b][cb][k]);
	TLatex eAndSignif;
	eAndSignif.SetNDC();
	eAndSignif.SetTextFont(52);
	eAndSignif.SetTextSize(0.038);
	//eAndSignif.DrawLatex(0.7,0.29,s_chi2);
	if(!bkgOnly) eAndSignif.DrawLatex(0.7,0.26,eff4);
	if(!bkgOnly) eAndSignif.DrawLatex(0.7,0.33,eff5);
	if(!bkgOnly) eAndSignif.DrawLatex(0.7,0.4,eff3);
	eAndSignif.DrawLatex(0.7,bkgOnly?0.33:0.47,eff2);
	if(!bkgOnly) eAndSignif.DrawLatex(0.7,0.54,eff1);

	//DRAW CMS Preliminary
	TLatex CMStag;
	CMStag.SetNDC();
	CMStag.SetTextFont(42);
	CMStag.SetTextSize(0.042);
	CMStag.DrawLatex(0.13,0.85,"#font[61]{CMS "+(TString)(ispp?Form("pp} (%.0f pb^{-1}, 5.02 TeV)",L_pp):Form("PbPb} (%.2f nb^{-1}, 5.02 TeV)",L_PbPb*1e3)));
	CMStag.DrawLatex(0.13,0.80,"#font[52]{Preliminary}");

	//DRAW PULL histos
	gStyle->SetBarWidth(0.8);
	if(k==0) c7->cd(2);
	else c6->cd(_nChan(ispp)+k);
	h_BcM_pull[b][cb][k]->Draw("Bsame");

	//SAVE canvases
	//For inclusive plot
	if(k==0){
	  c7->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_noBDTcut_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+metafitSyst+systExt+".pdf");
	  c7->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_noBDTcut_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+metafitSyst+systExt+".png");}
    
	// //For each BDT bin
	// else{    
	//   TCanvas* ctmp = new TCanvas("ctmp", "Bc mass in a BDT bin", 1500, 1500);
	//   TPad* ctmp1 = (TPad*)(c6->cd(k))->Clone();
	//   ctmp1->SetPad(0.,0.2,1.,1.);
	//   ctmp->cd(); 
	//   ctmp1->DrawClone();
	//   TPad* ctmp2 = (TPad*)(c6->cd(k+_nChan(ispp)))->Clone();
	//   ctmp2->SetPad(0.,0.,1.,0.2);  //necessary to recreate shape of canvas
	//   ctmp->cd(); 
	//   ctmp2->DrawClone();
      
	//   ctmp->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_BDTbin"+(TString)(to_string(k))+"_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+metafitSyst+systExt+".pdf");
	//   ctmp->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_BDTbin"+(TString)(to_string(k))+"_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+metafitSyst+systExt+".png");
	//   delete ctmp,ctmp1,ctmp2;
	// }
    
	//Canvas with all BDT bins
	if(k==_nChan(ispp)) {
	  c6->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_allBDTbins_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+metafitSyst+systExt+".pdf");
	  c6->SaveAs("figs/Bc_mass_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_allBDTbins_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+metafitSyst+systExt+".png");}

      } //end loop on BDT bins

      //********************************************************
      //DRAWING signal and background efficiencies, and the signal significances, versus BDT lower cut value
      //********************************************************
      if(bkgOnly) return;

      TCanvas* c5 = new TCanvas("c5", "Efficiencies and significance", 2000, 1500);
      effSig[b][cb]->SetLineColor(kRed);
      effSig[b][cb]->SetLineWidth(3);
      effSig[b][cb]->GetYaxis()->SetRangeUser(0,1);
      effSig[b][cb]->GetXaxis()->SetTitle("BDT");
      effSig[b][cb]->Draw("hist");
      rejBkg[b][cb]->SetLineWidth(3);
      rejBkg[b][cb]->Draw("histsame");
      c5->Update();

      //scale sigSignif to the pad coordinates
      Float_t rightmax = 1.1*sigSignif[b][cb]->GetMaximum();
      Float_t scale = gPad->GetUymax()/rightmax;
      sigSignif[b][cb]->SetLineWidth(3);
      sigSignif[b][cb]->SetLineColor(kGreen+3);
      sigSignif[b][cb]->Scale(scale);
      sigSignif[b][cb]->Draw("histsame");
      sigSignifAsim[b][cb]->SetLineWidth(3);
      sigSignifAsim[b][cb]->SetLineColor(kGreen);
      sigSignifAsim[b][cb]->Scale(scale);
      sigSignifAsim[b][cb]->Draw("histsame");

      //draw an axis on the right side
      TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
				gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
      axis->SetTitle("S/#sqrt{S+B}");
      axis->SetTitleColor(kGreen+2);
      axis->SetLineColor(kGreen+2);
      axis->SetLabelColor(kGreen+2);
      axis->SetLabelSize(0.03);
      axis->SetTextFont(52);
      axis->Draw();

      auto legend2 = new TLegend(0.2,0.65,0.48,0.9);
      legend2->SetBorderSize(0);
      legend2->SetFillStyle(0);
      legend2->SetTextSize(0.031);
      legend2->AddEntry(effSig[b][cb], "signal efficiency" );
      legend2->AddEntry(rejBkg[b][cb], "background rejection" );
      legend2->AddEntry(sigSignif[b][cb], "signal significance" );
      legend2->AddEntry(sigSignifAsim[b][cb], "signal significance (Asimov)" );
      legend2->Draw("same");

      if(!bkgOnly){
	c5->SaveAs("figs/efficiencySigBkgAndSignificance_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_vs_BDTbin"+metafitSyst+systExt+".png");
	c5->SaveAs("figs/efficiencySigBkgAndSignificance_"+(TString)(secondStep?"2ndStep_":"")+(TString)(ispp?"pp":"PbPb")+(TString)((b==0 && cb==0)?"_integrated":("_KinBin"+(TString)to_string(b)))+(TString)((cb==0)?"":("_CentBin"+(TString)to_string(cb)))+"_vs_BDTbin"+metafitSyst+systExt+".pdf");}

      if(b<_NanaBins) delete c5,c6,c7;
    } //end loop on cent bins
  } //end loop on pT bins
  
}

