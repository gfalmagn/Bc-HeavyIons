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
#include "../BDT/Definitions.h"

float quadSum(float f1, float f2, float f3=0, float f4=0){
  return sqrt(f1*f1+f2*f2+f3*f3+f4*f4);
}

void plotCombineOutput(bool ispp=true, bool bkgOnly=false, bool BDTuncorrFromM=false){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  int ntrees = 9;

  bool ignoreBin2 = false;
  bool ignoreBin1 = ignoreBin2 || false;
  bool separateCombinatorics = false; // in non-prompt Jpsi MC, separate the events where muW and Jpsi are not from same decay, and gather them with Prompt Jpsi MC to make a combinatorics sample
  bool useFlipJpsi = ispp;
  bool flipJpsiSameSide = false; // whether to keep only events with flipJpsi angle on same |eta| side
  bool bToJpsiOnly = false;//ispp;
  int nbin = 38;
  int nParam = 6; //number of parameters from the template fit, !!!hard-coded

  vector<float> BDTcut = (BDTuncorrFromM?_corrBDTcuts:_BDTcuts)(ispp); //vector of BDT cut values put into array
  //  if(ignore1stBin) BDTcut.erase(0);
  //int nchan = BDTcut.size() -1;
  float BDTcut_l[_nChan(ispp)+2]; 
  BDTcut_l[0] = -1;
  for(int k=0;k<=_nChan(ispp);k++){
    BDTcut_l[k+1] = BDTcut[k];
  }

  vector<vector<TH1F*> > h_BcM(ntrees, vector<TH1F*>(_nChan(ispp)+1));
  vector<TH1F*> h_BcM_exp(vector<TH1F*>(_nChan(ispp)+1)); //expected
  vector<TGraph*> h_BcM_pull(vector<TGraph*>(_nChan(ispp)+1));
  //vector<vector<TH1F*> > h_BDT(ntrees, vector<TH1F*>(_nChan(ispp)+1));
  // vector<vector<TH1F*> > h_BcPt(ntrees, vector<TH1F*>(_nChan(ispp),new TH1F( "Bc_Pt", "Bc transverse momentum", 50, 0, 30 )));
  vector<vector<TH1F*> > h_BcM_prefit(ntrees, vector<TH1F*>(_nChan(ispp)+1,new TH1F( "Bc_M_prefit", "Bc mass", _nbinM(ispp), _Mbinning(ispp) )));
  vector<vector<TH1F*> > h_BcM_postfit(ntrees, vector<TH1F*>(_nChan(ispp)+1,new TH1F( "Bc_M_postfit", "Bc mass", _nbinM(ispp), _Mbinning(ispp) )));
  vector<vector<TH1F*> > h_BcM_DrawPostfit(ntrees, vector<TH1F*>(_nChan(ispp)+1,new TH1F( "Bc_M_DrawPostfit", "Bc mass", _nbinM(ispp), _Mbinning(ispp) )));
  // vector<vector<TH1F*> > h_QQM(ntrees, vector<TH1F*>(_nChan(ispp),new TH1F( "QQ_M", "Jpsi mass", nbinM, 2,4 )));
  // vector<vector<TH1F*> > h_QQ2M(ntrees, vector<TH1F*>(_nChan(ispp),new TH1F( "QQ2_M", "Jpsi mass", nbinM, 2,4 )));
  // vector<vector<TH1F*> > h_QQ3M(ntrees, vector<TH1F*>(_nChan(ispp),new TH1F( "QQ3_M", "Jpsi mass", nbinM, 2,4 )));
  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi"};
  TString prettyName[] = {"WrongSign","J/#psi sidebands","High mass control","signal region","MC signal",
  			  "MC NonPromptJpsi","MC PromptJpsi","dimuon+track (misID)","flipped J/#psi"};
  if(separateCombinatorics) {
    prettyName[5] = "NonPromptJpsi (muW from same decay)";
    prettyName[6] = "J/#psi-#mu combinatorics";}

  TString procNameDef[] = {"Wrongsign","FakeJpsi","TrueJpsi","data_obs","BcSig",(ispp?"JpsiMC":"NPJpsi"),"PromptJpsi","dimuTrk","flipJpsi"};
  vector<vector<TString>> procName{{"Wrongsign"},{"FakeJpsi"},{"TrueJpsi"},{"data_obs"},{"BcSig"},
				  {"NonPromptJpsi","bToJpsi"}, {"PromptJpsi", "JpsiMC_uncorr"},
				  {"dimuTrk"}, {"flipJpsi","flipJpsiSameSide" }};
  vector<int> systIdx{0,0,0,0,0,
      (int)bToJpsiOnly,//nonprompt
      0,//prompt
      0,//dimuTrk
      (int)flipJpsiSameSide,//flipJpsi
      };

  bool drawShape[] = {true,true,false,true,true,true,//nonprompt
		      !ispp,//prompt
		      false,//dimutrk
		      useFlipJpsi//flipJpsi
  }; //mostly for drawing WrongSign on top
  bool usedForFit[] = {false,true,false,true,true,true,//nonprompt
		      !ispp,//prompt
		      false,//dimutrk
		      useFlipJpsi//flipJpsi
  };
  int nPiled = 0; for(int j=0;j<ntrees;j++) nPiled+= usedForFit[j];
  nPiled -= 1;
  int idxPiled[] = {1,5,8,4};//{1,5,8,4};
  if(!ispp) idxPiled[2] = 6; //PromptJpsi instead of flipJpsi

  TH1F* effSig = new TH1F( "effSig", "signal efficiency", _nChan(ispp)+1, BDTcut_l );
  TH1F* rejBkg = new TH1F( "rejBkg", "background rejection", _nChan(ispp)+1, BDTcut_l );
  TH1F* sigSignif = new TH1F( "sigSignif", "signal significance S/sqrt(S+B)", _nChan(ispp)+1, BDTcut_l );
  TH1F* sigSignifAsim = new TH1F( "sigSignifAsim", "signal significance (Asimov improved) S/sqrt(S+B)", _nChan(ispp)+1, BDTcut_l );

  //********************************************************
  //Extract BcM HISTOGRAMS used as combine input
  //********************************************************
  auto histFile = TFile::Open("InputForCombine_"+(TString)(BDTuncorrFromM?"BDTuncorrFromM_":"")+(TString)(ispp?"pp":"PbPb")+".root");
  for(int i=0; i<ntrees; i++){
    if(!drawShape[i]) continue;
    if(i==7) continue; //no dimutrk for now
    for(int k=0;k<=_nChan(ispp);k++){
      int kk = (k==0)?1:k;
      //h_BDT[i][k] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/"+procName[i][systIdx[i]]+"/BDTv");//new TH1F( "BDT_"+treeName[i]+"_BDTbin"+(TString)(to_string(k+1)), "BDTvalue "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp)[k])+","+(TString)(_BDTcut_s(ispp)[k+1])+"]", 50, BDTcut[0], BDTcut[_nChan(ispp)] );
      h_BcM[i][k] = (TH1F*)histFile->Get("BDT"+(TString)(to_string(kk))+"/"+procName[i][systIdx[i]]+"/BcM");//new TH1F( "Bc_M_"+treeName[i]+"_BDTbin"+(TString)(to_string(k+1)), "B_{c} mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp)[k])+","+(TString)(_BDTcut_s(ispp)[k+1])+"]", _nbinM(ispp), _Mbinning(ispp) );
      // h_BcPt[i][k] = new TH1F( "Bc_Pt_"+treeName[i]+"_BDTbin"+(TString)(to_string(k+1)), "B_{c} mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp)[k])+","+(TString)(_BDTcut_s(ispp)[k+1])+"]", 50, 0, 30 );
      // h_QQM[i][k] = new TH1F( "QQ_M_"+treeName[i]+"_BDTbin"+(TString)(to_string(k+1)), "J/#psi mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp)[k])+","+(TString)(_BDTcut_s(ispp)[k+1])+"]", nbinM, 2,4 );
      // h_QQ2M[i][k] = new TH1F( "QQ2_M_"+treeName[i]+"_BDTbin"+(TString)(to_string(k+1)), "QQ2 mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp)[k])+","+(TString)(_BDTcut_s(ispp)[k+1])+"]", nbinM, 2,4 );
      // h_QQ3M[i][k] = new TH1F( "QQ3_M_"+treeName[i]+"_BDTbin"+(TString)(to_string(k+1)), "QQ3 mass "+prettyName[i]+", BDT #in ["+(TString)(_BDTcut_s(ispp)[k])+","+(TString)(_BDTcut_s(ispp)[k+1])+"]", nbinM, 2,4 );
      
      if(k>1) {
	//h_BDT[i][0]->Add(h_BDT[i][k]);
	h_BcM[i][0]->Add(h_BcM[i][k]);}
    }
  }
  
  //********************************************************
  //Record PREFIT NORMALISATIONS
  //********************************************************
  vector<vector<float>> yields_prefit(ntrees, vector<float>(_nChan(ispp)+1));
  for(int i=0; i<ntrees; i++){
    if(i==7) continue;
    if(!drawShape[i]) continue;
    for(int k=0;k<=_nChan(ispp);k++){
      if(k==0) yields_prefit[i][k] = 0;
      else{ 
	yields_prefit[i][k] = h_BcM[i][k]->Integral();
	yields_prefit[i][0] += yields_prefit[i][k]; }
    }
  }

  //********************************************************
  //DRAWING Bc mass for all trees and all BDT bins, PREFIT
  //********************************************************
  gStyle->SetOptStat(0);

  Color_t cols[] = {kCyan, kMagenta+1, kGreen, kRed-2, kBlue, kOrange-7, kOrange, kGreen+4, kGreen+1};

  TCanvas* c3 = new TCanvas("c3", "Bc mass, in BDT bins, prefit", _nChan(ispp)*1500, 1500);
  c3->Divide(_nChan(ispp),1);
  TCanvas* c4 = new TCanvas("c4", "Bc mass, no BDT cuts, prefit", 1500, 1500);

  cout<<"\n----------------------\n Drawing PREFIT mass plots"<<endl;
  for(int k=0;k<=_nChan(ispp);k++){
    if(k==0) {
      c4->cd();
      c4->SetRightMargin(0.04);}
    else {
      c3->cd(k);
      c3->cd(k)->SetRightMargin(0.04);}

    cout<<"\n **** For BDTval in ["<<_BDTcut_s(ispp,BDTuncorrFromM)[max(k-1,0)]<<", "<<_BDTcut_s(ispp,BDTuncorrFromM)[(k==0)?_nChan(ispp):k]<<"]"<<endl;
    // cout<<"background efficiency = "<<rejBkg->GetBinContent(k+1)<<endl;
    // cout<<"MC signal efficiency = "<<effSig->GetBinContent(k+1)<<endl;
    // cout<<"Number of expected (MC) signal events = "<<nSig[k]<<endl;
    // cout<<"Number of events in signal region = "<<nData[k]<<endl;
    // cout<<"Theoretical significance S/sqrt(S+B) = "<<sigSignif->GetBinContent(k+1)<<endl;
    // cout<<"Theoretical significance (Asimov improved) S/sqrt(S+B) = "<<sigSignifAsim->GetBinContent(k+1)<<endl;
    
    //Formating all histos
    for(int i=0; i<ntrees; i++){
      if(i==7) continue;
      if(!drawShape[i]) continue;
      h_BcM[i][k]->SetFillStyle(usedForFit[i]?1001:0);
      h_BcM[i][k]->SetFillColor(cols[i]);
      h_BcM[i][k]->SetLineColor(((i==4)?kBlack:(cols[i]+2)));
      h_BcM[i][k]->SetLineWidth(3);
    }

    //Building the PILED (summed) histograms
    for(int i=0; i<nPiled; i++){
      int iold = idxPiled[i];
      h_BcM_prefit[i][k] = (TH1F*)h_BcM[iold][k]->Clone("Bc_M_prefit_"+treeName[iold]);
      
      if(i>0){
	h_BcM_prefit[i][k]->Add(h_BcM_prefit[i-1][k]);
      }

      h_BcM_prefit[i][k]->SetFillStyle(3002);
      h_BcM_prefit[i][k]->SetFillColor(cols[iold]);
      h_BcM_prefit[i][k]->SetLineColor(((iold==4)?kBlack:(cols[iold]+2)));
    }

    //DRAW data signal region
    //cout<<"--- Drawing: "<<prettyName[3]<<endl;
    h_BcM[3][k]->SetTitle("B_{c} candidates mass"+(TString)((k==0)?", no BDT cuts":(" with BDT #in "+(TString)((k==1)?"]-#infty":("["+_BDTcut_s(ispp,BDTuncorrFromM)[k-1]))+","+(TString)((k==_nChan(ispp))?"+#infty[":(_BDTcut_s(ispp,BDTuncorrFromM)[k]+"]"))) ));
    h_BcM[3][k]->GetXaxis()->SetTitle("trimuon mass [GeV]");
    h_BcM[3][k]->GetYaxis()->SetTitle("candidates/("+(TString)Form("%.2f",(m_Bc-3.3)/_nbinMSR(ispp)) +" GeV)");
    h_BcM[3][k]->GetYaxis()->SetRangeUser(0, 1.1*((h_BcM[3][k]->GetMaximum() > h_BcM_prefit[nPiled-1][k]->GetMaximum())?h_BcM[3][k]->GetMaximum():h_BcM_prefit[nPiled-1][k]->GetMaximum())  );
    h_BcM[3][k]->GetYaxis()->SetTitleOffset(1.2);
    h_BcM[3][k]->Draw("E");

    //DRAW PILED histos
    auto legend = new TLegend(0.63,0.63,0.96,0.9);
    for(int i=nPiled-1; i>=0; i--){
      int iold = idxPiled[i];
      //cout<<"--- Drawing piled backgrounds and expected signal: "<<prettyName[iold]<<endl;
      h_BcM_prefit[i][k]->Draw("histsame");      
      legend->AddEntry(h_BcM_prefit[i][k],  prettyName[iold], "lf");
    }

    //DRAW data signal region again
    //cout<<"--- Drawing again the data: "<<prettyName[3]<<endl;
    h_BcM[3][k]->Draw("Esame");
    legend->AddEntry(h_BcM[3][k], prettyName[3] ,"l");
    
    //DRAW remaining NON-PILED histos
    for(int i=0; i<ntrees; i++){
      if(drawShape[i] && !usedForFit[i]){
	//cout<<"--- Drawing other (non-piled) background: "<<prettyName[i]<<endl;
	h_BcM[i][k]->Draw("histsame");      
	legend->AddEntry(h_BcM[i][k],  prettyName[i], "l");
      }
    }

    legend->SetTextSize(0.034);
    legend->Draw("same");

    //DRAW CMS Preliminary
    TLatex CMStag;
    CMStag.SetNDC();
    CMStag.SetTextFont(42);
    CMStag.SetTextSize(0.038);
    CMStag.DrawLatex(0.13,0.85,"#font[61]{CMS "+(TString)(ispp?"pp":"PbPb")+"}");
    CMStag.DrawLatex(0.13,0.80,"#font[52]{Work in progress}");

    //Save canvas
    ((k==0)?c4:(c3->cd(k)))->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+(TString)((k==0)?"_noBDTcut":"_BDTbin"+(TString)(to_string(k)))+"_PREFIT.pdf");
    ((k==0)?c4:(c3->cd(k)))->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+(TString)((k==0)?"_noBDTcut":"_BDTbin"+(TString)(to_string(k)))+"_PREFIT.png");
    if(k==_nChan(ispp)) {
      c3->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_allBDTbins_PREFIT.png");
      c3->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_allBDTbins_PREFIT.pdf");}

  }

  //********************************************************
  //Extract results (POSTFIT NORMALISATIONS & SHAPES) from combine output
  //********************************************************
  TString normFileName = "./CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+(TString)(bToJpsiOnly?"bToJpsi":"NonPromptJpsi")+(TString)(useFlipJpsi?(flipJpsiSameSide?"_flipJpsiSameSide":"_flipJpsi"):(ispp?"":"_PromptJpsi"))+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ignoreBin2?"_ignoreBin1-2":(ignoreBin1?"_ignoreBin1":""))+(TString)(ispp?"_pp":"_PbPb")+".root";
    cout<<"Extract normalisations from file "<<normFileName<<endl;
  auto normFile = TFile::Open(normFileName);
  RooArgSet *Yields = (RooArgSet*)normFile->Get("norm_fit_"+(TString)(bkgOnly?"b":"s"));
  vector<vector<float>> yields(ntrees, vector<float>(_nChan(ispp)+1));
  for(int i=0; i<ntrees; i++){
    if(i==7) continue;
    if(!drawShape[i]) continue;
    for(int k=0;k<=_nChan(ispp);k++){

      //Clone trimuon mass histos prefit into postfit. Processes not used in fitting will not change histos postfit
      h_BcM_postfit[i][k] = (TH1F*)h_BcM[i][k]->Clone("Bc_M_postfit_"+treeName[i]); //Do this to recover proper X axis
      h_BcM_postfit[i][k]->SetDirectory(0);

      if(i==3 || !usedForFit[i])// || (k==1 && ignoreBin1) || (k==2 && ignoreBin2)) 
	{yields[i][k] = yields_prefit[i][k];}
      else{

	if(k==0) yields[i][k] = 0;
	else{ 
	  //fetch postfit yields
	  yields[i][k] = Yields->getRealValue("BDT"+(TString)(to_string(k))+"/"+procNameDef[i]);
	  yields[i][0] += yields[i][k]; 
	}

	//fetch postfit shapes
	TH1F* h_tmp = (k==0)?NULL:((TH1F*) normFile->Get("shapes_fit_"+(TString)(bkgOnly?"b":"s")+"/BDT"+(TString)(to_string(k))+"/"+procNameDef[i]));
	for(int bin=0;bin<=h_BcM_postfit[i][k]->GetNbinsX();bin++){
	  h_BcM_postfit[i][k]->SetBinContent(bin, (k==0)?0:(h_tmp->GetBinContent(bin)));
	  h_BcM_postfit[i][k]->SetBinError(bin, (k==0)?0:(h_tmp->GetBinError(bin)));
	}
	delete h_tmp;
      }

      //inclusive on BDT bins
      if(k>0) h_BcM_postfit[i][0]->Add(h_BcM_postfit[i][k]);
	
    }
    if(i!=3 && usedForFit[i]) {
      if(ignoreBin1){
	yields[i][1] *= yields[i][ignoreBin2?3:2]/yields_prefit[i][ignoreBin2?3:2];
	yields[i][0] += yields[i][1];}
      if(ignoreBin2){
	yields[i][2] *= yields[i][3]/yields_prefit[i][3];
	yields[i][0] += yields[i][2];}
    }
    for(int k=0;k<=_nChan(ispp);k++) cout<<"yield tree "<<i<<" BDTbin "<<k<<" = "<<yields[i][k] <<"    ratio to prefit = "<<yields[i][k]/yields_prefit[i][k]<<endl;
    
  }
  normFile->Close();

  // //********************************************************
  // //build inclusive histo (and SCALE histos to normalisations from combine)
  // //********************************************************
  // for(int i=0; i<ntrees; i++){
  //   if(!drawShape[i]) continue;
  //   for(int k=0;k<=_nChan(ispp);k++){
  //     if(i!=3 && usedForFit[i]){//no yield scaling for data, nor for the shapes not used in the fit
  //     	h_BDT[i][k]->Scale( yields[i][k] / h_BDT[i][k]->Integral());
  //     	h_BcM[i][k]->Scale( yields[i][k] / h_BcM[i][k]->Integral());}
  //     if(k==1) h_BDT[i][0] = h_BDT[i][k];
  //     if(k>1) {
  //     	h_BDT[i][0]->Add(h_BDT[i][k]);
  //     	h_BcM[i][0]->Add(h_BcM[i][k]);}
  //   }    
  // }
 
  //********************************************************
  //Measure SIGNAL SIGNIFICANCE and yields
  //********************************************************
  vector<float> nBkg = vector<float>(_nChan(ispp)+1); //initialized to zeros
  vector<float> nSig = vector<float>(_nChan(ispp)+1);
  vector<float> nData = vector<float>(_nChan(ispp)+1);
  for(int k=0;k<=_nChan(ispp);k++){
    nData[k] = yields[3][k];
    nSig[k] = yields[4][k];
    nBkg[k] = 0;
    for(int i=0; i<ntrees; i++){
      if(!usedForFit[i] || i==3 || i==4) continue; //only yields of histos used in combine fit
      nBkg[k] += yields[i][k];
    }
  }

  //Fill histos of efficiencies and significance
  for(int k=0;k<=_nChan(ispp);k++){    
    effSig->SetBinContent(k+1, nSig[k]/nSig[0]);
    rejBkg->SetBinContent(k+1, 1- nBkg[k]/nBkg[0]);
    sigSignif->SetBinContent(k+1, (nSig[k]+nBkg[k]>0)?( 
						       nSig[k]/sqrt(nBkg[k]+nSig[k]) ):0 );
    sigSignifAsim->SetBinContent(k+1, (nSig[k]+nBkg[k]>0 && ((2*nSig[k]+nBkg[k])*log(1+nSig[k]/(nSig[k]+nBkg[k])) > nSig[k]) )?(   
																sqrt(2*( (2*nSig[k]+nBkg[k])*log(1+nSig[k]/(nSig[k]+nBkg[k])) - nSig[k]))   ):0);
  }
  
  //********************************************************
  //Make PULL graphs, and measure chi2
  //********************************************************
  float pullMax = 3.8;
  vector<int> ndf = vector<int>(_nChan(ispp)+1 , _nbinM(ispp)-nParam);
  vector<float> chi2 = vector<float>(_nChan(ispp)+1, 0);
  const int n = h_BcM[1][1]->GetNbinsX();
  double x[_nChan(ispp)+1][n];
  double y[_nChan(ispp)+1][n];
  for(int k=0;k<=_nChan(ispp);k++){
    bool seen = false;
    for(int i=0; i<ntrees; i++){
      if(!usedForFit[i] || i==3) continue;
      if(!seen) {h_BcM_exp[k] = (TH1F*) h_BcM_postfit[i][k]->Clone();
	seen=true;}
      else h_BcM_exp[k]->Add(h_BcM_postfit[i][k]); 
    }

    for(int b=1; b<=n; b++){
      x[k][b-1] = h_BcM_exp[k]->GetXaxis()->GetBinCenter(b);

      float binc = h_BcM[3][k]->GetBinContent(b);
      float bine = h_BcM[3][k]->GetBinError(b);
      float binc2 = h_BcM_exp[k]->GetBinContent(b);
      y[k][b-1] = (bine<=0)?0:( (binc-binc2)/bine );
      if(y[k][b-1]>=pullMax) y[k][b-1] = pullMax*0.99999; //for plotting
      if(y[k][b-1]<=-pullMax) y[k][b-1] = -pullMax*0.99999; //for plotting
      chi2[k] += pow(y[k][b-1],2);
      if(y[k][b-1]==0) ndf[k] -= 1;
      //cout<<b<<" "<<x[k][b-1]<<" "<<y[k][b-1]<<endl;
    }
    cout<<"chi2/ndf data vs expectation = "<<chi2[k]/ndf[k]<<endl;//" "<<chi2[k]<<" "<<ndf[k]<<endl;

    h_BcM_pull[k] = new TGraph(n,x[k],y[k]);
    h_BcM_pull[k]->SetMarkerColor(kBlack);    
    h_BcM_pull[k]->SetFillColor(kBlack);    
    h_BcM_pull[k]->SetMarkerSize(1.);    
  }

  //dummy histo to draw axis in the pull graph's pad
  float padratio = 5.5, padcorr=0.9;
  TH1F* dummy = new TH1F("dummy", "", _nbinM(ispp), _Mbinning(ispp));
  dummy->GetXaxis()->SetTickLength(padcorr*padratio * h_BcM_exp[0]->GetXaxis()->GetTickLength());
  dummy->GetXaxis()->SetLabelSize(padcorr*padratio * h_BcM_exp[0]->GetLabelSize("X"));
  dummy->GetYaxis()->SetLabelSize(padcorr*padratio * h_BcM_exp[0]->GetLabelSize("Y"));
  dummy->GetYaxis()->SetRangeUser(-pullMax,pullMax);  
  dummy->GetYaxis()->SetNdivisions(7);
  //  dummy->GetXaxis()->SetRangeUser(h_BcM_exp[0]->GetBinLowEdge(1), h_BcM_exp[0]->GetBinLowEdge(n+1));
  dummy->GetXaxis()->SetTitle("trimuon mass [GeV]");    
  dummy->GetYaxis()->SetTitle("pull");    
  dummy->GetYaxis()->SetTitleSize(padcorr*padratio * h_BcM_exp[0]->GetXaxis()->GetTitleSize()); 
  dummy->GetXaxis()->SetTitleSize(padcorr*padratio * h_BcM_exp[0]->GetXaxis()->GetTitleSize()); 
  dummy->GetYaxis()->SetTitleOffset(0.19); 
  dummy->GetXaxis()->SetLabelOffset(0.03); 

  //********************************************************
  //DRAWING Bc mass for all trees and all BDT bins, POSTFIT
  //********************************************************
  TCanvas* c6 = new TCanvas("c6", "Bc mass, in BDT bins, postfit", _nChan(ispp)*1500, 1500);
  c6->Divide(_nChan(ispp),2);
  TCanvas* c7 = new TCanvas("c7", "Bc mass, no BDT cuts, postfit", 1500, 1500);
  c7->Divide(1,2);

  cout<<"\n----------------------\n Drawing POSTFIT mass plots"<<endl;
  for(int k=0;k<=_nChan(ispp);k++){
    if(k==0) {  
      c7->cd(2)->SetPad(0.,0.,1.,1/padratio);//lower pad for pull
      gPad->SetTickx(2);
      gPad->SetGridy();
      gPad->SetMargin(0.1,0.04, 0.38,0.);//left,right,bottom,top
      dummy->DrawClone();
      c7->cd(1)->SetPad(0.,1/padratio,1.,1.); //upper pad for trimuon mass distro
      gPad->SetMargin(0.1,0.04, 0.,0.1);
      c7->cd(1);
    }
    else {
      c6->cd(_nChan(ispp)+k)->SetPad((k-1)/(float)_nChan(ispp),0.,k/(float)_nChan(ispp),1/(padratio+1));//lower pad for pull
      gPad->SetTickx(2);
      gPad->SetGridy();
      gPad->SetMargin(0.1,0.04, 0.38,0.);//left,right,bottom,top
      dummy->DrawClone();
      c6->cd(k)->SetPad((k-1)/(float)_nChan(ispp),1/(padratio+1),k/(float)_nChan(ispp),1.);//upper pad for trimuon mass distro
      gPad->SetMargin(0.1,0.04, 0.,0.1);
      c6->cd(k);
    }

    cout<<"\n **** For BDTval in ["<<_BDTcut_s(ispp,BDTuncorrFromM)[max(k-1,0)]<<", "<<_BDTcut_s(ispp,BDTuncorrFromM)[(k==0)?_nChan(ispp):k]<<"]"<<endl;
    cout<<"background rejection = "<<rejBkg->GetBinContent(k+1)<<endl;
    if(!bkgOnly) cout<<"MC signal efficiency = "<<effSig->GetBinContent(k+1)<<endl;
    if(!bkgOnly) cout<<"Number of expected (MC postfit norm) signal events = "<<nSig[k]<<endl;
    cout<<"Number of events in signal region = "<<nData[k]<<endl;
    if(!bkgOnly) cout<<"Theoretical significance S/sqrt(S+B) = "<<sigSignif->GetBinContent(k+1)<<endl;
    if(!bkgOnly) cout<<"Theoretical significance (Asimov improved) S/sqrt(S+B) = "<<sigSignifAsim->GetBinContent(k+1)<<endl;

    auto legend = new TLegend(0.63,0.63,0.96,0.9);

    //DRAW data signal region
    h_BcM[3][k]->Draw("E");

    //Building the PILED (summed) histograms
    for(int i=0; i<nPiled; i++){
      int iold = idxPiled[i];
      h_BcM_DrawPostfit[i][k] = (TH1F*)h_BcM_postfit[iold][k]->Clone("Bc_M_DrawPostfit_"+treeName[iold]);    

      if(i>0){
	h_BcM_DrawPostfit[i][k]->Add(h_BcM_DrawPostfit[i-1][k]);
      }

      h_BcM_DrawPostfit[i][k]->SetFillStyle(3002);
      h_BcM_DrawPostfit[i][k]->SetFillColor(cols[iold]);
      h_BcM_DrawPostfit[i][k]->SetLineColor(((iold==4)?kBlack:(cols[iold]+2)));
    }

    //DRAW PILED histos
    for(int i=nPiled-1; i>=0; i--){
      int iold = idxPiled[i];
      //cout<<"--- Drawing piled backgrounds and expected signal: "<<prettyName[iold]<<endl;
      h_BcM_DrawPostfit[i][k]->Draw("histsame");      
      legend->AddEntry(h_BcM_DrawPostfit[i][k],  prettyName[iold], "lf");
    }

    //DRAW data signal region again
    h_BcM[3][k]->GetYaxis()->SetRangeUser(0, 1.1*((h_BcM[3][k]->GetMaximum() > h_BcM_DrawPostfit[nPiled-1][k]->GetMaximum())?h_BcM[3][k]->GetMaximum():h_BcM_DrawPostfit[nPiled-1][k]->GetMaximum())  );
    //cout<<"--- Drawing signal region: "<<prettyName[3]<<endl;
    h_BcM[3][k]->Draw("Esame");
    legend->AddEntry(h_BcM[3][k], prettyName[3] ,"l");

    //DRAW remaining NON-PILED histos
    for(int i=0; i<ntrees; i++){
      if(drawShape[i] && !usedForFit[i]){
	//cout<<"--- Drawing other (non-piled) background: "<<prettyName[i]<<endl;
	h_BcM[i][k]->Draw("histsame");      
	legend->AddEntry(h_BcM[i][k],  prettyName[i], "l");
      }
    }

    legend->SetTextSize(0.034);
    legend->Draw("same");

    // //ADDING normalization error to MC expectation + b->Jpsi MC
    // float MCnormRelErr = 0.34/2.54;
    // float MCbnormRelErr = quadSum(0.031,0.167, 2.596*0.023) / 2.596; //last error, from the normalization fit, must be checked
    // float MCpromptnormRelErr = quadSum(0.12,0.59, 10.04*0.012) / 10.04; //last error, from the normalization fit, must be checked
    // for (int i=0;i<nbinM;i++){
    //   h_BcM_DrawPostfit[nPiledTrees-1][k]->SetBinError(i, quadSum( MCpromptnormRelErr*h_BcM[6][k]->GetBinContent(i) , MCbnormRelErr*h_BcM[5][k]->GetBinContent(i) , MCnormRelErr*h_BcM[4][k]->GetBinContent(i) , h_BcM_DrawPostfit[nPiledTrees-1][k]->GetBinError(i) ) );
    // }
       
    //DRAW values of efficiencies and significance
    char eff1[100]; sprintf(eff1,"#varepsilon_{signal} = %.3f",effSig->GetBinContent(k+1));
    char eff2[100]; sprintf(eff2,"#varepsilon_{background} = %.3f",1-rejBkg->GetBinContent(k+1));
    char eff3[100]; sprintf(eff3,"S/#sqrt{S+B} = %.1f",sigSignif->GetBinContent(k+1));
    char eff4[100]; sprintf(eff4,"N_{signal}^{postfit} = %.0f",nSig[k]);
    char s_chi2[100]; sprintf(s_chi2,"#chi^{2}/ndf = %.2f",chi2[k]/ndf[k]);
    TLatex eAndSignif;
    eAndSignif.SetNDC();
    eAndSignif.SetTextFont(52);
    eAndSignif.SetTextSize(0.038);
    eAndSignif.DrawLatex(0.7,0.29,s_chi2);
    if(!bkgOnly) eAndSignif.DrawLatex(0.7,0.36,eff4);
    if(!bkgOnly) eAndSignif.DrawLatex(0.7,0.43,eff3);
    eAndSignif.DrawLatex(0.7,bkgOnly?0.36:0.5,eff2);
    if(!bkgOnly) eAndSignif.DrawLatex(0.7,0.57,eff1);

    //DRAW CMS Preliminary
    TLatex CMStag;
    CMStag.SetNDC();
    CMStag.SetTextFont(42);
    CMStag.SetTextSize(0.038);
    CMStag.DrawLatex(0.13,0.85,"#font[61]{CMS "+(TString)(ispp?"pp 2017":"PbPb 2018")+"}");
    CMStag.DrawLatex(0.13,0.80,"#font[52]{Work in progress}");

    //DRAW PULL histos
    if(k==0) c7->cd(2);
    else c6->cd(_nChan(ispp)+k);
    h_BcM_pull[k]->Draw("Bsame");


    //SAVE canvases
    //For inclusive plot
    if(k==0){
      c7->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_noBDTcut_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+".pdf");
      c7->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_noBDTcut_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+".png");}
    
    //For each BDT bin
    else{    
      TCanvas* ctmp = new TCanvas("ctmp", "Bc mass in a BDT bin", 1500, 1500);
      TPad* ctmp1 = (TPad*)(c6->cd(k))->Clone();
      ctmp1->SetPad(0.,0.2,1.,1.);
      ctmp->cd(); 
      ctmp1->DrawClone();
      TPad* ctmp2 = (TPad*)(c6->cd(k+_nChan(ispp)))->Clone();
      ctmp2->SetPad(0.,0.,1.,0.2);  //necessary to recreate shape of canvas
      ctmp->cd(); 
      ctmp2->DrawClone();
      
      ctmp->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_BDTbin"+(TString)(to_string(k))+"_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+".pdf");
      ctmp->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_BDTbin"+(TString)(to_string(k))+"_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+".png");
      delete ctmp,ctmp1,ctmp2;
    }
    
    //Canvas with all BDT bins
    if(k==_nChan(ispp)) {
      c6->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_allBDTbins_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+".pdf");
      c6->SaveAs("figs/Bc_mass_"+(TString)(ispp?"pp":"PbPb")+"_allBDTbins_POSTFIT"+(TString)(bkgOnly?"_bkgOnly":"")+".png");}

  }

  //********************************************************
  //DRAWING signal and background efficiencies, and the signal significances, versus BDT lower cut value
  //********************************************************
  if(bkgOnly) return;

  TCanvas* c5 = new TCanvas("c5", "Efficiencies and significance", 2000, 1500);
  effSig->SetLineColor(kRed);
  effSig->SetLineWidth(3);
  effSig->GetYaxis()->SetRangeUser(0,1);
  effSig->GetXaxis()->SetTitle("BDT");
  effSig->Draw("hist");
  rejBkg->SetLineWidth(3);
  rejBkg->Draw("histsame");
  c5->Update();

  //scale sigSignif to the pad coordinates
  Float_t rightmax = 1.1*sigSignif->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;
  sigSignif->SetLineWidth(3);
  sigSignif->SetLineColor(kGreen+3);
  sigSignif->Scale(scale);
  sigSignif->Draw("histsame");
  sigSignifAsim->SetLineWidth(3);
  sigSignifAsim->SetLineColor(kGreen);
  sigSignifAsim->Scale(scale);
  sigSignifAsim->Draw("histsame");

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

  auto legend2 = new TLegend(0.11,0.65,0.4,0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->SetTextSize(0.031);
  legend2->AddEntry(effSig, "signal efficiency" );
  legend2->AddEntry(rejBkg, "background rejection" );
  legend2->AddEntry(sigSignif, "signal significance" );
  legend2->AddEntry(sigSignifAsim, "signal significance (Asimov)" );
  legend2->Draw("same");

  if(!bkgOnly){
    c5->SaveAs("figs/efficiencySigBkgAndSignificance_"+(TString)(ispp?"pp":"PbPb")+"_vs_BDTbin.png");
    c5->SaveAs("figs/efficiencySigBkgAndSignificance_"+(TString)(ispp?"pp":"PbPb")+"_vs_BDTbin.pdf");}
  
}

