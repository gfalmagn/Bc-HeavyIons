#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLatex.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

vector<float> DetermineCuts(bool ispp=true, bool BDTuncorrFromM=false, int kinBin=0, bool withTM=false, bool firstTime=false){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  auto fullFile = TFile::Open("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  TTree* T = (TTree*)fullFile->Get("signal_MC");

  //Initialization of variables
  float Bc_M;
  float Bc_Pt;
  float Bc_Y;
  float weight;
  float w_simple2;
  float BDT;
  T->SetBranchAddress("Bc_M", &Bc_M);
  T->SetBranchAddress("Bc_Pt", &Bc_Pt);
  T->SetBranchAddress("Bc_Y", &Bc_Y);
  T->SetBranchAddress("BDT", &BDT);
  if(firstTime)
    T->SetBranchAddress("w_simple2", &w_simple2);
  else
    T->SetBranchAddress("weight", &weight);

  //*******************************************
  //Fetch BDT correction = f(M) to be subtracted from BDT, to uncorrelate it from mass
  TFile* f_BDTuncorrel = new TFile("../templateFit/BDTuncorrFromM_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  TH1F* h_correctBDT = BDTuncorrFromM?( (TH1F*) f_BDTuncorrel->Get("avBDTvsM_bkg") ):NULL; //make the BDT of the expected background (postfit) uncorrelated with mass

  TH1F* h_BDT = new TH1F("bdt","bdt",300,-1.5,1.5);

  //BEGIN event loop
  for(int j=0; j<T->GetEntries(); j++){
    T->GetEntry(j);

    float bdtcorr = BDTuncorrFromM?( h_correctBDT->GetBinContent(h_correctBDT->FindBin(Bc_M)) ):0;    

    //Keep events from the wanted analysis bin
    if(kinBin>0 && !( Bc_Pt>_BcPtmin[kinBin] && Bc_Pt<_BcPtmax[kinBin] && fabs(Bc_Y)>_BcYmin[kinBin] && fabs(Bc_Y)<_BcYmax[kinBin])) continue;
    if(kinBin==0){//integrated bin
      bool inFidCuts = Bc_Pt>_BcPtmin[0] && Bc_Pt<_BcPtmax[0] && fabs(Bc_Y)>_BcYmin[0] && fabs(Bc_Y)<_BcYmax[0];
      inFidCuts = inFidCuts || (Bc_Pt>_BcPtmin[1] && Bc_Pt<_BcPtmax[1] && fabs(Bc_Y)>_BcYmin[1] && fabs(Bc_Y)<_BcYmax[1]);
      if(!inFidCuts) continue;
    }

    h_BDT->Fill(BDT-bdtcorr, firstTime?w_simple2:weight);

  }
  //END event loop
  
  vector<float> res(0);//size 0
  float ntot = h_BDT->Integral();
  float ncur = 0;
  float efflim1 = _withTM?0.25:0.2;
  float efflim2 = _withTM?0.65:0.6;

  for(int b=1; b<=h_BDT->GetNbinsX();b++){
    if(h_BDT->GetBinContent(b)>0 && res.size()==0) res.push_back(h_BDT->GetBinLowEdge(b));
    ncur += h_BDT->GetBinContent(b);
    if(_withTM && ncur/ntot>0.05 && res.size()==1) res.push_back(h_BDT->GetBinLowEdge(b));
    if(ncur/ntot>efflim1 && res.size()==(_withTM?2:1)) res.push_back(h_BDT->GetBinLowEdge(b));
    if(ncur/ntot>efflim2 && res.size()==(_withTM?3:2)) res.push_back(h_BDT->GetBinLowEdge(b));
    if(ncur/ntot>0.9995 && res.size()==(_withTM?4:3)) res.push_back(h_BDT->GetBinLowEdge(b+3));
  }
  if(ncur/ntot<=0.9995 && res.size()==(_withTM?4:3)) res.push_back(h_BDT->GetBinLowEdge(h_BDT->GetNbinsX()));

  return res;
}

void DetermineBDTcuts(bool firstTime=false){

  ofstream outfile;
  outfile.open("../helpers/Cuts_BDT_preliminary.h");
  outfile<<"#include \"Definitions.h\"\n\n";
  outfile<< "std::vector<float> _BDTcuts(bool ispp, int kinBin=0, bool BDTuncorrFromM=false){\n";
  for(bool type : { true, false }){
    outfile<< "  "<<(type?"if(ispp){":"else{") <<"\n";
    for(bool withtm : { true, false }){
      outfile<< "    "<<(withtm?"if(_withTM){":"else{") <<"\n";
      for(bool BDTuncorr : { true, false }){
	outfile<< "      "<<(BDTuncorr?"if(BDTuncorrFromM){":"else{") <<"\n";
	for(int kinb=0;kinb<=_NanaBins;kinb++){
	  outfile<< "        if(kinBin=="<<kinb<<"){\n";
	  vector<float> res = DetermineCuts(type, BDTuncorr, kinb, withtm, firstTime);
	  outfile<< "          return std::vector<float>{";
	  for(int ix=0;ix<res.size();ix++){
	    outfile<< Form("%.2f",res[ix]);
	    outfile<< ((ix==res.size()-1)?"};}\n":",");
	  }
	}
	outfile<< "      }\n";
      }
      outfile<< "    }\n";
    }
    outfile<< "  }\n";
  }

  outfile<< "  return std::vector<float>{};\n}\n";
  outfile.close();
  
}
