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

vector<float> DetermineCuts(bool ispp=true, bool BDTuncorrFromM=false, int kinBin=0, int centBin=0, bool withTM=false, int varyBDTbin = 0, bool secondStep=false, bool firstTime=false){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  if(firstTime) BDTuncorrFromM=false;

  auto fullFile = TFile::Open("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  TTree* T = (TTree*)fullFile->Get("signal_MC");

  //Initialization of variables
  float Bc_M;
  float Bc_Pt;
  float Bc_Y;
  float weight;
  float w_simple2;
  float BDT;
  int hiBin;
  T->SetBranchAddress("Bc_M", &Bc_M);
  T->SetBranchAddress("Bc_Pt", &Bc_Pt);
  T->SetBranchAddress("Bc_Y", &Bc_Y);
  T->SetBranchAddress(secondStep?"BDT2":"BDT", &BDT);
  if(firstTime)
    T->SetBranchAddress("w_simple2", &w_simple2);
  else
    T->SetBranchAddress(secondStep?"weight2":"weight", &weight);
  if(!ispp) T->SetBranchAddress("Centrality", &hiBin);

  //*******************************************
  //Fetch BDT correction = f(M) to be subtracted from BDT, to uncorrelate it from mass
  TFile* f_BDTuncorrel = new TFile("../BDT/BDTuncorrFromM_"+(TString)(ispp?"pp":"PbPb")+".root","READ");
  TString s_Bin = (centBin==0)?((kinBin==0)?"_integrated":("_KinBin"+(TString)to_string(kinBin))):("_CentBin"+(TString)to_string(centBin));
  TH1F* h_correctBDT = BDTuncorrFromM?( (TH1F*) f_BDTuncorrel->Get("avBDTvsM_bkg"+s_Bin+(TString)(secondStep?"_2ndStep":"")) ):NULL; //make the BDT of the expected background (postfit) uncorrelated with mass
  TH1F* h_correctrmsBDT = BDTuncorrFromM?( (TH1F*) f_BDTuncorrel->Get("rmsBDTvsM_bkg"+s_Bin+(TString)(secondStep?"_2ndStep":"")) ):NULL;
  TH1F* h_BDT = new TH1F("bdt","bdt",300,BDTuncorrFromM?(-3.5):(-1.5),BDTuncorrFromM?5.5:1.5);

  //BEGIN event loop
  for(int j=0; j<T->GetEntries(); j++){
    T->GetEntry(j);

    float bdtcorr = BDTuncorrFromM?( h_correctBDT->Interpolate(Bc_M) ):0.;
    float bdtcorrrms = BDTuncorrFromM?( h_correctrmsBDT->Interpolate(Bc_M) ):1.;

    //Keep events from the wanted analysis bin
    if(!inFidCuts(kinBin,Bc_Pt,Bc_Y)) continue;
    if(!ispp && ((float)hiBin < 2*_Centmin[centBin] || (float)hiBin >= 2*_Centmax[centBin])) continue;

    h_BDT->Fill((BDT-bdtcorr)/bdtcorrrms, firstTime?w_simple2:weight);

  }
  //END event loop
  
  vector<float> res(0);//size 0
  float ntot = h_BDT->Integral();
  float ncur = 0;
  float efflim1 = _withTM?0.25:(0.25+varyBDTbin*0.05);
  float efflim2 = _withTM?0.65:(0.65+varyBDTbin*0.1);

  for(int b=1; b<=h_BDT->GetNbinsX();b++){
    //    if(h_BDT->GetBinContent(b)>0 && res.size()==0) res.push_back(h_BDT->GetBinLowEdge(max(b-3,1)));
    ncur += h_BDT->GetBinContent(b);
    if(ncur/ntot>0.001 && res.size()==0) res.push_back(h_BDT->GetBinLowEdge(b));
    if(_withTM && ncur/ntot>0.05 && res.size()==1) res.push_back(h_BDT->GetBinLowEdge(b));
    if(ncur/ntot>efflim1 && res.size()==(_withTM?2:1)) res.push_back(h_BDT->GetBinLowEdge(b));
    if(ncur/ntot>efflim2 && res.size()==(_withTM?3:2)) res.push_back(h_BDT->GetBinLowEdge(b));
    if(ncur/ntot>0.99999 && res.size()==(_withTM?4:3)) res.push_back(h_BDT->GetBinLowEdge(b+2));
  }
  if(ncur/ntot<=0.99999 && res.size()==(_withTM?4:3)) res.push_back(h_BDT->GetBinLowEdge(h_BDT->GetNbinsX()));

  return res;
}

void DetermineBDTcuts(bool secondStep=false, bool firstTime=false, bool firstStepForBDTuncorr=false){

  ofstream outfile;
  outfile.open("../helpers/Cuts_BDT_preliminary.h");
  outfile<<"#include \"Definitions.h\"\n\n";
  outfile<< "std::vector<float> _BDTcuts(bool ispp, int kinBin=0, int centBin=0, bool secondStep=false, bool BDTuncorrFromM=false, int varyBDTbin=0){\n";
  for(bool secStep : { true, false }){
    if(!secondStep && secStep) continue;
    if(secondStep) outfile<< "  "<<(secStep?"if(secondStep){":"else{") <<"\n";
    for(bool type : { true, false }){
      outfile<< "  "<<(type?"if(ispp){":"else{") <<"\n";
      for(bool withtm : { true, false }){
	outfile<< "    "<<(withtm?"if(_withTM){":"else{") <<"\n";
	for(bool BDTuncorr : { true, false }){
	  outfile<< "      "<<(BDTuncorr?"if(BDTuncorrFromM){":"else{") <<"\n";
	  for(int kinb=0;kinb<=_NanaBins;kinb++){
	    outfile<< "        if(kinBin=="<<kinb<<"){\n";
	    for(int centb=0;centb<=_NanaBins;centb++){
	      if(centb>0 && (type || kinb>0)) continue; //centrality bin only in PbPb and for integrated pT bin
	      if(!type && kinb==0) outfile<< "          if(centBin=="<<centb<<"){\n";
	      for(int varyBinning=-1;varyBinning<=1;varyBinning++){
		if(withtm && varyBinning!=0) continue; //vary BDT binning only for !_withtm
		vector<float> res = DetermineCuts(type, BDTuncorr, kinb, centb, withtm, varyBinning, BDTuncorr?(secStep && !firstStepForBDTuncorr):secStep, firstTime);
		if(!withtm) outfile<< "            if(varyBDTbin=="<<varyBinning<<"){\n";
		outfile<< "              return std::vector<float>{";
		for(int ix=0;ix<res.size();ix++){
		  outfile<< Form("%.2f",res[ix]);
		  outfile<< ((ix==res.size()-1)?"};":",");
		}
		if(!withtm) outfile<< "}\n";
	      }
	      if(withtm) {outfile<< "}\n";} else {outfile<< "          }\n";}
	    }
	    if(!type && kinb==0) outfile<< "        }\n";
	  }
	  outfile<< "      }\n";
	}
	outfile<< "    }\n";
      }
      outfile<< "  }\n";
    }
    if(secondStep) outfile<< "  }\n";
  }

  outfile<< "  return std::vector<float>{};\n}\n";
  outfile.close();
  
}
