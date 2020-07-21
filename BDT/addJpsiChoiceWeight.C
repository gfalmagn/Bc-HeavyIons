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
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"

void addJpsiChoiceW(bool ispp=true, bool useBDTbins=false, vector<float> BDTcuts = vector<float>()){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  //need to start from untouched saved file
  //  gSystem->Exec("cp BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+"_copy"+(TString)(useBDTbins?"step2":"")+".root BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root");
  //if(useBDTbins) System->Exec("cp BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+"_copystep2.root");
  auto fullFile = TFile::Open("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+".root","UPDATE");
  //  fullFile->Cp("BDT_InputTree_"+(TString)(ispp?"pp":"PbPb")+"_copystep"+(TString)(useBDTbins?"2":"0")+".root");
  int ntrees = 9;

  float muW_eta[ntrees];
  float mumi_eta[ntrees];
  float mupl_eta[ntrees];
  float Bc_Pt[ntrees];
  float Bc_Y[ntrees];
  float QQ_M[ntrees];
  float QQ2_M[ntrees];
  float QQ2_dca[ntrees];
  float QQ2_VtxProb[ntrees];
  float w_simple[ntrees];
  float w_unblind[ntrees];
  float weight[ntrees];
  float BDT[ntrees];

  TString treeName[] = {"bkgWRONGSIGN","bkgBCMASS","bkgTRUEJPSI","sigRegion","signal_MC","bToJpsi_MC","PromptJpsi_MC","dimuonTrk","flipJpsi"};
  TString prettyName[] = {"WRONGSIGN","J/Psi sidebands","High mass control","signal region","MC signal expectation",
			  "MC NonPromptJpsi","MC PromptJpsi","dimuon+track (misID)","flipped J/Psi"};
  vector<TTree*> T;
  for(int itree=0;itree<ntrees;itree++){
    T.push_back((TTree*)fullFile->Get(treeName[itree]));
  }

  //***************** Jpsi mass histos and binning
  vector<vector<TH1F*>> h_QQM;
  vector<vector<TH1F*>> h_QQM_tight;
  int nbinsPeak = 10, nbinsPeakT = 8;
  float binsM[7+nbinsPeak]; 
  binsM[0] = m_Jpsi-0.4;
  binsM[1] = m_Jpsi-JpeakLo-JpeakLoBuf-JpeakWid/2;
  binsM[2] = m_Jpsi-JpeakLo-JpeakLoBuf;
  for(int i=0; i<=nbinsPeak;i++) binsM[3+i] = m_Jpsi-JpeakLo+i*JpeakWid/nbinsPeak;
  binsM[nbinsPeak+4] = m_Jpsi+JpeakHi+JpeakHiBuf;
  binsM[nbinsPeak+5] = m_Jpsi+JpeakHi+JpeakHiBuf+JpeakWid/2;
  binsM[nbinsPeak+6] = m_Jpsi+0.4;

  float binsMT[7+nbinsPeakT];
  binsMT[0] = m_Jpsi-0.4;
  binsMT[1] = m_Jpsi-JpeakLoT-JpeakLoBuf-JpeakWidT/2;
  binsMT[2] = m_Jpsi-JpeakLoT-JpeakLoBuf;
  for(int i=0; i<=nbinsPeakT;i++) binsMT[3+i] = m_Jpsi-JpeakLoT+i*JpeakWidT/nbinsPeakT;
  binsMT[nbinsPeakT+4] = m_Jpsi+JpeakHiT+JpeakHiBuf;
  binsMT[nbinsPeakT+5] = m_Jpsi+JpeakHiT+JpeakHiBuf+JpeakWidT/2;
  binsMT[nbinsPeakT+6] = m_Jpsi+0.4;

  int nBDTb = BDTcuts.size()-1;
  if(!useBDTbins) nBDTb = 0;
  for(int icut=0;icut<=nBDTb;icut++){
    h_QQM.push_back(vector<TH1F*>());
    h_QQM_tight.push_back(vector<TH1F*>());
    for(int itree=0;itree<ntrees;itree++){
      TString cutstr = "BDT" + ((icut==0)?"all":("bin"+(TString)std::to_string(icut)));
      h_QQM[icut].push_back(new TH1F( "QQ_M_"+treeName[itree]+cutstr, "Jpsi mass", nbinsPeak+6, binsM ));
      h_QQM_tight[icut].push_back(new TH1F( "QQ_M_tight_"+treeName[itree]+cutstr, "Jpsi mass", nbinsPeakT+6, binsMT ));
    }
  }

  vector<TH1F*> JpsiM; //goes in output file
  vector<TH1F*> JpsiM_tight; //goes in output file

  //***************** Extract needed branches
  for(int iT=0; iT<(int)T.size(); iT++){
    T[iT]->SetBranchAddress("muW_eta", &muW_eta[iT]);
    T[iT]->SetBranchAddress("mumi_eta", &mumi_eta[iT]);
    T[iT]->SetBranchAddress("mupl_eta", &mupl_eta[iT]);
    T[iT]->SetBranchAddress("Bc_Pt", &Bc_Pt[iT]);
    T[iT]->SetBranchAddress("Bc_Y", &Bc_Y[iT]);
    T[iT]->SetBranchAddress("QQ_M", &QQ_M[iT]);
    T[iT]->SetBranchAddress("QQ2_M", &QQ2_M[iT]);
    T[iT]->SetBranchAddress("QQ2_dca", &QQ2_dca[iT]);
    T[iT]->SetBranchAddress("QQ2_VtxProb", &QQ2_VtxProb[iT]);
    T[iT]->SetBranchAddress("w_simple", &w_simple[iT]);
    if(!ispp && iT==3) T[iT]->SetBranchAddress("w_unblind", &w_unblind[iT]);
    if(useBDTbins) T[iT]->SetBranchAddress("BDT", &BDT[iT]);

    // T[iT]->SetBranchStatus("*",0); //disable all branches
    // T[iT]->SetBranchStatus("muW_eta",1);
    // T[iT]->SetBranchStatus("mumi_eta",1);
    // T[iT]->SetBranchStatus("mupl_eta",1);
    // T[iT]->SetBranchStatus("QQ_M",1);
    // T[iT]->SetBranchStatus("QQ2_M",1);
    // T[iT]->SetBranchStatus("QQ2_dca",1);
    // T[iT]->SetBranchStatus("QQ2_VtxProb",1);
    // T[iT]->SetBranchStatus("w_simple",1);
  }

  //*******************************************
  //Fill the QQ_M histograms
  for(int iT=1; iT<(int)T.size(); iT++){
    if(iT!=1 && iT!=3) continue;//consider ONLY DATA unambiguous events    //if(iT==0 || iT==7) continue; //forget WRONGSIGN and dimuon+track and MCs
    for(int j=0; j<T[iT]->GetEntries(); j++){//T[iT]->GetEntries()

      T[iT]->GetEntry(j);
      if(!(Bc_Pt[iT]>_BcPtmin[0] && Bc_Pt[iT]<_BcPtmax[0] && fabs(Bc_Y[iT])<_BcYmax[0] && (Bc_Pt[iT]>_BcPtmin[1] || fabs(Bc_Y[iT])>_BcYmin[2]) )) continue; //fiducial cuts

      float maxEta = max(fabs(muW_eta[iT]),max(fabs(mumi_eta[iT]),fabs(mupl_eta[iT])));
      int kbin = nBDTb; //kbin 0 is for all BDT values
      for(int k=0;k<nBDTb;k++){
	if(BDT[iT]>BDTcuts[k] && BDT[iT]<BDTcuts[k+1]) kbin = k+1;
      }

      if(( inJpsiMassRange(QQ_M[iT], maxEta<1.5) && !(inJpsiMassSB(QQ2_M[iT], maxEta<1.5)
						      && QQ2_VtxProb[iT]>_QQvtxProb_cut && QQ2_dca[iT]<_QQdca_cut && QQ2_dca[iT]>0) )
	 || (inJpsiMassSB(QQ_M[iT], maxEta<1.5) && !(inJpsiMassRange(QQ2_M[iT], maxEta<1.5)
						     && QQ2_VtxProb[iT]>_QQvtxProb_cut && QQ2_dca[iT]<_QQdca_cut && QQ2_dca[iT]>0) )
	 ){
	float w = (!ispp && iT==3)?w_unblind[iT]:w_simple[iT];
	if(maxEta<1.5) {
	  float norm = inJpsiMassRange(QQ_M[iT],true)?1:((float)nbinsPeakT/2.);
	  if(useBDTbins) h_QQM_tight[kbin][iT]->Fill(QQ_M[iT], fabs(w)/norm );//not considering the weight of Jpsi SB subtraction
	  h_QQM_tight[0][iT]->Fill(QQ_M[iT], fabs(w_simple[iT])/norm );
	}
	else {float norm = inJpsiMassRange(QQ_M[iT],false)?1:((float)nbinsPeak/2.);
	  if(useBDTbins) h_QQM[kbin][iT]->Fill(QQ_M[iT], fabs(w)/norm );//not considering the weight of Jpsi SB subtraction
	  h_QQM[0][iT]->Fill(QQ_M[iT], fabs(w_simple[iT])/norm );
	} 
      }
    }
  }

  for(int k=0;k<=nBDTb;k++){
    h_QQM[k][3]->Add(h_QQM[k][1]);
    h_QQM_tight[k][3]->Add(h_QQM_tight[k][1]);
  }


  //*******************************************
  //forget histos with too few entries
  for(int k=0;k<=nBDTb;k++){

    //*** loose SB
    if(h_QQM[k][3]->GetEntries() >100){
      JpsiM.push_back((TH1F*)h_QQM[k][3]->Clone("JpsiM_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));
    }
    else{
      if(h_QQM[0][3]->GetEntries() >100){
	//if the histo from this BDT bin is not good enough, take the histo for integrated
	JpsiM.push_back((TH1F*)h_QQM[0][3]->Clone("JpsiM_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));
      } else{ 
	//if nothing is good enough, a bit arbitrary here
	JpsiM.push_back((TH1F*)h_QQM[0][3]->Clone("JpsiM_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));
	for(int bin=1;bin<=JpsiM[k]->GetNbinsX();bin++){
	  if(inJpsiMassRange(JpsiM[k]->GetBinCenter(bin),false)){
	    JpsiM[k]->SetBinContent(bin,0.8+min(k,4)*0.05);}
	  else if(inJpsiMassSB(JpsiM[k]->GetBinCenter(bin),false)){
	    JpsiM[k]->SetBinContent(bin,0.2-min(k,4)*0.05);}
	  else{
	    JpsiM[k]->SetBinContent(bin,0.);}
	}
      }
    }

    //*** tight SB
    if(h_QQM_tight[k][3]->GetEntries() >100){
      JpsiM_tight.push_back((TH1F*)h_QQM_tight[k][3]->Clone("JpsiM_tight_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));
    }
    else if(h_QQM[k][3]->GetEntries() >100){
      //first try if the Jpsi proba from the loose SB histo is ok
      JpsiM_tight.push_back((TH1F*)h_QQM[k][3]->Clone("JpsiM_tight_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));
    } 
    else if(h_QQM_tight[0][3]->GetEntries() >100){//if not, then take the integrated BDT bin histos for tight SB
      JpsiM_tight.push_back((TH1F*)h_QQM_tight[0][3]->Clone("JpsiM_tight_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));	    
    }
    else if(h_QQM[0][3]->GetEntries() >100){//finally, try loose SB integrated BDT
      JpsiM_tight.push_back((TH1F*)h_QQM[0][3]->Clone("JpsiM_tight_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));
    }
    else { 
      //if nothing is good enough, a bit arbitrary here
      JpsiM_tight.push_back((TH1F*)h_QQM_tight[0][3]->Clone("JpsiM_tight_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)) ) ));
      for(int bin=1;bin<=JpsiM_tight[k]->GetNbinsX();bin++){
	if(inJpsiMassRange(JpsiM_tight[k]->GetBinCenter(bin),false)){
	  JpsiM_tight[k]->SetBinContent(bin,0.8+min(k,4)*0.05);}
	else if(inJpsiMassSB(JpsiM_tight[k]->GetBinCenter(bin),false)){
	  JpsiM_tight[k]->SetBinContent(bin,0.2-min(k,4)*0.05);}
	else{
	  JpsiM_tight[k]->SetBinContent(bin,0.);}
      }
    }

  }

  //*******************************************
  //Fill the weight branch, with weights containing the weight for Jpsi choice
  vector<TBranch*> b_weight;

  for(int iT=0; iT<(int)T.size(); iT++){
    std::cout << "--- Processing: " << T[iT]->GetEntries() << " events of tree "<< treeName[iT] << std::endl;

    b_weight.push_back( T[iT]->Branch(((useBDTbins)?"weight":"w_simple2"),&weight[iT],((useBDTbins)?"weight/F":"w_simple2/F") ) );

    //BEGIN event loop on the analyzed tree
    for(int j=0; j<T[iT]->GetEntries(); j++){//T[iT]->GetEntries()

      T[iT]->GetEntry(j);

      weight[iT] = w_simple[iT];
      if(!ispp && iT==3 && BDT[iT]<(_withTM?-0.2:-0.3)) weight[iT] = w_unblind[iT]; //unblind data events in the low-BDT CR

      if(iT!=0 && iT!=7) { //forget WRONGSIGN and dimuon+track 

	float maxEta = max(fabs(muW_eta[iT]),max(fabs(mumi_eta[iT]),fabs(mupl_eta[iT])));
	int kbin = nBDTb; //kbin 0 is for all BDT values
	for(int k=0;k<nBDTb;k++){
	  if(BDT[iT]>BDTcuts[k] && BDT[iT]<BDTcuts[k+1]) kbin = k+1;
	}

	if((inJpsiMassRange(QQ_M[iT], maxEta<1.5) && inJpsiMassSB(QQ2_M[iT], maxEta<1.5)
	    && QQ2_VtxProb[iT]>_QQvtxProb_cut && QQ2_dca[iT]<_QQdca_cut && QQ2_dca[iT]>0)
	   || (inJpsiMassSB(QQ_M[iT], maxEta<1.5) && inJpsiMassRange(QQ2_M[iT], maxEta<1.5)
	       && QQ2_VtxProb[iT]>_QQvtxProb_cut && QQ2_dca[iT]<_QQdca_cut && QQ2_dca[iT]>0)
	   ){
	   
	  //consider ONLY DATA unambiguous events
	  //*** loose SB
	  float binc_QQ1 = JpsiM[kbin]->GetBinContent(JpsiM[kbin]->FindBin(QQ_M[iT]));
	  float binc_QQ2 = JpsiM[kbin]->GetBinContent(JpsiM[kbin]->FindBin(QQ2_M[iT]));

	  //*** tight SB
	  if(maxEta<1.5){
	    binc_QQ1 = JpsiM_tight[kbin]->GetBinContent(JpsiM_tight[kbin]->FindBin(QQ_M[iT]));
	    binc_QQ2 = JpsiM_tight[kbin]->GetBinContent(JpsiM_tight[kbin]->FindBin(QQ2_M[iT]));
	  }

	  if (binc_QQ1==0) weight[iT] = 0;
	  else weight[iT] *= binc_QQ1 / (binc_QQ1+binc_QQ2);
	  //if(iT==1) cout<<"bin content of hist#"<<ihist <<" at value QQM,QQ2M = "<<QQ_M[iT]<<" "<<QQ2_M[iT]<<" = "<<binc_QQ1<<" "<<binc_QQ2<<endl;
	}
      }

      //if(w_simple[iT]!=weight[iT] || iT==1) cout<<w_simple[iT]<<" "<<weight[iT]<<endl;
      b_weight[iT]->Fill();
    }
    //END event loop

    //T[iT]->Print(); 
    T[iT]->Write("",TObject::kOverwrite); //overwrite, or two versions of the trees are saved 
  } 
  //END loop on trees
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  TLegend *leg1 = new TLegend(0.7,0.7,0.9,0.9);
  //  c1->Divide(3,2);
  int ic=1;
  for(int iT=1; iT<(int)T.size(); iT++){
    if(iT!=3) continue;
    //if(iT<2 || iT==7) continue; //forget WRONGSIGN and dimuon+track and MCs
    //c1->cd(ic);
    h_QQM[useBDTbins?1:0][iT]->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    h_QQM[useBDTbins?1:0][iT]->SetTitle("Dimuon mass of unambiguous candidates");
    if(!useBDTbins) h_QQM[0][iT]->Draw();
    else{
      for(int k=1; k<=nBDTb; k++){
	h_QQM[k][iT]->SetLineColor(k);
	h_QQM[k][iT]->Draw((k==1)?"":"same");
	leg1->AddEntry(h_QQM[k][iT],"BDT bin"+(TString)to_string(k));
      }
      leg1->Draw("same");
    }
    ic+=1;
  }
  c1->SaveAs("figs/JpsiMassUnambiguousCands_looseRange_"+(TString)(ispp?"pp":"PbPb")+(TString)(useBDTbins?"_inBDTbins":"")+".pdf");

  TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
  TLegend *leg2 = new TLegend(0.7,0.7,0.9,0.9);
  //c2->Divide(3,2);
  int ic2=1;
  for(int iT=1; iT<(int)T.size(); iT++){
    if(iT!=3) continue;
    //if(iT<2 || iT==7) continue; //forget WRONGSIGN and dimuon+track and MCs
    //c2->cd(ic2);
    h_QQM_tight[useBDTbins?1:0][iT]->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    h_QQM_tight[useBDTbins?1:0][iT]->SetTitle("Dimuon mass of unambiguous candidates");
    if(!useBDTbins) h_QQM_tight[0][iT]->Draw();
    else{
      for(int k=1; k<=nBDTb; k++){
	h_QQM_tight[k][iT]->SetLineColor(k);
	h_QQM_tight[k][iT]->Draw((k==1)?"":"same"); 
	leg2->AddEntry(h_QQM_tight[k][iT],"BDT bin"+(TString)to_string(k));
      }
      leg1->Draw("same");
    }
    ic2+=1;
  }
  c2->SaveAs("figs/JpsiMassUnambiguousCands_tightRange_"+(TString)(ispp?"pp":"PbPb")+(TString)(useBDTbins?"_inBDTbins":"")+".pdf");

  //*******************************************
  //output file for data Jpsi mass distro
  if(useBDTbins){
    TFile *outfile = TFile::Open("JpsiMassDistr.root","UPDATE");
    for(int k=0;k<=nBDTb;k++){
      JpsiM[k]->Write("JpsiMass_data_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)))+(TString)(ispp?"_pp":"_PbPb"));
      JpsiM_tight[k]->Write("JpsiMassCentralEta_data_"+(TString)((k==0)?"allBDTbins":("BDTbin"+(TString)to_string(k)))+(TString)(ispp?"_pp":"_PbPb"));
    }
  }
  
}

void addJpsiChoiceWeight(bool ispp=true, bool useBDTbins=false, vector<float> BDTcuts_ = vector<float>()){
  if (useBDTbins && BDTcuts_.size()<2){
    BDTcuts_ = _BDTcuts(ispp,0); //we don't bother to make this dependent on the kinematic bin
  }
  
  addJpsiChoiceW(ispp,useBDTbins,BDTcuts_);

}
