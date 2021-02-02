#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"
#include "SgMuonAcceptanceCuts.h"
#include "../helpers/AccEff2DBinning.h"

void BuildAcceptanceMap(bool runAEtoys=true, bool withTM = false){

  //**************************************************************  
  //Grab the variations of the pT bias of MC, from first step r1 and r2
  vector<vector<TH1F*> > bias = vector<vector<TH1F*> >(2);
  if(runAEtoys){
    TFile *BiasFile = TFile::Open("../twoSteps/pTBiases.root","READ");
    for(int col=0;col<2;col++){
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	bias[col].push_back((TH1F*)BiasFile->Get("pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)));
	//cout<<"col variation value_at_11 = "<<col<<" "<<v<<" "<<bias[col][v]->Eval(11)<<endl;
      }
    }
  }

  //**************************************************************  
  //Create Tree 
  TFile *file = TFile::Open("/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/acceptance/BcToJpsiMuNu_BCVEGPY_PYTHIA8_pp5TeV_RunIIpp5Spring18DR-00093_acceptance_14092020_ONIATREE_3640k.root");///data_CMS/cms/falmagne/tuples/Bc_GenOnly/BcToJpsiMuNu_BCVEGPY_PYTHIA8_GenOnly_11052020_ONIATREE.root
  TTree* T = (TTree*)file->Get("hionia/myTree");
  int nentries = T->GetEntries();
  std::cout<<"nevents = "<<nentries<<"\n";
  //  Trees[0]->Print();

  //**************************************************************  
  //Get input branches

  UInt_t eventNb;
  T->SetBranchAddress("eventNb", &eventNb);
  Short_t Gen_Bc_size;
  T->SetBranchAddress("Gen_Bc_size", &Gen_Bc_size);
  Short_t Gen_mu_size;
  T->SetBranchAddress("Gen_mu_size", &Gen_mu_size);
  Short_t Gen_Bc_muW_idx[10];
  T->SetBranchAddress("Gen_Bc_muW_idx", &Gen_Bc_muW_idx);
  Short_t Gen_Bc_QQ_idx[10];
  T->SetBranchAddress("Gen_Bc_QQ_idx", &Gen_Bc_QQ_idx);
  Short_t Gen_QQ_mumi_idx[10];
  T->SetBranchAddress("Gen_QQ_mumi_idx", &Gen_QQ_mumi_idx);
  Short_t Gen_QQ_mupl_idx[10];
  T->SetBranchAddress("Gen_QQ_mupl_idx", &Gen_QQ_mupl_idx);

  TClonesArray *Gen_3mu_4mom = new TClonesArray(); TBranch *b_Gen_3mu_4mom;
  T->SetBranchAddress("Gen_3mu_4mom", &Gen_3mu_4mom, &b_Gen_3mu_4mom);
  TClonesArray *Gen_Bc_4mom = new TClonesArray(); TBranch *b_Gen_Bc_4mom;
  T->SetBranchAddress("Gen_Bc_4mom", &Gen_Bc_4mom, &b_Gen_Bc_4mom);
  TClonesArray *Gen_QQ_4mom = new TClonesArray(); TBranch *b_Gen_QQ_4mom;
  T->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  TClonesArray *Gen_mu_4mom = new TClonesArray(); TBranch *b_Gen_mu_4mom;
  T->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  //**************************************************************
  //histos
  const int nPtBin = 21, nYBin = 17;
  float PtBins[nPtBin+1], absYBins[nYBin+1];
  for(int i=0;i<14;i++) PtBins[i] = (float)i;
  for(int i=14;i<=17;i++) PtBins[i] = 14+ 2*(i-14);
  PtBins[18] = 24; PtBins[19] = 29; PtBins[20] = 35; PtBins[21] = 50;
  //PtBins[23] = 50;
  //  PtBins[24] = 70;
  for(int i=0;i<=13;i++) absYBins[i] = i*0.15;
  absYBins[14] = 2.05; absYBins[15] = 2.2; absYBins[16] = 2.3; absYBins[17] = 2.4;

  TH2F* h_all = new TH2F("h_all","h_all",nYBin,absYBins,nPtBin,PtBins);
  TH2F* h_fid = new TH2F("h_fid","h_fid",nYBin,absYBins,nPtBin,PtBins); //Bc's in the fiducial cuts
  TH2F* h_acc = new TH2F("h_acc","h_acc",nYBin,absYBins,nPtBin,PtBins); //Bc's passing fiducial cuts and single muon acceptances
  TH1F* hPt_all = new TH1F("hPt_all","hPt_all",50,0,50);
  TH1F* hPt_fid = new TH1F("hPt_fid","hPt_fid",50,0,50);
  TH1F* hPt_acc = new TH1F("hPt_acc","hPt_acc",50,0,50);
  TH1F* hY_all = new TH1F("hY_all","hY_all",30,0,2.4);
  TH1F* hY_fid = new TH1F("hY_fid","hY_fid",30,0,2.4);
  TH1F* hY_acc = new TH1F("hY_acc","hY_acc",30,0,2.4);
  
  TH2Poly *hp =_hp();
  TH2Poly *hp_all = (TH2Poly*) hp->Clone("hp_all");
  TH2Poly *hp_acc = (TH2Poly*) hp->Clone("hp_acc");
  vector<float> gen_oneBinned(_NanaBins+1,0);
  vector<float> accepted_oneBinned(_NanaBins+1,0);
  vector<float> acc_oneBinned(_NanaBins+1,1);
  vector<vector<vector<float> > > gen_oneB_biased(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) ));
  vector<vector<vector<float> > > accepted_oneB_biased(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) ));
  vector<vector<vector<float> > > acc_oneB_biased(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 1) ));

  float nall=0;
  //**************************************************************
  //event loop
  for(int j=0; j<nentries; j++){
    Gen_3mu_4mom->Clear();
    Gen_Bc_4mom->Clear();
    Gen_QQ_4mom->Clear();
    Gen_mu_4mom->Clear();

    if(j%100000==0){ cout<<"Scanned "<<100.*(double)j/nentries<<"% of tree"<<endl; }

    T->GetEntry(j);

    for(int iBc=0;iBc<Gen_Bc_size;iBc++){
      TLorentzVector *genBc = (TLorentzVector*) Gen_Bc_4mom->At(iBc);
      TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom->At(iBc);

      int QQidx = Gen_Bc_QQ_idx[iBc];
      TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom->At(QQidx);
      TLorentzVector *genBc_muW = (TLorentzVector*) Gen_mu_4mom->At(Gen_Bc_muW_idx[iBc]);
      TLorentzVector *genBc_mumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[QQidx]);
      TLorentzVector *genBc_mupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[QQidx]);

      h_all->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt());
      hp_all->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt());
      hPt_all->Fill(gen3mu->Pt());
      hY_all->Fill(fabs(gen3mu->Rapidity()));
      nall +=1;

      if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM)){ 
	h_acc->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt());
	hp_acc->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt());
      }

      if(fabs(gen3mu->Rapidity())<_BcYmax[0] ){ //Y cut for pt distro
	hPt_fid->Fill(gen3mu->Pt());
	if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM))
	  hPt_acc->Fill(gen3mu->Pt());
      }
      if(gen3mu->Pt()>_BcPtmin[0]){ //pt cut for Y distro
	hY_fid->Fill(fabs(gen3mu->Rapidity()));	
	if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM))
	  hY_acc->Fill(fabs(gen3mu->Rapidity()));
      }

      for(int b=1;b<=_NanaBins;b++){
	if(fabs(gen3mu->Rapidity())>_BcYmin[b] && fabs(gen3mu->Rapidity())<_BcYmax[b] && gen3mu->Pt()>_BcPtmin[b] && gen3mu->Pt()<_BcPtmax[b] ){
	  h_fid->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt());
	  gen_oneBinned[0] += 1;
	  gen_oneBinned[b] += 1;

	  //biased MC
	  if(runAEtoys){
	    for(int col=0;col<2;col++){
	      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
		gen_oneB_biased[col][0][v] += getBias( bias[col][v] , gen3mu->Pt());
		gen_oneB_biased[col][b][v] += getBias( bias[col][v] , gen3mu->Pt());
	      }
	    }
	  }

	  if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM)){
	    accepted_oneBinned[0] += 1;
	    accepted_oneBinned[b] += 1;	  

	    //biased MC
	    if(runAEtoys){
	      for(int col=0;col<2;col++){
		for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
		  accepted_oneB_biased[col][0][v] += getBias( bias[col][v] , gen3mu->Pt());
		  accepted_oneB_biased[col][b][v] += getBias( bias[col][v] , gen3mu->Pt());
		}
	      }
	    }
	  }//end is accepted

	}
      }//end loop on bins

    }//end Bc loop
  } //end event loop

  for(int b=0;b<=_NanaBins;b++){
    acc_oneBinned[b] = accepted_oneBinned[b]/gen_oneBinned[b];
    if(runAEtoys){
      for(int col=0;col<2;col++){
	for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	  acc_oneB_biased[col][b][v] = accepted_oneB_biased[col][b][v]/gen_oneB_biased[col][b][v];	
	}
      }
    }
  }

  //**************************************************************
  //Lines for fiducial cuts
  TLine *line1 = new TLine(0,11,1.3,11);
  TLine *line2 = new TLine(1.3,6,1.3,11);
  TLine *line3 = new TLine(1.3,6,2.3,6);
  TLine *line4 = new TLine(2.3,6,2.3,50);
  line1->SetLineWidth(4);  line1->SetLineColor(kBlack);
  line2->SetLineWidth(4);  line2->SetLineColor(kBlack);
  line3->SetLineWidth(4);  line3->SetLineColor(kBlack);
  line4->SetLineWidth(4);  line4->SetLineColor(kBlack);

  //**************************************************************
  //Draw histos
  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(50);

  TCanvas *c1 = new TCanvas("c1","c1",3900,1300);
  c1->Divide(3,1);

  TH1F* h_acceptance = (TH1F*)h_acc->Clone("h_acceptance"); 
  c1->cd(1);
  gPad->SetLogz();
  h_acc->GetZaxis()->SetRangeUser(1,100);
  h_acc->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  h_acc->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  h_acc->SetTitle("Accepted B_{c}'s");
  h_acc->Draw("COLZ");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  // c1->cd(2);
  // h_acceptance->Divide(h_all); //h_fid
  // h_acceptance->GetZaxis()->SetRangeUser(0,1);
  // h_acceptance->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  // h_acceptance->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  // h_acceptance->SetTitle("Acceptance");
  // h_acceptance->Draw("COLZ");
  // line1->Draw("same");
  // line2->Draw("same");
  // line3->Draw("same");
  // line4->Draw("same");

  c1->cd(2);
  hPt_acc->Divide(hPt_fid);
  hPt_acc->GetXaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  hPt_acc->SetTitle("Acceptance (|Y^{vis}(B_{c})|<2.3)");
  hPt_acc->SetLineWidth(2);
  hPt_acc->Draw();

  c1->cd(3);
  hY_acc->Divide(hY_fid);
  hY_acc->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  hY_acc->SetTitle("Acceptance (p_{T}^{vis}(B_{c})>4GeV)");
  hY_acc->SetLineWidth(2);
  hY_acc->Draw();

  c1->SaveAs("figs/AcceptanceMap_regularBins"+(TString)(withTM?"_withTrackerMu":"")+".pdf");
  c1->SaveAs("figs/AcceptanceMap_regularBins"+(TString)(withTM?"_withTrackerMu":"")+".png");

  //**************************************************************
  //Draw TH2Poly

  TCanvas *c2 = new TCanvas("c2","c2",3000,1500);
  c2->Divide(2,1);

  TH2Poly* hp_acceptance = (TH2Poly*)hp_acc->Clone("hp_acceptance"); 
  c2->cd(1);
  gPad->SetLogz();
  hp_acc->GetZaxis()->SetRangeUser(1,220);
  hp_acc->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  hp_acc->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  hp_acc->SetTitle("Accepted B_{c}'s");
  hp_acc->Draw("COLZ0");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  c2->cd(2);
  gPad->SetLogz();
  hp_acceptance->Divide(hp_all);
  hp_acceptance->GetZaxis()->SetRangeUser(1e-4,1);
  hp_acceptance->GetXaxis()->SetTitle("|Y^{vis}(B_{c})|");
  hp_acceptance->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
  hp_acceptance->SetTitle("Acceptance");
  hp_acceptance->Draw("COLZ0");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  c2->SaveAs("figs/AcceptanceMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+".pdf");
  c2->SaveAs("figs/AcceptanceMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+".png");

  //**************************************************************
  //output file
  TFile out_file("acceptanceMap.root","RECREATE");
  hp_acceptance->Write();
  out_file.WriteObject(&acc_oneBinned,"acceptance_oneBinned");
  if(runAEtoys) out_file.WriteObject(&acc_oneB_biased,"acceptance_oneBinned_biased");

  out_file.Close();
  
  cout<<"nall, n_fiducial, n_accepted, acceptance = "<<nall<<" "<<gen_oneBinned[0]<<" "<<accepted_oneBinned[0]<<" "<<accepted_oneBinned[0]/gen_oneBinned[0]<<endl;

}
