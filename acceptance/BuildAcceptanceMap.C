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
#include "../helpers/SgMuonAcceptanceCuts.h"
#include "../helpers/AccEff2DBinning.h"

void BuildAcceptanceMap(bool runAEtoys=true, bool secondStep=false, bool runMCclos=false, bool withTM = false){

  if(runAEtoys && !secondStep) cout<<"!!!!! WARNING, you should set runAEtoys=true only with secondStep !"<<endl;

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  //**************************************************************  
  //Grab the variations of the pT bias of MC, from first step r1 and r2
  vector<vector<TH1F*> > bias = vector<vector<TH1F*> >(2);
  vector<vector<TH1F*> > bias_2ndStep = vector<vector<TH1F*> >(2);
  TFile *BiasFile = TFile::Open("../twoSteps/pTBiases.root","READ");
  if(secondStep || runAEtoys){
    for(int col=0;col<2;col++){
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	bias[col].push_back((TH1F*)BiasFile->Get("pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)));
	if(runAEtoys)
	  bias_2ndStep[col].push_back((TH1F*)BiasFile->Get("pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)+(TString)(secondStep?"_2ndStep":"")));
      }
    }
  }

  vector<TH1F*> biasMCclos = vector<TH1F*>();
  vector<TH1F*> biasMCclos_2ndStep = vector<TH1F*>();
  if(secondStep || runAEtoys){
    for(int t=0;t<_nMCclos;t++){
      biasMCclos.push_back(runMCclos?((TH1F*)BiasFile->Get("pTbias_PbPb_MCclosure_toy"+(TString)to_string(t))):NULL);
      if(runAEtoys)
	biasMCclos_2ndStep.push_back((secondStep && runMCclos)?((TH1F*)BiasFile->Get("pTbias_PbPb_MCclosure_toy"+(TString)to_string(t)+"_2ndStep")):NULL);
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
  const int nPtBin = 20, nYBin = 17;
  float PtBins[nPtBin+1], absYBins[nYBin+1];
  for(int i=0;i<14;i++) PtBins[i] = (float)i;
  for(int i=14;i<=17;i++) PtBins[i] = 14+ 2*(i-14);
  PtBins[18] = 24; PtBins[19] = 29; PtBins[20] = 35;
  //PtBins[21] = 50;
  //  PtBins[22] = 70;
  for(int i=0;i<=13;i++) absYBins[i] = i*0.15;
  absYBins[14] = 2.05; absYBins[15] = 2.2; absYBins[16] = 2.3; absYBins[17] = 2.4;

  vector<TH2F*> h_all,h_fid,h_acc;
  vector<TH1F*> hPt_all,hPt_fid,hPt_acc,hY_all,hY_fid,hY_acc;
  for(int col=0;col<2;col++){
    h_all.push_back(new TH2F("h_all"+(TString)to_string(col),"h_all",nYBin,absYBins,nPtBin,PtBins));
    h_fid.push_back(new TH2F("h_fid"+(TString)to_string(col),"h_fid",nYBin,absYBins,nPtBin,PtBins)); //Bc's in the fiducial cuts
    h_acc.push_back(new TH2F("h_acc"+(TString)to_string(col),"h_acc",nYBin,absYBins,nPtBin,PtBins)); //Bc's passing fiducial cuts and single muon acceptances
    hPt_all.push_back(new TH1F("hPt_all"+(TString)to_string(col),"hPt_all",35,0,35));
    hPt_fid.push_back(new TH1F("hPt_fid"+(TString)to_string(col),"hPt_fid",35,0,35));
    hPt_acc.push_back(new TH1F("hPt_acc"+(TString)to_string(col),"hPt_acc",35,0,35));
    hY_all.push_back(new TH1F("hY_all"+(TString)to_string(col),"hY_all",30,0,2.4));
    hY_fid.push_back(new TH1F("hY_fid"+(TString)to_string(col),"hY_fid",30,0,2.4));
    hY_acc.push_back(new TH1F("hY_acc"+(TString)to_string(col),"hY_acc",30,0,2.4));
  }
  
  TH2Poly* hp = _hp();
  vector<TH2Poly*> hp_all, hp_acc;
  for(int col=0;col<2;col++){
    hp_all.push_back((TH2Poly*) hp->Clone("hp_all"+(TString)to_string(col)));
    hp_acc.push_back((TH2Poly*) hp->Clone("hp_acc"+(TString)to_string(col)));
  }
  vector<vector<float> > gen_oneBinned(2, vector<float>(_NanaBins+1,0));
  vector<vector<float> > accepted_oneBinned(2, vector<float>(_NanaBins+1,0));
  vector<vector<float> > acc_oneBinned(2, vector<float>(_NanaBins+1,1));
  vector<vector<vector<float> > > gen_oneB_biased(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) ));
  vector<vector<vector<float> > > accepted_oneB_biased(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) ));
  vector<vector<vector<float> > > acc_oneB_biased(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 1) ));
  vector<vector<float> > gen_MCclos(_nMCclos, vector<float>(_NanaBins+1,0));
  vector<vector<float> > accepted_MCclos(_nMCclos, vector<float>(_NanaBins+1,0));
  vector<vector<float> > acc_MCclos(_nMCclos, vector<float>(_NanaBins+1,1));

  float baseW[] = {_scaleMCsigAcc[true] , _scaleMCsigAcc[false]};

  vector<float> nall(2,0);
  //**************************************************************
  //event loop
  for(int j=0; j<nentries; j++){
    Gen_3mu_4mom->Clear();
    Gen_Bc_4mom->Clear();
    Gen_QQ_4mom->Clear();
    Gen_mu_4mom->Clear();

    if(j%100000==0){ cout<<"Scanned "<<100.*(double)j/nentries<<"% of tree"<<endl; }

    T->GetEntry(j);
    //    cout<<"j="<<j<<endl;
    for(int iBc=0;iBc<Gen_Bc_size;iBc++){
      TLorentzVector *genBc = (TLorentzVector*) Gen_Bc_4mom->At(iBc);
      TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom->At(iBc);

      int QQidx = Gen_Bc_QQ_idx[iBc];
      TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom->At(QQidx);
      TLorentzVector *genBc_muW = (TLorentzVector*) Gen_mu_4mom->At(Gen_Bc_muW_idx[iBc]);
      TLorentzVector *genBc_mumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[QQidx]);
      TLorentzVector *genBc_mupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[QQidx]);
      
      for(int col=0;col<2;col++){
	//if(col>0 && !secondStep && !runAEtoys) continue;
	float w = baseW[col];
	w *= (secondStep?( getBias(bias[col][_nomMethVar],gen3mu->Pt()) ):1);
	float wadd = (secondStep && runAEtoys)?( getBias(bias_2ndStep[col][_nomMethVar],gen3mu->Pt()) ):1;
	w *= wadd;

	//for MC closure test
	vector<float> wt = vector<float>(_nMCclos, baseW[col]);
	if(runMCclos){
	  for(int t=0;t<_nMCclos;t++){
	    if(!secondStep && !runAEtoys) wt[t] *= MCclosurePTw(gen3mu->Pt(),t);
	    wt[t] *= secondStep?( getBias(biasMCclos[t],gen3mu->Pt()) ):1;
	    wt[t] *= (secondStep && runAEtoys)?( getBias(biasMCclos_2ndStep[t],gen3mu->Pt()) ):1;
	  }
	}

	h_all[col]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt() , w);
	hp_all[col]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt() , w);
	hPt_all[col]->Fill(gen3mu->Pt() , w);
	hY_all[col]->Fill(fabs(gen3mu->Rapidity()) , w);
	nall[col] +=1;

	if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM)){ 
	  h_acc[col]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt() , w);
	  hp_acc[col]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt() , w);
	}

	if(fabs(gen3mu->Rapidity())<_BcYmax[0] ){ //Y cut for pt distro
	  hPt_fid[col]->Fill(gen3mu->Pt() , w);
	  if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM))
	    hPt_acc[col]->Fill(gen3mu->Pt() , w);
	}
	if(gen3mu->Pt()>_BcPtmin[0]){ //pt cut for Y distro
	  hY_fid[col]->Fill(fabs(gen3mu->Rapidity()) , w);
	  if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM))
	    hY_acc[col]->Fill(fabs(gen3mu->Rapidity()) , w);
	}

	for(int b=1;b<=_NanaBins;b++){
	  if(inFidCuts(b,gen3mu->Pt(),gen3mu->Rapidity())){

	    h_fid[col]->Fill(fabs(gen3mu->Rapidity()),gen3mu->Pt() , w);
	    gen_oneBinned[col][0] += w;
	    gen_oneBinned[col][b] += w;

	    //MC closure
	    if(runMCclos && col==1){
	      for(int t=0;t<_nMCclos;t++){
		gen_MCclos[t][0] += wt[t];
		gen_MCclos[t][b] += wt[t];
	      }
	    }

	    //biased MC
	    if(runAEtoys && !runMCclos){
	      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
		gen_oneB_biased[col][0][v] += w * getBias( (secondStep?bias_2ndStep:bias)[col][v] , gen3mu->Pt()) / wadd;
		gen_oneB_biased[col][b][v] += w * getBias( (secondStep?bias_2ndStep:bias)[col][v] , gen3mu->Pt()) / wadd;
	      }
	    }

	    if(InAcc(*genBc_muW,*genBc_mumi,*genBc_mupl,withTM)){
	      accepted_oneBinned[col][0] += w;
	      accepted_oneBinned[col][b] += w;	  

	      //MC closure
	      if(runMCclos && col==1){
		for(int t=0;t<_nMCclos;t++){
		  accepted_MCclos[t][0] += wt[t];
		  accepted_MCclos[t][b] += wt[t];
		}
	      }

	      //biased MC
	      if(runAEtoys && !runMCclos){
		for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
		  accepted_oneB_biased[col][0][v] += w * getBias( (secondStep?bias_2ndStep:bias)[col][v] , gen3mu->Pt()) / wadd;
		  accepted_oneB_biased[col][b][v] += w * getBias( (secondStep?bias_2ndStep:bias)[col][v] , gen3mu->Pt()) / wadd;
		}
	      }
	    }//end is accepted

	  }
	}//end loop on bins
      }

    }//end Bc loop
  } //end event loop

  for(int col=0;col<2;col++){
    for(int b=0;b<=_NanaBins;b++){
      acc_oneBinned[col][b] = accepted_oneBinned[col][b]/gen_oneBinned[col][b];
      if(runAEtoys && !runMCclos){
	for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	  acc_oneB_biased[col][b][v] = accepted_oneB_biased[col][b][v]/gen_oneB_biased[col][b][v];	
	}

      }

      //MC closure
      if(runMCclos && col==1){
	for(int t=0;t<_nMCclos;t++){
	  acc_MCclos[t][b] = accepted_MCclos[t][b]/gen_MCclos[t][b]; //record acc_MCclos and gen_MCclos(=Ncorr)
	}
      }
      //if(!secondStep) acc_oneBinned[1][b] = acc_oneBinned[0][b];
    }
  }

  vector<TH1F*> h_acceptance;
  vector<TH2Poly*> hp_acceptance;
  for(int col=0;col<2;col++){
    //if(col>0 && !secondStep) continue;

    //**************************************************************
    //Lines for fiducial cuts
    TLine *line1 = new TLine(0,11,1.3,11);
    TLine *line2 = new TLine(1.3,6,1.3,11);
    TLine *line3 = new TLine(1.3,6,2.3,6);
    TLine *line4 = new TLine(2.3,6,2.3,35);
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

    h_acceptance.push_back( (TH1F*)h_acc[col]->Clone("h_acceptance"+(TString)to_string(col)) ); 
    c1->cd(1);
    gPad->SetLogz();
    h_acc[col]->GetZaxis()->SetRangeUser(1,100);
    h_acc[col]->GetXaxis()->SetTitle("|y^{vis}(B_{c})|");
    h_acc[col]->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
    h_acc[col]->SetTitle("Accepted B_{c}'s");
    h_acc[col]->Draw("COLZ");
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");

    // c1->cd(2);
    // h_acceptance->Divide(h_all); //h_fid
    // h_acceptance->GetZaxis()->SetRangeUser(0,1);
    // h_acceptance->GetXaxis()->SetTitle("|y^{vis}(B_{c})|");
    // h_acceptance->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
    // h_acceptance->SetTitle("Acceptance");
    // h_acceptance->Draw("COLZ");
    // line1->Draw("same");
    // line2->Draw("same");
    // line3->Draw("same");
    // line4->Draw("same");

    c1->cd(2);
    hPt_acc[col]->Divide(hPt_fid[col]);
    hPt_acc[col]->GetXaxis()->SetTitle("p_{T}^{vis}(B_{c})");
    hPt_acc[col]->SetTitle("Acceptance (|y^{vis}(B_{c})|<2.3)");
    hPt_acc[col]->SetLineWidth(2);
    hPt_acc[col]->Draw();

    c1->cd(3);
    hY_acc[col]->Divide(hY_fid[col]);
    hY_acc[col]->GetXaxis()->SetTitle("|y^{vis}(B_{c})|");
    hY_acc[col]->SetTitle("Acceptance (p_{T}^{vis}(B_{c})>4GeV)");
    hY_acc[col]->SetLineWidth(2);
    hY_acc[col]->Draw();
    
    if(!runMCclos){
      c1->SaveAs("figs/AcceptanceMap_regularBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep_":"_2ndStep_")+(TString)((col==0)?"pp":"PbPb")):"")+".pdf");
      c1->SaveAs("figs/AcceptanceMap_regularBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep_":"_2ndStep_")+(TString)((col==0)?"pp":"PbPb")):"")+".png");
    }

    //**************************************************************
    //Draw TH2Poly

    TCanvas *c2 = new TCanvas("c2","c2",3000,1500);
    c2->Divide(2,1);

    hp_acceptance.push_back( (TH2Poly*)hp_acc[col]->Clone("hp_acceptance"+(TString)to_string(col)) ); 
    c2->cd(1);
    gPad->SetLogz();
    hp_acc[col]->GetZaxis()->SetRangeUser(1,220);
    hp_acc[col]->GetXaxis()->SetTitle("|y^{vis}(B_{c})|");
    hp_acc[col]->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
    hp_acc[col]->SetTitle("Accepted B_{c}'s");
    hp_acc[col]->Draw("COLZ0");
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");

    c2->cd(2);
    gPad->SetLogz();
    hp_acceptance[col]->Divide(hp_all[col]);
    hp_acceptance[col]->GetZaxis()->SetRangeUser(1e-4,1);
    hp_acceptance[col]->GetXaxis()->SetTitle("|y^{vis}(B_{c})|");
    hp_acceptance[col]->GetYaxis()->SetTitle("p_{T}^{vis}(B_{c})");
    hp_acceptance[col]->SetTitle("Acceptance");
    hp_acceptance[col]->Draw("COLZ0");
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");

    if(!runMCclos){
      c2->SaveAs("figs/AcceptanceMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep_":"_2ndStep_")+(TString)((col==0)?"pp":"PbPb")):"")+".pdf");
      c2->SaveAs("figs/AcceptanceMap_tunedBins"+(TString)(withTM?"_withTrackerMu":"")+(TString)(secondStep?((TString)(runAEtoys?"_3rdStep_":"_2ndStep_")+(TString)((col==0)?"pp":"PbPb")):"")+".png");
    }

    cout<<"pp/PbPb, bin, nall, n_fiducial, n_accepted, acceptance = "<<col<<" all "<<nall[col]<<" "<<gen_oneBinned[col][0]<<" "<<accepted_oneBinned[col][0]<<" "<<accepted_oneBinned[col][0]/gen_oneBinned[col][0]<<endl;
    cout<<"pp/PbPb, bin, nall, n_fiducial, n_accepted, acceptance = "<<col<<" 1 "<<nall[col]<<" "<<gen_oneBinned[col][1]<<" "<<accepted_oneBinned[col][1]<<" "<<accepted_oneBinned[col][1]/gen_oneBinned[col][1]<<endl;
    cout<<"pp/PbPb, bin, nall, n_fiducial, n_accepted, acceptance = "<<col<<" 2 "<<nall[col]<<" "<<gen_oneBinned[col][2]<<" "<<accepted_oneBinned[col][2]<<" "<<accepted_oneBinned[col][2]/gen_oneBinned[col][2]<<endl;
  }
  if(!secondStep) hp_acceptance.push_back(hp_acceptance[0]);

  //**************************************************************
  //output file
  TFile out_file("acceptanceMap.root","UPDATE");
  if(!runMCclos){
    hp_acceptance[0]->Write("hp_acceptance_pp"+(TString)(secondStep?(TString)(runAEtoys?"_3rdStep":"_2ndStep"):""));
    hp_acceptance[1]->Write("hp_acceptance_PbPb"+(TString)(secondStep?(TString)(runAEtoys?"_3rdStep":"_2ndStep"):""));
    out_file.WriteObject(&acc_oneBinned,"acceptance_oneBinned"+(TString)(secondStep?(TString)(runAEtoys?"_3rdStep":"_2ndStep"):""));
    if(runAEtoys) out_file.WriteObject(&acc_oneB_biased,"acceptance_oneBinned_biased"+(TString)(secondStep?"_2ndStep":""));
  }else{
    out_file.WriteObject(&acc_MCclos,"acceptance_MCclosure"+(TString)(secondStep?(TString)(runAEtoys?"_3rdStep":"_2ndStep"):""));
    out_file.WriteObject(&gen_MCclos,"NcorrGen_MCclosure"+(TString)(secondStep?(TString)(runAEtoys?"_3rdStep":"_2ndStep"):""));
  }
  out_file.Close();
  

}
