#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "../PbPb18/Utilities/EVENTUTILS.h"

Double_t fctexpfitratio(Double_t *x, Double_t *par)
{
  float pt =x[0];
  float tau = par[0];
  float ptlim = par[1];
  float C = par[2];
  float pt0 = TMath::Log(C)+ ptlim/tau;

  if(pt>ptlim) return C;
  else return TMath::Exp(pt0 - pt/tau);
  //SetParameter(3.8,11,0.15)
}

Double_t fctfitratio(Double_t *x, Double_t *par)
{
  float pt = x[0];
  float n = par[0];
  float ptlim = par[1];
  float C = par[2];
  float A = C * pow(ptlim, n); //guarantees continuity

  if(pt>ptlim) return C;
  else return A*pow(pt, -n);
}

void GetNormalization(bool ispp = true, bool isPrompt = false){

  auto h_test = new TH1D();
  h_test->SetDefaultSumw2(true);

  //****** Open the tree and make it scan branches one by one (SetMakeClass, to study one branch at a time) 
  TFile *fileMC = TFile::Open(ispp?(isPrompt?
			      "/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/PromptJpsi/Oniatree_MC_trimuons_PromptJpsi_ptHatMinCombined_05082019.root"
				    :"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/NonPromptJpsi/MConiatree/crab_BJPsiMM_TuneCUETP8M1_5p02TeV_pythia8_05082019_wLambdabFor10_ptHatMinCombined_ONIATREE.root"):
			      (isPrompt?"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/PromptJpsi/Jpsi_pThat-2_TuneCP5_HydjetDrumMB_HINPbPbAutumn18DR_trimuons_oniatree_25012021.root" //Jpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_trimuons_oniatree_09012020.root
			       :"/data_CMS/cms/falmagne/tuples/PbPb18/Bc/TripleMu/MC/NonPromptJpsi/BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_trimuons_oniatree_09012020.root"),"READ");
  TTree *TMC = (TTree*)fileMC->Get("hionia/myTree");
  int neventsMC = (int)TMC->GetEntries();
  std::cout<<"nevents MC = "<<neventsMC<<"\n";
  TMC->SetMakeClass(1);

  //****** Gen info 
  Short_t genQQ_size;
  TBranch *b_genQQ_size = TMC->GetBranch("Gen_QQ_size");
  b_genQQ_size->SetAddress(&genQQ_size);

  Short_t genmu_size;
  TBranch *b_genmu_size = TMC->GetBranch("Gen_mu_size");
  b_genmu_size->SetAddress(&genmu_size);

  // int Reco_QQ_whichGen[100];
  // TBranch *b_Reco_QQ_whichGen = TMC->GetBranch("Reco_QQ_whichGen");
  // b_Reco_QQ_whichGen->SetAddress(&Reco_QQ_whichGen);

  float pthatweight;
  TBranch *b_pthatweight = TMC->GetBranch("pthatweight");
  if(ispp)
    b_pthatweight->SetAddress(&pthatweight);
  
  float Gen_weight;
  TBranch *b_Gen_weight = TMC->GetBranch("Gen_weight");
  if(!ispp)
    b_Gen_weight->SetAddress(&Gen_weight);

  int Centrality;
  TBranch *b_Centrality = TMC->GetBranch("Centrality");
  if(!ispp)
    b_Centrality->SetAddress(&Centrality);

  TClonesArray *Gen_QQ_4mom = new TClonesArray();
  TBranch *b_Gen_QQ_4mom; 
  TMC->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);

  TClonesArray *Gen_mu_4mom = new TClonesArray();
  TBranch *b_Gen_mu_4mom; 
  TMC->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  Short_t Gen_QQ_mumi_idx[10];
  TBranch *b_Gen_QQ_mumi_idx = TMC->GetBranch("Gen_QQ_mumi_idx");
  b_Gen_QQ_mumi_idx->SetAddress(&Gen_QQ_mumi_idx);

  Short_t Gen_QQ_mupl_idx[10];
  TBranch *b_Gen_QQ_mupl_idx = TMC->GetBranch("Gen_QQ_mupl_idx");
  b_Gen_QQ_mupl_idx->SetAddress(&Gen_QQ_mupl_idx);

  float NcollW;

  //***** Some histograms & get the data XS
  TFile *fileXS = TFile::Open(isPrompt?"HEPData-ins1644903-v1-Table_3.root":"HEPData-ins1644903-v1-Table_4.root");
  TGraphErrors *dataXS = (TGraphErrors*) fileXS->Get((TString)(isPrompt?"Table 3":"Table 4")+(TString)(ispp?"/Graph1D_y1":"/Graph1D_y2"));

  int nbins = dataXS->GetN();
  double *Xcenter = dataXS->GetX();
  double Xerror[nbins];
  for(int bin=0;bin<nbins;bin++){
    Xerror[bin] = dataXS->GetErrorX(bin);}
  double binLimits[nbins+1];
  
  binLimits[0] = Xcenter[0] - Xerror[0];
  for(int bin=1;bin<nbins+1;bin++){
    binLimits[bin] = Xcenter[bin-1] + Xerror[bin-1];
  }
  TH1F *h_QQ_Pt = new TH1F("h_QQ_Pt",((TString)(isPrompt?"P":"Non-p")+"rompt J/#psi MC cross-section, "+(ispp?"pp":"PbPb")),nbins,binLimits);//nbins,binLimits

  //***************** Loop on events (MC, background)
  for(int i=0;i<neventsMC; i++){//neventsMC

    if(i%500000==0){ std::cout<<"Looping over MC events. Part of tree scanned : "<< 100*(double)i / (double)neventsMC<<" % \n";}
    
    //***** Initialize pointers for arrays
    Gen_QQ_4mom->Clear();
    Gen_mu_4mom->Clear();

    //***** Choose needed branches
    b_genQQ_size->GetEntry(i);
    //    b_genmu_size->GetEntry(i);
    b_Gen_QQ_4mom->GetEntry(i);
    b_Gen_mu_4mom->GetEntry(i);
    b_Gen_QQ_mumi_idx->GetEntry(i);
    b_Gen_QQ_mupl_idx->GetEntry(i);
    if(ispp){
      b_pthatweight->GetEntry(i);}
    else{
      b_Gen_weight->GetEntry(i);
      b_Centrality->GetEntry(i);
    }

    NcollW = ispp?1:findNcoll((Centrality<0 || Centrality>199)?0:Centrality);
    for (int QQNb=0; QQNb<(int)genQQ_size; QQNb++){
      TLorentzVector *QQ = (TLorentzVector*) Gen_QQ_4mom->At(QQNb);    
      //      if(!(Gen_QQ_mumi_idx[QQNb]>-1 && Gen_QQ_mumi_idx[QQNb]<(int)genmu_size && Gen_QQ_mupl_idx[QQNb]>-1 && Gen_QQ_mupl_idx[QQNb]<(int)genmu_size)) continue;
      TLorentzVector *mumi = (TLorentzVector*) Gen_mu_4mom->At((int)Gen_QQ_mumi_idx[QQNb]);    
      TLorentzVector *mupl = (TLorentzVector*) Gen_mu_4mom->At((int)Gen_QQ_mupl_idx[QQNb]);    
      if(fabs(QQ->Rapidity())<2.4 && fabs(mumi->Eta())<2.5 && fabs(mupl->Eta())<2.5 && mumi->Pt()>0.5  && mupl->Pt()>0.5){ //cuts in the gen MC
	h_QQ_Pt->Fill(QQ->Pt(), NcollW*((ispp)?pthatweight:Gen_weight) ); //(isPrompt || (!ispp))?pthatweight:1
	//	if(QQ->Pt()>6.5 && QQ->Pt()<7.5) cout<<"centrality, NcollWeight, Gen_weight, NcollW*GenW = "<<Centrality<<" "<<NcollW<<" "<<Gen_weight<<" "<<NcollW*((ispp)?pthatweight:Gen_weight)<<endl; 
      }
    }
  }

  //Fetch 'acceptance' corresponding to cuts on gen Jpsi MC
  TFile* facc = TFile::Open((TString)(isPrompt?"Prompt":"NonPrompt")+"Jpsi_acceptance.root","update");
  TH1F* JpsiAcc = (TH1F*) facc->Get("h_QQ_PtAcc");
  for(int bin=1;bin<nbins+1;bin++){
    cout<<"acceptance for bin "<<bin<<" = "<<JpsiAcc->GetBinContent(bin)<<endl;
  }

  for(int bin=1;bin<nbins+1;bin++){
    cout<<h_QQ_Pt->GetBinContent(bin)<<" "<<binLimits[bin] - binLimits[bin-1]<<endl;
    h_QQ_Pt->SetBinContent(bin,h_QQ_Pt->GetBinContent(bin) / ( JpsiAcc->GetBinContent(bin)*( binLimits[bin] - binLimits[bin-1]))); //In units of pb/GeV
  }
  h_QQ_Pt->Scale(0.06/1000); //in units of nb/GeV, multiplied by BF

  //Prepare graph for data/MC ratio
  double yr[nbins];
  double yrerr[nbins];
  for(int bin=0;bin<nbins;bin++){
    yr[bin] = dataXS->GetY()[bin] / h_QQ_Pt->GetBinContent(bin+1) ;
    yrerr[bin] = fabs(yr[bin]) * sqrt( pow(dataXS->GetErrorY(bin)/dataXS->GetY()[bin],2) + pow(h_QQ_Pt->GetBinError(bin+1)/h_QQ_Pt->GetBinContent(bin+1),2) );
    cout<<"bin "<<bin<<" Xcenter "<<Xcenter[bin]<<" ratio ="<<yr[bin]<<endl;
  }
  TGraphErrors *dataMCratio = new TGraphErrors(nbins,Xcenter,yr,Xerror,yrerr);

  //Draw everything
  gStyle->SetOptStat(0);
  
  TCanvas * c1 = new TCanvas("c1","XS(pt)",2000,2000);
  c1->Divide(1,2);
  c1->SetLeftMargin(0.11);
  c1->cd(1)->SetBottomMargin(0.12);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetLogx();
  h_QQ_Pt->SetLineColor(kBlue);
  h_QQ_Pt->SetLineWidth(2);
  h_QQ_Pt->GetXaxis()->SetTitle("Gen J/#psi p_{T} [GeV]");
  h_QQ_Pt->GetYaxis()->SetTitle("BF(J/#psi#rightarrow#mu#mu) #times d#sigma/dp_{T} [nb/GeV]");
  h_QQ_Pt->SetMinimum(0.4*dataXS->GetY()[nbins-1]);
  h_QQ_Pt->SetMaximum(2.5*h_QQ_Pt->GetMaximum());
  h_QQ_Pt->GetXaxis()->SetMoreLogLabels();
  h_QQ_Pt->GetXaxis()->SetTitleSize(0.051);
  h_QQ_Pt->GetYaxis()->SetTitleSize(0.051);
  h_QQ_Pt->Draw();
  dataXS->SetMarkerStyle(20);
  dataXS->SetMarkerSize(2.5);
  dataXS->SetMarkerColor(kRed);
  dataXS->SetLineWidth(2);
  dataXS->SetLineColor(kRed);
  dataXS->GetXaxis()->SetLimits(binLimits[0],binLimits[nbins]);
  dataXS->Draw("Psame");

  
  TLegend lgd(0.65,0.68,0.9,0.9);
  lgd.SetFillStyle(0);
  lgd.SetBorderSize(0);
  lgd.SetTextSize(0.05);
  lgd.AddEntry(h_QQ_Pt, "#sigma#times BF MC", "l");
  lgd.AddEntry(dataXS,"#sigma#times BF data","p");
  lgd.DrawClone("same");

  c1->cd(2)->SetLogy();
  c1->cd(2)->SetLogx();
  c1->cd(2)->SetBottomMargin(0.12);
  dataMCratio->SetTitle("cross section data/MC ratio");
  dataMCratio->GetYaxis()->SetTitle("#sigma_{data}/#sigma_{MC}");
  dataMCratio->GetXaxis()->SetTitle("Gen J/#psi p_{T} [GeV]");
  dataMCratio->GetXaxis()->SetLimits(binLimits[0],binLimits[nbins]);
  dataMCratio->SetMarkerStyle(20);
  dataMCratio->SetMarkerSize(2.5);
  dataMCratio->SetMarkerColor(kBlack);
  dataMCratio->SetLineWidth(2);
  dataMCratio->SetLineColor(kBlack);
  dataMCratio->GetXaxis()->SetMoreLogLabels();
  dataMCratio->GetYaxis()->SetMoreLogLabels();
  dataMCratio->GetXaxis()->SetTitleSize(0.051);
  dataMCratio->GetYaxis()->SetTitleSize(0.051);
  dataMCratio->Draw("AP");

  //Fit the ratio with exponential/constant function
  TFile* fout = TFile::Open((TString)(isPrompt?"Prompt":"NonPrompt")+"JpsiXS_scalefactor.root","update");
  TF1 *fratio = new TF1("ffitratio",fctfitratio,1.0,binLimits[nbins],3);
  // fratio->SetParameters(3.5,11,0.15); //for exponential fit function
  // fratio->SetParLimits(0,0.2,9);
  // fratio->SetParLimits(1,7,20);
  // fratio->SetParLimits(2,0,1.1);
  // fratio->SetParNames("tau","ptlim","constant");  
  
  if(isPrompt){
    if(ispp){ 
      fratio->SetParameters(0.6,20,0.4);
      fratio->SetParLimits(0,0.5,0.8);
      fratio->SetParLimits(1,19,24);
      fratio->SetParLimits(2,0.33,0.52);
    }else{
      fratio->SetParameters(0.5,13,0.0001);
      fratio->SetParLimits(0,0.2,1.4);
      fratio->SetParLimits(1,10,20);
      fratio->SetParLimits(2,0.,0.001);
    }
  }
  else{
    if(ispp){
      fratio->SetParameters(0.,20.,0.04);
      //fratio->FixParameter(0,1.05);
      //fratio->FixParameter(1,1.5);
      //fratio->SetParLimits(2,0,0.06);
    }
    else{
      fratio->SetParameters(0.,20,0.04);
      fratio->FixParameter(0,1.05);
      fratio->FixParameter(1,1.5);
      //fratio->SetParLimits(2,0,0.06);
    }
  }
  fratio->SetParNames("n","ptlim","constant");  
  fratio->Draw("same");

  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-5);
  dataMCratio->Fit("ffitratio","V","",binLimits[0],binLimits[nbins]);
  fratio->Draw("same");
  fratio->Write("SF"+(TString)(ispp?"pp":"PbPb"));
  fout->Close();

  c1->SaveAs((TString)(isPrompt?"Prompt":"NonPrompt")+"JpsiXS_dataMCcomp_"+(TString)(ispp?"pp":"PbPb")+".pdf");
}
