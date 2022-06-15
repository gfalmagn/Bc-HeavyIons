#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFormula.h"
#include "TStyle.h"
#include "../../helpers/SgMuonAcceptanceCuts.h"

typedef std::map<std::string, std::string> StrToStr;

double min(double a, double b){
  if(a<b) {return a;}
  else{return b;}
}

double GetZVal(TLorentzVector* mu, TH2D *h){
  int b = h->FindBin(mu->Eta() , fabs(mu->Pt()) ) ;
  return h->GetBinContent(b);
}

bool IsAcceptedMuon(TLorentzVector* mu, int ST, int SelType){
  bool res = false;
  double feta = fabs(mu->Eta());
  double pt = mu->Pt();

  if(SelType==0){
    //0 is: harsh global muons
    if((ST&2)>0 && feta<2.4){
      if(feta>=2.1 && pt>1.8){res = true;}
      else if(feta<1.2 && pt>3.5 ){res = true;}
      else if(feta>=1.2 && feta<2.1 && (pt > (5.77-1.8*feta)) ){res = true;}
    }
  }
  
  else if(SelType==1){
    //1 is: Loose, global muons
    if((ST&2)>0 && feta<2.4){
      if(feta<=1.0 && pt>3.3){res = true;}
      if(feta<=1.35 && feta>1.0 && (pt>6.73-3.43*feta) ){res = true;}
      if(feta>1.35 && (pt > (3.52-1.05*feta)) ){res = true;}
    }
  }
  
  else if(SelType==2){
    //2 is: Loose, tracker muons
    if((ST&8)>0 && feta<2.4){
      if(feta<=0.8 && pt>3.3){res = true;}
      if(feta<=1.5 && feta>0.8 && (pt>5.81-3.14*feta) ){res = true;}
      if(feta>1.5 && (pt > (1.89-0.526*feta)) ){res = true;}
    }
  }

  else{
    std::cout<<"Wrong muon selection type: accepted int's are 0 (harsh global), 1 (loose global), and 2 (loose tracker)"<<std::endl;}

  return res;
}

void drawSimpleSgMuon(){

  bool genOnly = false;

  auto h_test = new TH1D();
  h_test->SetDefaultSumw2(true);

  //****** Open the tree and make it scan branches one by one (SetMakeClass, to study one branch at a time) 
  TFile *file = TFile::Open("/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/acceptance/BcToJpsiMuNu_BCVEGPY_PYTHIA8_pp5TeV_RunIIpp5Spring18DR-00093_acceptance_14092020_ONIATREE_3640k.root","READ");//"/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/Oniatree_MC_Bc_trimuons_21112019.root");
  TTree *T = (TTree*)file->Get("hionia/myTree");
  int nevents= (int)T->GetEntries();
  std::cout<<"nevents = "<<nevents<<"\n";
  T->SetMakeClass(1);


  //****** Define adress for all branches
  Short_t genmu_size;
  TBranch *b_genmu_size = T->GetBranch("Gen_mu_size");
  b_genmu_size->SetAddress(&genmu_size);

  Short_t genJpsi_size;
  TBranch *b_genJpsi_size = T->GetBranch("Gen_QQ_size");
  b_genJpsi_size->SetAddress(&genJpsi_size);

  TClonesArray *Gen_QQ_4mom = new TClonesArray();
  TBranch *b_Gen_QQ_4mom; 
  T->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);

  Short_t Gen_QQ_mumi_idx[10];
  TBranch *b_Gen_QQ_mumi_idx = T->GetBranch("Gen_QQ_mumi_idx");
  b_Gen_QQ_mumi_idx->SetAddress(&Gen_QQ_mumi_idx);

  Short_t Gen_QQ_mupl_idx[10];
  TBranch *b_Gen_QQ_mupl_idx = T->GetBranch("Gen_QQ_mupl_idx");
  b_Gen_QQ_mupl_idx->SetAddress(&Gen_QQ_mupl_idx); 

  Short_t Gen_QQ_whichRec[10];
  TBranch *b_Gen_QQ_whichRec = T->GetBranch("Gen_QQ_whichRec");
  b_Gen_QQ_whichRec->SetAddress(&Gen_QQ_whichRec); 

  //***** Gen Bc
  Short_t genBc_size;
  TBranch *b_genBc_size = T->GetBranch("Gen_Bc_size");
  b_genBc_size->SetAddress(&genBc_size);

  Short_t Gen_Bc_muW_idx[10];
  TBranch *b_Gen_Bc_muW_idx = T->GetBranch("Gen_Bc_muW_idx");
  b_Gen_Bc_muW_idx->SetAddress(&Gen_Bc_muW_idx);

  Short_t Gen_Bc_QQ_idx[10];
  TBranch *b_Gen_Bc_QQ_idx = T->GetBranch("Gen_Bc_QQ_idx");
  b_Gen_Bc_QQ_idx->SetAddress(&Gen_Bc_QQ_idx);

  //***** Jpsi and its daughter muons
  Short_t recJpsi_size;
  TBranch *b_recJpsi_size = T->GetBranch("Reco_QQ_size");
  b_recJpsi_size->SetAddress(&recJpsi_size);

  TClonesArray *Reco_QQ_4mom = new TClonesArray();
  TBranch *b_Reco_QQ_4mom; 
  T->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);

  Short_t Reco_QQ_mumi_idx[100];
  TBranch *b_Reco_QQ_mumi_idx = T->GetBranch("Reco_QQ_mumi_idx");
  b_Reco_QQ_mumi_idx->SetAddress(&Reco_QQ_mumi_idx);

  Short_t Reco_QQ_mupl_idx[100];
  TBranch *b_Reco_QQ_mupl_idx = T->GetBranch("Reco_QQ_mupl_idx");
  b_Reco_QQ_mupl_idx->SetAddress(&Reco_QQ_mupl_idx);

  //***** All muons
  Short_t recmu_size;
  TBranch *b_recmu_size = T->GetBranch("Reco_mu_size");
  b_recmu_size->SetAddress(&recmu_size);

  TClonesArray *Reco_mu_4mom = new TClonesArray();
  TBranch *b_Reco_mu_4mom; 
  T->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);

  Short_t Reco_mu_whichGen[100];
  TBranch *b_Reco_mu_whichGen = T->GetBranch("Reco_mu_whichGen");
  b_Reco_mu_whichGen->SetAddress(&Reco_mu_whichGen);

  int Reco_mu_nPixWMea[100];
  TBranch *b_Reco_mu_nPixWMea = T->GetBranch("Reco_mu_nPixWMea");
  b_Reco_mu_nPixWMea->SetAddress(&Reco_mu_nPixWMea);

  bool Reco_mu_highPurity[100];
  TBranch *b_Reco_mu_highPurity = T->GetBranch("Reco_mu_highPurity");
  b_Reco_mu_highPurity->SetAddress(&Reco_mu_highPurity);

  int Reco_mu_nTrkWMea[100];
  TBranch *b_Reco_mu_nTrkWMea = T->GetBranch("Reco_mu_nTrkWMea");
  b_Reco_mu_nTrkWMea->SetAddress(&Reco_mu_nTrkWMea);

  ULong64_t Reco_mu_trig[100];
  TBranch *b_Reco_mu_trig = T->GetBranch("Reco_mu_trig");
  b_Reco_mu_trig->SetAddress(&Reco_mu_trig);

  Short_t Gen_mu_whichRec[100];
  TBranch *b_Gen_mu_whichRec = T->GetBranch("Gen_mu_whichRec");
  b_Gen_mu_whichRec->SetAddress(&Gen_mu_whichRec);

  TClonesArray *Gen_mu_4mom = new TClonesArray();
  TBranch *b_Gen_mu_4mom; 
  T->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  TClonesArray *Gen_3mu_4mom = new TClonesArray();
  TBranch *b_Gen_3mu_4mom; 
  T->SetBranchAddress("Gen_3mu_4mom", &Gen_3mu_4mom, &b_Gen_3mu_4mom);

  int recmu_SelType[100];
  TBranch *b_recmu_SelType = T->GetBranch("Reco_mu_SelectionType");
  b_recmu_SelType->SetAddress(&recmu_SelType);

  float Reco_mu_dxy[100];
  TBranch *b_Reco_mu_dxy = T->GetBranch("Reco_mu_dxy");
  b_Reco_mu_dxy->SetAddress(&Reco_mu_dxy);

  float Reco_mu_dz[100];
  TBranch *b_Reco_mu_dz = T->GetBranch("Reco_mu_dz");
  b_Reco_mu_dz->SetAddress(&Reco_mu_dz);
  
  //***** Some histograms
  TH1D *h_softestMu = new TH1D("h_softestMu","pT of the softest muon from the Gen Bc",80,0,5);
  TH1D *h_hardestMu = new TH1D("h_hardestMu","pT of the hardest muon from the Gen Bc",100,0,12);
  TH2D *h2D_softestMu = new TH2D("h2D_softestMu","pT of the softest muon from the Gen Bc",30,0,3.0,50,0,5);
  TH2D *h2D_otherMu = new TH2D("h2D_otherMu","pT of the softest muon from the Gen Bc",30,0,3.0,50,0,5);
  TH2D *h2D_otherMu2 = new TH2D("h2D_otherMu2","pT of the softest muon from the Gen Bc",30,0,3.0,50,0,5);
  TH1D *h_softestMu2 = new TH1D("h_softestMu2","pT of the softest muon from the Gen Bc",80,0,5);
  //  TH1D *h_softestMu3 = new TH1D("h_softestMu3","pT of the softest muon from the Gen Bc",80,0,5);
  TH1D *h_otherMu = new TH1D("h_otherMu","pT of the other 2 muons from the Gen Bc",80,0,5);
  TH1D *h_otherMu2 = new TH1D("h_otherMu2","pT of the other 2 muons from the Gen Bc",80,0,5);
  TH2D *h2D_genmu = new TH2D("h2D_genmu","Generated muons from B_{c};Gen muon |#eta|;Gen muon p_{T} [GeV]",39,0,2.6,60,0,6);

  int nall=0,npass=0,npass2=0,npass3=0;
  //***** Loop on events
  for(int i=0;i<nevents;i++){//nevents

    if(i%100000==0){ std::cout<<"Part of tree scanned : "<< 100*(double)i / (double)nevents<<" % \n";}

    //***** Initialize pointers for arrays
    //Reco_QQ_4mom->Clear();
    Reco_mu_4mom->Clear();
    //Gen_QQ_4mom->Clear();
    Gen_mu_4mom->Clear();
    Gen_3mu_4mom->Clear();
    
    //***** Choose needed branches
    b_genBc_size->GetEntry(i);
    //b_Gen_QQ_whichRec->GetEntry(i);
    b_Gen_QQ_mupl_idx->GetEntry(i);
    b_Gen_QQ_mumi_idx->GetEntry(i);
    b_Gen_Bc_QQ_idx->GetEntry(i);
    b_Gen_Bc_muW_idx->GetEntry(i);
    //    b_Reco_QQ_4mom->GetEntry(i);
    //    b_Gen_QQ_4mom->GetEntry(i);

    b_genmu_size->GetEntry(i);
    b_Gen_mu_4mom->GetEntry(i);
    b_Gen_mu_whichRec->GetEntry(i);

    if(!genOnly){
      //    b_recmu_size->GetEntry(i);
      //if(recmu_size>0){
      b_Reco_mu_4mom->GetEntry(i);
      //    b_Reco_mu_trig->GetEntry(i);
      b_Reco_mu_nPixWMea->GetEntry(i);
      b_Reco_mu_nTrkWMea->GetEntry(i);
      b_Reco_mu_dxy->GetEntry(i);
      b_Reco_mu_dz->GetEntry(i);
      //b_Reco_mu_whichGen->GetEntry(i);
      //    b_Reco_mu_highPurity->GetEntry(i);
      b_recmu_SelType->GetEntry(i);  
      //}
    }

    for(int BcNb=0;BcNb<genBc_size;BcNb++){
      TLorentzVector *genmuW = (TLorentzVector*) Gen_mu_4mom->At(Gen_Bc_muW_idx[BcNb]);
      Short_t mumiidx = Gen_QQ_mumi_idx[Gen_Bc_QQ_idx[BcNb]];
      Short_t muplidx = Gen_QQ_mupl_idx[Gen_Bc_QQ_idx[BcNb]];
      TLorentzVector *genmumi = (TLorentzVector*) Gen_mu_4mom->At(mumiidx);
      TLorentzVector *genmupl = (TLorentzVector*) Gen_mu_4mom->At(muplidx);
      // if(Gen_Bc_muW_idx[BcNb]>=genmu_size || muplidx>=genmu_size || mumiidx>=genmu_size){
      // 	cout<<"wrong index for muon from gen Bc. muWidx, mumiidx, muplidx, Genmusize = "<<Gen_Bc_muW_idx[BcNb]<<" "<<muplidx<<" "<<mumiidx<<" "<<genmu_size<<endl;
      // 	continue;}

      h2D_genmu->Fill(fabs(genmuW->Eta()) , genmuW->Pt());  
      h2D_genmu->Fill(fabs(genmumi->Eta()) , genmumi->Pt());  
      h2D_genmu->Fill(fabs(genmupl->Eta()) , genmupl->Pt());  

      if(!genOnly){

	//hardest muon
	Short_t muhardidx = Gen_Bc_muW_idx[BcNb]; Short_t muh1idx = mumiidx; Short_t muh2idx = muplidx;
	if(genmumi->Pt()>genmuW->Pt() && genmumi->Pt()>genmupl->Pt()) {
	  muhardidx = mumiidx;
	  muh1idx = Gen_Bc_muW_idx[BcNb];}
	if(genmupl->Pt()>genmuW->Pt() && genmupl->Pt()>genmumi->Pt()) {
	  muhardidx = muplidx;
	  muh2idx = Gen_Bc_muW_idx[BcNb];}

	Short_t recmuh1idx = Gen_mu_whichRec[muh1idx];
	Short_t recmuh2idx = Gen_mu_whichRec[muh2idx];      
	Short_t recmuhardidx = Gen_mu_whichRec[muhardidx];      

	if(recmuh1idx>-1 && recmuh2idx>-1){
	  if((recmu_SelType[muh1idx]&2)>0 && (recmu_SelType[muh2idx]&2)>0 &&
	     (recmu_SelType[muh1idx]&8)>0 && (recmu_SelType[muh2idx]&8)>0 &&
	     Reco_mu_nPixWMea[recmuh1idx]>0 && Reco_mu_nPixWMea[recmuh2idx]>0 &&
	     Reco_mu_nTrkWMea[recmuh1idx]>5 && Reco_mu_nTrkWMea[recmuh2idx]>5 &&
	     fabs(Reco_mu_dxy[recmuh1idx])<0.3 && fabs(Reco_mu_dxy[recmuh2idx])<0.3 &&
	     fabs(Reco_mu_dz[recmuh1idx])<20 && fabs(Reco_mu_dz[recmuh2idx])<20 &&
	     looseAcc(((TLorentzVector*) Reco_mu_4mom->At(recmuh1idx))->Pt(),((TLorentzVector*) Reco_mu_4mom->At(recmuh1idx))->Eta()) && looseAcc(((TLorentzVector*) Reco_mu_4mom->At(recmuh2idx))->Pt(),((TLorentzVector*) Reco_mu_4mom->At(recmuh2idx))->Eta()) &&
	     (tightAcc(((TLorentzVector*) Reco_mu_4mom->At(recmuh1idx))->Pt(),((TLorentzVector*) Reco_mu_4mom->At(recmuh1idx))->Eta()) || tightAcc(((TLorentzVector*) Reco_mu_4mom->At(recmuh2idx))->Pt(),((TLorentzVector*) Reco_mu_4mom->At(recmuh2idx))->Eta())) &&
	     fabs(((TLorentzVector*) Gen_mu_4mom->At(muhardidx))->Eta())<2.4
	     ){
	    h_hardestMu->Fill(((TLorentzVector*) Gen_mu_4mom->At(muhardidx))->Pt());
	  }

	}

	//softest muon
	Short_t musoftidx = Gen_Bc_muW_idx[BcNb]; Short_t mu1idx = mumiidx; Short_t mu2idx = muplidx;
	if(genmumi->Pt()<genmuW->Pt() && genmumi->Pt()<genmupl->Pt()) {
	  musoftidx = mumiidx;
	  mu1idx = Gen_Bc_muW_idx[BcNb];}
	if(genmupl->Pt()<genmuW->Pt() && genmupl->Pt()<genmumi->Pt()) {
	  musoftidx = muplidx;
	  mu2idx = Gen_Bc_muW_idx[BcNb];}

	nall+=1;
	Short_t recmu1idx = Gen_mu_whichRec[mu1idx];
	Short_t recmu2idx = Gen_mu_whichRec[mu2idx];      
	Short_t recmusoftidx = Gen_mu_whichRec[musoftidx];      
	if(recmu1idx==-1 || recmu2idx==-1) continue;

	if((recmu_SelType[mu1idx]&2)>0 && (recmu_SelType[mu2idx]&2)>0 && 
	   (recmu_SelType[mu1idx]&8)>0 && (recmu_SelType[mu2idx]&8)>0 && 
	   Reco_mu_nPixWMea[recmu1idx]>0 && Reco_mu_nPixWMea[recmu2idx]>0 && 
	   Reco_mu_nTrkWMea[recmu1idx]>5 && Reco_mu_nTrkWMea[recmu2idx]>5 && 
	   fabs(Reco_mu_dxy[recmu1idx])<0.3 && fabs(Reco_mu_dxy[recmu2idx])<0.3 &&
	   fabs(Reco_mu_dz[recmu1idx])<20 && fabs(Reco_mu_dz[recmu2idx])<20
	   ){
	  npass+=1;
	  h_softestMu->Fill(((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Pt());
	  if(tightAcc(((TLorentzVector*) Gen_mu_4mom->At(mu1idx))->Pt(),((TLorentzVector*) Gen_mu_4mom->At(mu1idx))->Eta()) && tightAcc(((TLorentzVector*) Gen_mu_4mom->At(mu2idx))->Pt(),((TLorentzVector*) Gen_mu_4mom->At(mu2idx))->Eta()) ){
	    h_softestMu2->Fill(((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Pt());
	    npass2+=1;
	  }

	  if(tightAcc(((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Pt(),((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Eta())) npass3+=1;
	  //h_otherMu->Fill(((TLorentzVector*) Gen_mu_4mom->At(mu1idx))->Pt());
	  // h2D_softestMu->Fill(fabs(((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Eta()) , ((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Pt());
	  // if(recmusoftidx>-1)
	  // 	h2D_otherMu->Fill(fabs(((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Eta()) , ((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Pt());
	  // else
	  // 	h2D_otherMu2->Fill(fabs(((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Eta()) , ((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Pt());
	  //	h_otherMu2->Fill(((TLorentzVector*) Gen_mu_4mom->At(mu2idx))->Pt() - ((TLorentzVector*) Gen_mu_4mom->At(musoftidx))->Pt());
	  //h_otherMu->Fill(((TLorentzVector*) Gen_mu_4mom->At(mu2idx))->Pt());
	}
      }
      
    }
  }

  cout<<"Bc's number: all, pass 2 hybrid-soft, pass 2 HS+TightAcc, pass 3rd muon in tightAcc = "<<nall<<" "<<npass<<" "<<npass2<<" "<<npass3<<endl;


  //****************************************************************
  //Drawing all
  //****************************************************************

  gStyle->SetOptStat(0);

  if(!genOnly){
    TCanvas * c2 = new TCanvas("c2","c2",1500,1500);  
    c2->SetLeftMargin(0.12);
    c2->SetTopMargin(0.04);
    h_softestMu->SetTitle(";p_{T}(#mu) [GeV];counts");
    h_softestMu->SetLineWidth(3);
    h_softestMu->Draw("E");
    h_softestMu2->Scale(1*h_softestMu->Integral()/h_softestMu2->Integral());
    h_softestMu2->SetLineWidth(3);
    h_softestMu2->SetLineColor(kRed);
    //    h_softestMu2->Draw("Esame");

    // h_softestMu3->Scale(0.8*h_softestMu->Integral()/h_softestMu3->Integral());
    // h_softestMu3->SetLineWidth(3);
    // h_softestMu3->SetLineColor(kGreen+3);
    // h_softestMu3->Draw("Esame");
    c2->SaveAs("pT_softestMuon_otherMuonsAreHybridSoft.pdf");

    TCanvas * c3 = new TCanvas("c3","c3",1500,1500);  
    c3->SetLeftMargin(0.12);
    c3->SetTopMargin(0.04);
    h_hardestMu->SetTitle(";p_{T}(#mu) [GeV];counts");
    h_hardestMu->SetLineWidth(3);
    h_hardestMu->Draw("E");
    c3->SaveAs("pT_hardestMuon_otherMuonsAreHybridSoftAndInLooseAcc_oneMuInTightAcc.pdf");
  }
    
  gStyle->SetPalette(kBlueRedYellow);

  TCanvas * c3 = new TCanvas("c3","c3",1500,1500);  
  c3->SetRightMargin(0.15);
  c3->SetLeftMargin(0.08);
  c3->SetBottomMargin(0.08);
  //  h2D_genmu->SetTitleOffset(0.1);
  h2D_genmu->Draw("COLZ");
  c3->SaveAs("PtEtaMap_generatedMuonsFromBc_noCuts.pdf");

}
