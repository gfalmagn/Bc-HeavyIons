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
#include "../../PbPb18/Utilities/EVENTUTILS.h"

typedef std::map<std::string, std::string> StrToStr;

double min(double a, double b){
  if(a<b) {return a;}
  else{return b;}
}

double GetZVal(TLorentzVector* mu, TH2D *h){
  int b = h->FindBin(mu->Eta() , fabs(mu->Pt()) ) ;
  return h->GetBinContent(b);
}

bool looseAcc(float pt, float eta, bool TM = false){
  if(TM) return (fabs(eta) < 2.4 &&
                 ((fabs(eta) < 1.1 && pt >= 3.3) ||
                  (1.1 <= fabs(eta) && fabs(eta) < 1.3 && pt >= 13.2-9.*fabs(eta) ) ||
                  (1.3 <= fabs(eta) && pt >= 0.8 && pt >= 3.02-1.17*fabs(eta) )));
  else return (fabs(eta) < 2.4 &&
               ((fabs(eta) < 0.3 && pt >= 3.4) ||
                (fabs(eta) > 0.3 && fabs(eta) < 1.1 && pt >= 3.3) ||
                (fabs(eta) > 1.1 && fabs(eta) < 1.55 && pt >= 2.1 && pt >= 7.7-4.0*fabs(eta) ) ||
                (fabs(eta) > 1.55 && pt >= 1.2 && pt >= 4.25-1.39*fabs(eta)) ));
}

bool tightAcc(float pt, float eta){
  return (fabs(eta) < 2.4 &&
	  ((fabs(eta) < 1.2 && pt >= 3.5) ||
	   (1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 5.47-1.89*fabs(eta)) ||
	   (2.1 <= fabs(eta) && pt >= 1.5)));
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

void get_acceptance(bool ispp = true){

  bool withTrig = false;
  bool L1trig = false;
  bool doTnP = withTrig;
  bool doJpsiEff = false;
  bool doJpsiAcc = false;
  bool addNcollW = true && (!ispp);

  auto h_test = new TH1D();
  h_test->SetDefaultSumw2(true);

  //****** Open the tree and make it scan branches one by one (SetMakeClass, to study one branch at a time) 
  TFile *file = TFile::Open(ispp?"/data_CMS/cms/falmagne/tuples/pp17/dimuons/MC/PromptJpsi/JPsiMM_TuneCUETP8M1_5p02TeV_pythia8_ptHatMin2_dimuons_ONIATREE_27122019.root":"/data_CMS/cms/falmagne/tuples/PbPb18/PromptJpsi/MC/JPsi_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8_dimuons_oniatree_19052020.root","READ"); //"/home/llr/cms/falmagne/production/OniaDev/pp17/CMSSW_9_4_12/src/Oniatree_MC_Bc_trimuons_22012019.root"
  TTree *T = (TTree*)file->Get("hionia/myTree");
  int nevents= (int)T->GetEntries();
  std::cout<<"nevents = "<<nevents<<"\n";
  T->SetMakeClass(1);


  //****** Define adress for all branches
  float Gen_weight; 
  TBranch *b_Gen_weight = T->GetBranch("Gen_weight");
  b_Gen_weight->SetAddress(&Gen_weight);

  int Centrality; 
  TBranch *b_Centrality = ispp?NULL:T->GetBranch("Centrality");
  if(!ispp) b_Centrality->SetAddress(&Centrality);

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

  int recmu_SelType[100];
  TBranch *b_recmu_SelType = T->GetBranch("Reco_mu_SelectionType");
  b_recmu_SelType->SetAddress(&recmu_SelType);

  float Reco_mu_dxy[100];
  TBranch *b_Reco_mu_dxy = T->GetBranch("Reco_mu_dxy");
  b_Reco_mu_dxy->SetAddress(&Reco_mu_dxy);

  float Reco_mu_dz[100];
  TBranch *b_Reco_mu_dz = T->GetBranch("Reco_mu_dz");
  b_Reco_mu_dz->SetAddress(&Reco_mu_dz);
  
  //****** Some parameters
  // double ScaleToXS = 208*208 * 460 * 2.54e-3 * 0.668 / (double)nevents; // A^2 * Lumi_PbPb[mub-1] * (XS_Bc_pp * BF((J/psi -> mu mu) mu nu))[mub] * (XS(5.02 TeV) / XS(7 TeV))
  // std::cout<<"Scaling the number of generated events by : "<<ScaleToXS<<std::endl;
  TString MuTypeForEfficiency = withTrig?"globaltracker":"tracker"; //"all", "global", "nonglobal", "tracker","globaltracker"
  std::string QQSelectionAcc = "all";
  int trigbit =  ispp?8:(L1trig?32768:131072) ;  //2^15,16,17 = L1,L2,L3 step //2^3 for L1DoubleMu0 in pp

  //***** Some histograms
  int AEnbins = 27;
  double acceffBins[AEnbins+1]; 
  for(int l=0;l<19;l++) acceffBins[l] = l*0.5;
  for(int l=19;l<25;l++) acceffBins[l] = 9. + (l-19);
  acceffBins[25] = 16.;  acceffBins[26] = 19;  acceffBins[27] = 25;

  TH1D *h_JpsiNb = new TH1D("JpsiNb","Jpsi number",2,0,2);
  TH1D *h_recJpsi_Pt = new TH1D("recJpsiPt","transverse momentum of reconstructed J/#psi;p_{T}(J/#psi) [GeV]",30,2,20);
  TH1D *h_rec2Jpsi_Pt = new TH1D("rec2JpsiPt","transverse momentum of reconstructed J/#psi + new cuts;p_{T}(J/#psi) [GeV]",30,2,20);
  TH2D *h_recJpsi_YPt = new TH2D("recJpsiYPt","Reconstructed J/#psi;|Rapidity|;p_{T} [GeV]",13,0,2.6,AEnbins,acceffBins);
  TH1D *h_genJpsi_Pt = new TH1D("genJpsiPt","transverse momentum of observable J/#psi;p_{T}(J/#psi) [GeV]",30,2,20);
  TH1D *h_genJpsi_Y = new TH1D("genJpsiY","rapidity of observable J/#psi;|Rapidity|",50,0,3);
  TH2D *h_genJpsi_YPt = new TH2D("genJpsiYPt","All generated J/#psi;|Rapidity|;p_{T} [GeV]",13,0,2.6,AEnbins,acceffBins);
  int muptbins = doTnP?60:90;
  int muetabins = doTnP?39:52;
  double ptlow = doTnP?0.03:0; double pthigh = doTnP?6.03:6;
  TH2D *h_genmu_EtaPt = new TH2D("genmuEtaPt","All generated #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_recmu_EtaPt = new TH2D("recmuEtaPt","All reconstructed #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_genmu_CentPt = new TH2D("genmuCentPt","All generated #mu;Centrality;p_{T} [GeV]",withTrig?30:50,0,100,(int)(0.7*muptbins),0.8,7);
  TH2D *h_recmu_CentPt = new TH2D("recmuCentPt","All reconstructed #mu;Centrality;p_{T} [GeV]",withTrig?30:50,0,100,(int)(0.7*muptbins),0.8,7);
  TH2D *h_genmu_EtaPt020 = new TH2D("genmuEtaPt0-20","All generated #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_recmu_EtaPt020 = new TH2D("recmuEtaPt0-20","All reconstructed #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_genmu_EtaPt2050 = new TH2D("genmuEtaPt20-50","All generated #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_recmu_EtaPt2050 = new TH2D("recmuEtaPt20-50","All reconstructed #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_genmu_EtaPt50100 = new TH2D("genmuEtaPt50-100","All generated #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);
  TH2D *h_recmu_EtaPt50100 = new TH2D("recmuEtaPt50-100","All reconstructed #mu;|#eta|;p_{T} [GeV]",muetabins,0,2.6,muptbins,ptlow,pthigh);

  int NQQ_recoAccepted =0;
  int NQQ_accepted =0;
  int NQQ_reco =0;
  int genmuNb_goodTag = 0;
  int recmuNb_goodTag = 0;

  //***** Loop on events
  for(int i=0;i<nevents;i++){//nevents

    if(i%100000==0){ std::cout<<"Part of tree scanned : "<< 100*(double)i / (double)nevents<<" % \n";}

    //***** Initialize pointers for arrays
    //Reco_QQ_4mom->Clear();
    Reco_mu_4mom->Clear();
    Gen_QQ_4mom->Clear();
    Gen_mu_4mom->Clear();

    if(!ispp && addNcollW) b_Centrality->GetEntry(i);
    
    //***** Choose needed branches
    if(doTnP){
      b_genJpsi_size->GetEntry(i);
      b_Gen_QQ_whichRec->GetEntry(i);
      b_Gen_QQ_mupl_idx->GetEntry(i);
      b_Gen_QQ_mumi_idx->GetEntry(i);
      //    b_Reco_QQ_4mom->GetEntry(i);
      b_Gen_QQ_4mom->GetEntry(i);
    }

    b_genmu_size->GetEntry(i);
    if(genmu_size>0){
      b_Gen_mu_4mom->GetEntry(i);
      b_Gen_mu_whichRec->GetEntry(i);}

    b_recmu_size->GetEntry(i);
    if(recmu_size>0){
      b_Reco_mu_4mom->GetEntry(i);
      b_Reco_mu_trig->GetEntry(i);
      b_Reco_mu_nPixWMea->GetEntry(i);
      b_Reco_mu_nTrkWMea->GetEntry(i);
      b_Reco_mu_dxy->GetEntry(i);
      b_Reco_mu_dz->GetEntry(i);
      b_Reco_mu_whichGen->GetEntry(i);
      b_Reco_mu_highPurity->GetEntry(i);
      b_recmu_SelType->GetEntry(i);  
    }

    if(ispp) { Gen_weight = 1;
    }else{ b_Gen_weight->GetEntry(i);  }

    float wei = Gen_weight;
    if(addNcollW && !ispp) wei *= findNcoll(Centrality);

    //Only muons from a gen Jpsi
    if(doTnP || doJpsiEff){
      for (int genQQNb=0; genQQNb<genJpsi_size; genQQNb++){
    	TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom->At(genQQNb);    
    	//      h_genJpsi_YPt->Fill(fabs(genQQ->Rapidity()),genQQ->Pt(), wei);

    	//      int recQQIdx = Gen_QQ_whichRec[genQQNb];
    	int recmumiIdx = Gen_mu_whichRec[Gen_QQ_mumi_idx[genQQNb]];
    	int recmuplIdx = Gen_mu_whichRec[Gen_QQ_mupl_idx[genQQNb]];

    	int tagIdx, probeIdx, tagGenIdx, probeGenIdx;
    	//Randomize when the mumi or the mupl is the tag or the probe

    	if(doTnP){
    	  for (int k=0;k<2;k++){ //try both muons for the tag
    	    if(k==0){
    	      tagIdx = recmumiIdx; tagGenIdx = Gen_QQ_mumi_idx[genQQNb];
    	      probeIdx = recmuplIdx; probeGenIdx = Gen_QQ_mupl_idx[genQQNb];
    	    } else{
    	      tagIdx = recmuplIdx; tagGenIdx = Gen_QQ_mupl_idx[genQQNb];
    	      probeIdx = recmumiIdx; probeGenIdx = Gen_QQ_mumi_idx[genQQNb];
    	    }
	
    	    TLorentzVector *genmuTag = (TLorentzVector*) Gen_mu_4mom->At(tagGenIdx);
    	    if(tagIdx==-1) continue;
    	    //      TLorentzVector *recmuTag = (TLorentzVector*) Reco_mu_4mom->At(tagIdx);
    	    if(1==1
    	       && looseAcc(genmuTag->Pt(),genmuTag->Eta())
	       && (tightAcc(genmuTag->Pt(),genmuTag->Eta()) || !withTrig)
    	       // && ((recmu_SelType[tagIdx])&(int)pow(2,12))>0
    	       && Reco_mu_nPixWMea[tagIdx]>0
    	       && Reco_mu_nTrkWMea[tagIdx]>5
    	       && Reco_mu_dxy[tagIdx]<0.3
    	       && Reco_mu_dz[tagIdx]<20
    	       && (((Reco_mu_trig[tagIdx])&(ULong64_t)pow(2, (int)(ispp?16:23) ))>0 || ((Reco_mu_trig[tagIdx])&(ULong64_t)pow(2, (int)(ispp?25:32) ))>0) //2^23 = L3Mu3_NHitQ10 // 2^16 for pp L3Mu3 //+9 on trigger bit for version 2
	 
    	       && (//(MuTypeForEfficiency=="all")
    		   // || (MuTypeForEfficiency=="globaltracker" && 
		   (recmu_SelType[tagIdx]&2) > 0 && (recmu_SelType[tagIdx]&8) > 0
		   //)
    		   // || (MuTypeForEfficiency=="global" && (recmu_SelType[tagIdx]&2) > 0)
    		   // || (MuTypeForEfficiency=="nonglobal" && (recmu_SelType[tagIdx]&2) <= 0)
    		   // || (MuTypeForEfficiency=="tracker" && (recmu_SelType[tagIdx]&(ULong64_t)pow(2,12)) > 0 && Reco_mu_highPurity[tagIdx]) //TMOneStaTight
    		   )
    	       ){

	
    	      //Fill Pt,Eta for gen muons
    	      TLorentzVector *genmuProbe = (TLorentzVector*) Gen_mu_4mom->At(probeGenIdx);
    	      h_genmu_EtaPt->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei);
	      if(!ispp){
		if(looseAcc(genmuProbe->Pt(),genmuProbe->Eta())
		   && (tightAcc(genmuProbe->Pt(),genmuProbe->Eta()) || !withTrig)){
		  h_genmu_CentPt->Fill((float)Centrality/2,genmuProbe->Pt(), wei);}
		if(Centrality<2*20) h_genmu_EtaPt020->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei); 
		if(Centrality>2*20 && Centrality<2*50) h_genmu_EtaPt2050->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei);
		if(Centrality>2*50 && Centrality<2*100) h_genmu_EtaPt50100->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei);
	      }
    	      genmuNb_goodTag+=1;

    	      if(probeIdx==-1) continue;
    	      if(// ((recmu_SelType[probeIdx])&(int)pow(2,12))>0 &&
    		 Reco_mu_nPixWMea[probeIdx]>0
    		 && Reco_mu_nTrkWMea[probeIdx]>5
    		 && Reco_mu_dxy[probeIdx]<0.3
    		 && Reco_mu_dz[probeIdx]<20
    		 && ((MuTypeForEfficiency=="all")
    		     || (MuTypeForEfficiency=="globaltracker" && (recmu_SelType[probeIdx]&2) > 0 && (recmu_SelType[probeIdx]&8) > 0)
    		     || (MuTypeForEfficiency=="global" && (recmu_SelType[probeIdx]&2) > 0)
    		     || (MuTypeForEfficiency=="nonglobal" && (recmu_SelType[probeIdx]&2) <= 0)
    		     || (MuTypeForEfficiency=="tracker" && (recmu_SelType[probeIdx]&(ULong64_t)pow(2,12)) > 0 && Reco_mu_highPurity[probeIdx]) //TMOneStaTight
    		     )

    		 && (!withTrig || ((Reco_mu_trig[probeIdx])&(int)( trigbit ))>0)
    		 //	      && (genmu->P()>5.5 || fabs(genmu->Eta())<1)
    		 ){
    		//TLorentzVector *recmuProbe = (TLorentzVector*) Reco_mu_4mom->At(probeIdx);
    		h_recmu_EtaPt->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei);
		if(!ispp){
		  if(looseAcc(genmuProbe->Pt(),genmuProbe->Eta())
		     && (tightAcc(genmuProbe->Pt(),genmuProbe->Eta()) || !withTrig)){
		    h_recmu_CentPt->Fill((float)Centrality/2,genmuProbe->Pt(), wei);}
		  if(Centrality<2*20) h_recmu_EtaPt020->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei);
		  if(Centrality>2*20 && Centrality<2*50) h_recmu_EtaPt2050->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei);
		  if(Centrality>2*50 && Centrality<2*100) h_recmu_EtaPt50100->Fill(fabs(genmuProbe->Eta()),genmuProbe->Pt(), wei);
		}
    		recmuNb_goodTag+=1;
    	      }
    	    }

    	  }
    	}


    	if(doJpsiEff){
    	  TLorentzVector *genmupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[genQQNb]);
    	  TLorentzVector *genmumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[genQQNb]);
      
    	  if(looseAcc(genmupl->Pt(),genmupl->Eta()) && looseAcc(genmumi->Pt(),genmumi->Eta())
    	     ){
    	    h_genJpsi_YPt->Fill(genQQ->Rapidity(),genQQ->Pt(), wei);

    	    if (recmuplIdx==-1 || recmumiIdx==-1) continue;
    	    if( //((recmu_SelType[recmuplIdx])&(int)pow(2,12))>0 &&
    	       Reco_mu_nPixWMea[recmuplIdx]>0
    	       && Reco_mu_nTrkWMea[recmuplIdx]>5
    	       && Reco_mu_dxy[recmuplIdx]<0.3
    	       && Reco_mu_dz[recmuplIdx]<20
    	       && ((Reco_mu_trig[recmuplIdx])&(int)trigbit)>0
    	       && ((MuTypeForEfficiency=="all")
    		   || (MuTypeForEfficiency=="globaltracker" && (recmu_SelType[recmuplIdx]&2) > 0 && (recmu_SelType[recmuplIdx]&8) > 0)
    		   || (MuTypeForEfficiency=="global" && (recmu_SelType[recmuplIdx]&2) > 0)
    		   || (MuTypeForEfficiency=="nonglobal" && (recmu_SelType[recmuplIdx]&2) <= 0)
    		   || (MuTypeForEfficiency=="tracker" && (recmu_SelType[recmuplIdx]&(int)pow(2,12)) > 0 && Reco_mu_highPurity[recmuplIdx]) 
    		   )
    	       //	    && ((recmu_SelType[recmumiIdx])&(int)pow(2,12))>0
    	       && Reco_mu_nPixWMea[recmumiIdx]>0
    	       && Reco_mu_nTrkWMea[recmumiIdx]>5
    	       && Reco_mu_dxy[recmumiIdx]<0.3
    	       && Reco_mu_dz[recmumiIdx]<20
    	       && ((Reco_mu_trig[recmumiIdx])&(int)trigbit)>0
    	       && ((MuTypeForEfficiency=="all")
    		   || (MuTypeForEfficiency=="globaltracker" && (recmu_SelType[recmumiIdx]&2) > 0 && (recmu_SelType[recmumiIdx]&8) > 0)
    		   || (MuTypeForEfficiency=="global" && (recmu_SelType[recmumiIdx]&2) > 0)
    		   || (MuTypeForEfficiency=="nonglobal" && (recmu_SelType[recmumiIdx]&2) <= 0)
    		   || (MuTypeForEfficiency=="tracker" && (recmu_SelType[recmumiIdx]&(int)pow(2,12)) > 0 && Reco_mu_highPurity[recmumiIdx]) 
    		   )
    		){

    	      //      	h_genJpsi_Pt->Fill(genQQ->Pt(), wei);
    	      h_recJpsi_YPt->Fill(genQQ->Rapidity(),genQQ->Pt(), wei);
      
    	      // if(((Reco_mu_trig[recmuplIdx])&(int)trigbit)>0  && ((Reco_mu_trig[recmumiIdx])&(int)trigbit)>0){
    	      //   h_recJpsi_Pt->Fill(genQQ->Pt(), wei);}

    	      // if(tightAcc(genmupl->Pt(),genmupl->Eta()) && tightAcc(genmumi->Pt(),genmumi->Eta())){
    	      //   h_recJpsi_Pt->Fill(genQQ->Pt(), wei);}
    	      // if(looseAcc(genmupl->Pt(),genmupl->Eta()) && looseAcc(genmumi->Pt(),genmumi->Eta())){
    	      //   h_rec2Jpsi_Pt->Fill(genQQ->Pt(), wei);}

    	    }
    	  }
    	}

	


    	// bool AcceptanceAndEff = false;
    	// if(AcceptanceAndEff && recQQIdx >= 0){

    	// 	NQQ_reco+=1;
    	// 	int rec_mumi_idx = Reco_QQ_mumi_idx[recQQIdx];
    	// 	int rec_mupl_idx = Reco_QQ_mupl_idx[recQQIdx];
    	// 	if (rec_mumi_idx==-1 || rec_mupl_idx==-1){std::cout<<"Beware, one muon was not reconstructed, so the dimuon should not be reconstructed!"<<std::endl;}
    	// 	if (fabs(Reco_QQ_charge[recQQIdx])!=0){std::cout<<"Beware, QQ charge is wrong, so the dimuon should not be reconstructed!"<<std::endl;}

    	// 	if( IsAcceptedDimuon(rec_mumi_idx, rec_mupl_idx, Reco_mu_4mom, recmu_SelType, QQSelectionAcc, false
    	// 			      ) ){
	  
    	// 	  NQQ_recoAccepted+=1;
    	// 	  double w_QQ = 1;

    	// 	  h_genQQ_Pt->Fill(genQQ->Pt() , w_QQ);
    	// 	  h_genQQ_Y->Fill(fabs(genQQ->Rapidity()) , w_QQ);
    	// 	  h_genQQ_YPt_acc->Fill(fabs(genQQ->Rapidity()),genQQ->Pt() , w_QQ);

    	// 	}
    	// }
    	// if(!AcceptanceAndEff){

    	// 	if( IsAcceptedDimuon( Gen_QQ_mumi_idx[genJpsiNb], Gen_QQ_mupl_idx[genJpsiNb], Gen_mu_4mom, recmu_SelType, QQSelectionAcc, true
    	// 			      ) ){ //ignore SelectionType when checking only acceptance
	  
    	// 	  NQQ_accepted+=1;
    	// 	  double w_QQ = 1;

    	// 	  h_genQQ_Pt->Fill(genQQ->Pt() , w_QQ);
    	// 	  h_genQQ_Y->Fill(fabs(genQQ->Rapidity()) , w_QQ);
    	// 	  h_genQQ_YPt_acc->Fill(fabs(genQQ->Rapidity()),genQQ->Pt() , w_QQ);
    	// 	}
    	// }

      }
      ///////////////end Jpsi loop
    }

    if(!doTnP){ // without TnP method
      //***** efficiency of single muons
    
      // double wei = 1;
      for (int genmuNb=0; genmuNb<genmu_size; genmuNb++){
    	// cout<<genmu_size<<" genmu#"<<genmuNb<<" Gen_mu_whichRec[genmuNb] = "<<Gen_mu_whichRec[genmuNb]<<endl;
    	// cout<<Gen_mu_4mom<<endl;
    	// cout<<Reco_mu_4mom<<endl;
    	// cout<<((TLorentzVector*) Gen_mu_4mom->At(genmuNb))->Pt()<<endl;
    	TLorentzVector *genmu = (TLorentzVector*) Gen_mu_4mom->At(genmuNb);
    	h_genmu_EtaPt->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
	if(!ispp){
	  if(looseAcc(genmu->Pt(),genmu->Eta(), (MuTypeForEfficiency=="tracker") )
	     && (tightAcc(genmu->Pt(),genmu->Eta()) || !withTrig)){
	    h_genmu_CentPt->Fill((float)Centrality/2,genmu->Pt(), wei);}
	  if(Centrality<2*20) h_genmu_EtaPt020->Fill(fabs(genmu->Eta()),genmu->Pt(), wei); 
	  if(Centrality>2*20 && Centrality<2*50) h_genmu_EtaPt2050->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
	  if(Centrality>2*50 && Centrality<2*100) h_genmu_EtaPt50100->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
	}
    	int recmuNb = Gen_mu_whichRec[genmuNb];
	
    	if(recmuNb>-1){
    	  if( (MuTypeForEfficiency=="all")
    	      || (MuTypeForEfficiency=="globaltracker" && (recmu_SelType[recmuNb]&2) > 0 && (recmu_SelType[recmuNb]&8) > 0)
    	      || (MuTypeForEfficiency=="global" && (recmu_SelType[recmuNb]&2) > 0)
    	      || (MuTypeForEfficiency=="nonglobal" && (recmu_SelType[recmuNb]&2) <= 0)
    	      || (MuTypeForEfficiency=="tracker" && (recmu_SelType[recmuNb]&(int)pow(2,12)) > 0 && Reco_mu_highPurity[recmuNb])
    	      //|| (MuTypeForEfficiency=="tracker" && (recmu_SelType[recmuNb]&8) > 0)
    	      ){
    	    TLorentzVector *recmu = (TLorentzVector*) Reco_mu_4mom->At(recmuNb);
    	    if( //((recmu_SelType[recmuNb])&(int)pow(2,12))>0 &&
    	       Reco_mu_nPixWMea[recmuNb]>0
    	       && Reco_mu_nTrkWMea[recmuNb]>5
    	       && Reco_mu_dxy[recmuNb]<0.3
    	       && Reco_mu_dz[recmuNb]<20
    	       //&& ((Reco_mu_trig[recmuNb])&(int)pow(2,16))>0 //2^15,16,17 = L1,L2,L3 step
    		){
	      if(!ispp){
		if(looseAcc(genmu->Pt(),genmu->Eta(), (MuTypeForEfficiency=="tracker") )
		   && (tightAcc(genmu->Pt(),genmu->Eta()) || !withTrig)){
		  h_recmu_CentPt->Fill((float)Centrality/2,genmu->Pt(), wei);}
		if(Centrality<2*20) h_recmu_EtaPt020->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
		if(Centrality>2*20 && Centrality<2*50) h_recmu_EtaPt2050->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
		if(Centrality>2*50 && Centrality<2*100) h_recmu_EtaPt50100->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
	      }
    	      h_recmu_EtaPt->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
    	    }
    	  }
    	}
      }
    
       // for (int recmuNb=0; recmuNb<recmu_size; recmuNb++){
       // 	if(Reco_mu_whichGen[recmuNb]>-1
       // 	   //&& 
       // 	   ){
       // 	  if( (MuTypeForEfficiency=="all")
       // 	      || (MuTypeForEfficiency=="globaltracker" && (recmu_SelType[recmuNb]&2) > 0 && (recmu_SelType[recmuNb]&8) > 0)
       // 	      || (MuTypeForEfficiency=="global" && (recmu_SelType[recmuNb]&2) > 0)
       // 	      || (MuTypeForEfficiency=="nonglobal" && (recmu_SelType[recmuNb]&2) <= 0)
       // 	      || (MuTypeForEfficiency=="tracker" && (recmu_SelType[recmuNb]&8) > 0)
       // 	      ){
       // 	    TLorentzVector *genmu = (TLorentzVector*) Gen_mu_4mom->At(Reco_mu_whichGen[recmuNb]);
       // 	    TLorentzVector *recmu = (TLorentzVector*) Reco_mu_4mom->At(recmuNb);
       // 	    if( //((recmu_SelType[recmuNb])&(int)pow(2,12))>0 &&
       // 	       Reco_mu_nPixWMea[recmuNb]>0
       // 	       && Reco_mu_nTrkWMea[recmuNb]>5
       // 	       && Reco_mu_dxy[recmuNb]<0.3
       // 	       && Reco_mu_dz[recmuNb]<20
       // 	       //&& ((Reco_mu_trig[recmuNb])&(int)pow(2,16))>0 //2^15,16,17 = L1,L2,L3 step
       // 		){
       // 	      h_recmu_EtaPt->Fill(fabs(genmu->Eta()),genmu->Pt(), wei);
       // 	    }
       // 	  }
       // 	}
       // 	// if(recmu_size>genmu_size){
       // 	//   std::cout<<"Rec Px Py Pt : "<<recmu->Px()<<" "<<recmu->Py()<<"     "<<recmu->Pt()<<" GlobalMuon "<<((recmu_SelType[recmuNb]&2) > 0)<<std::endl;	  
       // 	// }
      //}

    }
    


    // if (genBc_size>1){
    //   std::cout<<genBc_whichJpsi[0]<<" "<<genBc_whichJpsi[1]<<std::endl;
    //   //T->Scan("","","",1,i);
    // }

  }

  //****************************************************************
  //Jpsi acceptance
  //****************************************************************

  //***** Some histograms
  TH2D *h_genJpsi_YPt_all = new TH2D("genJpsiYPtAccAll","Gen J/#psi;|Rapidity|;p_{T} [GeV]",13,0,2.6,AEnbins,acceffBins);
  TH2D *h_genJpsi_YPt_acc = new TH2D("genJpsiYPtAcc","Gen J/#psi with accepted muons;|Rapidity|;p_{T} [GeV]",13,0,2.6,AEnbins,acceffBins);

  if(doJpsiAcc){  
    //****** Open the tree and make it scan branches one by one (SetMakeClass, to study one branch at a time) 
    TFile *file = TFile::Open(ispp?"root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root":"root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root","READ");
    TTree *T = (TTree*)file->Get("hionia/myTree");
    int nevents= (int)T->GetEntries();
    std::cout<<"nevents GEN-only = "<<nevents<<"\n";
    T->SetMakeClass(1);

    //****** Define adress for all branches
    int genJpsi_size;
    TBranch *b_genJpsi_size = T->GetBranch("Gen_QQ_size");
    b_genJpsi_size->SetAddress(&genJpsi_size);

    TClonesArray *Gen_QQ_4mom = new TClonesArray();
    TBranch *b_Gen_QQ_4mom; 
    T->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);

    TClonesArray *Gen_QQ_mupl_4mom = new TClonesArray();
    TBranch *b_Gen_QQ_mupl_4mom; 
    T->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);

    TClonesArray *Gen_QQ_mumi_4mom = new TClonesArray();
    TBranch *b_Gen_QQ_mumi_4mom; 
    T->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);

    //***** Loop on events
    for(int i=0;i<nevents;i++){

      if(i%200000==0){ std::cout<<"Part of tree scanned : "<< 100*(double)i / (double)nevents<<" % \n";}

      //***** Initialize pointers for arrays
      Gen_QQ_4mom->Clear();
      Gen_QQ_mupl_4mom->Clear();
      Gen_QQ_mumi_4mom->Clear();
    
      //***** Choose needed branches
      b_genJpsi_size->GetEntry(i);
      b_Gen_QQ_mupl_4mom->GetEntry(i);
      b_Gen_QQ_mumi_4mom->GetEntry(i);
      b_Gen_QQ_4mom->GetEntry(i);

      for (int genQQNb=0; genQQNb<genJpsi_size; genQQNb++){
	TLorentzVector *genQQ = (TLorentzVector*) Gen_QQ_4mom->At(genQQNb);
	h_genJpsi_YPt_all->Fill(genQQ->Rapidity(),genQQ->Pt());

	TLorentzVector *genmupl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(genQQNb);
	TLorentzVector *genmumi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(genQQNb);
      
	if(looseAcc(genmupl->Pt(),genmupl->Eta()) && looseAcc(genmumi->Pt(),genmumi->Eta())
	   ){
	  h_genJpsi_YPt_acc->Fill(genQQ->Rapidity(),genQQ->Pt());
	}
      }
      
    }
  }

  //****************************************************************
  //Drawing all
  //****************************************************************
  
  // std::cout<<"Number of reconstructed QQ: "<<NQQ_reco<<std::endl;
  // std::cout<<"Number of reconstructed QQ with accepted muons: "<<NQQ_recoAccepted<<std::endl;
  // std::cout<<"Number of QQ with accepted muons: "<<NQQ_accepted<<std::endl;

  std::cout<<"Number of generated probe muons from Jpsi, with a passing tag = "<<genmuNb_goodTag<<std::endl;
  std::cout<<"Number of reconstructed(+cuts) probe muons from Jpsi, with a passing tag = "<<recmuNb_goodTag<<std::endl;
  
  //  gStyle->SetPalette(kRainBow);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.21);
  gStyle->SetStatH(0.16);
  gStyle->SetOptStat("nime");

  // Set Palette
  const Int_t NCont = 999;
  const Int_t NRGBs = 7;
  Double_t stops[NRGBs] = { 0.00, 0.0999, 0.1, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.99, 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.86, 0.9, 0.0, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.99, 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };
  // const Int_t NRGBs = 7;
  // Double_t stops[NRGBs] = { 0.00, 0.025, 0.1, 0.34, 0.61, 0.84, 1.00 };
  // Double_t red[NRGBs]   = { 0.96, 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
  // Double_t green[NRGBs] = { 0.85, 0.0, 0.0, 0.81, 1.00, 0.20, 0.00 };
  // Double_t blue[NRGBs]  = { 0.96, 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  // TCanvas * c2 = new TCanvas("c2","pt acceptance",800,700);
  // h_genQQ_Pt->Draw("B");
  // TCanvas * c3 = new TCanvas("c3","Y acceptance",800,700);
  // h_genQQ_Y->Draw("B");

  gStyle->SetOptStat(0);
  StrToStr leg;
  leg["2global1tracker"] = "2 global + 1 tracker";
  leg["2harsh1tracker"] = "2 harsh + 1 tracker";
  leg["3global"] = "3 global";
  leg["3tracker"] = "3 tracker";

  // TCanvas * c1 = new TCanvas("c1","acceptance",800,700);
  // c1->SetRightMargin(0.115);
  // TH2D *h_genQQ_YPt_accRatio = (TH2D*)h_genQQ_YPt_acc->Clone("genQQYPtAccRatio");
  // h_genQQ_YPt_accRatio->SetTitle(("Generated B_{c} acceptance x efficiency ("+leg[QQSelectionAcc]+" muons);Gen Y;Gen p_{T} [GeV]").c_str()); //x efficiency
  // h_genQQ_YPt_accRatio->Divide(h_genQQ_YPt);
  // h_genQQ_YPt_accRatio->GetZaxis()->SetRangeUser(0,0.25);
  // h_genQQ_YPt_accRatio->Draw("COLZ");
    
  //***** Lines to draw geometrical acceptance
  TLine *l1 = new TLine(0, 3.3, 1.1, 3.3);
  TLine *l2 = new TLine(1.1, 3.3, 1.3, 1.5);
  TLine *l3 = new TLine(1.3, 1.5, 1.9, 0.8);
  TLine *l3b = new TLine(1.9, 0.8, 2.4, 0.8);
  TLine *l4 = new TLine(2.4, 0.8, 2.4, 6);
  // if(MuTypeForEfficiency=="nonglobal" || MuTypeForEfficiency=="tracker"){
  // l1->SetX1(0); l1->SetY1(3.3); l1->SetX2(0.8); l1->SetY2(3.3);
  // l2->SetX1(0.8); l2->SetY1(3.3); l2->SetX2(1.5); l2->SetY2(1.1);
  // l3->SetX1(1.5); l3->SetY1(1.1); l3->SetX2(2.07); l3->SetY2(0.8);
  // l3b->SetX1(2.07); l3b->SetY1(0.8); l3b->SetX2(2.4); l3b->SetY2(0.8);
  // l4->SetX1(2.4); l4->SetY1(6); l4->SetX2(2.4); l4->SetY2(0.8);
  // }
  l1->SetLineColor(kBlue+1);
  l1->SetLineWidth(7);
  l2->SetLineColor(kBlue+1);
  l2->SetLineWidth(7);
  l3->SetLineColor(kBlue+1);
  l3->SetLineWidth(7);
  l3b->SetLineColor(kBlue+1);
  l3b->SetLineWidth(7);
  l4->SetLineColor(kBlue+1);
  l4->SetLineWidth(7);

  TLine *l13 = new TLine(0, 3.4, 0.3, 3.4);
  TLine *l13b = new TLine(0.3, 3.4, 0.3, 3.3);
  TLine *l13c = new TLine(0.3, 3.3, 1.1, 3.3);
  TLine *l14 = new TLine(1.1, 3.3, 1.4, 2.1);
  TLine *l15 = new TLine(1.4, 2.1, 1.55, 2.1);
  TLine *l16 = new TLine(1.55, 2.1, 2.2, 1.2);
  TLine *l17 = new TLine(2.2, 1.2, 2.4, 1.2);
  TLine *l18 = new TLine(2.4, 1.2, 2.4, 6);
  l13->SetLineColor(kGreen+1);
  l13->SetLineWidth(7);
  l13b->SetLineColor(kGreen+1);
  l13b->SetLineWidth(7);
  l13c->SetLineColor(kGreen+1);
  l13c->SetLineWidth(7);
  l14->SetLineColor(kGreen+1);
  l14->SetLineWidth(7);
  l15->SetLineColor(kGreen+1);
  l15->SetLineWidth(7);
  l16->SetLineColor(kGreen+1);
  l16->SetLineWidth(7);
  l17->SetLineColor(kGreen+1);
  l17->SetLineWidth(7);
  l18->SetLineColor(kGreen+1);
  l18->SetLineWidth(7);

  TLine *l9 = new TLine(0, 3.5, 1.2, 3.5);
  TLine *l9b = new TLine(1.2, 3.5, 1.2, 3.2);
  TLine *l10 = new TLine(1.2, 3.2, 2.1, 1.5);
  TLine *l11 = new TLine(2.1, 1.5, 2.4, 1.5);
  TLine *l12 = new TLine(2.4, 1.5, 2.4, 6);
  l9->SetLineColor(kMagenta+1);
  l9->SetLineWidth(7);
  l9b->SetLineColor(kMagenta+1);
  l9b->SetLineWidth(7);
  l10->SetLineColor(kMagenta+1);
  l10->SetLineWidth(7);
  l11->SetLineColor(kMagenta+1);
  l11->SetLineWidth(7);
  l12->SetLineColor(kMagenta+1);
  l12->SetLineWidth(7);

  bool drawHarsh=false;
  TLine *l5 = new TLine(0, 3.5, 1.2, 3.5);
  TLine *l6 = new TLine(1.2, 3.5, 2.1, 1.8);
  TLine *l7 = new TLine(2.1, 1.8, 2.4, 1.8);
  TLine *l8 = new TLine(2.4, 1.8, 2.4, 6);
  if(drawHarsh){
    l5->SetLineColor(kBlack);
    l5->SetLineWidth(7);
    l6->SetLineColor(kBlack);
    l6->SetLineWidth(7);
    l7->SetLineColor(kBlack);
    l7->SetLineWidth(7);
    l8->SetLineColor(kBlack);
    l8->SetLineWidth(7);
  }

  TCanvas * c5 = new TCanvas("c5","c5",800,700);
  //h_recmu_EtaPt->Scale(ScaleToXS);
  h_recmu_EtaPt->Draw("COLZ");
  // l1->DrawClone("same");
  // l2->DrawClone("same");
  // l3->DrawClone("same");
  // l4->DrawClone("same");
  // if(drawHarsh){
  //   l5->DrawClone("same");
  //   l6->DrawClone("same");
  //   l7->DrawClone("same");
  //   l8->DrawClone("same");
  //  }

  TCanvas * c6 = new TCanvas("c6","c6",800,700);
  //  h_genmu_EtaPt->Scale(ScaleToXS);
  h_genmu_EtaPt->Draw("COLZ");

  TCanvas * c4 = new TCanvas("c4","muon efficiency",1600,1400);
  c4->SetRightMargin(0.135);
  TH2D *h_muAccEff_EtaPt = (TH2D*)h_recmu_EtaPt->Clone("muAccEffEtaPt");
  h_muAccEff_EtaPt->SetTitle("(Reco+"+(TString)((MuTypeForEfficiency=="tracker" || MuTypeForEfficiency=="nonglobal")?"Soft":"HybSoft")+"ID"+(TString)(withTrig?"+"+(TString)(L1trig?"L1 ":"")+"Trigger":"")+")/Gen muons"+(TString)((ispp)?" (pp)":" (PbPb)")+";Gen |#eta|;Gen p_{T} [GeV];efficiency");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
  h_muAccEff_EtaPt->GetZaxis()->SetTitleOffset(1.1);
  h_muAccEff_EtaPt->Divide(h_genmu_EtaPt);
  h_muAccEff_EtaPt->GetZaxis()->SetRangeUser(0,1);

  // for (int i=1;i<h_muAccEff_EtaPt->GetNbinsX();i++){
  //   for (int j=1;j<h_muAccEff_EtaPt->GetNbinsY();j++){
  //     if(h_muAccEff_EtaPt->GetBinContent(i,j) < 0.1) h_muAccEff_EtaPt->SetBinContent(i,j,0);
  //   }
  // }

  h_muAccEff_EtaPt->Draw("COLZ");

  if(drawHarsh){
    l5->DrawClone("same");
    l6->DrawClone("same");
    l7->DrawClone("same");
    l8->DrawClone("same");
  }

  if(MuTypeForEfficiency!="globaltracker" && MuTypeForEfficiency!="global"){
    l1->DrawClone("same");
    l2->DrawClone("same");
    l3->DrawClone("same");
    l3b->DrawClone("same");
    l4->DrawClone("same");
  }
  // TF1 *pcut = new TF1("pcut","5.5 * TMath::ASin( 2* TMath::ATan( TMath::Exp(-x))) ",0,3); //the normalization is the P cut // (expo(2*x)-1)/(expo(2*x)+1)
  // pcut->SetLineColor(kYellow);
  // pcut->SetLineWidth(7);
  // pcut->Draw("same");

  if(MuTypeForEfficiency!="nonglobal" && MuTypeForEfficiency!="tracker"){
    l9->DrawClone("same");
    l9b->DrawClone("same");
    l10->DrawClone("same");
    l11->DrawClone("same");
    l12->DrawClone("same");
  }
  //} else if(!withTrig){
  l13->DrawClone("same");
  l13b->DrawClone("same");
  l13c->DrawClone("same");
  l14->DrawClone("same");
  l15->DrawClone("same");
  l16->DrawClone("same");
  l17->DrawClone("same");
  l18->DrawClone("same");
  //}

  if(doTnP){
    c4->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pp":"_PbPb")+".png");
    c4->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pp":"_PbPb")+".pdf");
  } else{
    c4->SaveAs("SingleMuAcceptance_ALLMUONS_"+(TString)((MuTypeForEfficiency=="tracker")?"SoftMuId":"HybridSoftMuId")+(TString)((ispp)?"_pp":"_PbPb")+".pdf");
    c4->SaveAs("SingleMuAcceptance_ALLMUONS_"+(TString)((MuTypeForEfficiency=="tracker")?"SoftMuId":"HybridSoftMuId")+(TString)((ispp)?"_pp":"_PbPb")+".png");
  }

  if(!ispp){
    TCanvas * c7 = new TCanvas("c7","muon efficiency",1600,1400);
    c7->SetRightMargin(0.135);
    TH2D *h_muAccEff_EtaPt020 = (TH2D*)h_recmu_EtaPt020->Clone("muAccEffEtaPt0-20");
    h_muAccEff_EtaPt020->SetTitle("(Reco+"+(TString)((MuTypeForEfficiency=="tracker" || MuTypeForEfficiency=="nonglobal")?"Soft":"HybSoft")+"ID"+(TString)(withTrig?"+"+(TString)(L1trig?"L1 ":"")+"Trigger":"")+")/Gen muons"+(TString)((ispp)?" (pp)":" (PbPb, cent [0,20%])")+";Gen |#eta|;Gen p_{T} [GeV];efficiency");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
    h_muAccEff_EtaPt020->GetZaxis()->SetTitleOffset(1.1);
    h_muAccEff_EtaPt020->Divide(h_genmu_EtaPt020);
    h_muAccEff_EtaPt020->GetZaxis()->SetRangeUser(0,1);
    h_muAccEff_EtaPt020->Draw("COLZ");
    if(doTnP){
      c7->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pp":"_PbPbCentrality020")+".pdf");
    } else{
      c7->SaveAs("SingleMuAcceptance_ALLMUONS_"+(TString)((MuTypeForEfficiency=="tracker")?"SoftMuId":"HybridSoftMuId")+(TString)((ispp)?"_pp":"_PbPbCentrality020")+".pdf");
    }

    TCanvas * c9 = new TCanvas("c9","muon efficiency",1600,1400);
    c9->SetRightMargin(0.135);
    TH2D *h_muAccEff_EtaPt2050 = (TH2D*)h_recmu_EtaPt2050->Clone("muAccEffEtaPt20-50");
    h_muAccEff_EtaPt2050->SetTitle("(Reco+"+(TString)((MuTypeForEfficiency=="tracker" || MuTypeForEfficiency=="nonglobal")?"Soft":"HybSoft")+"ID"+(TString)(withTrig?"+"+(TString)(L1trig?"L1 ":"")+"Trigger":"")+")/Gen muons"+(TString)((ispp)?" (pp)":" (PbPb, cent [20,50%])")+";Gen |#eta|;Gen p_{T} [GeV];efficiency");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
    h_muAccEff_EtaPt2050->GetZaxis()->SetTitleOffset(1.1);
    h_muAccEff_EtaPt2050->Divide(h_genmu_EtaPt2050);
    h_muAccEff_EtaPt2050->GetZaxis()->SetRangeUser(0,1);
    h_muAccEff_EtaPt2050->Draw("COLZ");
    if(doTnP){
      c9->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pp":"_PbPbCentrality2050")+".pdf");
    } else{
      c9->SaveAs("SingleMuAcceptance_ALLMUONS_"+(TString)((MuTypeForEfficiency=="tracker")?"SoftMuId":"HybridSoftMuId")+(TString)((ispp)?"_pp":"_PbPbCentrality2050")+".pdf");
    }

    TCanvas * c14 = new TCanvas("c14","muon efficiency",1600,1400);
    c14->SetRightMargin(0.135);
    TH2D *h_muAccEff_EtaPt50100 = (TH2D*)h_recmu_EtaPt50100->Clone("muAccEffEtaPt50-100");
    h_muAccEff_EtaPt50100->SetTitle("(Reco+"+(TString)((MuTypeForEfficiency=="tracker" || MuTypeForEfficiency=="nonglobal")?"Soft":"HybSoft")+"ID"+(TString)(withTrig?"+"+(TString)(L1trig?"L1 ":"")+"Trigger":"")+")/Gen muons"+(TString)((ispp)?" (pp)":" (PbPb, cent [50,100%])")+";Gen |#eta|;Gen p_{T} [GeV];efficiency");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
    h_muAccEff_EtaPt50100->GetZaxis()->SetTitleOffset(1.1);
    h_muAccEff_EtaPt50100->Divide(h_genmu_EtaPt50100);
    h_muAccEff_EtaPt50100->GetZaxis()->SetRangeUser(0,1);
    h_muAccEff_EtaPt50100->Draw("COLZ");
    if(doTnP){
      c14->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pp":"_PbPbCentrality50100")+".pdf");
    } else{
      c14->SaveAs("SingleMuAcceptance_ALLMUONS_"+(TString)((MuTypeForEfficiency=="tracker")?"SoftMuId":"HybridSoftMuId")+(TString)((ispp)?"_pp":"_PbPbCentrality50100")+".pdf");
    }

    TCanvas * c15 = new TCanvas("c15","muon efficiency",1600,1400);
    c15->SetRightMargin(0.135);
    TH2D *h_muAccEff_CentPt = (TH2D*)h_recmu_CentPt->Clone("muAccEffCentPt");
    h_muAccEff_CentPt->SetTitle("(Reco+"+(TString)((MuTypeForEfficiency=="tracker" || MuTypeForEfficiency=="nonglobal")?"Soft":"HybSoft")+"ID"+(TString)(withTrig?"+"+(TString)(L1trig?"L1 ":"")+"Trigger":"")+")/Gen muons"+(TString)((ispp)?" (pp)":" (PbPb)")+";Centrality;Gen p_{T} [GeV];efficiency");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
    h_muAccEff_CentPt->GetZaxis()->SetTitleOffset(1.1);
    h_muAccEff_CentPt->Divide(h_genmu_CentPt);
    h_muAccEff_CentPt->GetZaxis()->SetRangeUser(0,1);
    h_muAccEff_CentPt->Draw("COLZ");
    if(doTnP){
      c15->SaveAs("SingleMuAcceptance_"+(TString)(withTrig?(L1trig?"L1step":"L3step"):"HybridSoftMuId")+"_TagPassedL3Mu3"+(TString)((ispp)?"_pp":"_PbPbCentralityMap")+".pdf");
    } else{
      c15->SaveAs("SingleMuAcceptance_ALLMUONS_"+(TString)((MuTypeForEfficiency=="tracker")?"SoftMuId":"HybridSoftMuId")+(TString)((ispp)?"_pp":"_PbPbCentralityMap")+".pdf");
    }
  }

  // gStyle->SetOptStat(0);
  if(doJpsiEff){
    TCanvas * c8 = new TCanvas("c8","Jpsi Pt Y efficiency",1600,1400);
    c8->SetRightMargin(0.115);
    TH1D *h_JpsiEff_YPt = (TH1D*)h_recJpsi_YPt->Clone("JpsiEffYPt");
    h_JpsiEff_YPt->SetTitle("Efficiency of J/#psi's with accepted muons;Gen |Rapidity|; Gen J/#psi p_{T} [GeV]");// ("Reconstructed/Generated for muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
    h_JpsiEff_YPt->Divide(h_genJpsi_YPt);
    h_JpsiEff_YPt->GetZaxis()->SetRangeUser(0,1);
    h_JpsiEff_YPt->Draw("COLZ");
    c8->SaveAs("JpsiEfficiency_NewMuonKinCuts.png");
    c8->SaveAs("JpsiEfficiency_NewMuonKinCuts.pdf");

    if(doJpsiAcc){
      TCanvas * c12 = new TCanvas("c12","Jpsi Pt Y acceptance",1600,1400);
      c12->SetRightMargin(0.115);
      TH1D *h_JpsiAcc_YPt = (TH1D*)h_genJpsi_YPt_acc->Clone("JpsiAccYPt");
      h_JpsiAcc_YPt->SetTitle("Acceptance of generated J/#psi's;Gen |Rapidity|; Gen J/#psi p_{T} [GeV]");
      h_JpsiAcc_YPt->Divide(h_genJpsi_YPt_all);
      h_JpsiAcc_YPt->GetZaxis()->SetRangeUser(0,1);
      h_JpsiAcc_YPt->Draw("COLZ");
      c12->SaveAs("JpsiAcceptance_NewMuonKinCuts.png");
      c12->SaveAs("JpsiAcceptance_NewMuonKinCuts.pdf");

      TH1D *h_JpsiAccEff_YPt = (TH1D*)h_JpsiAcc_YPt->Clone("JpsiAccEffYPt");
      TCanvas * c13 = new TCanvas("c13","Jpsi Pt Y acceptance",1600,1400);
      h_JpsiAccEff_YPt->SetTitle("Acceptance #times Efficiency for generated J/#psi's;Gen |Rapidity|; Gen J/#psi p_{T} [GeV]");
      h_JpsiAccEff_YPt->Multiply(h_JpsiEff_YPt);
      h_JpsiAccEff_YPt->GetZaxis()->SetRangeUser(0,1);
      h_JpsiAccEff_YPt->Draw("COLZ");
      c13->SaveAs("JpsiAccEff_NewMuonKinCuts.png");
      c13->SaveAs("JpsiAccEff_NewMuonKinCuts.pdf");

    }
  }

  // TCanvas * c8 = new TCanvas("c8","Jpsi Pt efficiency",1600,1400);
  // c8->SetRightMargin(0.115);
  // TH1D *h_JpsiAccEff1_Pt = (TH1D*)h_recJpsi_Pt->Clone("JpsiAccEffPt");
  // h_JpsiAccEff1_Pt->SetTitle("Efficiency of good J/#psi's to pass dR3p5 trigger;Gen J/#psi p_{T} [GeV]");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
  // h_JpsiAccEff1_Pt->Divide(h_genJpsi_Pt);
  // h_JpsiAccEff1_Pt->GetZaxis()->SetRangeUser(0,1);
  // h_JpsiAccEff1_Pt->SetLineWidth(3);
  // h_JpsiAccEff1_Pt->Draw("hist");

  // TH1D *h_JpsiAccEff2_Pt = (TH1D*)h_rec2Jpsi_Pt->Clone("JpsiAccEff2Pt");
  // h_JpsiAccEff2_Pt->SetTitle("Efficiency of J/#psi passing all cuts + old/new muon acceptance;Gen J/#psi p_{T} [GeV]");// ("Reconstructed/Generated for "+MuTypeForEfficiency+" muons;Gen |#eta|;Gen p_{T} [GeV]"+"").c_str() 
  // h_JpsiAccEff2_Pt->Divide(h_genJpsi_Pt);
  // h_JpsiAccEff2_Pt->GetZaxis()->SetRangeUser(0,1);
  // h_JpsiAccEff2_Pt->SetLineWidth(3);
  // h_JpsiAccEff2_Pt->SetLineColor(kRed);
  // h_JpsiAccEff2_Pt->Draw("histsame");

  // auto legend = new TLegend(0.1,0.7,0.48,0.9);
  // legend->AddEntry(h_JpsiAccEff1_Pt,"old muon acceptance","l");
  // legend->AddEntry(h_JpsiAccEff2_Pt,"new muon acceptance","l");
  // legend->Draw();
  //  c8->SaveAs("JpsiEfficiency_dR3p5.png");

  // for(int bini=1;bini<=muptbins;bini++){
  //   for(int binj=1;binj<=muetabins;binj++){
  // 	if(h_muAccEff_EtaPt->GetBinContent(bini,binj)>1){
  // 	  std::cout<<"bin "<<bini<<", "<<binj<<" contains ratio >1"<<std::endl;
  // 	  std::cout<<"Number of reco muons in this bin : "<<h_recmu_EtaPt->GetBinContent(bini,binj)<<std::endl;
  // 	  std::cout<<"Number of gen muons in this bin : "<<h_genmu_EtaPt->GetBinContent(bini,binj)<<std::endl;
  // 	}
  //   }
  // }

  // TFile out_file("SingleMuonEfficiency_new.root","NEW");
  // h_muAccEff_EtaPt->Write();

}
