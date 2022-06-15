#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TRandom.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/hub.cpp"

void BCVEGPYppXS(){

  gStyle->SetOptStat(0);
  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);

  //Create Hub gathering various data 
  Hub H = Hub(true,true);
  H.SetpTbias();

  TFile *facc = TFile::Open(H.fAccName,"READ");
  TTree* T = (TTree*)facc->Get("hionia/myTree");
  int nentries = T->GetEntries();

  TBranch *b_Gen_3mu_4mom;
  TClonesArray *Gen_3mu_4mom = new TClonesArray();
  T->SetBranchAddress("Gen_3mu_4mom", &Gen_3mu_4mom, &b_Gen_3mu_4mom);
  
  TBranch *b_Gen_Bc_4mom;
  TClonesArray *Gen_Bc_4mom = new TClonesArray();
  T->SetBranchAddress("Gen_Bc_4mom", &Gen_Bc_4mom, &b_Gen_Bc_4mom);

  float npass[_NanaBins+1];
  float fpass[_NanaBins+1];
  float npassCor[_NanaBins+1];
  float ntotCor=0;
  float fpassCor[_NanaBins+1];

  for(int b=0;b<=_NanaBins;b++){
    npass[b] = 0;
    npassCor[b] = 0;
  }

  for(int j=0; j<nentries; j++){
    Gen_Bc_4mom->Clear();
    Gen_3mu_4mom->Clear();
    
    b_Gen_3mu_4mom->GetEntry(j);
    b_Gen_Bc_4mom->GetEntry(j);

    TLorentzVector *genBc = (TLorentzVector*) Gen_Bc_4mom->At(0); //assuming only 1 gen Bc per event, which is true with 4M events
    TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom->At(0);
    
    float wei = getBias( H.pTbias_step1[0][_nomMethVar], gen3mu->Pt()) * getBias( H.pTbias_step2[0][_nomMethVar], gen3mu->Pt()); //pT spectrum correction from our measurement
    ntotCor += wei;

    for(int b=0;b<=_NanaBins;b++){
      //count uncut events passing analysis bin cuts, with or w/o pT spectrum correction from our measurement
      if(inFidCuts(b, gen3mu->Pt(), gen3mu->Rapidity())){
	npass[b] += 1.;
	npassCor[b] += wei;
      }
    }

  }

  for(int b=0;b<=_NanaBins;b++){
    fpassCor[b] = npassCor[b]/ntotCor;
    fpass[b] = npass[b]/nentries;
  }

  cout<<"fraction of uncut MC in pT bin0, 1 and 2 = "<<fpass[0]<<" "<<fpass[1]<<" "<<fpass[2]<<endl;
  cout<<"fraction of uncut MC, with pT spectrum correction, in pT bin0, 1 and 2 = "<<fpassCor[0]<<" "<<fpassCor[1]<<" "<<fpassCor[2]<<endl;
  cout<<"BCVEGPY pT bin1/bin2 = "<<(fpass[1]/H.YBinWidth[1] / H.pTBinWidth[1]) / (fpass[2]/H.YBinWidth[2] / H.pTBinWidth[2])<<endl;

  //
  //Input from other measurements and BCVEGPY (BVP) total XS
  //

  float BVPXSmin = 2*79.9; //[120 nb] //only 1s0 and 3s1 states
  float BVPXSmax = 2*117.7; //[240 nb] //all ten states simulated in BCVEGPY, with 100% branching to Bc+X
  
  float BFJpsiMuMu = 0.0593; //PDG, low uncertainty
  vector<float> BFlept = {6.4, 1.4, 7.5, 1.9, 2.3, 2.7, 1.6, 1.7, 1.7, 1.9, 2.3, 2.2, 2.6, 2.5, 1.3, 1.4, 1.5, 1.9, 2.2, 2.91/0.469, 1.3/0.469, 2.1/0.469, 0.73/0.469, 1.3/0.469, 0.61/0.469, 1.1/0.469, 1.7/0.469, 2.0/0.469}; //[in %] //from reviews in https://arxiv.org/pdf/1910.13404.pdf Table 2, and https://journals-aps-org.ezproxy.cern.ch/prd/pdf/10.1103/PhysRevD.89.034008 + BF(hadr)/BF(lept) that has 14% uncertainty
  for(int ibf=0;ibf<BFlept.size();ibf++) BFlept[ibf] *= 1e-2;

  float XSBlhcb = 72000 * 0.688 / 0.245; //[202 microb] 
  //   XS(b-quark meson) in 2<eta<5 (LHCb https://arxiv.org/pdf/1612.05140.pdf)  , +-9.5%
  // * 7->5.02 TeV with BCVEGPY (confirmed by LHCb paper above), small error
  // / fraction of uncut BCVEGPY in 2<eta<5, unknown (large) error  
  float XStimesBFlhcb = XSBlhcb * BFJpsiMuMu * 5.0e-5; //[10.1 nb * BFJpsiMuMu] // Proba(b-quark -> Bc) * leptonic BF, from LHCb https://arxiv.org/pdf/1910.13404.pdf, +-5%
  
  float BCVPcorrTot = fpassCor[0]/fpass[0]; //about 1.5, to correct the fraction of BCVP that goes in the CMS phase space, with the measured pT spectrum
  vector<float> XStimesBFbcvp;
  for(int ibf=0;ibf<BFlept.size();ibf++){
    //consider here a 100% bf of Bc excited states to the ground state + correct for our measured pT distro
    XStimesBFbcvp.push_back( BVPXSmax * BCVPcorrTot * BFJpsiMuMu * BFlept[ibf] ); 
  }
  float XStimesBFlhcbCor = XStimesBFlhcb * BCVPcorrTot; //correct the fractions of BCVEGPY in bin 1 or 2 (to be applied later) for our measured pT distro

  double ptbins[] = {_BcPtmin[1],_BcPtmax[1],_BcPtmax[2]};
  TH2F *diffXStBFbcvp = new TH2F("diffXStBFbcvp","diffXStBFbcvp",_NanaBins,ptbins,100000,0.01,30); // [pb/GeV]

  for(int b=1;b<=_NanaBins;b++){
    for(int ibf=0;ibf<BFlept.size();ibf++){
      diffXStBFbcvp->Fill( _BcPtmin[b]+0.1, 
			   fpass[b] * XStimesBFbcvp[ibf] / H.YBinWidth[b] / H.pTBinWidth[b] * 1000); //go to pb/GeV //multiply by fraction of BCVEGPY in this bin
    }
  }

  diffXStBFbcvp->SetLineColor(kRed+3);

  // TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
  // c1->SetLogy();
  // diffXStBFbcvp->Draw("CANDLEX1");

  TFile *fout = new TFile("BCVEGPYppXS.root","recreate");
  diffXStBFbcvp->Write();

}
