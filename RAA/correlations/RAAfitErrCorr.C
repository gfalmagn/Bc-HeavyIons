#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "TRandom.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include <Math/GSLRndmEngines.h>
#include "../../helpers/Definitions.h"
#include "../../helpers/Cuts.h"

void RAAErrCorr(bool fromfit=true, bool fullerr=false, bool centDep=false){

  if(fromfit) fullerr=false;
  if(fullerr) fromfit=false;

  vector<float> *AEcorr12_pp;
  vector<float> *AEcorr12_PbPb;
  vector<vector<double> > *InvAccEff_pp;
  vector<vector<double> > *InvAccEff_PbPb;
  TFile *infile2 = new TFile("../../twoSteps/AccEffFrom2ndStepToys.root","READ");
  infile2->GetObject((TString)(fullerr?"CorrYields":"InvAccEff")+"FromCorrMC_LinearisedCorrelationFactor_pp_2ndStep", AEcorr12_pp);
  infile2->GetObject((TString)(fullerr?"CorrYields":"InvAccEff")+"FromCorrMC_LinearisedCorrelationFactor_PbPb"+(TString)(centDep?"_inCentBins":"")+"_2ndStep", AEcorr12_PbPb);
  infile2->GetObject((TString)(fullerr?"CorrYields":"InvAccEff")+"FromCorrMC_withSystErr_pp_2ndStep", InvAccEff_pp);
  infile2->GetObject((TString)(fullerr?"CorrYields":"InvAccEff")+"FromCorrMC_withSystErr_PbPb"+(TString)(centDep?"_inCentBins":"")+"_2ndStep", InvAccEff_PbPb);

  vector<float> *y_nom_pp;
  vector<vector<float> > *y_fitErr_pp;
  vector<float> *r1r2Corr_pp;
  vector<float> *y_nom_PbPb;
  vector<vector<float> > *y_fitErr_PbPb;
  vector<float> *r1r2Corr_PbPb;

  TFile *infile = new TFile("../../AccEffCorr/corrected_yields_3rdStep.root","READ");
  infile->GetObject("FinalCorrectedYield_pp", y_nom_pp);
  infile->GetObject("FinalCorrectedYield_fitError_pp", y_fitErr_pp);
  infile->GetObject("FinalCorrectedYield"+(TString)(centDep?"_centralityDep":"")+"_PbPb", y_nom_PbPb);
  infile->GetObject("FinalCorrectedYield_fitError"+(TString)(centDep?"_centralityDep":"")+"_PbPb", y_fitErr_PbPb);

  TFile *infile3 = new TFile("../../AccEffCorr/corrected_yields_2ndStep.root","READ");
  infile3->GetObject("r1r2Correlation_pp", r1r2Corr_pp);
  infile3->GetObject("r1r2Correlation"+(TString)(centDep?"_centralityDep":"")+"_PbPb", r1r2Corr_PbPb);

  float rho_pp = centDep?1.:( fromfit?((*r1r2Corr_pp)[0]):((*AEcorr12_pp)[0]) );
  float rho_PbPb = fromfit?((*r1r2Corr_PbPb)[0]):((*AEcorr12_PbPb)[0]);
  cout<<"corr factor (pT bins) pp = "<<rho_pp<<endl;
  cout<<"corr factor ("<<(TString)(centDep?"centrality":"pT")<<" bins) PbPb = "<<rho_PbPb<<endl;

  vector<float> yErr_pp;
  if(centDep){
    yErr_pp.push_back(fromfit?(*y_fitErr_pp)[0][0]:(
						    (float)(*InvAccEff_pp)[0][1] * (*y_nom_pp)[0] / (*InvAccEff_pp)[0][0])); //error of AccEff on y_nom
  } else {
    yErr_pp.push_back(fromfit?(*y_fitErr_pp)[1][0]:(
						    (float)(*InvAccEff_pp)[1][1] * (*y_nom_pp)[1] / (*InvAccEff_pp)[1][0])); //error of AccEff on y_nom
    yErr_pp.push_back(fromfit?(*y_fitErr_pp)[2][0]:(
						    (float)(*InvAccEff_pp)[2][1] * (*y_nom_pp)[2] / (*InvAccEff_pp)[2][0])); //error of AccEff on y_nom
  }

  vector<float> yErr_PbPb;
  yErr_PbPb.push_back(fromfit?(*y_fitErr_PbPb)[1][0]:(
						  (float)(*InvAccEff_PbPb)[1][1] * (*y_nom_PbPb)[1] / (*InvAccEff_PbPb)[1][0])); //error of AccEff on y_nom
  yErr_PbPb.push_back(fromfit?(*y_fitErr_PbPb)[2][0]:(
						  (float)(*InvAccEff_PbPb)[2][1] * (*y_nom_PbPb)[2] / (*InvAccEff_PbPb)[2][0])); //error of AccEff on y_nom
						    

  const int n = 200000;
  double X1_pp,X2_pp,X1_PbPb,X2_PbPb,X1_RAA,X2_RAA;
 
  // ROOT::Math::GSLRandomEngine rengine =  ROOT::Math::GSLRandomEngine();

  float normpp0 = 1.;//L_pp * (-_BcPtmin[0]+_BcPtmax[0]);
  float normPbPb0 = 1.;//NMB_PbPb * TAA_090 * (-_BcPtmin[0]+_BcPtmax[0]);
  float normpp1 = centDep?normpp0:1.;//L_pp * (-_BcPtmin[1]+_BcPtmax[1]);
  float normPbPb1 = 1.;//NMB_PbPb * TAA_090 * (-_BcPtmin[1]+_BcPtmax[1]);

  float norm2pp0 = L_pp * (-_BcPtmin[centDep?0:1]+_BcPtmax[centDep?0:1]) * (-_BcYmin[centDep?0:1]+_BcYmax[centDep?0:1]);
  float norm2PbPb0 = NMB_PbPb * TAA_090 * (-_BcPtmin[1]+_BcPtmax[1]) * (-_BcYmin[1]+_BcYmax[1]) * (_Centmax[0]/100 - _Centmin[0]/100);
  float norm2pp1 = L_pp * (-_BcPtmin[centDep?0:2]+_BcPtmax[centDep?0:2]) * (-_BcYmin[centDep?0:2]+_BcYmax[centDep?0:2]);
  float norm2PbPb1 = NMB_PbPb * TAA_090 * (-_BcPtmin[2]+_BcPtmax[2]) * (-_BcYmin[2]+_BcYmax[2]) * (_Centmax[0]/100 - _Centmin[0]/100);
  if(centDep){
    norm2PbPb0 = NMB_PbPb * TAA_020 * (-_BcPtmin[0]+_BcPtmax[0]) * (_BcYmax[0]-_BcYmin[0]) * (_Centmax[1]/100 - _Centmin[1]/100);
    norm2PbPb1 = NMB_PbPb * TAA_2090 * (-_BcPtmin[0]+_BcPtmax[0]) * (_BcYmax[0]-_BcYmin[0]) * (_Centmax[2]/100 - _Centmin[2]/100);
  }
  
  cout<<"pp nominal: "<<(*y_nom_pp)[centDep?0:1]/normpp0<<" "<<(*y_nom_pp)[centDep?0:2]/normpp1<<endl;
  cout<<"pp error: "<<yErr_pp[0]/normpp0<<" "<<yErr_pp[centDep?0:1]/normpp1<<endl;
  cout<<"PbPb nominal: "<<(*y_nom_PbPb)[1]/normPbPb0<<" "<<(*y_nom_PbPb)[2]/normPbPb1<<endl;
  cout<<"PbPb error: "<<yErr_PbPb[0]/normPbPb0<<" "<<yErr_PbPb[1]/normPbPb1<<endl;

  TH2D* h_pp = new TH2D("h_pp","pp",300,(*y_nom_pp)[centDep?0:1]/normpp0 - 4*yErr_pp[0]/normpp0,(*y_nom_pp)[centDep?0:1]/normpp0 + 4*yErr_pp[0]/normpp0,300,(*y_nom_pp)[centDep?0:2]/normpp1 - 4*yErr_pp[centDep?0:1]/normpp1,(*y_nom_pp)[centDep?0:2]/normpp1 + 4*yErr_pp[centDep?0:1]/normpp1);
  TH2D* h_PbPb = new TH2D("h_PbPb","PbPb",300,max((float)0.,(*y_nom_PbPb)[1]/normPbPb0 - 4*yErr_PbPb[0]/normPbPb0),(*y_nom_PbPb)[1]/normPbPb0 + 4*yErr_PbPb[0]/normPbPb0,300,max((float)0.,(*y_nom_PbPb)[2]/normPbPb1 - 4*yErr_PbPb[1]/normPbPb1),(*y_nom_PbPb)[2]/normPbPb1 + 4*yErr_PbPb[1]/normPbPb1);
  TH2D* h_RAA = new TH2D("h_RAA","RAA",300,0.,3.,300,0.,1.3);
  
  for(int i=0;i<n;i++){
    //    rengine.Gaussian2D((*y_fitErr_pp)[0],(*y_fitErr_pp)[1],rho_pp,&(x1_pp[i]),&(x2_pp[i]));
    //    rengine.Gaussian2D((*y_fitErr_PbPb)[0],(*y_fitErr_PbPb)[1],rho_PbPb,&(x1_PbPb[i]),&(x2_PbPb[i]));

    double z1 = gRandom->Gaus(),z2 = gRandom->Gaus(),z3 = gRandom->Gaus(),z4 = gRandom->Gaus();
    X1_pp = (yErr_pp[0] * z1 + (*y_nom_pp)[centDep?0:1])/normpp0; //correlate the X1 and X2 with rho_pp
    X2_pp = (yErr_pp[centDep?0:1] * ( rho_pp*z1+sqrt(1-pow(rho_pp,2))*z2 ) + (*y_nom_pp)[centDep?0:2])/normpp1;
    X1_PbPb = (yErr_PbPb[0] * z3 + (*y_nom_PbPb)[1])/normPbPb0;
    X2_PbPb = (yErr_PbPb[1] * ( rho_PbPb*z3+sqrt(1-pow(rho_PbPb,2))*z4 ) + (*y_nom_PbPb)[2])/normPbPb1;
    X1_PbPb = max(X1_PbPb,0.);
    X2_PbPb = max(X2_PbPb,0.);
    X1_RAA = X1_PbPb*norm2pp0/(X1_pp*norm2PbPb0) ;
    X2_RAA = X2_PbPb*norm2pp1/(X2_pp*norm2PbPb1);
    h_pp->Fill(X1_pp,X2_pp);
    h_PbPb->Fill(X1_PbPb,X2_PbPb);
    h_RAA->Fill(X1_RAA,X2_RAA);
    // cout<<"X1_pp z1 X1_PbPb z2 = "<<X1_pp<<" "<<z1<<" "<<X1_PbPb<<" "<<z2<<endl;
    // cout<<"X2_pp z3 X2_PbPb z4 = "<<X2_pp<<" "<<z3<<" "<<X2_PbPb<<" "<<z4<<endl;
  }

  cout<<"h_pp->GetCorrelationFactor() = "<<h_pp->GetCorrelationFactor()<<endl;
  cout<<"h_PbPb->GetCorrelationFactor() = "<<h_PbPb->GetCorrelationFactor()<<endl;
  cout<<"h_RAA->GetCorrelationFactor() = "<<h_RAA->GetCorrelationFactor()<<endl;
  vector<float> ppcorr;
  vector<float> PbPbcorr;
  vector<float> RAAcorr;
  ppcorr.push_back(rho_pp);
  PbPbcorr.push_back(rho_PbPb);
  RAAcorr.push_back(h_RAA->GetCorrelationFactor());

  TCanvas *c1 = new TCanvas("c1","c1",3000,1000);
  c1->Divide(3,1);
  c1->cd(1)->SetLeftMargin(0.15);
  h_pp->GetXaxis()->SetTitle("corrected yield (bin 1)");
  h_pp->GetYaxis()->SetTitle("corrected yield (bin 2)");
  h_pp->Draw("COLZ");
  c1->cd(2)->SetLeftMargin(0.15);
  h_PbPb->GetXaxis()->SetTitle("corrected yield (bin 1)");
  h_PbPb->GetYaxis()->SetTitle("corrected yield (bin 2)");
  h_PbPb->Draw("COLZ");
  c1->cd(3)->SetLeftMargin(0.15);
  h_RAA->GetXaxis()->SetTitle("R_{PbPb} (bin 1)");
  h_RAA->GetYaxis()->SetTitle("R_{PbPb} (bin 2)");
  h_RAA->Draw("COLZ");

  c1->SaveAs((TString)(fullerr?"Full":(fromfit?"Fit":"Acceff"))+"ErrorCorrelation"+(TString)(centDep?"_centralityBins":"")+".pdf");
  
  TFile * outf = new TFile("../../AccEffCorr/corrected_yields_3rdStep.root","UPDATE");
  if(!centDep) 
    outf->WriteObject(&ppcorr,(TString)(fromfit?"r1r2":"AcceffSyst_")+"Correlation_pp");
  cout<<"writing out "<<((TString)(fullerr?"FullError":(fromfit?"r1r2":"AcceffSyst_"))+"Correlation_PbPb"+(TString)(centDep?"_centralityBins":""))<<endl;
  outf->WriteObject(&PbPbcorr,(TString)(fullerr?"FullError":(fromfit?"r1r2":"AcceffSyst_"))+"Correlation_PbPb"+(TString)(centDep?"_centralityBins":""));
  outf->WriteObject(&RAAcorr,(TString)(fullerr?"FullError":(fromfit?"r1r2":"AcceffSyst_"))+"Correlation_RAA"+(TString)(centDep?"_centralityBins":""));
  outf->Close();

}

void RAAfitErrCorr(){
  cout<<"\ntotal\n"<<endl;
  RAAErrCorr(false,true,false);
  cout<<"\nfit\n"<<endl;
  RAAErrCorr(true,false,false);
  cout<<"\nacceptance*efficiency\n"<<endl;
  RAAErrCorr(false,false,false);
  cout<<"\nfit centrality-dep\n"<<endl;
  RAAErrCorr(true,false,true);
  cout<<"\nacceptance*efficiency centrality-dep\n"<<endl;
  RAAErrCorr(false,false,true);
}
