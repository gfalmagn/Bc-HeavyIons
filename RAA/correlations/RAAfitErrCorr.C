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

void RAAfitErrCorr(){

  vector<float> *y_nom_pp;
  vector<float> *y_fitErr_pp;
  vector<float> *r1r2Corr_pp;
  vector<float> *y_nom_PbPb;
  vector<float> *y_fitErr_PbPb;
  vector<float> *r1r2Corr_PbPb;

  TFile *infile = new TFile("../../AccEffCorr/corrected_yields.root","READ");
  infile->GetObject("FinalCorrectedYield_pp", y_nom_pp);
  infile->GetObject("FinalCorrectedYield_fitError_pp", y_fitErr_pp);
  infile->GetObject("r1r2Correlation_pp", r1r2Corr_pp);
  infile->GetObject("FinalCorrectedYield_PbPb", y_nom_PbPb);
  infile->GetObject("FinalCorrectedYield_fitError_PbPb", y_fitErr_PbPb);
  infile->GetObject("r1r2Correlation_PbPb", r1r2Corr_PbPb);

  float rho_pp = (*r1r2Corr_pp)[0];
  float rho_PbPb = (*r1r2Corr_PbPb)[0];
  cout<<"corr factor pp = "<<rho_pp<<endl;
  cout<<"corr factor PbPb = "<<rho_PbPb<<endl;

  const int n = 200000;
  double X1_pp,X2_pp,X1_PbPb,X2_PbPb,X1_RAA,X2_RAA;
 
  // ROOT::Math::GSLRandomEngine rengine =  ROOT::Math::GSLRandomEngine();

  float normpp0 = L_pp * (-_BcPtmin[0]+_BcPtmax[0]);
  float normPbPb0 = NMB_PbPb * TAA_090 * (-_BcPtmin[0]+_BcPtmax[0]);
  float normpp1 = L_pp * (-_BcPtmin[1]+_BcPtmax[1]);
  float normPbPb1 = NMB_PbPb * TAA_090 * (-_BcPtmin[1]+_BcPtmax[1]);
  
  cout<<"pp nominal: "<<(*y_nom_pp)[0]/normpp0<<" "<<(*y_nom_pp)[1]/normpp1<<endl;
  cout<<"pp erro: "<<(*y_fitErr_pp)[0]/normpp0<<" "<<(*y_fitErr_pp)[1]/normpp1<<endl;
  cout<<"PbPb nominal: "<<(*y_nom_PbPb)[0]/normPbPb0<<" "<<(*y_nom_PbPb)[1]/normPbPb1<<endl;
  cout<<"PbPb erro: "<<(*y_fitErr_PbPb)[0]/normPbPb0<<" "<<(*y_fitErr_PbPb)[1]/normPbPb1<<endl;

  TH2D* h_pp = new TH2D("h_pp","h_pp",300,(*y_nom_pp)[0]/normpp0 - 4*(*y_fitErr_pp)[0]/normpp0,(*y_nom_pp)[0]/normpp0 + 4*(*y_fitErr_pp)[0]/normpp0,300,(*y_nom_pp)[1]/normpp1 - 4*(*y_fitErr_pp)[1]/normpp1,(*y_nom_pp)[1]/normpp1 + 4*(*y_fitErr_pp)[1]/normpp1);
  TH2D* h_PbPb = new TH2D("h_PbPb","h_PbPb",300,(*y_nom_PbPb)[0]/normPbPb0 - 4*(*y_fitErr_PbPb)[0]/normPbPb0,(*y_nom_PbPb)[0]/normPbPb0 + 4*(*y_fitErr_PbPb)[0]/normPbPb0,300,(*y_nom_PbPb)[1]/normPbPb1 - 4*(*y_fitErr_PbPb)[1]/normPbPb1,(*y_nom_PbPb)[1]/normPbPb1 + 4*(*y_fitErr_PbPb)[1]/normPbPb1);
  TH2D* h_RAA = new TH2D("h_RAA","h_RAA",300,0.,3,300,0.2,1);
  
  for(int i=0;i<n;i++){
    //    rengine.Gaussian2D((*y_fitErr_pp)[0],(*y_fitErr_pp)[1],rho_pp,&(x1_pp[i]),&(x2_pp[i]));
    //    rengine.Gaussian2D((*y_fitErr_PbPb)[0],(*y_fitErr_PbPb)[1],rho_PbPb,&(x1_PbPb[i]),&(x2_PbPb[i]));

    double z1 = gRandom->Gaus(),z2 = gRandom->Gaus(),z3 = gRandom->Gaus(),z4 = gRandom->Gaus();

    X1_pp = ((*y_fitErr_pp)[0] * z1 + (*y_nom_pp)[0])/normpp0; //correlate the X1 and X2 with rho_pp
    X2_pp = ((*y_fitErr_pp)[1] * ( rho_pp*z1+sqrt(1-pow(rho_pp,2))*z2 ) + (*y_nom_pp)[1])/normpp1;
    X1_PbPb = ((*y_fitErr_PbPb)[0] * z3 + (*y_nom_PbPb)[0])/normPbPb0;
    X2_PbPb = ((*y_fitErr_PbPb)[1] * ( rho_PbPb*z3+sqrt(1-pow(rho_PbPb,2))*z4 ) + (*y_nom_PbPb)[1])/normPbPb1;
    X1_RAA = X1_PbPb/X1_pp;
    X2_RAA = X2_PbPb/X2_pp;
    h_pp->Fill(X1_pp,X2_pp);
    h_PbPb->Fill(X1_PbPb,X2_PbPb);
    h_RAA->Fill(X1_RAA,X2_RAA);
    // cout<<"X1_pp z1 X1_PbPb z2 = "<<X1_pp<<" "<<z1<<" "<<X1_PbPb<<" "<<z2<<endl;
    // cout<<"X2_pp z3 X2_PbPb z4 = "<<X2_pp<<" "<<z3<<" "<<X2_PbPb<<" "<<z4<<endl;
  }

  cout<<"h_pp->GetCorrelationFactor() = "<<h_pp->GetCorrelationFactor()<<endl;
  cout<<"h_PbPb->GetCorrelationFactor() = "<<h_PbPb->GetCorrelationFactor()<<endl;
  cout<<"h_RAA->GetCorrelationFactor() = "<<h_RAA->GetCorrelationFactor()<<endl;
  vector<float> RAAcorr;
  RAAcorr.push_back(h_RAA->GetCorrelationFactor());

  TCanvas *c1 = new TCanvas("c1","c1",3000,1000);
  c1->Divide(3,1);
  c1->cd(1);
  h_pp->GetXaxis()->SetTitle("corrected yield (bin 1)");
  h_pp->GetYaxis()->SetTitle("corrected yield (bin 2)");
  h_pp->Draw("COLZ");
  c1->cd(2);
  h_PbPb->GetXaxis()->SetTitle("corrected yield (bin 1)");
  h_PbPb->GetYaxis()->SetTitle("corrected yield (bin 2)");
  h_PbPb->Draw("COLZ");
  c1->cd(3);
  h_RAA->GetXaxis()->SetTitle("R_{PbPb} (bin 1)");
  h_RAA->GetYaxis()->SetTitle("R_{PbPb} (bin 2)");
  h_RAA->Draw("COLZ");

  c1->SaveAs("FitErrorCorrelation.pdf");
  
  TFile * outf = new TFile("../../AccEffCorr/corrected_yields.root","UPDATE");
  outf->WriteObject(&RAAcorr,"r1r2Correlation_RAA");

}
