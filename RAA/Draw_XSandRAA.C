#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

void Draw_XSandRAA(){

  vector<float> *y_nom_pp;
  vector<float> *y_acceffErr_pp;
  vector<float> *y_fitErr_pp;
  vector<float> *y_metafitRelErr_pp;
  vector<float> *r1r2Corr_pp;
  vector<float> *metafitErrCorr_pp;
  vector<float> *y_nom_PbPb;
  vector<float> *y_acceffErr_PbPb;
  vector<float> *y_fitErr_PbPb;
  vector<float> *y_metafitRelErr_PbPb;
  vector<float> *r1r2Corr_PbPb;
  vector<float> *metafitErrCorr_PbPb;
  vector<float> *y_metafitRelErr_RAA;
  vector<float> *metafitErrCorr_RAA;
  vector<float> *r1r2Corr_RAA;

  TFile *infile = new TFile("../AccEffCorr/corrected_yields.root","READ");
  infile->GetObject("FinalCorrectedYield_pp", y_nom_pp);
  infile->GetObject("FinalCorrectedYield_fitError_pp", y_fitErr_pp);
  infile->GetObject("FinalCorrectedYield_AccEffSystError_pp", y_acceffErr_pp);
  infile->GetObject("CorrectedYield_MetafitRelSystError_pp", y_metafitRelErr_pp); //or "CorrectedYield_MetafitRelSystErrorMaxDeviation_pp"
  infile->GetObject("CorrectedYield_MetafitSyst_LinearizedCorrelationMatrix_pp", metafitErrCorr_pp);
  infile->GetObject("r1r2Correlation_pp", r1r2Corr_pp);
  infile->GetObject("FinalCorrectedYield_PbPb", y_nom_PbPb);
  infile->GetObject("FinalCorrectedYield_fitError_PbPb", y_fitErr_PbPb);
  infile->GetObject("FinalCorrectedYield_AccEffSystError_PbPb", y_acceffErr_PbPb);
  infile->GetObject("CorrectedYield_MetafitRelSystError_PbPb", y_metafitRelErr_PbPb);
  infile->GetObject("CorrectedYield_MetafitSyst_LinearizedCorrelationMatrix_PbPb", metafitErrCorr_PbPb);
  infile->GetObject("r1r2Correlation_PbPb", r1r2Corr_PbPb);
  infile->GetObject("RAA_MetafitRelSystError", y_metafitRelErr_RAA);
  infile->GetObject("RAA_MetafitSyst_LinearizedCorrelationMatrix", metafitErrCorr_RAA);
  infile->GetObject("r1r2Correlation_RAA", r1r2Corr_RAA);

  const int nbins = _NanaBins;
  const int nYregions = 2;//hard-coded
  int nPts[nYregions+1] = {2,1,1};
  const int nMaxPts = 2;//hard-coded
  int nbins_test = 0;
  for(int b=1;b<=nYregions;b++) nbins_test += nPts[b];
  if(nbins_test!=nbins) cout<<"!!!!!!!!!!!! PROBLEM with number of analysis bins and/or number of Y regions ! Expected segfault"<<endl;

  double xErr[nYregions+1][nMaxPts], x[nYregions+1][nMaxPts], zero[nYregions+1][nMaxPts];
  double y_metafitErr_pp[nbins],y_metafitErr_PbPb[nbins],y_metafitErr_RAA[nbins];
  double yErrPart_pp[nbins],yErrPart_PbPb[nbins],yfitErr_pp[nbins],yfitErr_PbPb[nbins],yErrRelPart_RAA[nbins],yfitErr_RAA[nbins];
  double y_pp[nYregions+1][nMaxPts], yErr_pp[nYregions+1][nMaxPts], y_PbPb[nYregions+1][nMaxPts], yErr_PbPb[nYregions+1][nMaxPts], y_RAA[nYregions+1][nMaxPts], yErr_RAA[nYregions+1][nMaxPts];
  double corrtot[3];

  for(int b=0;b<=nYregions;b++){ //b==0 is the graph with all points, without Y discrimination
    for(int m=0;m<nPts[b];m++){

      int trueb = (b==0)?(m+1):b;
      x[b][m] = (_BcPtmin[trueb]+_BcPtmax[trueb])/2;
      xErr[b][m] = (-_BcPtmin[trueb]+_BcPtmax[trueb])/2;
      zero[b][m] = 0;

      float normpp = L_pp * (2*xErr[b][m]);
      y_pp[b][m] = (*y_nom_pp)[trueb-1] / normpp;
      y_metafitErr_pp[trueb-1] = y_pp[b][m] * (*y_metafitRelErr_pp)[trueb-1];
      yfitErr_pp[trueb-1] = (*y_fitErr_pp)[trueb-1] / normpp;
      yErrPart_pp[trueb-1] = sqrt(pow(yfitErr_pp[trueb-1] ,2) + pow((*y_acceffErr_pp)[trueb-1] / normpp ,2) ); //without metafit error
      yErr_pp[b][m] = sqrt(pow(y_metafitErr_pp[trueb-1] ,2) + pow(yErrPart_pp[trueb-1] ,2));

      float normPbPb = NMB_PbPb * TAA_090 * (2*xErr[b][m]); //Leq_PbPb if we choose the lumi option
      y_PbPb[b][m] = (*y_nom_PbPb)[trueb-1] / normPbPb;
      y_metafitErr_PbPb[trueb-1] = y_PbPb[b][m] * (*y_metafitRelErr_PbPb)[trueb-1];
      yfitErr_PbPb[trueb-1] = (*y_fitErr_PbPb)[trueb-1] / normPbPb;
      yErrPart_PbPb[trueb-1] = sqrt(pow(yfitErr_PbPb[trueb-1] ,2) + pow((*y_acceffErr_PbPb)[trueb-1] / normPbPb ,2) ); //without metafit error
      yErr_PbPb[b][m] = sqrt(pow(y_metafitErr_PbPb[trueb-1] ,2) + pow(yErrPart_PbPb[trueb-1] ,2));

      y_RAA[b][m] = y_PbPb[b][m]/y_pp[b][m];
      y_metafitErr_RAA[trueb-1] = y_RAA[b][m] * (*y_metafitRelErr_RAA)[trueb-1];
      yErrRelPart_RAA[trueb-1] = sqrt(pow(yErrPart_pp[trueb-1]/y_pp[b][m],2) + pow(yErrPart_PbPb[trueb-1]/y_PbPb[b][m],2)); //without metafit error
      yfitErr_RAA[trueb-1] = y_RAA[b][m] * sqrt(pow(yfitErr_pp[trueb-1]/y_pp[b][m],2) + pow(yfitErr_PbPb[trueb-1]/y_PbPb[b][m],2) );
      yErr_RAA[b][m] = y_RAA[b][m] * sqrt(pow((*y_metafitRelErr_RAA)[trueb-1],2) + pow(yErrRelPart_RAA[trueb-1],2));

      if(b==0){
	cout<<"point "<<m<<": y_pp, y_PbPb, y_RAA = "<<y_pp[b][m]<<" "<<y_PbPb[b][m]<<" "<<y_RAA[b][m]<<endl;
	cout<<"point "<<m<<": errors on y_pp, y_PbPb, y_RAA = "<<yErr_pp[b][m]<<" "<<yErr_PbPb[b][m]<<" "<<yErr_RAA[b][m]<<endl;
	cout<<"fit errors on pp, PbPb, RAA = "<<yfitErr_pp[trueb-1]<<" "<<yfitErr_PbPb[trueb-1]<<" "<<yfitErr_RAA[trueb-1]<<endl;
	cout<<"point "<<m<<": RAA errors from fit, metafit, and acceff = "<<yfitErr_RAA[trueb-1]<<" "<<y_metafitErr_RAA[trueb-1]<<" "<<y_RAA[b][m] * sqrt(pow(yErrRelPart_RAA[trueb-1],2) - pow(yfitErr_RAA[trueb-1]/y_RAA[b][m],2))<<endl;
      }

    }
  }

  corrtot[0] = ( (*r1r2Corr_pp)[0] * yfitErr_pp[0] * yfitErr_pp[1]  +  (*metafitErrCorr_pp)[0] * y_metafitErr_pp[0] * y_metafitErr_pp[1] )/( yErr_pp[0][0]*yErr_pp[0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  corrtot[1] = ( (*r1r2Corr_PbPb)[0] * yfitErr_PbPb[0] * yfitErr_PbPb[1]  +  (*metafitErrCorr_PbPb)[0] * y_metafitErr_PbPb[0] * y_metafitErr_PbPb[1] )/( yErr_PbPb[0][0]*yErr_PbPb[0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  corrtot[2] = ( (*r1r2Corr_RAA)[0] * yfitErr_RAA[0] * yfitErr_RAA[1]  +  (*metafitErrCorr_RAA)[0] * y_metafitErr_RAA[0] * y_metafitErr_RAA[1] )/( yErr_RAA[0][0]*yErr_RAA[0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  cout<<"correlation factor between 2 pT bins, for pp, PbPb, RAA = "<<corrtot[0]<<" "<<corrtot[1]<<" "<<corrtot[2]<<endl;

  vector<TGraphErrors*> g_pp, g_PbPb, g_RAA;
  for(int b=0;b<=nYregions;b++){ 
    g_pp.push_back( new TGraphErrors(nPts[b], x[b],y_pp[b], xErr[b],yErr_pp[b]) );
    g_PbPb.push_back( new TGraphErrors(nPts[b], x[b],y_PbPb[b], xErr[b],yErr_PbPb[b]) );
    g_RAA.push_back( new TGraphErrors(nPts[b], x[b],y_RAA[b], xErr[b],yErr_RAA[b]) );

    g_pp[b]->SetMarkerSize(5);
    g_pp[b]->SetMarkerStyle((b==1)?24:20);
    g_pp[b]->SetMarkerColor(kBlue+2);
    g_pp[b]->SetLineColor(kBlue);
    g_pp[b]->SetLineWidth(4);

    g_PbPb[b]->SetMarkerSize(5);
    g_PbPb[b]->SetMarkerStyle((b==1)?25:21);
    g_PbPb[b]->SetMarkerColor(kSpring-6);
    g_PbPb[b]->SetLineColor(kSpring-5);
    g_PbPb[b]->SetLineWidth(4);

    g_RAA[b]->SetMarkerSize(5);
    g_RAA[b]->SetMarkerStyle((b==1)?24:20);
    g_RAA[b]->SetMarkerColor(kAzure+4);
    g_RAA[b]->SetLineColor(kAzure+5);
    g_RAA[b]->SetLineWidth(3);
  }

  TCanvas *c1 = new TCanvas("c1","XS",2000,2000);
  c1->SetLeftMargin(0.14);
  c1->SetTopMargin(0.14);
  g_pp[0]->SetTitle("Cross-sections");
  g_pp[0]->GetYaxis()->SetRangeUser(0.5*y_PbPb[0][1],2*y_PbPb[0][0]);
  g_pp[0]->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);
  g_pp[0]->GetYaxis()->SetTitleOffset(1.4);
  g_pp[0]->GetXaxis()->SetTitle("p_{T}^{vis}(B_{c}) [GeV]");
  g_pp[0]->GetYaxis()->SetTitle("BF  #times  #left(#frac{d#sigma_{pp}}{dp_{T}^{vis}}    or    #frac{1}{N_{MB} T_{PbPb}} #frac{dN_{PbPb}^{corr}}{dp_{T}^{vis}}#right)  [pb/GeV]");
  g_pp[0]->SetMarkerColor(kWhite);
  g_pp[0]->SetLineColor(kWhite);
  g_pp[0]->Draw("AP");
  g_pp[1]->Draw("Psame");
  g_pp[2]->Draw("Psame");
  g_PbPb[1]->Draw("Psame");
  g_PbPb[2]->Draw("Psame");

  vector<TLegend*> leg;
  for(int b=0;b<nYregions;b++){
    leg.push_back(new TLegend((b==0)?0.48:0.73,0.65,(b==0)?0.67:0.9,0.8));
    leg[b]->AddEntry(g_pp[b+1],(b==0)?"pp":"pp");
    leg[b]->AddEntry(g_PbPb[b+1],(b==0)?"PbPb":"PbPb");
    leg[b]->SetTextSize(0.039);
    leg[b]->SetBorderSize(0);
    leg[b]->Draw("same");
  }
  TLatex legEntry;
  legEntry.SetNDC();
  legEntry.SetTextFont(42);
  legEntry.SetTextSize(0.035);
  legEntry.DrawLatex(0.45,0.8,"1.3<|Y^{vis}|<2.3");
  legEntry.DrawLatex(0.71,0.8,"0<|Y^{vis}|<2.3");
  legEntry.DrawLatex(0.65,0.52,Form("#rho_{1-2}^{pp} = %.2f",corrtot[0]));
  legEntry.DrawLatex(0.65,0.45,Form("#rho_{1-2}^{PbPb} = %.2f",corrtot[1]));

  gPad->SetLogy();

  //DRAW CMS Preliminary                                                                                                                                                                                                 
  TLatex CMStag;
  CMStag.SetNDC();
  CMStag.SetTextFont(42);
  CMStag.SetTextSize(0.035);
  CMStag.DrawLatex(0.1,0.87,"#font[61]{CMS} #font[52]{Work in progress}");
  CMStag.DrawLatex(0.46,0.87, Form("pp (%.0f pb^{-1}), PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));

  c1->SaveAs("CrossSections.pdf");

  TCanvas *c2 = new TCanvas("c2","RAA",2000,2000);
  //c2->SetLeftMargin(0.14);
  c2->SetTopMargin(0.11);
  g_RAA[0]->SetTitle("R_{PbPb}");
  g_RAA[0]->SetMarkerColor(kWhite);
  g_RAA[0]->SetLineColor(kWhite);
  g_RAA[0]->GetYaxis()->SetTitleOffset(1.2);
  g_RAA[0]->GetYaxis()->SetRangeUser(0,2.5);
  g_RAA[0]->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);
  g_RAA[0]->GetXaxis()->SetTitle("p_{T}^{vis}(B_{c}) [GeV]");
  g_RAA[0]->GetYaxis()->SetTitle("R_{PbPb}(B_{c})");
  g_RAA[0]->Draw("AP");
  g_RAA[1]->Draw("Psame");
  g_RAA[2]->Draw("Psame");

  TLegend *leg2 = new TLegend(0.5,0.65,0.8,0.82);
  leg2->AddEntry(g_RAA[1], "1.3<|Y^{vis}|<2.3");
  leg2->AddEntry(g_RAA[2], "0<|Y^{vis}|<2.3");
  leg2->SetTextSize(0.039);
  leg2->SetBorderSize(0);
  leg2->Draw("same");

  legEntry.DrawLatex(0.65,0.57,Form("#rho_{1-2} = %.2f",corrtot[2]));

  CMStag.DrawLatex(0.1,0.9,"#font[61]{CMS} #font[52]{Work in progress}");
  CMStag.DrawLatex(0.46,0.9, Form("pp (%.0f pb^{-1}), PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));

  TLine *one = new TLine(_BcPtmin[0],1,_BcPtmax[0],1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->Draw("same");

  c2->SaveAs("RAA.pdf");

}
