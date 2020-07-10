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
  vector<float> *y_systErr_pp;
  vector<float> *y_fitErr_pp;
  vector<float> *r1r2Corr_pp;
  vector<float> *y_nom_PbPb;
  vector<float> *y_systErr_PbPb;
  vector<float> *y_fitErr_PbPb;
  vector<float> *r1r2Corr_PbPb;

  TFile *infile = new TFile("../AccEffCorr/corrected_yields.root","READ");
  infile->GetObject("FinalCorrectedYield_pp", y_nom_pp);
  infile->GetObject("FinalCorrectedYield_fitError_pp", y_fitErr_pp);
  infile->GetObject("FinalCorrectedYield_systError_pp", y_systErr_pp);
  infile->GetObject("r1r2Correlation_pp", r1r2Corr_pp);
  infile->GetObject("FinalCorrectedYield_PbPb", y_nom_PbPb);
  infile->GetObject("FinalCorrectedYield_fitError_PbPb", y_fitErr_PbPb);
  infile->GetObject("FinalCorrectedYield_systError_PbPb", y_systErr_PbPb);
  infile->GetObject("r1r2Correlation_PbPb", r1r2Corr_PbPb);

  const int nbins = _NanaBins;
  const int nYregions = 2;//hard-coded
  int nPts[nYregions+1] = {2,1,1};
  const int nMaxPts = 2;//hard-coded
  int nbins_test = 0;
  for(int b=1;b<=nYregions;b++) nbins_test += nPts[b];
  if(nbins_test!=nbins) cout<<"!!!!!!!!!!!! PROBLEM with number of analysis bins and/or number of Y regions ! Expected segfault"<<endl;

  double xErr[nYregions+1][nMaxPts], x[nYregions+1][nMaxPts], zero[nYregions+1][nMaxPts];
  double y_pp[nYregions+1][nMaxPts], yErr_pp[nYregions+1][nMaxPts], y_PbPb[nYregions+1][nMaxPts], yErr_PbPb[nYregions+1][nMaxPts], y_RAA[nYregions+1][nMaxPts], yErr_RAA[nYregions+1][nMaxPts];

  for(int b=0;b<=nYregions;b++){ //b==0 is the graph with all points, without Y discrimination
    for(int m=0;m<nPts[b];m++){
      cout<<"analysis bin #"<<b<<endl;

      int trueb = (b==0)?(m+1):b;
      x[b][m] = (_BcPtmin[trueb]+_BcPtmax[trueb])/2;
      xErr[b][m] = (-_BcPtmin[trueb]+_BcPtmax[trueb])/2;
      zero[b][m] = 0;

      float normpp = L_pp * (2*xErr[b][m]);
      y_pp[b][m] = (*y_nom_pp)[b] / normpp;
      yErr_pp[b][m] = sqrt(pow((*y_fitErr_pp)[b] ,2) + pow((*y_systErr_pp)[b] ,2) ) / normpp;

      float normPbPb = NMB_PbPb * TAA_090 * (2*xErr[b][m]); //Leq_PbPb if we choose the lumi option
      y_PbPb[b][m] = (*y_nom_PbPb)[b] / normPbPb;
      yErr_PbPb[b][m] = sqrt(pow((*y_fitErr_PbPb)[b] ,2) + pow((*y_systErr_PbPb)[b] ,2) ) / normPbPb;

      y_RAA[b][m] = y_PbPb[b][m]/y_pp[b][m];
      yErr_RAA[b][m] = sqrt(pow(yErr_pp[b][m],2) + pow(yErr_PbPb[b][m],2));

    }
  }

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
  g_pp[0]->SetTitle(":p_{T}^{vis}(B_{c}) [GeV]:#frac{d#sigma_{pp}}{dp_{T}^{vis}} or #frac{1}{N_{MB} T_{PbPb}}#frac{dN_{PbPb}^{corr}}{dp_{T}^{vis}}");
  g_pp[0]->Draw("AP");
  g_pp[1]->Draw("Psame");
  g_pp[2]->Draw("Psame");
  g_PbPb[1]->Draw("Psame");
  g_PbPb[2]->Draw("Psame");

  vector<TLegend*> leg;
  for(int b=0;b<nYregions;b++){
    leg.push_back(new TLegend((b==0)?0.7:0.8,0.7,(b==0)?0.8:0.9,0.85));
    leg[b]->AddEntry(g_pp[b+1],(b==0)?"pp":" ");
    leg[b]->AddEntry(g_PbPb[b+1],(b==0)?"PbPb":" ");
    leg[b]->SetTextSize(0.039);
    leg[b]->SetBorderSize(0);
    leg[b]->Draw("same");
  }
  TLatex legEntry;
  legEntry.SetNDC();
  legEntry.SetTextFont(42);
  legEntry.SetTextSize(0.035);
  legEntry.DrawLatex(0.7,0.85,"1.3<|Y^{vis}|<2.3");
  legEntry.DrawLatex(0.8,0.85,"0<|Y^{vis}|<2.3");

  gPad->SetLogy();

  //DRAW CMS Preliminary                                                                                                                                                                                                 
  TLatex CMStag;
  CMStag.SetNDC();
  CMStag.SetTextFont(42);
  CMStag.SetTextSize(0.035);
  CMStag.DrawLatex(0.1,0.9,"#font[61]{CMS } #font[52]{Work in progress}");
  CMStag.DrawLatex(0.7,0.9, Form("pp (%.0f pb^{-1}), PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));

  c1->SaveAs("CrossSections.pdf");

  TCanvas *c2 = new TCanvas("c2","RAA",2000,2000);
  g_RAA[0]->SetTitle(":p_{T}^{vis}(B_{c}) [GeV]:R_{PbPb}(B_{c})");
  g_RAA[0]->Draw("AP");
  g_RAA[1]->Draw("Psame");
  g_RAA[2]->Draw("Psame");

  TLegend *leg2 = new TLegend(0.6,0.6,0.9,0.9);
  leg2->AddEntry(g_RAA[1], "1.3<|Y^{vis}|<2.3");
  leg2->AddEntry(g_RAA[2], "0<|Y^{vis}|<2.3");
  leg2->SetTextSize(0.039);
  leg2->SetBorderSize(0);
  leg2->Draw("same");

  CMStag.DrawLatex(0.1,0.9,"#font[61]{CMS } #font[52]{Work in progress}");
  CMStag.DrawLatex(0.7,0.9, Form("pp (%.0f pb^{-1}), PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));

  TLine *one = new TLine(_BcPtmin[0],1,_BcPtmax[0],1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->Draw("same");

  c2->SaveAs("RAA.pdf");

}
