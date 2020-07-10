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

void Draw_corrYields(bool ispp=true){

  vector<vector<vector<float> > > *Yields_prefit;
  vector<vector<vector<float> > > *Yields_postfit;
  vector<vector<float> > *corrYield;
  vector<vector<float> > *corrYield_BDT23;
  vector<vector<float> > *corrYield_BDT3;
  vector<vector<float> > *corrYield_MC;
  vector<vector<float> > *corrYield_MCv2;
  vector<float> *rsig_relerr;
  vector<float> *eff_oneBinned;
  vector<float> *acc_oneBinned;
  vector<float> *r1r2Corr;

  TFile *infile3 = new TFile("../acceptance/acceptanceMap.root","READ");
  infile3->GetObject("acceptance_oneBinned", acc_oneBinned);

  TFile *infile2 = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
  infile2->GetObject("efficiency_oneBinned"+(TString)(ispp?"_pp":"_PbPb"), eff_oneBinned);

  bool BDTuncorrFromM = false;
  TFile *infile = new TFile("corrected_yields.root","READ");
  infile->GetObject("Yields_prefit"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), Yields_prefit);
  infile->GetObject("Yields_postfit"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), Yields_postfit);
  infile->GetObject("corrYield"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield);
  infile->GetObject("corrYield_BDTeff23"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_BDT23);
  infile->GetObject("corrYield_BDTeff3"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_BDT3);
  infile->GetObject("corrYield_MC"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_MC);
  infile->GetObject("corrYield_MCv2"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_MCv2);
  infile->GetObject("rsig_relerr"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), rsig_relerr);
  infile->GetObject("r1r2Correlation"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), r1r2Corr);

  const int nbins = _NanaBins;
  double xErr[nbins], yErr[nbins], yErrMC[nbins], zero[nbins], x[nbins];
  double y_nominal[nbins], y_no1stBDTbin[nbins], y_BDTeff23[nbins], y_no1st2ndBDTbin[nbins], y_BDTeff3[nbins], y_no1stBDTbin_BDTeffWAccEff[nbins], y_no1st2ndBDTbin_BDTeffWAccEff[nbins], y_oneBinned[nbins], y_MC[nbins], y_MCv2[nbins];
  for(int b=1;b<=nbins;b++){
    cout<<"analysis bin #"<<b<<endl;
    cout<<"one-binned acceptance, efficiency = "<<(*acc_oneBinned)[b]<<" "<<(*eff_oneBinned)[b]<<endl;

    float BDTeffWAccEff_bins23 = ((*corrYield_MCv2)[b][2] + (*corrYield_MCv2)[b][3]) / (*corrYield_MCv2)[b][0];
    float BDTeffWAccEff_bin3 = (*corrYield_MCv2)[b][3] / (*corrYield_MCv2)[b][0];
    float BDTeff_bins23 = ((*Yields_postfit)[0][b][2] + (*Yields_postfit)[0][b][3]) / (*Yields_postfit)[0][b][0];
    float BDTeff_bin3 = (*Yields_postfit)[0][b][3] / (*Yields_postfit)[0][b][0];
    cout<<"BDT efficiency for bins 2-3 and for bin 3 = "<<BDTeff_bins23<<" "<<BDTeff_bin3<<endl;
    cout<<"BDT efficiency for bins 2-3 and for bin 3 (WAccEff) = "<<BDTeffWAccEff_bins23<<" "<<BDTeffWAccEff_bin3<<endl;

    x[b-1] = (_BcPtmin[b]+_BcPtmax[b])/2;
    xErr[b-1] = (-_BcPtmin[b]+_BcPtmax[b])/2;
    zero[b-1] = 0;
    y_nominal[b-1] = (*corrYield)[b][0];
    yErr[b-1] = (*rsig_relerr)[b] * y_nominal[b-1];
    y_BDTeff23[b-1] = (*corrYield_BDT23)[b][2] + (*corrYield_BDT23)[b][3];
    y_BDTeff3[b-1] = (*corrYield_BDT3)[b][3];
    y_MC[b-1] = (*corrYield_MC)[b][0];
    yErrMC[b-1] = (*rsig_relerr)[b] * y_MC[b-1];
    y_MCv2[b-1] = (*corrYield_MCv2)[b][0];
    //add prefit yields for MC in infile, to correct by prefit BDT efficiency of being in 3rd bin
    y_no1stBDTbin_BDTeffWAccEff[b-1] = ((*corrYield)[b][2] + (*corrYield)[b][3]) / BDTeffWAccEff_bins23;
    y_no1st2ndBDTbin_BDTeffWAccEff[b-1] = (*corrYield)[b][3] / BDTeffWAccEff_bin3;
    y_no1stBDTbin[b-1] = ((*corrYield)[b][2] + (*corrYield)[b][3]) / BDTeff_bins23;
    y_no1st2ndBDTbin[b-1] = (*corrYield)[b][3] / BDTeff_bin3;
    y_oneBinned[b-1] = (*Yields_postfit)[0][b][0] / ((*acc_oneBinned)[b] * (*eff_oneBinned)[b]);
    cout<<"b, y_nominal, yErr = "<<b<<" "<<y_nominal[b-1]<<" "<<yErr[b-1]<<endl;
    cout<<"y_BDTeff23[b-1], corrYield_BDT23[b][2], corrYield_BDT23)[b][3] = "<<y_BDTeff23[b-1]<<" "<<(*corrYield_BDT23)[b][2] <<" "<< (*corrYield_BDT23)[b][3]<<endl; 
    cout<<"y_BDTeff3[b-1] = "<<y_BDTeff3[b-1]<<endl; 
  }

  TGraphErrors *g_nominal = new TGraphErrors(nbins, x,y_nominal, xErr,zero);
  TGraphErrors *g_oneBinned = new TGraphErrors(nbins, x,y_oneBinned, xErr,zero);
  TGraphErrors *g_BDTeff23 = new TGraphErrors(nbins, x,y_BDTeff23, xErr,zero);
  TGraphErrors *g_BDTeff3 = new TGraphErrors(nbins, x,y_BDTeff3, xErr,zero);
  TGraphErrors *g_no1stBDTbin = new TGraphErrors(nbins, x,y_no1stBDTbin, xErr,zero);
  TGraphErrors *g_no1st2ndBDTbin = new TGraphErrors(nbins, x,y_no1st2ndBDTbin, xErr,zero);
  TGraphErrors *g_no1stBDTbin_BDTeffWAccEff = new TGraphErrors(nbins, x,y_no1stBDTbin_BDTeffWAccEff, xErr,zero);
  TGraphErrors *g_no1st2ndBDTbin_BDTeffWAccEff = new TGraphErrors(nbins, x,y_no1st2ndBDTbin_BDTeffWAccEff, xErr,zero);
  TGraphErrors *g_MC = new TGraphErrors(nbins, x,y_MC, xErr,yErrMC);
  TGraphErrors *g_MCv2 = new TGraphErrors(nbins, x,y_MCv2, xErr,zero);

  TCanvas *c1 = new TCanvas("c1","c1",2000,2000);
  g_nominal->SetMarkerSize(3);
  g_nominal->SetMarkerColor(kBlack);
  g_nominal->SetMarkerStyle(20);
  g_oneBinned->SetMarkerSize(3);
  g_oneBinned->SetMarkerColor(kMagenta);
  g_oneBinned->SetMarkerStyle(34);
  g_no1stBDTbin->SetMarkerSize(3);
  g_no1stBDTbin->SetMarkerColor(kBlue);
  g_no1stBDTbin->SetMarkerStyle(21);
  g_BDTeff23->SetMarkerSize(3);
  g_BDTeff23->SetMarkerColor(kBlue+2);
  g_BDTeff23->SetMarkerStyle(25);
  // g_no1stBDTbin_BDTeffWAccEff->SetMarkerSize(3);
  // g_no1stBDTbin_BDTeffWAccEff->SetMarkerColor(kBlue+2);
  // g_no1stBDTbin_BDTeffWAccEff->SetMarkerStyle(25);
  g_no1st2ndBDTbin->SetMarkerSize(4);
  g_no1st2ndBDTbin->SetMarkerColor(kCyan);
  g_no1st2ndBDTbin->SetMarkerStyle(29); //stars need larger size
  g_BDTeff3->SetMarkerSize(4);
  g_BDTeff3->SetMarkerColor(kCyan-5);
  g_BDTeff3->SetMarkerStyle(30); //stars need larger size
  // g_no1st2ndBDTbin_BDTeffWAccEff->SetMarkerSize(4);
  // g_no1st2ndBDTbin_BDTeffWAccEff->SetMarkerColor(kCyan-5);
  // g_no1st2ndBDTbin_BDTeffWAccEff->SetMarkerStyle(30); //stars need larger size
  g_MC->SetMarkerSize(3);
  g_MC->SetMarkerColor(kGreen);
  g_MC->SetMarkerStyle(22);
  g_MCv2->SetMarkerSize(3);
  g_MCv2->SetMarkerColor(kGreen+2);
  g_MCv2->SetMarkerStyle(23);

  g_MC->GetYaxis()->SetRangeUser(0.5*y_MC[nbins-1],2*y_MC[0]);//+3*yErr[0]);
  g_MC->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);//+3*yErr[0]);
  g_MC->SetTitle("Systematic variations of corrected yields");
  g_MC->GetYaxis()->SetTitle("corrected yield");
  g_MC->GetXaxis()->SetTitle("p_{T}(#mu#mu#mu) [GeV]");

  TLegend *leg = new TLegend(0.5,0.6,0.9,0.9);
  if(ispp) leg->AddEntry(g_nominal,"nominal");
  leg->AddEntry(g_oneBinned,"one-binned AccEff");
  leg->AddEntry(g_MC,"MC v1: method crosscheck");
  leg->AddEntry(g_MCv2,"MC v2: effect of AccEff maps (close to one-binned)");
  leg->AddEntry(g_no1stBDTbin,"no 1st BDT bin");
  leg->AddEntry(g_BDTeff23,"no 1st BDT bin, BDT eff map");
  leg->AddEntry(g_no1st2ndBDTbin,"no 1st & 2nd BDT bin");
  leg->AddEntry(g_BDTeff3,"no 1st & 2nd BDT bin, BDT eff map");
  leg->SetTextSize(0.034);
  leg->SetBorderSize(0);

  g_MC->Draw("AP");
  if(ispp) g_nominal->Draw("Psame");
  g_oneBinned->Draw("Psame");
  g_no1stBDTbin->Draw("Psame");
  g_no1st2ndBDTbin->Draw("Psame");
  g_BDTeff23->Draw("Psame");
  g_BDTeff3->Draw("Psame");
  g_MCv2->Draw("Psame");
  leg->Draw("same");
  gPad->SetLogy();

  //DRAW CMS Preliminary                                                                                                                                                                                                                  
  TLatex CMStag;
  CMStag.SetNDC();
  CMStag.SetTextFont(42);
  CMStag.SetTextSize(0.035);
  CMStag.DrawLatex(0.19,0.85,"#font[61]{CMS "+(TString)(ispp?"pp":"PbPb")+"}");
  CMStag.DrawLatex(0.19,0.80,"#font[52]{Work in progress}");

  c1->SetLeftMargin(0.15);
  c1->SaveAs("correctedYields_AccEffSystVariations_"+(TString)(ispp?"pp":"PbPb")+".pdf");

  //******************************************************
  //Record the nominal result + systematic error
  //******************************************************
  vector<float> y_nom(nbins);
  vector<float> y_systErr(nbins);
  vector<float> y_fitErr(nbins);

  for(int b=0;b<nbins;b++){
    if(ispp) {
      y_nom[b] = (y_nominal[b] + y_oneBinned[b] + y_BDTeff23[b])/3;
      //systErr = max deviation from average of the 3 methods
      y_systErr[b] = fabs(y_nom[b] - y_nominal[b]);
      y_systErr[b] = max(y_systErr[b] , (float)fabs(y_nom[b] - y_oneBinned[b]) );
      y_systErr[b] = max(y_systErr[b] , (float)fabs(y_nom[b] - y_BDTeff23[b]) );
    }
    else{
      y_nom[b] = (y_oneBinned[b] + y_BDTeff23[b])/2;
      //systErr = max deviation from average of the 2 methods
      y_systErr[b] = fabs(y_nom[b] - y_BDTeff3[b]);
      y_systErr[b] = max(y_systErr[b] , (float)fabs(y_nom[b] - y_oneBinned[b]) );
    }

    y_fitErr[b] = (*rsig_relerr)[b+1] * y_nom[b];
  }

  TFile * outf = new TFile("corrected_yields.root","UPDATE");
  outf->WriteObject(&y_nom,"FinalCorrectedYield"+(TString)(ispp?"_pp":"_PbPb"));
  outf->WriteObject(&y_systErr,"FinalCorrectedYield_systError"+(TString)(ispp?"_pp":"_PbPb"));
  outf->WriteObject(&y_fitErr,"FinalCorrectedYield_fitError"+(TString)(ispp?"_pp":"_PbPb"));
  outf->WriteObject(r1r2Corr,"r1r2Correlation"+(TString)(ispp?"_pp":"_PbPb"));

}
