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

void corrYields(bool ispp=true, bool secondStep=true, bool thirdStep=false, bool centDep=false){

  vector<vector<vector<float> > > *Yields_postfit;
  vector<vector<float> > *corrYield;
  vector<vector<float> > *corrYield_BDT23;
  vector<vector<float> > *corrYield_BDT3;
  vector<vector<float> > *corrYieldErr;
  vector<vector<float> > *corrYieldErr_BDT23;
  vector<vector<float> > *corrYieldErr_BDT3;
  vector<vector<float> > *corrYield_MC;
  vector<vector<float> > *corrYield_MCv2;
  vector<vector<float> > *rsig_relerr;
  vector<vector<float> > *eff_oneBinned;
  vector<vector<float> > *acc_oneBinned;
  vector<vector<float> > *eff_oneBinned_2ndStep;
  vector<vector<float> > *acc_oneBinned_2ndStep;
  vector<vector<double> > *InvAccEff;
  

  TFile *infile3 = new TFile("../acceptance/acceptanceMap.root","READ");
  TFile *infile2 = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
  infile3->GetObject("acceptance_oneBinned", acc_oneBinned);
  infile2->GetObject("efficiency_oneBinned"+(TString)((centDep && !ispp)?"_centDep":"")+(TString)(ispp?"_pp":"_PbPb"), eff_oneBinned);
  if(secondStep){
    infile3->GetObject("acceptance_oneBinned_2ndStep", acc_oneBinned_2ndStep);
    infile2->GetObject("efficiency_oneBinned"+(TString)((centDep && !ispp)?"_centDep":"")+(TString)(ispp?"_pp":"_PbPb")+(TString)"_2ndStep", eff_oneBinned_2ndStep);
  }

  TFile *infile4 = (TFile*)(thirdStep?(new TFile("../twoSteps/AccEffFrom2ndStepToys.root","READ")):NULL);
  if(thirdStep) infile4->GetObject("InvAccEffFromCorrMC_withSystErr_"+(TString)(ispp?"pp":"PbPb")+(TString)(centDep?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""), InvAccEff);

  bool BDTuncorrFromM = false;
  TFile *infile = new TFile("corrected_yields"+(TString)(secondStep?"_2ndStep":"")+".root","READ");
  infile->GetObject("Yields_postfit"+(TString)((centDep && !ispp)?"_centralityDep":"")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), Yields_postfit);
  infile->GetObject("corrYield"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield);
  infile->GetObject("corrYield_BDTeff23"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_BDT23);
  infile->GetObject("corrYield_BDTeff3"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_BDT3);
  infile->GetObject("corrYieldErr"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYieldErr);
  infile->GetObject("corrYieldErr_BDTeff23"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYieldErr_BDT23);
  infile->GetObject("corrYieldErr_BDTeff3"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYieldErr_BDT3);
  infile->GetObject("corrYield_MC"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_MC);
  infile->GetObject("corrYield_MCv2"+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), corrYield_MCv2);
  infile->GetObject("rsig_relerr"+(TString)((centDep && !ispp)?"_centralityDep":"")+(TString)(BDTuncorrFromM?"_BDTuncorrFromM":"")+(TString)(ispp?"_pp":"_PbPb"), rsig_relerr);

  const int nbins = _NanaBins+1;
  double xErr[nbins], yErr[nbins], yErrMC[nbins], zero[nbins], x[nbins];
  double x_EByE[nbins], x_no1stBDTbin[nbins], x_no1st2ndBDTbin[nbins], x_2steps[nbins], x_3steps[nbins];
  double yErr_EByE[nbins], yErr_no1stBDTbin[nbins], yErr_no1st2ndBDTbin[nbins];
  double y_2steps[nbins], yErr_2steps[nbins];
  double y_3steps[nbins], yErr_3steps[nbins];
  double y_EByE[nbins], y_no1stBDTbin[nbins], y_BDTeff23[nbins], y_no1st2ndBDTbin[nbins], y_BDTeff3[nbins], y_no1stBDTbin_BDTeffWAccEff[nbins], y_no1st2ndBDTbin_BDTeffWAccEff[nbins], y_oneBinned[nbins], y_MC[nbins], y_MCv2[nbins];
  for(int b=0;b<nbins;b++){
    cout<<"analysis bin #"<<b<<endl;
    cout<<"one-binned acceptance, efficiency = "<<(*acc_oneBinned)[1-(int)ispp][centDep?0:b]<<" "<<(*eff_oneBinned)[b][0]<<endl;

    float BDTeffWAccEff_bins23 = ((*corrYield_MCv2)[b][2] + (*corrYield_MCv2)[b][3]) / (*corrYield_MCv2)[b][0];
    float BDTeffWAccEff_bin3 = (*corrYield_MCv2)[b][3] / (*corrYield_MCv2)[b][0];
    float BDTeff_bins23 = ((*Yields_postfit)[0][b][2] + (*Yields_postfit)[0][b][3]) / (*Yields_postfit)[0][b][0];
    float BDTeff_bin3 = (*Yields_postfit)[0][b][3] / (*Yields_postfit)[0][b][0];
    cout<<"BDT efficiency for bins 2-3 and for bin 3 = "<<BDTeff_bins23<<" "<<BDTeff_bin3<<endl;
    cout<<"BDT efficiency for bins 2-3 and for bin 3 (WAccEff) = "<<BDTeffWAccEff_bins23<<" "<<BDTeffWAccEff_bin3<<endl;

    x[b] = (_BcPtmin[b]+_BcPtmax[b])/2;
    xErr[b] = (-_BcPtmin[b]+_BcPtmax[b])/2;
    x_EByE[b] = x[b] + 0.06*((b==1)?3:1)*xErr[b];
    x_no1stBDTbin[b] = x[b] + 0.12*((b==1)?3:1)*xErr[b];
    x_no1st2ndBDTbin[b] = x[b] + 0.16*((b==1)?3:1)*xErr[b];
    x_2steps[b] = x[b] + 0.20*((b==1)?3:1)*xErr[b];
    x_3steps[b] = x[b] + 0.24*((b==1)?3:1)*xErr[b];
    zero[b] = 0;
    y_EByE[b] = (*corrYield)[b][0];
    y_BDTeff23[b] = (*corrYield_BDT23)[b][2] + (*corrYield_BDT23)[b][3];
    y_BDTeff3[b] = (*corrYield_BDT3)[b][3];
    y_MC[b] = (*corrYield_MC)[b][0];
    yErrMC[b] = (*rsig_relerr)[b][0] * y_MC[b];
    y_MCv2[b] = (*corrYield_MCv2)[b][0];
    //add prefit yields for MC in infile, to correct by prefit BDT efficiency of being in 3rd bin
    y_no1stBDTbin_BDTeffWAccEff[b] = ((*corrYield)[b][2] + (*corrYield)[b][3]) / BDTeffWAccEff_bins23;
    y_no1st2ndBDTbin_BDTeffWAccEff[b] = (*corrYield)[b][3] / BDTeffWAccEff_bin3;
    y_no1stBDTbin[b] = ((*corrYield)[b][2] + (*corrYield)[b][3]) / BDTeff_bins23;
    y_no1st2ndBDTbin[b] = (*corrYield)[b][3] / BDTeff_bin3;
    y_oneBinned[b] = (*Yields_postfit)[0][b][0] / ((*acc_oneBinned)[1-(int)ispp][centDep?0:b] * (*eff_oneBinned)[b][0]);
    yErr[b] = (*rsig_relerr)[b][0] * y_oneBinned[b];
    if(secondStep){
      y_2steps[b] = (*Yields_postfit)[0][b][0]  / ((*acc_oneBinned_2ndStep)[1-(int)ispp][centDep?0:b] * (*eff_oneBinned_2ndStep)[b][0]);
      yErr_2steps[b] = (*rsig_relerr)[b][0] * y_2steps[b];
    }
    if(thirdStep){
      y_3steps[b] = (*Yields_postfit)[0][b][0] * (*InvAccEff)[b][0];
      yErr_3steps[b] = (*Yields_postfit)[0][b][0] * (*InvAccEff)[b][1];//sqrt(pow((*Yields_postfit)[0][b][0] * (*InvAccEff)[b][1], 2) + pow( (*rsig_relerr)[b][0] * y_2steps[b] ,2));
    }
    yErr_EByE[b] = (*corrYieldErr)[b][0];
    yErr_no1stBDTbin[b] = (*corrYieldErr_BDT23)[b][0];
    yErr_no1st2ndBDTbin[b] = (*corrYieldErr_BDT3)[b][0];
    cout<<"b, y_EByE, yErr, y_MC = "<<b<<" "<<y_EByE[b]<<" "<<yErr[b]<<" "<<y_MC[b]<<endl;
    cout<<"Yields_postfit, y_oneBinned, InvAccEff (if 3rd step), y_2steps, yErr_2steps = "<<(*Yields_postfit)[0][b][0]<<" "<<y_oneBinned[b]<<" "<<(thirdStep?((*InvAccEff)[b][0]):0)<<" "<<y_2steps[b]<<" "<<yErr_2steps[b]<<endl;
    // cout<<"y_BDTeff23[b], corrYield_BDT23[b][2], corrYield_BDT23)[b][3] = "<<y_BDTeff23[b]<<" "<<(*corrYield_BDT23)[b][2] <<" "<< (*corrYield_BDT23)[b][3]<<endl; 
    // cout<<"y_BDTeff3[b] = "<<y_BDTeff3[b]<<endl; 
  }

  if(!centDep){
    TGraphErrors *g_EByE = new TGraphErrors(nbins, x_EByE+1,y_EByE+1,zero+1,yErr_EByE+1);
    TGraphErrors *g_2steps = new TGraphErrors(nbins, x_2steps+1,y_2steps+1,zero+1,yErr_2steps+1);
    TGraphErrors *g_3steps = new TGraphErrors(nbins, x_3steps+1,y_3steps+1,zero+1,yErr_3steps+1);
    TGraphErrors *g_oneBinned = new TGraphErrors(nbins, x+1,y_oneBinned+1, xErr+1,yErr+1);
    TGraphErrors *g_BDTeff23 = new TGraphErrors(nbins, x_no1stBDTbin+1,y_BDTeff23+1, zero+1,zero+1);
    TGraphErrors *g_BDTeff3 = new TGraphErrors(nbins, x_no1st2ndBDTbin+1,y_BDTeff3+1, zero+1,zero+1);
    TGraphErrors *g_no1stBDTbin = new TGraphErrors(nbins, x_no1stBDTbin+1,y_no1stBDTbin+1, zero+1,yErr_no1stBDTbin+1);
    TGraphErrors *g_no1st2ndBDTbin = new TGraphErrors(nbins, x_no1st2ndBDTbin+1,y_no1st2ndBDTbin+1, zero+1,yErr_no1st2ndBDTbin+1);
    TGraphErrors *g_no1stBDTbin_BDTeffWAccEff = new TGraphErrors(nbins, x+1,y_no1stBDTbin_BDTeffWAccEff+1, xErr+1,zero+1);
    TGraphErrors *g_no1st2ndBDTbin_BDTeffWAccEff = new TGraphErrors(nbins, x+1,y_no1st2ndBDTbin_BDTeffWAccEff+1, xErr+1,zero+1);
    TGraphErrors *g_MC = new TGraphErrors(nbins, x+1,y_MC+1, xErr+1,zero+1);
    TGraphErrors *g_MCv2 = new TGraphErrors(nbins, x+1,y_MCv2+1, xErr+1,zero+1);

    TCanvas *c1 = new TCanvas("c1","c1",2000,2000);
    g_2steps->SetMarkerSize(4);
    g_2steps->SetMarkerColor(kOrange+9);
    g_2steps->SetMarkerStyle(33);
    g_3steps->SetMarkerSize(4);
    g_3steps->SetMarkerColor(kViolet);
    g_3steps->SetMarkerStyle(21);
    g_EByE->SetMarkerSize(3);
    g_EByE->SetMarkerColor(kBlack);
    g_EByE->SetMarkerStyle(20);
    g_oneBinned->SetMarkerSize(3);
    g_oneBinned->SetMarkerColor(kMagenta);
    g_oneBinned->SetMarkerStyle(34);
    // g_no1stBDTbin->SetMarkerSize(3);
    // g_no1stBDTbin->SetMarkerColor(kBlue);
    // g_no1stBDTbin->SetMarkerStyle(22);
    g_BDTeff23->SetMarkerSize(3);
    g_BDTeff23->SetMarkerColor(kBlue+2);
    g_BDTeff23->SetMarkerStyle(22);
    // g_no1stBDTbin_BDTeffWAccEff->SetMarkerSize(3);
    // g_no1stBDTbin_BDTeffWAccEff->SetMarkerColor(kBlue+2);
    // g_no1stBDTbin_BDTeffWAccEff->SetMarkerStyle(25);
    // g_no1st2ndBDTbin->SetMarkerSize(4);
    // g_no1st2ndBDTbin->SetMarkerColor(kCyan);
    // g_no1st2ndBDTbin->SetMarkerStyle(23);
    g_BDTeff3->SetMarkerSize(4);
    g_BDTeff3->SetMarkerColor(kCyan-5);
    g_BDTeff3->SetMarkerStyle(23); //stars need larger size
    // g_no1st2ndBDTbin_BDTeffWAccEff->SetMarkerSize(4);
    // g_no1st2ndBDTbin_BDTeffWAccEff->SetMarkerColor(kCyan-5);
    // g_no1st2ndBDTbin_BDTeffWAccEff->SetMarkerStyle(29); //stars need larger size
    g_MC->SetMarkerSize(3);
    g_MC->SetMarkerColor(kGreen);
    g_MC->SetMarkerStyle(22);
    g_MCv2->SetMarkerSize(3);
    g_MCv2->SetMarkerColor(kGreen+2);
    g_MCv2->SetMarkerStyle(23);

    g_MC->GetYaxis()->SetRangeUser((ispp?0.7:0.5)*y_oneBinned[nbins-1],(ispp?1.2:2.)*y_oneBinned[1]);
    g_MC->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);
    g_MC->SetTitle("Variations of acceptance and efficiency method");
    g_MC->GetYaxis()->SetTitle("corrected yield");
    g_MC->GetXaxis()->SetTitle("p_{T}(#mu#mu#mu) [GeV]");

    TLegend *leg = new TLegend(0.45,0.6,0.9,0.9);
    if(thirdStep) leg->AddEntry(g_3steps,"2-steps (nominal)");
    if(secondStep) leg->AddEntry(g_2steps,"1^{st}-step corrected MC");
    leg->AddEntry(g_oneBinned,"one-binned AccEff","lpe");
    if(ispp) leg->AddEntry(g_EByE,"full event-by-event");
    //leg->AddEntry(g_MCv2,"MC event-by-event v2 crosscheck");
    //  leg->AddEntry(g_no1stBDTbin,"no BDT bin 1");
    leg->AddEntry(g_BDTeff23,"no BDT bin 1");//, BDT eff map
    //leg->AddEntry(g_no1st2ndBDTbin,"no BDT bin 1-2");
    leg->AddEntry(g_BDTeff3,"no BDT bin 1-2");//, BDT eff map
    leg->AddEntry(g_MC,"MC event-by-event crosscheck");
    leg->SetTextSize(0.034);
    leg->SetBorderSize(0);

    g_MC->Draw("AP");
    if(ispp) g_EByE->Draw("Psame");
    g_oneBinned->Draw("Psame");
    if(secondStep) g_2steps->Draw("Psame");
    if(thirdStep) g_3steps->Draw("Psame");
    // g_no1stBDTbin->Draw("Psame");
    // g_no1st2ndBDTbin->Draw("Psame");
    g_BDTeff23->Draw("Psame");
    g_BDTeff3->Draw("Psame");
    //  g_MCv2->Draw("Psame");
    leg->Draw("same");
    if(!ispp) gPad->SetLogy();

    //DRAW CMS Preliminary                                                                                                                                                                                                                  
    TLatex CMStag;
    CMStag.SetNDC();
    CMStag.SetTextFont(42);
    CMStag.SetTextSize(0.035);
    CMStag.DrawLatex(0.19,0.85,"#font[61]{CMS "+(TString)(ispp?"pp":"PbPb")+"}");
    CMStag.DrawLatex(0.19,0.80,"#font[52]{Work in progress}");

    c1->SetLeftMargin(0.15);
    c1->SaveAs("figs/correctedYields_AccEffSystVariations_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+".pdf");
  }

  //******************************************************
  //Record the nominal result
  //******************************************************
  vector<float> y_nom(nbins);
  vector<vector<float> > y_fitErr(nbins, vector<float>(3));
  vector<float> y_TnPrelErr(nbins);

  for(int b=0;b<nbins;b++){ //nbins=2+1 (with integrated at position 0)
    y_nom[b] =  thirdStep?y_3steps[b]:(secondStep?y_2steps[b]:y_oneBinned[b]);

    cout<<"b, y  = "<<b<<" "<<y_nom[b]<<endl;

    y_fitErr[b][1] = (*rsig_relerr)[b][1] * y_nom[b];
    y_fitErr[b][2] = (*rsig_relerr)[b][2] * y_nom[b];
    y_fitErr[b][0] = (y_fitErr[b][1] + y_fitErr[b][2])/2;

    y_TnPrelErr[b] = (*eff_oneBinned)[b][3] / (*eff_oneBinned)[b][0];
    //cout<<"fit err, fit relative err, tnp relative err = "<< y_fitErr[b][0]<<" "<<(*rsig_relerr)[b][1]<<" "<<y_TnPrelErr[b] <<endl;
  }

  TFile * outf = new TFile("corrected_yields"+(TString)(thirdStep?"_3rdStep":(secondStep?"_2ndStep":""))+".root", "UPDATE");
  outf->WriteObject(&y_nom,"FinalCorrectedYield"+(TString)(centDep?"_centralityDep":"")+(TString)(ispp?"_pp":"_PbPb"));
  outf->WriteObject(&y_fitErr,"FinalCorrectedYield_fitError"+(TString)(centDep?"_centralityDep":"")+(TString)(ispp?"_pp":"_PbPb"));
  outf->WriteObject(&y_TnPrelErr,"TagAndProbe_relError"+(TString)(centDep?"_centralityDep":"")+(TString)(ispp?"_pp":"_PbPb"));
  outf->Close();

}

void Draw_corrYields(bool secondStep=true, bool thirdStep=false){
  cout<<"\n******** pp pT dependence *******\n"<<endl;
  corrYields(true, secondStep, thirdStep, false);
  cout<<"\n******** PbPb pT dependence *******\n"<<endl;
  corrYields(false, secondStep, thirdStep, false);
  cout<<"\n******** PbPb centrality dependence *******\n"<<endl;
  corrYields(false, secondStep, thirdStep, true);
}
