#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../helpers/Tools.h"
using namespace ROOT;

//!!! THIS MACRO NEEDS ROOT6.14+, FOR RDATAFRAME, AND NO CMSSW ENVIRONMENT !
void drawChecks(bool ispp, float GOFdataval, string seed);

void drawFitChecks(){
  cout<<"\n______________ PP _______________________\n"<<endl;
  drawChecks(true,138.4,"2235");
  cout<<"\n______________ PBPB _______________________\n"<<endl;
  drawChecks(false,56.4,"2235");
}

void drawChecks(bool ispp=true, float GOFdataval=100, string seed="2235"){

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);
  string s_ispp = (ispp?"pp":"PbPb");

  //********************************
  //FIT RESULTS
  auto filename3 = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+s_ispp+"_2bins.root";  
  auto fFitRes = new TFile((TString)filename3, "read");
  
  //postfit signal yields
  RooArgSet *Yields = (RooArgSet*)fFitRes->Get("norm_fit_s");
  double pf_yields[_NanaBins][_nChan(ispp)+1];
  pf_yields[0][0] = 0;
  for(int b=1;b<=_NanaBins;b++){
    pf_yields[b][0] = 0;
    for(int k=1;k<=_nChan(ispp);k++){      
      pf_yields[b][k] = Yields->getRealValue("BDT"+(TString)(to_string(k))+"Kin"+(TString)(to_string(b))+"/BcSig");
      pf_yields[b][0] += pf_yields[b][k];
    }
    pf_yields[0][0] += pf_yields[b][0];
  }

  //postfit r1 and r2 with errors
  RooFitResult *fitres = (RooFitResult*)fFitRes->Get("fit_s");
  RooArgSet fitargs = (RooArgSet)fitres->floatParsFinal();
  RooRealVar *r1 = (RooRealVar*)(((RooArgSet)fitres->floatParsFinal())["r1"]).Clone();
  RooRealVar *r2 = (RooRealVar*)(((RooArgSet)fitres->floatParsFinal())["r2"]).Clone();
  double r1errlo = fabs(r1->getErrorLo()), r1errhi = fabs(r1->getErrorHi());
  double r2errlo = fabs(r2->getErrorLo()), r2errhi = fabs(r2->getErrorHi());
  cout<<"1st bin signal normalisation = "<<r1->getValV()<<" + "<<r1errhi<<" - "<<r1errlo<<endl;
  cout<<"2nd bin signal normalisation = "<<r2->getValV()<<" + "<<r2errhi<<" - "<<r2errlo<<endl;
  double r1relerr = (r1errlo+r1errhi)/(2*r1->getValV());
  double r2relerr = (r2errlo+r2errhi)/(2*r2->getValV());
  cout<<"relative errors on r1 and r2 = "<<r1relerr<<" "<<r2relerr<<endl;
  cout<<"Post-fit signal yields in bin 1 and 2 = "<<pf_yields[1][0]<<" +- "<<pf_yields[1][0]*r1relerr<<" and "<<pf_yields[2][0]<<" +- "<<pf_yields[2][0]*r2relerr<<endl;
  double rho = fitres->correlation("r1","r2");
  cout<<"correlation factor between r1 and r2 = "<<rho<<endl;
  cout<<endl;

  //********************************
  //NLL SCAN
  gStyle->SetOptStat(0);

  auto filename2 = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/higgsCombine_"+s_ispp+".MultiDimFit.mH120.root";
  cout<<"NLL SCAN : read file "<<filename2<<endl;
  ROOT::RDataFrame d_toys("limit", filename2); //filename cannot handle TString pfffrr....

  //NLL scan
  auto d2_toys = d_toys.Define("chi2", "2*deltaNLL");
  auto NLL = d2_toys.Filter("chi2<25").Profile2D({"nll","",25,static_cast<double>(ispp?0.45:0.),static_cast<double>(ispp?1.25:2.8),
	                                                   25,static_cast<double>(ispp?0.8:0.),static_cast<double>(ispp?1.6:1.1)},
                                                  "r1","r2","chi2");
  auto NLLv2 = d2_toys.Profile2D(      {"nll","",25,static_cast<double>(ispp?0.5:0.),static_cast<double>(ispp?1.3:2.7),
	                         25,static_cast<double>(ispp?0.7:0.),static_cast<double>(ispp?1.5:1.2)},
    "r1","r2","chi2");

  //Evaluate significance by extrapolating NLL to (r1=0,r2=0)
  double x1 = (NLLv2->GetXaxis())->GetBinCenter(1); double y1 = (NLLv2->GetYaxis())->GetBinCenter(1); double z1 = NLLv2->GetBinContent(1,1);
  double x2 = (NLLv2->GetXaxis())->GetBinCenter(2); double y2 = (NLLv2->GetYaxis())->GetBinCenter(2); double z2 = NLLv2->GetBinContent(2,2);
  //cout<<"1st bin from 0, center x, center y, nll value = "<<x1<<" "<<y1<<" "<<z1<<endl;
  //cout<<"2nd bin from 0, center x, center y, nll value = "<<x2<<" "<<y2<<" "<<z2<<endl;
  double nllAt0v1 = z1 + (z1-z2) * x1/(x2-x1);
  double nllAt0v2 = z1 + (z1-z2) * y1/(y2-y1); //  f(x0) = f(x1) + (x0-x1) * (f(x2)-f(x1)) / x2-x1
  double nllAt0 = (nllAt0v1+nllAt0v2)/2;
  cout<<"nll at (r1=0,r2=0) v1,v2,average = "<<nllAt0v1<<" "<<nllAt0v2<<" "<<nllAt0<<endl;
  cout<<"significance of non-zero signal (assuming gaussian likelihood) = "<<sqrt(nllAt0)<<endl;

  //Best fit graph
  float xBF[] = {(*d_toys.Filter("quantileExpected == -1").Take<float>("r1"))[0]}; //dereference, then transform the vector<double> into array for TGraph
  float yBF[] = {(*d_toys.Filter("quantileExpected == -1").Take<float>("r2"))[0]};
  auto bestFit = new TGraph(1,xBF,yBF);
  bestFit->SetMarkerSize(5);
  bestFit->SetMarkerStyle(34);

  //Draw NLL scan
  TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
  c2->SetRightMargin(0.16);
  c2->SetTopMargin(0.03);

  NLL->GetXaxis()->SetTitle("r1");
  NLL->GetYaxis()->SetTitle("r2");
  NLL->GetZaxis()->SetTitle("2 #times #DeltaNLL");
  NLL->SetMaximum(25);
  NLL->SetContour(25);

  NLL->Draw("COLZ");
  bestFit->Draw("p same");

  //moving z color palette
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)NLL->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.85);
  palette->SetX2NDC(0.9);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.97);
  gPad->Modified();
  gPad->Update();

  c2->SaveAs("figs/NLLscan_"+(TString)s_ispp+".pdf");

  cout<<endl;

  //********************************
  //GOODNESS OF FIT

  gStyle->SetOptStat("emr");
  gStyle->SetStatY(0.97);
  gStyle->SetStatH(0.2);
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.3);

  auto filename = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/higgsCombine_"+s_ispp+".GoodnessOfFit.mH120."+seed+".root";  
  cout<<"GOODNESS OF FIT : read file "<<filename<<endl;
  ROOT::RDataFrame d_gfit("limit", filename); //filename cannot handle TString pfffrr....
  auto testStat = d_gfit.Histo1D({"","",40,static_cast<double>(ispp?70:23),static_cast<double>(ispp?168:83)},"limit");

  testStat->SetLineWidth(3);
  testStat->GetXaxis()->SetTitle("test statistic");
  testStat->GetYaxis()->SetTitle("number of toys");

  TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
  c1->SetTopMargin(0.03);
  c1->SetRightMargin(0.05);
  testStat->Draw("hist");
  float maxpad = testStat->GetMaximum();

  TLatex t1(0.9*GOFdataval,.52*maxpad,"data value");
  t1.SetTextSize(0.05);
  t1.Draw("same");

  TArrow *ar = new TArrow(GOFdataval,0.03*maxpad,GOFdataval,0.5*maxpad,0.05,"<");
  ar->SetLineWidth(3);
  float maxX = testStat->GetXaxis()->GetXmax();
  if(GOFdataval>maxX) ar->DrawArrow(0.75*maxX,0.2*maxpad, 0.995*maxX,0.2*maxpad);
  else ar->Draw();

  double integral1 = testStat->Integral(testStat->FindBin(GOFdataval), testStat->GetNbinsX());
  double integral2 = testStat->Integral(testStat->FindBin(GOFdataval)+1, testStat->GetNbinsX());
  //  cout<<"p-value(prob of toys having testStat>(dataValue_"+s_ispp+"="<<GOFdataval<<") ) = "<< integral1 / testStat->Integral() <<endl;
  //cout<<"or "<< integral2 / testStat->Integral() <<endl;
  cout<<"p-value (average) (prob of toys having testStat>(dataValue_"+s_ispp+"="<<GOFdataval<<") ) = "<< (integral1+integral2)/(2*testStat->Integral())<<endl;
  c1->SaveAs("figs/GOF_"+(TString)s_ispp+".pdf");

  cout<<endl;

  //********************************
  //SIGNAL TOYS

  auto filename4 = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+s_ispp+"_2bins_wSigToys.root";
  cout<<"SIGNAL TOYS : read file "<<filename4<<endl;
  ROOT::RDataFrame d_sigtoys("tree_fit_sb", filename4); //filename cannot handle TString pfffrr....
  
  auto r1r2 = d_sigtoys.Histo2D({"",";r1;r2;",25,static_cast<double>(ispp?0.6:0.),static_cast<double>(ispp?1.1:2.3),25,static_cast<double>(ispp?0.95:0.15),static_cast<double>(ispp?1.45:0.8)},"r1","r2");
  auto r1err = d_sigtoys.Histo1D({"",";r1 uncertainty;n toys",30,static_cast<double>(ispp?0.068:0.22),static_cast<double>(ispp?0.084:0.43)},"r1Err");
  auto r2err = d_sigtoys.Histo1D({"",";r2 uncertainty;n toys",30,static_cast<double>(ispp?0.052:0.065),static_cast<double>(ispp?0.063:0.11)},"r2Err");
				 
  TCanvas *c3 = new TCanvas("c3","c3",4500,1500);
  c3->Divide(3,1);

  c3->cd(2)->SetRightMargin(0.05);
  c3->cd(2)->SetTopMargin(0.03);
  r1err->SetLineWidth(2);
  r1err->Draw("hist");
  
  c3->cd(3)->SetRightMargin(0.05);
  c3->cd(3)->SetTopMargin(0.03);
  r2err->SetLineWidth(2);
  r2err->Draw("hist");

  gStyle->SetStatH(0.15);

  c3->cd(1)->SetRightMargin(0.05);
  c3->cd(1)->SetTopMargin(0.03);
  c3->cd(1)->SetLeftMargin(0.12);
  r1r2->SetFillColor(kRed);
  r1r2->Draw("BOX");

  cout<<"r1-r2 correlation factor from the signal toys = "<<r1r2->GetCorrelationFactor();

  c3->SaveAs("figs/SignalToys_POIvalAndErrors_"+(TString)s_ispp+".pdf");

  cout<<endl;
  
}
