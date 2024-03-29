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
#include "../helpers/hub.cpp"
#include "../helpers/Cuts_BDT.h"
using namespace ROOT;

//!!! THIS MACRO NEEDS ROOT6.14+, FOR RDATAFRAME, AND NO CMSSW ENVIRONMENT !
void drawChecks(bool ispp, bool secondStep, float GOFdataval, string seed);

void drawFitChecks(bool secondStep=true){
  cout<<"\n______________ PP _______________________\n"<<endl;
  drawChecks(true,secondStep,100.9,"2235");
  cout<<"\n______________ PBPB _______________________\n"<<endl;
  drawChecks(false,secondStep,47.9,"2235");
}

void drawChecks(bool ispp=true, bool secondStep=true, float GOFdataval=100, string seed="2235"){

  Hub H = Hub(secondStep,secondStep);
  // H.SetFit(true);
  H.SetMetafit();
  // H.SetAccEffFinal(true);
  // H.SetAccEff();
  // H.SetTnP();
  // H.SetBcTau();
  // H.SetLumi();
  // H.SetxLW();
  // H.SetpTbias();
  // H.SetFullErr();
  // H.ScaleByLumi();

  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);
  string s_ispp = (ispp?"pp":"PbPb");

  //********************************
  //FIT RESULTS
  auto filename3 = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+s_ispp+"_2bins"+(TString)(secondStep?"_2ndStep":"")+".root";  
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

  //NLL at 0,0
  float NLL0;
  double p0,fitsignif;
  if(!ispp){
    auto filename5 = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/higgsCombine_PbPb_NLL0_2ndStep.MultiDimFit.mH120.root";
    cout<<"NLL at (r1=0,r2=0) : read file "<<filename5<<endl;
    auto f5 = new TFile(filename5,"read");
    TTree* T = (TTree*)f5->Get("limit"); 

    T->SetBranchAddress("deltaNLL",&NLL0);
    T->GetEntry(1);
    cout<<"2*NLL at (r1=0,r2=0) = "<<2*NLL0<<endl;
    p0 = TMath::Prob(2*NLL0,2)/2;
    fitsignif = ROOT::Math::gaussian_quantile(1-p0,1);//sigma=1
    cout<<"associated 1-sided pvalue and significance of non-zero signal (assuming gaussian likelihood) (fit only) = "<<p0<<" "<<fitsignif<<endl;
  }

  //NLL Scan
  auto filename2 = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/higgsCombine_"+s_ispp+".MultiDimFit.mH120.root";
  cout<<"NLL SCAN : read file "<<filename2<<endl;
  ROOT::RDataFrame d_toys("limit", filename2); //filename cannot handle TString pfffrr....

  //Best fit graph
  float xBF[] = {(*d_toys.Filter("quantileExpected == -1").Take<float>("r1"))[0]}; //dereference, then transform the vector<double> into array for TGraph
  float yBF[] = {(*d_toys.Filter("quantileExpected == -1").Take<float>("r2"))[0]};
  float one[] = {1.};
  auto bestFit = new TGraph(1,one,one);//xBF,yBF);
  bestFit->SetMarkerSize(5);
  bestFit->SetMarkerStyle(34);

  //NLL scan
  auto d2_toys = d_toys.Define("chi2", "2*deltaNLL").Define("scaledr1",(const char*)TString::Format("r1/%.4f",xBF[0])).Define("scaledr2",(const char*)TString::Format("r2/%.4f",yBF[0]));
  auto NLL = d2_toys.Filter("chi2<28.7").Profile2D({"nll","",30,static_cast<double>(ispp?0.75:0.)/xBF[0],static_cast<double>(ispp?1.65:3.)/xBF[0],
	                                                   30,static_cast<double>(ispp?0.6:0.15)/yBF[0],static_cast<double>(ispp?1.5:2.3)/yBF[0]},
                                                  "scaledr1","scaledr2","chi2");


  //Extract the contribution of metafit error to the significance
  vector<double> metafitRErrLo = {H.metafit.PbPb_pt.RelErrLo[1], H.metafit.PbPb_pt.RelErrLo[2]};
  double metafitCorr = H.metafit.PbPb_pt.Corr[0];
  double metafitNLL = (1/(1-pow(metafitCorr,2))) * ( 1/pow(metafitRErrLo[0],2) + 1/pow(metafitRErrLo[1],2) - 2*metafitCorr/(metafitRErrLo[0]*metafitRErrLo[1]) ); //2*NLL(r1,r2=0,0) //actually just -2log(Gaussian2D)
  double metafitP0 = TMath::Prob(metafitNLL,2)/2;
  if(!ispp){
    cout<<"metafit relative uncertainty (low) for PbPb, pT bin 1,2 = "<<metafitRErrLo[0]<<" "<<metafitRErrLo[1]<<endl;
    cout<<"        bin-to-bin correlation = "<<metafitCorr<<endl;
    cout<<"metafit error gives (2*NLL) = "<<metafitNLL<<endl;

    cout<<"p-value of fit and from metafit = "<<p0<<" "<<metafitP0<<endl;
    double p0Comb = 1-(1-p0)*(1-metafitP0);
    cout<<"combined p-value = "<<p0Comb<<endl;
    cout<<"significance of non-zero signal (assuming gaussian likelihood) including metafit uncertainty from pT dependence = "<<ROOT::Math::gaussian_quantile(1-p0Comb,1)<<endl;
    //    cout<<"significance of non-zero signal (assuming gaussian likelihood) including metafit uncertainty from pT dependence = "<<1/sqrt(1/nllAt0+1/metafitNLL)<<endl;
  }

  gStyle->SetPalette(kLightTemperature);

  //Draw NLL scan
  TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
  c2->SetRightMargin(0.17);
  c2->SetTopMargin(0.03);

  NLL->GetXaxis()->SetTitle("r1/r1_{best}");
  NLL->GetYaxis()->SetTitle("r2/r2_{best}");
  if(!ispp) NLL->GetYaxis()->SetRangeUser(0,2.3);
  NLL->GetZaxis()->SetTitle("2 #times #DeltaNLL   ");
  NLL->GetZaxis()->SetTitleOffset(1.4);
  NLL->SetMaximum(28.7);
  NLL->SetContour(100);

  NLL->Draw("COLZ");
  bestFit->Draw("p same");

  //moving z color palette
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)NLL->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.84);
  palette->SetX2NDC(0.89);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.97);
  gPad->Modified();
  gPad->Update();

  TLatex CMStag;
  CMStag.SetNDC();
  CMStag.SetTextFont(42);
  CMStag.SetTextSize(0.045);
  CMStag.SetTextAlign(13);
  CMStag.DrawLatex(0.135,0.955,"#font[61]{CMS}");
  //  CMStag.DrawLatex(0.135,0.91,"#font[52]{Preliminary}");

  TLatex tsigma;
  tsigma.SetNDC();
  tsigma.SetTextFont(42);
  tsigma.SetTextSize(0.035);
  tsigma.SetTextAlign(12);
  tsigma.DrawLatex(0.92,0.1+(2.3/28.7)*(0.97-0.1) +0.015,"#font[61]{1#sigma}");
  tsigma.DrawLatex(0.92,0.1+(2.3/28.7)*(0.97-0.1) -0.019,"#font[61]{#scale[0.8]{(68%)}}");
  tsigma.DrawLatex(0.92,0.1+(6.2/28.7)*(0.97-0.1) +0.015,"#font[61]{2#sigma}");
  tsigma.DrawLatex(0.92,0.1+(6.2/28.7)*(0.97-0.1) -0.019,"#font[61]{#scale[0.8]{(95%)}}");
  tsigma.DrawLatex(0.92,0.1+(11.8/28.7)*(0.97-0.1), "#font[61]{3#sigma}");
  tsigma.DrawLatex(0.92,0.1+(19.3/28.7)*(0.97-0.1), "#font[61]{ 4#sigma}");
  tsigma.DrawLatex(0.92,0.97,"#font[61]{5#sigma}");

  TLine *hor = new TLine;
  hor->SetLineWidth(3);
  hor->SetLineStyle(2);
  hor->DrawLineNDC(0.84,0.1+(2.3/28.7)*(0.97-0.1),0.91,0.1+(2.3/28.7)*(0.97-0.1));
  hor->DrawLineNDC(0.84,0.1+(6.2/28.7)*(0.97-0.1),0.91,0.1+(6.2/28.7)*(0.97-0.1));
  hor->DrawLineNDC(0.84,0.1+(11.8/28.7)*(0.97-0.1),0.91,0.1+(11.8/28.7)*(0.97-0.1));
  hor->DrawLineNDC(0.84,0.1+(19.3/28.7)*(0.97-0.1),0.91,0.1+(19.3/28.7)*(0.97-0.1));
  hor->DrawLineNDC(0.84,0.97,0.91,0.97);

  c2->SaveAs("figs/NLLscan_"+(TString)s_ispp+(TString)(secondStep?"_2ndStep":"")+".pdf");

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
  auto testStat = d_gfit.Histo1D({"","",40,static_cast<double>(ispp?40:24),static_cast<double>(ispp?150:89)},"limit");

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
  c1->SaveAs("figs/GOF_"+(TString)s_ispp+(TString)(secondStep?"_2ndStep":"")+".pdf");

  cout<<endl;

  //********************************
  //SIGNAL TOYS

  auto filename4 = "/home/llr/cms/falmagne/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/fitDiagnostics_"+s_ispp+"_2bins_wSigToys.root";
  cout<<"SIGNAL TOYS : read file "<<filename4<<endl;
  ROOT::RDataFrame d_sigtoys("tree_fit_sb", filename4); //filename cannot handle TString pfffrr....
  
  auto r1r2 = d_sigtoys
    .Define("scaledr1","r1/1.").Define("scaledr2","r2/1.")
    .Histo2D({"",";r1/r1_{true};r2/r2_{true};",25,static_cast<double>(ispp?0.7:0.),static_cast<double>(ispp?1.3:2.5),25,static_cast<double>(ispp?0.8:0.5),static_cast<double>(ispp?1.2:1.5)},"scaledr1","scaledr2");
  auto r1err = d_sigtoys
    .Define("scaledr1Err","r1Err/1.")
    .Histo1D({"",";(r1 uncertainty)/r1_{true};n toys",30,static_cast<double>(ispp?0.092:0.26),static_cast<double>(ispp?0.12:0.55)},"scaledr1Err");
  auto r2err = d_sigtoys
    .Define("scaledr2Err","r2Err/1.")
    .Histo1D({"",";(r2 uncertainty)/r2_{true};n toys",30,static_cast<double>(ispp?0.046:0.14),static_cast<double>(ispp?0.06:0.23)},"scaledr2Err");
				 
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

  c3->SaveAs("figs/SignalToys_POIvalAndErrors_"+(TString)s_ispp+(TString)(secondStep?"_2ndStep":"")+".pdf");

  cout<<endl;
  
}
