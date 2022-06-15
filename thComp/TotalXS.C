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

void TotalXS(){

  gStyle->SetOptStat(0);
  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);
  TVirtualFitter::SetMaxIterations( 10000 );//default is 5000
  gRandom = new TRandom3(234); //some seed give one or two doubly-failed fits

  TF1* P2Log = new TF1("P2Log", "[0]*(x**( [1]+[2]*TMath::Log(x) ))",0,50);
  P2Log->SetParameter(1, 1.88);
  P2Log->SetParameter(0, 77.91/2);
  P2Log->SetParameter(2, -1.239);

  cout<<"Integral of P2Log = "<<P2Log->Integral(0,50)<<endl;

  TF1* Kaplan = new TF1("Kaplan", "[0] / (1+ (x/[1])**2 )**[2]",0,50);
  // Kaplan->SetParameter(1, 7.72);
  // Kaplan->SetParameter(0,178.5/2);
  // Kaplan->SetParameter(2, 3.1);
  Kaplan->SetParameter(1, 8.84);
  Kaplan->SetParameter(0,61.3/2);
  Kaplan->SetParameter(2, 2.96);

  cout<<"Integral of Kaplan = "<<Kaplan->Integral(0,50)<<endl;

}
