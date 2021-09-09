#include <TGraphAsymmErrors.h>
#include <TH1D.h>

namespace Bs
{
  TGraphAsymmErrors* grStat = new TGraphAsymmErrors(2);
  void RAA_markers(){
    grStat->Draw("psame");    
  }

  void RAA(bool noXerr=false)
  {
    //=========Macro generated from canvas: c2/RAA
    //=========  (Tue May 11 00:30:31 2021) by ROOT version 6.22/03
   
    grStat->SetName("Bs_graphStat");
    grStat->SetTitle("");

    grStat->SetLineColor(kRed+2);
    grStat->SetLineWidth(3);
    grStat->SetMarkerColor(kRed+2);
    grStat->SetMarkerStyle(34);
    grStat->SetMarkerSize(1.5);
    grStat->SetFillColor(kRed-9);
    grStat->SetFillColorAlpha(kRed-9,0.5);

    grStat->SetPoint(0, 11., 1.51);
    grStat->SetPointError(0, 4., 4., 0.61, 0.61);
    //    grStat->SetPointEY(0, 0.61, 0.61);
    grStat->SetPoint(1, 32.5, 0.87);
    grStat->SetPointError(1, 17.5, 17.5, 0.30, 0.30);
    //grStat->SetPointEY(1, 0.30, 0.30);
   
    TGraphAsymmErrors* grSyst = new TGraphAsymmErrors(2);
    grSyst->SetName("Bs_graphSyst");
    grSyst->SetTitle("");

    grSyst->SetLineWidth(0);
    grSyst->SetMarkerColor(kRed+2);
    grSyst->SetMarkerStyle(34);
    grSyst->SetMarkerSize(1.5);
    grSyst->SetFillColor(kRed-9);
    grSyst->SetFillColorAlpha(kRed-9,0.5);

    grSyst->SetPoint(0, 11., 1.51);
    grSyst->SetPointError(0, 4.,4., 0.50, 0.50);
    //    grSyst->SetPointEY(0, 0.50, 0.50);
    grSyst->SetPoint(1, 32.5, 0.87);
    grSyst->SetPointError(1, 17.5, 17.5, 0.17, 0.17);
    //grSyst->SetPointEY(1, 0.17, 0.17);

    if(noXerr){
      SetEx(grStat,0);
    }
    
    grSyst->Draw("2same");    

  }

}
