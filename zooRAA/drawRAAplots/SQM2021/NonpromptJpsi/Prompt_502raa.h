// https://www.hepdata.net/record/ins1644903

#include <TGraphErrors.h>

namespace pjpsi5
{
  TGraphErrors* gre5TeV;
  TGraphErrors* gre5TeV2;
  TGraphErrors* gre5TeVstat;
  TGraphErrors* gre5TeV2stat;

  Int_t cifillJ;

  void expBeautyCMS_marker()
  {
    gre5TeV2stat->Draw("pe same");
    gre5TeVstat->Draw("pe same");
  }

  void expBeautyCMS()
  {
    //=========Macro generated from canvas: c2/c2
    //=========  (Wed Feb  1 19:16:45 2017) by ROOT version5.34/36

    Int_t cimarker = kOrange+2;
    Int_t cifill = kOrange+1;

    cifillJ = cifill;

    Double_t bin_3_syst_fx3011[3] = {
      3.75,
      5,
      6};
    Double_t bin_3_syst_fy3011[3] = {
      0.553,
      0.365,
      0.328};
    Double_t bin_3_syst_felx3011[3] = {
      3.75*0.08,
      5*0.08,
      6*0.08};
    Double_t bin_3_syst_fely3011[3] = {
      0.113,
      0.05,
      0.036};

    gre5TeV2 = new TGraphErrors(3,bin_3_syst_fx3011,bin_3_syst_fy3011,bin_3_syst_felx3011,bin_3_syst_fely3011);
    gre5TeV2->SetName("Graph1");
    gre5TeV2->SetTitle("Graph1");
    gre5TeV2->SetFillColor(cifill);
    gre5TeV2->SetFillColorAlpha(cifill, 0.5);
    gre5TeV2->SetMarkerColor(cimarker);
    gre5TeV2->SetMarkerStyle(34);
    gre5TeV2->SetMarkerSize(1.7);
    gre5TeV2->SetLineColor(0);

    gre5TeV2->Draw("2 same");
  
    Double_t bin_3_fx3012[3] = {
      3.75,
      5,
      6};
    Double_t bin_3_fy3012[3] = {
      0.553,
      0.365,
      0.328};
    Double_t bin_3_felx3012[3] = {
      0.75,
      0.5,
      0.5};
    Double_t bin_3_fely3012[3] = {
      0.047,
      0.024,
      0.022};

    gre5TeV2stat = new TGraphErrors(3,bin_3_fx3012,bin_3_fy3012,bin_3_felx3012,bin_3_fely3012);
    gre5TeV2stat->SetName("Graph1stat");
    gre5TeV2stat->SetTitle("Graph1stat");
    gre5TeV2stat->SetMarkerColor(cimarker);
    gre5TeV2stat->SetLineColor(cimarker);
    gre5TeV2stat->SetMarkerStyle(34);
    gre5TeV2stat->SetMarkerSize(1.7);

    // gre5TeV2stat->Draw("p same");
  
    Double_t bin_0_syst_fx3013[12] = {
      7,
      8,
      9,
      10.25,
      12,
      14,
      16.25,
      18.75,
      22.5,
      27.5,
      32.5,
      42.5};
    Double_t bin_0_syst_fy3013[12] = {
      0.361,
      0.348,
      0.334,
      0.33,
      0.339,
      0.36,
      0.36,
      0.364,
      0.425,
      0.521,
      0.495,
      0.521};
    Double_t bin_0_syst_felx3013[12] = {
      7*0.08,
      8*0.08,
      9*0.08,
      10.25*0.08,
      12*0.08,
      14*0.08,
      16.25*0.08,
      18.75*0.08,
      22.5*0.08,
      27.5*0.08,
      32.5*0.08,
      42.5*0.08};
    Double_t bin_0_syst_fely3013[12] = {
      0.035,
      0.027,
      0.022,
      0.019,
      0.019,
      0.019,
      0.019,
      0.023,
      0.023,
      0.035,
      0.048,
      0.047};
    gre5TeV = new TGraphErrors(12,bin_0_syst_fx3013,bin_0_syst_fy3013,bin_0_syst_felx3013,bin_0_syst_fely3013);
    gre5TeV->SetName("Graph0");
    gre5TeV->SetTitle("Graph0");
  
    gre5TeV->SetFillColor(cifill);
    gre5TeV->SetFillColorAlpha(cifill, 0.5);
    gre5TeV->SetMarkerColor(cimarker);
    gre5TeV->SetMarkerStyle(29);
    gre5TeV->SetMarkerSize(2);
    gre5TeV->SetLineColor(0);

    gre5TeV->Draw("2 same");
  
    Double_t bin_0_fx3014[12] = {
      7,
      8,
      9,
      10.25,
      12,
      14,
      16.25,
      18.75,
      22.5,
      27.5,
      32.5,
      42.5};
    Double_t bin_0_fy3014[12] = {
      0.361,
      0.348,
      0.334,
      0.33,
      0.339,
      0.36,
      0.36,
      0.364,
      0.425,
      0.521,
      0.495,
      0.521};
    Double_t bin_0_felx3014[12] = {
      0.5,
      0.5,
      0.5,
      0.75,
      1,
      1,
      1.25,
      1.25,
      2.5,
      2.5,
      2.5,
      7.5};
    Double_t bin_0_fely3014[12] = {
      0.012,
      0.01,
      0.009,
      0.008,
      0.009,
      0.011,
      0.013,
      0.017,
      0.021,
      0.039,
      0.055,
      0.068};

    gre5TeVstat = new TGraphErrors(12,bin_0_fx3014,bin_0_fy3014,bin_0_felx3014,bin_0_fely3014);
    gre5TeVstat->SetName("Graph0stat");
    gre5TeVstat->SetTitle("Graph0stat");
    gre5TeVstat->SetLineColor(cimarker);
    gre5TeVstat->SetMarkerColor(cimarker);
    gre5TeVstat->SetMarkerStyle(29);
    gre5TeVstat->SetMarkerSize(2);

    // gre5TeVstat->Draw("p same");
  
  }

}
