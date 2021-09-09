#include <TGraphAsymmErrors.h>
#include <TH1D.h>

namespace bplus
{
  TGraphAsymmErrors* grae;
  TH1D *hNuclearModification2;

  void canvasRAAPbPb_0_100_ThmRAA_marker(bool noXerr=false)
  {
    hNuclearModification2->Draw(noXerr?"peX0 same":"pe same");
  }

  void canvasRAAPbPb_0_100_ThmRAA()
  {
    //=========Macro generated from canvas: canvasRAA/canvasRAA
    //=========  (Mon Apr 10 17:52:34 2017) by ROOT version6.02/10

    Int_t ci;      // for color index setting
    TColor *color; // for color definition with alpha

    Double_t gNuclearModification_fx3002[5] = {
      8.5,
      12.5,
      17.5,
      25,
      40};
    Double_t gNuclearModification_fy3002[5] = {
      0.3458558,
      0.4480066,
      0.4402087,
      0.615181,
      0.3486562};
    Double_t gNuclearModification_felx3002[5] = {
      1.5,
      2.5,
      2.5,
      5,
      10};
    Double_t gNuclearModification_fely3002[5] = {
      0.06570926,
      0.07712284,
      0.07409171,
      0.1016717,
      0.05734825};
    Double_t gNuclearModification_fehx3002[5] = {
      1.5,
      2.5,
      2.5,
      5,
      10};
    Double_t gNuclearModification_fehy3002[5] = {
      0.06570926,
      0.07712284,
      0.07409171,
      0.1016717,
      0.05734825};
   
    grae = new TGraphAsymmErrors(5,gNuclearModification_fx3002,gNuclearModification_fy3002,gNuclearModification_felx3002,gNuclearModification_fehx3002,gNuclearModification_fely3002,gNuclearModification_fehy3002);
    grae->SetName("gBNuclearModification");
    // grae->SetTitle("Graph");
    grae->SetFillColor(kViolet-8);//kAzure+7
    grae->SetFillColorAlpha(kViolet-8, 0.5);//kAzure+7

    ci = TColor::GetColor("#0033cc");
    grae->SetMarkerColor(kViolet-7);//ci
    grae->SetMarkerStyle(33); // 22
    grae->SetMarkerSize(1.6); // 1.2
    grae->SetLineColor(kViolet-7);//ci
    grae->SetLineWidth(0);
   
    TH1F* Graph_gNuclearModification3002 = new TH1F("Graph_gNuclearModification3002","Graph",100,2.7,54.3);
    Graph_gNuclearModification3002->SetMinimum(0.2364759);
    Graph_gNuclearModification3002->SetMaximum(0.7605233);
    Graph_gNuclearModification3002->SetDirectory(0);
    Graph_gNuclearModification3002->SetStats(0);

    ci = TColor::GetColor("#000099");
    Graph_gNuclearModification3002->SetLineColor(kViolet-7);//ci
    Graph_gNuclearModification3002->SetMarkerStyle(20);
    Graph_gNuclearModification3002->GetXaxis()->SetLabelFont(42);
    Graph_gNuclearModification3002->GetXaxis()->SetLabelSize(0.035);
    Graph_gNuclearModification3002->GetXaxis()->SetTitleSize(0.035);
    Graph_gNuclearModification3002->GetXaxis()->SetTitleFont(42);
    Graph_gNuclearModification3002->GetYaxis()->SetLabelFont(42);
    Graph_gNuclearModification3002->GetYaxis()->SetLabelSize(0.035);
    Graph_gNuclearModification3002->GetYaxis()->SetTitleSize(0.035);
    Graph_gNuclearModification3002->GetYaxis()->SetTitleFont(42);
    Graph_gNuclearModification3002->GetZaxis()->SetLabelFont(42);
    Graph_gNuclearModification3002->GetZaxis()->SetLabelSize(0.035);
    Graph_gNuclearModification3002->GetZaxis()->SetTitleSize(0.035);
    Graph_gNuclearModification3002->GetZaxis()->SetTitleFont(42);
    grae->SetHistogram(Graph_gNuclearModification3002);
   
    grae->Draw("2 same");
    Double_t xAxis1[6] = {7, 10, 15, 20, 30, 50}; 
   
    hNuclearModification2 = new TH1D("hNuclearModification2","",5, xAxis1);
    hNuclearModification2->SetBinContent(1,0.3458558);
    hNuclearModification2->SetBinContent(2,0.4480066);
    hNuclearModification2->SetBinContent(3,0.4402087);
    hNuclearModification2->SetBinContent(4,0.615181);
    hNuclearModification2->SetBinContent(5,0.3486562);
    hNuclearModification2->SetBinError(1,0.1086955);
    hNuclearModification2->SetBinError(2,0.07447295);
    hNuclearModification2->SetBinError(3,0.07507186);
    hNuclearModification2->SetBinError(4,0.09192405);
    hNuclearModification2->SetBinError(5,0.1120147);
    hNuclearModification2->SetEntries(109.8058);

    ci = TColor::GetColor("#0033cc");
    hNuclearModification2->SetLineColor(kViolet-7);//ci
    hNuclearModification2->SetLineWidth(3);

    ci = TColor::GetColor("#0033cc");
    hNuclearModification2->SetMarkerColor(kViolet-7);//ci
    hNuclearModification2->SetMarkerStyle(33); // 22
    hNuclearModification2->SetMarkerSize(1.6); // 1.2
    hNuclearModification2->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
    hNuclearModification2->GetXaxis()->SetLabelFont(42);
    hNuclearModification2->GetXaxis()->SetLabelSize(0.035);
    hNuclearModification2->GetXaxis()->SetTitleSize(0.035);
    hNuclearModification2->GetXaxis()->SetTitleFont(42);
    hNuclearModification2->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
    hNuclearModification2->GetYaxis()->SetLabelFont(42);
    hNuclearModification2->GetYaxis()->SetLabelSize(0.035);
    hNuclearModification2->GetYaxis()->SetTitleSize(0.035);
    hNuclearModification2->GetYaxis()->SetTitleFont(42);
    hNuclearModification2->GetZaxis()->SetLabelFont(42);
    hNuclearModification2->GetZaxis()->SetLabelSize(0.035);
    hNuclearModification2->GetZaxis()->SetTitleSize(0.035);
    hNuclearModification2->GetZaxis()->SetTitleFont(42);
    // hNuclearModification2->Draw("pe same");

  }
}
