#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH1D.h>

namespace npd
{
  TGraphAsymmErrors* grae2;
  TH1D *hNuclearModification2;

  void canvasRAAPbPb_0_100_BtoDRAA_marker()
  {
    hNuclearModification2->Draw("pe same");  
  }

  void canvasRAAPbPb_0_100_BtoDRAA()
  {
    //=========Macro generated from canvas: canvasRAA/canvasRAA
    //=========  (Mon Apr 10 17:52:34 2017) by ROOT version6.02/10

    Int_t ci;      // for color index setting
    TColor *color; // for color definition with alpha

    Double_t gNuclearModification_fx3002[12] = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 16, 30, 50, 80};
    Double_t gNuclearModification_fy3002[12] = {0.757746, 0.566197, 0.532394, 0.569014, 0.512676, 0.44507, 0.504225, 0.552113, 0.538028, 0.552113, 0.605634, 0.757746};
    Double_t gNuclearModification_felx3002[12] = {0.2, 0.28, 0.36, 0.44, 0.52, 0.6, 0.72, 0.88, 1.28, 2.4, 4, 6.4};
    Double_t gNuclearModification_fehx3002[12] = {0.2, 0.28, 0.36, 0.44, 0.52, 0.6, 0.72, 0.88, 1.28, 2.4, 4, 6.4};
    Double_t gNuclearModification_fely3002[12] = {0.166197, 0.115493, 0.104225, 0.112676, 0.0985915, 0.084507, 0.101408, 0.112676, 0.121127, 0.146479, 0.177465, 0.28169};
    Double_t gNuclearModification_fehy3002[12] = {0.166197, 0.115493, 0.104225, 0.112676, 0.0985915, 0.084507, 0.101408, 0.112676, 0.121127, 0.146479, 0.177465, 0.28169};

    grae2 = new TGraphAsymmErrors(12,gNuclearModification_fx3002,gNuclearModification_fy3002,gNuclearModification_felx3002,gNuclearModification_fehx3002,gNuclearModification_fely3002,gNuclearModification_fehy3002);
    grae2->SetName("gBNuclearModification2");
    // grae2->SetTitle("Graph");
    grae2->SetFillColor(kViolet-8);
    grae2->SetFillColorAlpha(kViolet-8, 0.4);

    ci = kViolet-7;
    grae2->SetMarkerColor(ci);
    grae2->SetMarkerStyle(47); // 22
    grae2->SetMarkerSize(1.4); // 1.2
    grae2->SetLineColor(0);
    grae2->SetLineWidth(0);
   
    TH1F* Graph_gNuclearModification3002 = new TH1F("Graph_gNuclearModification3002","Graph",100,2.7,54.3);
    Graph_gNuclearModification3002->SetMinimum(0.2364759);
    Graph_gNuclearModification3002->SetMaximum(0.7605233);
    Graph_gNuclearModification3002->SetDirectory(0);
    Graph_gNuclearModification3002->SetStats(0);

    ci = TColor::GetColor("#000099");
    Graph_gNuclearModification3002->SetLineColor(ci);
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
    grae2->SetHistogram(Graph_gNuclearModification3002);
   
    grae2->Draw("2 same");
    Double_t xAxis1[13] = {2., 3., 4., 5., 6., 7., 8., 10., 12., 20., 40., 60., 100.};
   
    hNuclearModification2 = new TH1D("hNuclearModification2a","",12, xAxis1);
    hNuclearModification2->SetBinContent(1,0.757746);
    hNuclearModification2->SetBinContent(2,0.566197);
    hNuclearModification2->SetBinContent(3,0.532394);
    hNuclearModification2->SetBinContent(4,0.569014);
    hNuclearModification2->SetBinContent(5,0.512676);
    hNuclearModification2->SetBinContent(6,0.44507);
    hNuclearModification2->SetBinContent(7,0.504225);
    hNuclearModification2->SetBinContent(8,0.552113);
    hNuclearModification2->SetBinContent(9,0.538028);
    hNuclearModification2->SetBinContent(10,0.552113);
    hNuclearModification2->SetBinContent(11,0.605634);
    hNuclearModification2->SetBinContent(12,0.757746);
    hNuclearModification2->SetBinError(1,0.129577);
    hNuclearModification2->SetBinError(2,0.0676056);
    hNuclearModification2->SetBinError(3,0.0507042);
    hNuclearModification2->SetBinError(4,0.0591549);
    hNuclearModification2->SetBinError(5,0.0619718);
    hNuclearModification2->SetBinError(6,0.0535211);
    hNuclearModification2->SetBinError(7,0.0647887);
    hNuclearModification2->SetBinError(8,0.107042);
    hNuclearModification2->SetBinError(9,0.0985915);
    hNuclearModification2->SetBinError(10,0.0985915);
    hNuclearModification2->SetBinError(11,0.138028);
    hNuclearModification2->SetBinError(12,0.312676);

    ci = TColor::GetColor("#0033cc");
    hNuclearModification2->SetLineColor(ci);
    hNuclearModification2->SetLineWidth(1);

    ci = TColor::GetColor("#0033cc");
    hNuclearModification2->SetMarkerColor(ci);
    hNuclearModification2->SetMarkerStyle(47); // 22
    hNuclearModification2->SetMarkerSize(1.4); // 1.2
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
