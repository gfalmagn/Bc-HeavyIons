#include <TGraphAsymmErrors.h>
#include <TGraphMultiErrors.h>
#include <TH1D.h>

namespace Bc
{
  void RAA()
  {
    //=========Macro generated from canvas: c2/RAA
    //=========  (Tue May 11 00:30:31 2021) by ROOT version 6.22/03

    float transp = 0.65;

    TGraphMultiErrors* tgme = new TGraphMultiErrors(2, 4);
    tgme->SetName("Graph0");
    tgme->SetTitle("");

    Int_t ci;      // for color index setting
    TColor *color; // for color definition with alpha
    ci = TColor::GetColor("#8cd1e0");
    tgme->SetFillColor(ci);
    tgme->SetLineColor(0);
    tgme->SetLineWidth(0);
    tgme->SetMarkerColor(0);
    tgme->SetMarkerStyle(20);
    tgme->SetMarkerSize(1.4);

    ci = TColor::GetColor("#8cd1e0");
    tgme->GetAttFill(0)->SetFillColor(ci);
    //tgme->GetAttFill(0)->SetFillStyle(3001);
    tgme->GetAttFill(0)->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(0)->SetLineColor(ci);
    tgme->GetAttLine(0)->SetLineWidth(3);

    ci = TColor::GetColor("#8cd1e0");
    tgme->GetAttFill(1)->SetFillColor(ci);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(1)->SetLineColor(ci);
    tgme->GetAttLine(1)->SetLineWidth(3);
    tgme->GetAttFill(2)->SetFillColor(19);
    tgme->GetAttFill(2)->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(2)->SetLineColor(ci);
    tgme->GetAttLine(2)->SetLineStyle(2);
    tgme->GetAttLine(2)->SetLineWidth(3);
    tgme->GetAttFill(3)->SetFillColor(19);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(3)->SetLineColor(ci);
    tgme->GetAttLine(3)->SetLineWidth(3);
    tgme->SetPoint(0, 8.059, 1.665756);
    tgme->SetPointEX(0, 2.059, 2.941);
    tgme->SetPointEY(0, 0, 0.5781341, 0.9120825);
    tgme->SetPointEY(0, 1, 0.5259643, 0.851301);
    tgme->SetPointEY(0, 2, 0.5781341, 0.9120825);
    tgme->SetPointEY(0, 3, 0.5259643, 0.851301);
    tgme->SetPoint(1, 17.165, 0.6961609);
    tgme->SetPointEX(1, 6.165001, 17.835);
    tgme->SetPointEY(1, 0, 0.1213511, 0.1403348);
    tgme->SetPointEY(1, 1, 0.1104005, 0.1309829);
    tgme->SetPointEY(1, 2, 0.1213511, 0.1403348);
    tgme->SetPointEY(1, 3, 0.1104005, 0.1309829);
   
    TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","",100,0,35.95);
    Graph_Graph01->SetMinimum(0);
    Graph_Graph01->SetMaximum(3.2);
    Graph_Graph01->SetDirectory(0);
    Graph_Graph01->SetStats(0);

    ci = TColor::GetColor("#000099");
    Graph_Graph01->SetLineColor(ci);
    Graph_Graph01->GetXaxis()->SetTitle("p_{T}^{#mu#mu#mu} [GeV]");
    Graph_Graph01->GetXaxis()->SetLabelFont(42);
    Graph_Graph01->GetXaxis()->SetTitleOffset(1.2);
    Graph_Graph01->GetXaxis()->SetTitleFont(42);
    Graph_Graph01->GetYaxis()->SetTitle("R_{PbPb}(B_{c})");
    Graph_Graph01->GetYaxis()->SetLabelFont(42);
    Graph_Graph01->GetYaxis()->SetTitleOffset(1.3);
    Graph_Graph01->GetYaxis()->SetTitleFont(42);
    Graph_Graph01->GetZaxis()->SetLabelFont(42);
    Graph_Graph01->GetZaxis()->SetTitleOffset(1);
    Graph_Graph01->GetZaxis()->SetTitleFont(42);
    tgme->SetHistogram(Graph_Graph01);
   
    tgme->Draw("p");
   
    tgme = new TGraphMultiErrors(1, 4);
    tgme->SetName("Graph1");
    tgme->SetTitle("Graph");

    ci = TColor::GetColor("#8cd1e0");
    tgme->SetFillColor(ci);
    tgme->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->SetLineColor(ci);

    ci = TColor::GetColor("#003366");
    tgme->SetMarkerColor(ci);
    tgme->SetMarkerStyle(20);
    tgme->SetMarkerSize(1.4);

    ci = TColor::GetColor("#8cd1e0");
    tgme->GetAttFill(0)->SetFillColor(ci);
    //tgme->GetAttFill(0)->SetFillStyle(3001);
    tgme->GetAttFill(0)->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(0)->SetLineColor(ci);
    tgme->GetAttLine(0)->SetLineWidth(3);

    ci = TColor::GetColor("#8cd1e0");
    tgme->GetAttFill(1)->SetFillColor(ci);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(1)->SetLineColor(ci);
    tgme->GetAttLine(1)->SetLineWidth(3);
    tgme->GetAttFill(2)->SetFillColor(19);
    tgme->GetAttFill(2)->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(2)->SetLineColor(ci);
    tgme->GetAttLine(2)->SetLineStyle(2);
    tgme->GetAttLine(2)->SetLineWidth(3);
    tgme->GetAttFill(3)->SetFillColor(19);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(3)->SetLineColor(ci);
    tgme->GetAttLine(3)->SetLineWidth(3);
    tgme->SetPoint(0, 8.059, 1.665756);
    tgme->SetPointEX(0, 2.059, 2.941);
    tgme->SetPointEY(0, 0, 0.5781341, 0.9120825);
    tgme->SetPointEY(0, 1, 0.5259643, 0.851301);
    tgme->SetPointEY(0, 2, 0.5781341, 0.9120825);
    tgme->SetPointEY(0, 3, 0.5259643, 0.851301);

    //remove partial error here
    tgme->SetPointEY(0, 1, 0., 0.); 
    tgme->SetPointEY(0, 3, 0., 0.); 
   
    TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","Graph",100,5.5,11.5);
    Graph_Graph12->SetMinimum(0.8767821);
    Graph_Graph12->SetMaximum(2.659318);
    Graph_Graph12->SetDirectory(0);
    Graph_Graph12->SetStats(0);

    ci = TColor::GetColor("#000099");
    Graph_Graph12->SetLineColor(ci);
    Graph_Graph12->GetXaxis()->SetLabelFont(42);
    Graph_Graph12->GetXaxis()->SetTitleOffset(1);
    Graph_Graph12->GetXaxis()->SetTitleFont(42);
    Graph_Graph12->GetYaxis()->SetLabelFont(42);
    Graph_Graph12->GetYaxis()->SetTitleFont(42);
    Graph_Graph12->GetZaxis()->SetLabelFont(42);
    Graph_Graph12->GetZaxis()->SetTitleOffset(1);
    Graph_Graph12->GetZaxis()->SetTitleFont(42);
    tgme->SetHistogram(Graph_Graph12);
   
    tgme->Draw("px s; 2; 2; p s=0; p s=0");
   
    tgme->SetLineWidth(0);











    tgme = new TGraphMultiErrors(1, 4);
    tgme->SetName("Graph2");
    tgme->SetTitle("Graph");

    ci = TColor::GetColor("#8cd1e0");
    tgme->SetFillColor(ci);
    tgme->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->SetLineColor(ci);

    ci = TColor::GetColor("#003366");
    tgme->SetMarkerColor(ci);
    tgme->SetMarkerStyle(89);
    tgme->SetMarkerSize(1.4);

    ci = TColor::GetColor("#8cd1e0");
    tgme->GetAttFill(0)->SetFillColor(ci);
    //tgme->GetAttFill(0)->SetFillStyle(3001);
    tgme->GetAttFill(0)->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(0)->SetLineColor(ci);
    tgme->GetAttLine(0)->SetLineWidth(3);

    ci = TColor::GetColor("#8cd1e0");
    tgme->GetAttFill(1)->SetFillColor(ci);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(1)->SetLineColor(ci);
    tgme->GetAttLine(1)->SetLineWidth(3);
    tgme->GetAttFill(2)->SetFillColor(19);
    tgme->GetAttFill(2)->SetFillColorAlpha(ci,transp);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(2)->SetLineColor(ci);
    tgme->GetAttLine(2)->SetLineStyle(2);
    tgme->GetAttLine(2)->SetLineWidth(3);
    tgme->GetAttFill(3)->SetFillColor(19);

    ci = TColor::GetColor("#3399cc");
    tgme->GetAttLine(3)->SetLineColor(ci);
    tgme->GetAttLine(3)->SetLineWidth(3);
    tgme->SetPoint(0, 17.165, 0.6961609);
    tgme->SetPointEX(0, 6.165001, 17.835);
    tgme->SetPointEY(0, 0, 0.1213511, 0.1403348);
    tgme->SetPointEY(0, 1, 0.1104005, 0.1309829);
    tgme->SetPointEY(0, 2, 0.1213511, 0.1403348);
    tgme->SetPointEY(0, 3, 0.1104005, 0.1309829);

    //remove partial error here
    tgme->SetPointEY(0, 1, 0., 0.); 
    tgme->SetPointEY(0, 3, 0., 0.); 
   
    TH1F *Graph_Graph23 = new TH1F("Graph_Graph23","Graph",100,8.6,37.4);
    Graph_Graph23->SetMinimum(0.5405611);
    Graph_Graph23->SetMaximum(0.8545857);
    Graph_Graph23->SetDirectory(0);
    Graph_Graph23->SetStats(0);

    ci = TColor::GetColor("#000099");
    Graph_Graph23->SetLineColor(ci);
    Graph_Graph23->GetXaxis()->SetLabelFont(42);
    Graph_Graph23->GetXaxis()->SetTitleOffset(1);
    Graph_Graph23->GetXaxis()->SetTitleFont(42);
    Graph_Graph23->GetYaxis()->SetLabelFont(42);
    Graph_Graph23->GetYaxis()->SetTitleFont(42);
    Graph_Graph23->GetZaxis()->SetLabelFont(42);
    Graph_Graph23->GetZaxis()->SetTitleOffset(1);
    Graph_Graph23->GetZaxis()->SetTitleFont(42);
    tgme->SetHistogram(Graph_Graph23);
   
    tgme->Draw("px s ; 2; 2; p s=0; p s=0");

    tgme->SetLineWidth(0);
   
    // TLegend *leg = new TLegend(0.55,0.62,0.94,0.78,NULL,"brNDC");
    // leg->SetBorderSize(0);
    // leg->SetTextSize(0.039);
    // leg->SetLineColor(1);
    // leg->SetLineStyle(1);
    // leg->SetLineWidth(1);
    // leg->SetFillColor(0);
    // leg->SetFillStyle(1001);
    // TLegendEntry *entry=leg->AddEntry("Graph1","1.3<|y^{#mu#mu#mu}|<2.3","lpf");

    // ci = TColor::GetColor("#8cd1e0");
    // entry->SetFillColor(ci);
    // entry->SetFillStyle(1001);

    // ci = TColor::GetColor("#3399cc");
    // entry->SetLineColor(ci);
    // entry->SetLineStyle(1);
    // entry->SetLineWidth(3);

    // ci = TColor::GetColor("#003366");
    // entry->SetMarkerColor(ci);
    // entry->SetMarkerStyle(20);
    // entry->SetMarkerSize(2);
    // entry->SetTextFont(42);
    // entry=leg->AddEntry("Graph2","0<|y^{#mu#mu#mu}|<2.3","lpf");

    // ci = TColor::GetColor("#8cd1e0");
    // entry->SetFillColor(ci);
    // entry->SetFillStyle(1001);

    // ci = TColor::GetColor("#3399cc");
    // entry->SetLineColor(ci);
    // entry->SetLineStyle(1);
    // entry->SetLineWidth(3);

    // ci = TColor::GetColor("#003366");
    // entry->SetMarkerColor(ci);
    // entry->SetMarkerStyle(89);
    // entry->SetMarkerSize(2);
    // entry->SetTextFont(42);
    // leg->Draw();
    // TLatex *   tex = new TLatex(0.35,0.93,"pp (302 pb^{-1}) + PbPb (1.61 nb^{-1}), 5.02 TeV");
    // tex->SetNDC();
    // tex->SetTextFont(42);
    // tex->SetTextSize(0.035);
    // tex->SetLineWidth(2);
    // tex->Draw();
    // tex = new TLatex(0.47,0.87,"#font[61]{B_{c}^{#pm} #rightarrow (J/#psi #rightarrow #mu^{-} #mu^{+}) #mu^{#pm} #nu_{#mu} }");
    // tex->SetNDC();
    // tex->SetTextFont(42);
    // tex->SetTextSize(0.043);
    // tex->SetLineWidth(2);
    // tex->Draw();
    // tex = new TLatex(327.9,0.68,"Centrality 0-20%");
    // tex->SetTextAlign(21);
    // tex->SetTextFont(42);
    // tex->SetTextSize(0.034);
    // tex->SetLineWidth(2);
    // tex->Draw();
    // tex = new TLatex(73.5,0.68,"Centrality 20-90%");
    // tex->SetTextAlign(21);
    // tex->SetTextFont(42);
    // tex->SetTextSize(0.034);
    // tex->SetLineWidth(2);
    // tex->Draw();
  }

}
