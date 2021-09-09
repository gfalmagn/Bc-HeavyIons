#include <TColor.h>
#include <TBox.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>

#ifndef _CHAR_LEG_GR_
#define _CHAR_LEG_GR_
TGraphAsymmErrors* gTrackPt_leg;
#endif
TH1D* RAA_0_100584__1__1__1;

void RAA_0_100_marker(bool noXerr=false)
{
  RAA_0_100584__1__1__1->Draw(noXerr?"EX0same":"same");  
}
void RAA_0_100()
{
  //=========Macro generated from canvas: RAA/RAA
  //=========  (Sat Oct 22 22:09:39 2016) by ROOT version5.34/34
  /*
    TCanvas *RAA = new TCanvas("RAA", "RAA",8,31,800,676);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    RAA->Range(-0.9791072,-0.3,2.895341,2.2);
    RAA->SetFillColor(0);
    RAA->SetBorderMode(0);
    RAA->SetBorderSize(2);
    RAA->SetLogx();
    RAA->SetTickx(1);
    RAA->SetTicky(1);
    RAA->SetLeftMargin(0.15);
    RAA->SetRightMargin(0.04);
    RAA->SetTopMargin(0.08);
    RAA->SetBottomMargin(0.12);
    RAA->SetFrameFillStyle(0);
    RAA->SetFrameBorderMode(0);
    RAA->SetFrameFillStyle(0);
    RAA->SetFrameBorderMode(0);
  */

  Int_t ci = TColor::GetColor("#ffcc00");

  gTrackPt_leg = new TGraphAsymmErrors(1);
  gTrackPt_leg->SetMarkerStyle(20);
  gTrackPt_leg->SetMarkerSize(1.0); // 0.8
  gTrackPt_leg->SetMarkerColor(kRed);
  gTrackPt_leg->SetLineColor(0);
  gTrackPt_leg->SetFillColor(ci);

  Double_t xAxis1[38] = {0.4, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24, 28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400}; 
   
  RAA_0_100584__1__1__1 = new TH1D("RAA_0_100584__1__1__1","",37, xAxis1);
  RAA_0_100584__1__1__1->SetBinContent(1,-1);
  RAA_0_100584__1__1__1->SetBinContent(2,-1);
  /*
    RAA_0_100584__1__1__1->SetBinContent(3,0.3605047);
    RAA_0_100584__1__1__1->SetBinContent(4,0.3850582);
    RAA_0_100584__1__1__1->SetBinContent(5,0.4188398);
  */
  RAA_0_100584__1__1__1->SetBinContent(3,-1);
  RAA_0_100584__1__1__1->SetBinContent(4,-1);
  RAA_0_100584__1__1__1->SetBinContent(5,-1);
  RAA_0_100584__1__1__1->SetBinContent(6,0.4411318);
  RAA_0_100584__1__1__1->SetBinContent(7,0.4564881);
  RAA_0_100584__1__1__1->SetBinContent(8,0.4798668);
  RAA_0_100584__1__1__1->SetBinContent(9,0.4985427);
  RAA_0_100584__1__1__1->SetBinContent(10,0.516094);
  RAA_0_100584__1__1__1->SetBinContent(11,0.5160761);
  RAA_0_100584__1__1__1->SetBinContent(12,0.5157593);
  RAA_0_100584__1__1__1->SetBinContent(13,0.5090351);
  RAA_0_100584__1__1__1->SetBinContent(14,0.4718224);
  RAA_0_100584__1__1__1->SetBinContent(15,0.381288);
  RAA_0_100584__1__1__1->SetBinContent(16,0.3016595);
  RAA_0_100584__1__1__1->SetBinContent(17,0.2568378);
  RAA_0_100584__1__1__1->SetBinContent(18,0.2362394);
  RAA_0_100584__1__1__1->SetBinContent(19,0.2313929);
  RAA_0_100584__1__1__1->SetBinContent(20,0.2397292);
  RAA_0_100584__1__1__1->SetBinContent(21,0.2795414);
  RAA_0_100584__1__1__1->SetBinContent(22,0.3109757);
  RAA_0_100584__1__1__1->SetBinContent(23,0.359047);
  RAA_0_100584__1__1__1->SetBinContent(24,0.4103232);
  RAA_0_100584__1__1__1->SetBinContent(25,0.4535739);
  RAA_0_100584__1__1__1->SetBinContent(26,0.4871136);
  RAA_0_100584__1__1__1->SetBinContent(27,0.5390386);
  RAA_0_100584__1__1__1->SetBinContent(28,0.6064001);
  RAA_0_100584__1__1__1->SetBinContent(29,0.6583816);
  RAA_0_100584__1__1__1->SetBinContent(30,0.6914895);
  RAA_0_100584__1__1__1->SetBinContent(31,0.7510186);
  RAA_0_100584__1__1__1->SetBinContent(32,0.8027461);
  RAA_0_100584__1__1__1->SetBinContent(33,0.8064291);
  RAA_0_100584__1__1__1->SetBinContent(34,0.8144861);
  RAA_0_100584__1__1__1->SetBinContent(35,0.8835603);
  RAA_0_100584__1__1__1->SetBinContent(36,0.9077175);
  RAA_0_100584__1__1__1->SetBinContent(37,0.8767255);
  RAA_0_100584__1__1__1->SetBinError(1,1.279903e-005);
  RAA_0_100584__1__1__1->SetBinError(2,1.387861e-005);
  RAA_0_100584__1__1__1->SetBinError(3,1.439234e-005);
  RAA_0_100584__1__1__1->SetBinError(4,1.673993e-005);
  RAA_0_100584__1__1__1->SetBinError(5,1.996438e-005);
  RAA_0_100584__1__1__1->SetBinError(6,2.296277e-005);
  RAA_0_100584__1__1__1->SetBinError(7,2.652961e-005);
  RAA_0_100584__1__1__1->SetBinError(8,2.303589e-005);
  RAA_0_100584__1__1__1->SetBinError(9,2.956636e-005);
  RAA_0_100584__1__1__1->SetBinError(10,3.708897e-005);
  RAA_0_100584__1__1__1->SetBinError(11,4.504344e-005);
  RAA_0_100584__1__1__1->SetBinError(12,5.385318e-005);
  RAA_0_100584__1__1__1->SetBinError(13,6.332594e-005);
  RAA_0_100584__1__1__1->SetBinError(14,4.278537e-005);
  RAA_0_100584__1__1__1->SetBinError(15,6.266848e-005);
  RAA_0_100584__1__1__1->SetBinError(16,8.299302e-005);
  RAA_0_100584__1__1__1->SetBinError(17,0.0001099248);
  RAA_0_100584__1__1__1->SetBinError(18,0.0001477243);
  RAA_0_100584__1__1__1->SetBinError(19,0.0002002684);
  RAA_0_100584__1__1__1->SetBinError(20,0.0001962782);
  RAA_0_100584__1__1__1->SetBinError(21,0.00043599);
  RAA_0_100584__1__1__1->SetBinError(22,0.0007955915);
  RAA_0_100584__1__1__1->SetBinError(23,0.001121086);
  RAA_0_100584__1__1__1->SetBinError(24,0.002349183);
  RAA_0_100584__1__1__1->SetBinError(25,0.003934258);
  RAA_0_100584__1__1__1->SetBinError(26,0.005431548);
  RAA_0_100584__1__1__1->SetBinError(27,0.007593713);
  RAA_0_100584__1__1__1->SetBinError(28,0.008532134);
  RAA_0_100584__1__1__1->SetBinError(29,0.006423681);
  RAA_0_100584__1__1__1->SetBinError(30,0.009668077);
  RAA_0_100584__1__1__1->SetBinError(31,0.009013466);
  RAA_0_100584__1__1__1->SetBinError(32,0.01230905);
  RAA_0_100584__1__1__1->SetBinError(33,0.01985833);
  RAA_0_100584__1__1__1->SetBinError(34,0.02497965);
  RAA_0_100584__1__1__1->SetBinError(35,0.04145952);
  RAA_0_100584__1__1__1->SetBinError(36,0.0459072);
  RAA_0_100584__1__1__1->SetBinError(37,0.152905);
  RAA_0_100584__1__1__1->SetMinimum(0);
  RAA_0_100584__1__1__1->SetMaximum(2);
  RAA_0_100584__1__1__1->SetEntries(11078.9);
  RAA_0_100584__1__1__1->SetDirectory(0);

  RAA_0_100584__1__1__1->SetLineWidth(2);
  RAA_0_100584__1__1__1->SetLineColor(kRed);
  RAA_0_100584__1__1__1->SetMarkerStyle(20);
  RAA_0_100584__1__1__1->SetMarkerSize(1.0); // 0.8
  RAA_0_100584__1__1__1->SetMarkerColor(kRed);

  /*
    RAA_0_100584__1__1__1->SetFillColor(ci);
    ci = TColor::GetColor("#000099");
    RAA_0_100584__1__1__1->SetLineColor(ci);
    RAA_0_100584__1__1__1->SetLineWidth(2);
    RAA_0_100584__1__1__1->SetMarkerStyle(20);
    RAA_0_100584__1__1__1->SetMarkerSize(1.3);
    RAA_0_100584__1__1__1->GetXaxis()->SetTitle("p_{T} (GeV)");
    RAA_0_100584__1__1__1->GetXaxis()->SetRange(1,38);
    RAA_0_100584__1__1__1->GetXaxis()->CenterTitle(true);
    RAA_0_100584__1__1__1->GetXaxis()->SetLabelFont(42);
    RAA_0_100584__1__1__1->GetXaxis()->SetLabelOffset(-0.002);
    RAA_0_100584__1__1__1->GetXaxis()->SetLabelSize(0.045);
    RAA_0_100584__1__1__1->GetXaxis()->SetTitleSize(0.05);
    RAA_0_100584__1__1__1->GetXaxis()->SetTitleFont(42);
    RAA_0_100584__1__1__1->GetYaxis()->SetTitle("Nuclear modification factor");
    RAA_0_100584__1__1__1->GetYaxis()->CenterTitle(true);
    RAA_0_100584__1__1__1->GetYaxis()->SetNdivisions(508);
    RAA_0_100584__1__1__1->GetYaxis()->SetLabelFont(42);
    RAA_0_100584__1__1__1->GetYaxis()->SetLabelSize(0.045);
    RAA_0_100584__1__1__1->GetYaxis()->SetTitleSize(0.05);
    RAA_0_100584__1__1__1->GetYaxis()->SetTitleFont(42);
    RAA_0_100584__1__1__1->GetZaxis()->SetLabelFont(42);
    RAA_0_100584__1__1__1->GetZaxis()->SetLabelSize(0.05);
    RAA_0_100584__1__1__1->GetZaxis()->SetTitleSize(0.09);
    RAA_0_100584__1__1__1->GetZaxis()->SetTitleFont(42);
    RAA_0_100584__1__1__1->Draw("same");
  */

  TBox *box = new TBox(0.7,0.3322697,0.8,0.3887397);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  //box->Draw();
  box = new TBox(0.8,0.3547214,0.9,0.4153951);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  //box->Draw();
  box = new TBox(0.9,0.3856854,1,0.4519942);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  //box->Draw();
  box = new TBox(1,0.4035821,1.1,0.4786815);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(1.1,0.4178419,1.2,0.4951344);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(1.2,0.4336979,1.4,0.5260356);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(1.4,0.4471445,1.6,0.5499409);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(1.6,0.454559,1.8,0.5776292);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(1.8,0.4545544,2,0.5775976);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(2,0.4458492,2.2,0.5856694);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(2.2,0.4400381,2.4,0.5780321);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(2.4,0.4004238,3.2,0.5432211);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(3.2,0.3194642,4,0.443112);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(4,0.2555583,4.8,0.3477606);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(4.8,0.2236467,5.6,0.290029);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(5.6,0.2091373,6.4,0.2633415);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(6.4,0.2084873,7.2,0.2542985);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(7.2,0.218131,9.6,0.2613273);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(9.6,0.254349,12,0.3047338);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(12,0.279773,14.4,0.3421784);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(14.4,0.3230158,19.2,0.395078);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(19.2,0.3691498,24,0.4514966);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(24,0.4081148,28.8,0.4990331);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(28.8,0.4382444,35.2,0.5359827);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(35.2,0.4849575,41.6,0.5931196);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(41.6,0.54671,48,0.6660902);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(48,0.5934783,60.8,0.7232848);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(60.8,0.6229389,73.6,0.7600401);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(73.6,0.6765564,86.4,0.8254808);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(86.4,0.7227269,103.6,0.8827653);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(103.6,0.7255117,120.8,0.8873466);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(120.8,0.732652,140,0.8963203);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(140,0.7948418,165,0.9722788);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(165,0.8158651,250,0.9995696);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();
  box = new TBox(250,0.7867229,400,0.9667282);

  ci = TColor::GetColor("#ffcc00");
  box->SetFillColor(ci);
  box->Draw();

  // RAA_0_100584__1__1__1->Draw("same");

  /*
    TLatex *   tex = new TLatex(0.96,0.936,"27.4 pb^{-1} (pp) + 35 nb^{-1} (pPb) + 404 #mub^{-1} (PbPb) 5.02 TeV");
    tex->SetNDC();
    tex->SetTextAlign(31);
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.82645,0.892,"CMS");
    tex->SetNDC();
    tex->SetTextAlign(13);
    tex->SetTextFont(61);
    tex->SetTextSize(0.06);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.5,0.2,"0-100%");
    tex->SetTextAlign(13);
    tex->SetTextFont(61);
    tex->SetLineWidth(2);
    tex->Draw();
    box = new TBox(0.5,0.5741027,0.6,0.6452332);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(0.6,0.5998874,0.7,0.6741247);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(0.7,0.6431771,0.8,0.7226311);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(0.8,0.6575168,0.9,0.7429621);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(0.9,0.7007984,1,0.7918376);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(1,0.7199725,1.1,0.8150175);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(1.1,0.7479496,1.2,0.8466735);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(1.2,0.7628813,1.4,0.881854);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(1.4,0.8133276,1.6,0.9401723);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(1.6,0.8240483,1.8,0.9784846);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(1.8,0.8612759,2,1.022693);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(2,0.863603,2.2,1.053784);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(2.2,0.8809571,2.4,1.074963);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(2.4,0.905089,3.2,1.127269);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(3.2,0.9327053,4,1.185313);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(4,0.951561,4.8,1.201849);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(4.8,0.9722057,5.6,1.202725);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(5.6,0.9779246,6.4,1.181506);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(6.4,0.9821993,7.2,1.173107);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(7.2,1.000774,9.6,1.176233);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(9.6,0.9807633,12,1.149596);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(12,0.9964006,14.4,1.168115);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(14.4,1.001861,19.2,1.197871);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(19.2,1.041055,24,1.244671);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(24,1.053542,28.8,1.259239);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(28.8,1.03686,35.2,1.239524);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(35.2,1.045464,41.6,1.249769);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(41.6,1.093577,48,1.307943);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(48,1.095519,60.8,1.31062);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(60.8,1.092843,73.6,1.308307);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(73.6,1.090087,86.4,1.305008);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(86.4,1.089054,103.6,1.30414);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(103.6,1.076576,120.8,1.29108);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(120.8,0,140,0);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(140,0,165,0);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(165,0,250,0);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(250,0,400,0);

    ci = TColor::GetColor("#00944c");
    box->SetFillColor(ci);

    ci = TColor::GetColor("#00944c");
    box->SetLineColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(0.55,0.9553,0.65,1.0447);

    ci = TColor::GetColor("#8fbc8f");
    box->SetFillColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(0.65,0.911,0.77,1.089);

    ci = TColor::GetColor("#9999ff");
    box->SetFillColor(ci);
    box->SetLineWidth(0);
    box->Draw();
    box = new TBox(0.77,0.977,0.9,1.023);
    box->SetFillColor(17);
    box->SetLineWidth(0);
    box->Draw();
    TLine *line = new TLine(0.4,1,550,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
    Double_t xAxis2[38] = {0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24, 28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400}; 
   
    TH1D *pPb_data = new TH1D("pPb_data","",37, xAxis2);
    pPb_data->SetBinContent(1,0.609864);
    pPb_data->SetBinContent(2,0.6372112);
    pPb_data->SetBinContent(3,0.6831244);
    pPb_data->SetBinContent(4,0.7004547);
    pPb_data->SetBinContent(5,0.7465476);
    pPb_data->SetBinContent(6,0.767272);
    pPb_data->SetBinContent(7,0.7970799);
    pPb_data->SetBinContent(8,0.8221631);
    pPb_data->SetBinContent(9,0.8765319);
    pPb_data->SetBinContent(10,0.9010772);
    pPb_data->SetBinContent(11,0.9417868);
    pPb_data->SetBinContent(12,0.9588383);
    pPb_data->SetBinContent(13,0.978108);
    pPb_data->SetBinContent(14,1.016318);
    pPb_data->SetBinContent(15,1.059064);
    pPb_data->SetBinContent(16,1.076763);
    pPb_data->SetBinContent(17,1.08753);
    pPb_data->SetBinContent(18,1.080052);
    pPb_data->SetBinContent(19,1.077315);
    pPb_data->SetBinContent(20,1.083601);
    pPb_data->SetBinContent(21,1.059035);
    pPb_data->SetBinContent(22,1.076021);
    pPb_data->SetBinContent(23,1.093905);
    pPb_data->SetBinContent(24,1.136667);
    pPb_data->SetBinContent(25,1.150111);
    pPb_data->SetBinContent(26,1.132018);
    pPb_data->SetBinContent(27,1.14139);
    pPb_data->SetBinContent(28,1.194263);
    pPb_data->SetBinContent(29,1.196569);
    pPb_data->SetBinContent(30,1.194113);
    pPb_data->SetBinContent(31,1.191102);
    pPb_data->SetBinContent(32,1.190166);
    pPb_data->SetBinContent(33,1.177516);
    pPb_data->SetBinError(1,0.002782783);
    pPb_data->SetBinError(2,0.002901004);
    pPb_data->SetBinError(3,0.003108645);
    pPb_data->SetBinError(4,0.003177324);
    pPb_data->SetBinError(5,0.003372345);
    pPb_data->SetBinError(6,0.003479202);
    pPb_data->SetBinError(7,0.003563534);
    pPb_data->SetBinError(8,0.002745707);
    pPb_data->SetBinError(9,0.002849567);
    pPb_data->SetBinError(10,0.002942639);
    pPb_data->SetBinError(11,0.003001263);
    pPb_data->SetBinError(12,0.003069924);
    pPb_data->SetBinError(13,0.003257273);
    pPb_data->SetBinError(14,0.001788977);
    pPb_data->SetBinError(15,0.001954066);
    pPb_data->SetBinError(16,0.001841783);
    pPb_data->SetBinError(17,0.002316191);
    pPb_data->SetBinError(18,0.002484718);
    pPb_data->SetBinError(19,0.003397278);
    pPb_data->SetBinError(20,0.0122483);
    pPb_data->SetBinError(21,0.01373172);
    pPb_data->SetBinError(22,0.01781856);
    pPb_data->SetBinError(23,0.03796211);
    pPb_data->SetBinError(24,0.008864955);
    pPb_data->SetBinError(25,0.01105857);
    pPb_data->SetBinError(26,0.01207979);
    pPb_data->SetBinError(27,0.01533868);
    pPb_data->SetBinError(28,0.01539577);
    pPb_data->SetBinError(29,0.01411351);
    pPb_data->SetBinError(30,0.02025792);
    pPb_data->SetBinError(31,0.02360509);
    pPb_data->SetBinError(32,0.03653756);
    pPb_data->SetBinError(33,0.05304864);
    pPb_data->SetMinimum(0.5);
    pPb_data->SetMaximum(1.5);
    pPb_data->SetEntries(130699.2);

    ci = TColor::GetColor("#00944c");
    pPb_data->SetFillColor(ci);
    pPb_data->SetLineStyle(0);
    pPb_data->SetLineWidth(2);
    pPb_data->SetMarkerStyle(21);
    pPb_data->GetXaxis()->SetTitle("pt");
    pPb_data->GetXaxis()->SetLabelFont(42);
    pPb_data->GetXaxis()->SetLabelOffset(0.007);
    pPb_data->GetXaxis()->SetLabelSize(0.05);
    pPb_data->GetXaxis()->SetTitleSize(0.06);
    pPb_data->GetXaxis()->SetTitleOffset(0.9);
    pPb_data->GetXaxis()->SetTitleFont(42);
    pPb_data->GetYaxis()->SetTitle("RpPb");
    pPb_data->GetYaxis()->SetLabelFont(42);
    pPb_data->GetYaxis()->SetLabelOffset(0.007);
    pPb_data->GetYaxis()->SetLabelSize(0.05);
    pPb_data->GetYaxis()->SetTitleSize(0.06);
    pPb_data->GetYaxis()->SetTitleOffset(1.25);
    pPb_data->GetYaxis()->SetTitleFont(42);
    pPb_data->GetZaxis()->SetLabelFont(42);
    pPb_data->GetZaxis()->SetLabelOffset(0.007);
    pPb_data->GetZaxis()->SetLabelSize(0.05);
    pPb_data->GetZaxis()->SetTitleSize(0.06);
    pPb_data->GetZaxis()->SetTitleFont(42);
    //pPb_data->Draw("same");
    Double_t xAxis3[38] = {0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24, 28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400}; 
   
    TH1D *RAA_0_100584__2__2__2 = new TH1D("RAA_0_100584__2__2__2","",37, xAxis3);
    RAA_0_100584__2__2__2->SetBinContent(1,-1);
    RAA_0_100584__2__2__2->SetBinContent(2,-1);
    RAA_0_100584__2__2__2->SetBinContent(3,0.3605047);
    RAA_0_100584__2__2__2->SetBinContent(4,0.3850582);
    RAA_0_100584__2__2__2->SetBinContent(5,0.4188398);
    RAA_0_100584__2__2__2->SetBinContent(6,0.4411318);
    RAA_0_100584__2__2__2->SetBinContent(7,0.4564881);
    RAA_0_100584__2__2__2->SetBinContent(8,0.4798668);
    RAA_0_100584__2__2__2->SetBinContent(9,0.4985427);
    RAA_0_100584__2__2__2->SetBinContent(10,0.516094);
    RAA_0_100584__2__2__2->SetBinContent(11,0.5160761);
    RAA_0_100584__2__2__2->SetBinContent(12,0.5157593);
    RAA_0_100584__2__2__2->SetBinContent(13,0.5090351);
    RAA_0_100584__2__2__2->SetBinContent(14,0.4718224);
    RAA_0_100584__2__2__2->SetBinContent(15,0.381288);
    RAA_0_100584__2__2__2->SetBinContent(16,0.3016595);
    RAA_0_100584__2__2__2->SetBinContent(17,0.2568378);
    RAA_0_100584__2__2__2->SetBinContent(18,0.2362394);
    RAA_0_100584__2__2__2->SetBinContent(19,0.2313929);
    RAA_0_100584__2__2__2->SetBinContent(20,0.2397292);
    RAA_0_100584__2__2__2->SetBinContent(21,0.2795414);
    RAA_0_100584__2__2__2->SetBinContent(22,0.3109757);
    RAA_0_100584__2__2__2->SetBinContent(23,0.359047);
    RAA_0_100584__2__2__2->SetBinContent(24,0.4103232);
    RAA_0_100584__2__2__2->SetBinContent(25,0.4535739);
    RAA_0_100584__2__2__2->SetBinContent(26,0.4871136);
    RAA_0_100584__2__2__2->SetBinContent(27,0.5390386);
    RAA_0_100584__2__2__2->SetBinContent(28,0.6064001);
    RAA_0_100584__2__2__2->SetBinContent(29,0.6583816);
    RAA_0_100584__2__2__2->SetBinContent(30,0.6914895);
    RAA_0_100584__2__2__2->SetBinContent(31,0.7510186);
    RAA_0_100584__2__2__2->SetBinContent(32,0.8027461);
    RAA_0_100584__2__2__2->SetBinContent(33,0.8064291);
    RAA_0_100584__2__2__2->SetBinContent(34,0.8144861);
    RAA_0_100584__2__2__2->SetBinContent(35,0.8835603);
    RAA_0_100584__2__2__2->SetBinContent(36,0.9077175);
    RAA_0_100584__2__2__2->SetBinContent(37,0.8767255);
    RAA_0_100584__2__2__2->SetBinError(1,1.279903e-005);
    RAA_0_100584__2__2__2->SetBinError(2,1.387861e-005);
    RAA_0_100584__2__2__2->SetBinError(3,1.439234e-005);
    RAA_0_100584__2__2__2->SetBinError(4,1.673993e-005);
    RAA_0_100584__2__2__2->SetBinError(5,1.996438e-005);
    RAA_0_100584__2__2__2->SetBinError(6,2.296277e-005);
    RAA_0_100584__2__2__2->SetBinError(7,2.652961e-005);
    RAA_0_100584__2__2__2->SetBinError(8,2.303589e-005);
    RAA_0_100584__2__2__2->SetBinError(9,2.956636e-005);
    RAA_0_100584__2__2__2->SetBinError(10,3.708897e-005);
    RAA_0_100584__2__2__2->SetBinError(11,4.504344e-005);
    RAA_0_100584__2__2__2->SetBinError(12,5.385318e-005);
    RAA_0_100584__2__2__2->SetBinError(13,6.332594e-005);
    RAA_0_100584__2__2__2->SetBinError(14,4.278537e-005);
    RAA_0_100584__2__2__2->SetBinError(15,6.266848e-005);
    RAA_0_100584__2__2__2->SetBinError(16,8.299302e-005);
    RAA_0_100584__2__2__2->SetBinError(17,0.0001099248);
    RAA_0_100584__2__2__2->SetBinError(18,0.0001477243);
    RAA_0_100584__2__2__2->SetBinError(19,0.0002002684);
    RAA_0_100584__2__2__2->SetBinError(20,0.0001962782);
    RAA_0_100584__2__2__2->SetBinError(21,0.00043599);
    RAA_0_100584__2__2__2->SetBinError(22,0.0007955915);
    RAA_0_100584__2__2__2->SetBinError(23,0.001121086);
    RAA_0_100584__2__2__2->SetBinError(24,0.002349183);
    RAA_0_100584__2__2__2->SetBinError(25,0.003934258);
    RAA_0_100584__2__2__2->SetBinError(26,0.005431548);
    RAA_0_100584__2__2__2->SetBinError(27,0.007593713);
    RAA_0_100584__2__2__2->SetBinError(28,0.008532134);
    RAA_0_100584__2__2__2->SetBinError(29,0.006423681);
    RAA_0_100584__2__2__2->SetBinError(30,0.009668077);
    RAA_0_100584__2__2__2->SetBinError(31,0.009013466);
    RAA_0_100584__2__2__2->SetBinError(32,0.01230905);
    RAA_0_100584__2__2__2->SetBinError(33,0.01985833);
    RAA_0_100584__2__2__2->SetBinError(34,0.02497965);
    RAA_0_100584__2__2__2->SetBinError(35,0.04145952);
    RAA_0_100584__2__2__2->SetBinError(36,0.0459072);
    RAA_0_100584__2__2__2->SetBinError(37,0.152905);
    RAA_0_100584__2__2__2->SetMinimum(0);
    RAA_0_100584__2__2__2->SetMaximum(1.8);
    RAA_0_100584__2__2__2->SetEntries(11078.9);
    RAA_0_100584__2__2__2->SetDirectory(0);

    ci = TColor::GetColor("#ffcc00");
    RAA_0_100584__2__2__2->SetFillColor(ci);
    RAA_0_100584__2__2__2->SetLineWidth(2);
    RAA_0_100584__2__2__2->SetMarkerStyle(20);
    RAA_0_100584__2__2__2->SetMarkerSize(1.3);
    RAA_0_100584__2__2__2->GetXaxis()->SetTitle("p_{T} (GeV)");
    RAA_0_100584__2__2__2->GetXaxis()->SetRange(1,37);
    RAA_0_100584__2__2__2->GetXaxis()->CenterTitle(true);
    RAA_0_100584__2__2__2->GetXaxis()->SetLabelFont(42);
    RAA_0_100584__2__2__2->GetXaxis()->SetLabelOffset(-0.002);
    RAA_0_100584__2__2__2->GetXaxis()->SetLabelSize(0.045);
    RAA_0_100584__2__2__2->GetXaxis()->SetTitleSize(0.05);
    RAA_0_100584__2__2__2->GetXaxis()->SetTitleFont(42);
    RAA_0_100584__2__2__2->GetYaxis()->SetTitle("R_{AA}, R_{pA}");
    RAA_0_100584__2__2__2->GetYaxis()->CenterTitle(true);
    RAA_0_100584__2__2__2->GetYaxis()->SetNdivisions(508);
    RAA_0_100584__2__2__2->GetYaxis()->SetLabelFont(42);
    RAA_0_100584__2__2__2->GetYaxis()->SetLabelSize(0.045);
    RAA_0_100584__2__2__2->GetYaxis()->SetTitleSize(0.05);
    RAA_0_100584__2__2__2->GetYaxis()->SetTitleFont(42);
    RAA_0_100584__2__2__2->GetZaxis()->SetLabelFont(42);
    RAA_0_100584__2__2__2->GetZaxis()->SetLabelSize(0.05);
    RAA_0_100584__2__2__2->GetZaxis()->SetTitleSize(0.09);
    RAA_0_100584__2__2__2->GetZaxis()->SetTitleFont(42);
    RAA_0_100584__2__2__2->Draw("same");
   
    TLegend *leg = new TLegend(0.1729825,0.7558462,0.5038095,0.8558462,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04607692);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    TLegendEntry *entry=leg->AddEntry("pPb_data","R_{pA}","pf");

    ci = TColor::GetColor("#00944c");
    entry->SetFillColor(ci);
    entry->SetFillStyle(1001);
    entry->SetLineColor(1);
    entry->SetLineWidth(2);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1);
    entry->SetTextFont(62);
    entry=leg->AddEntry("NULL","R_{AA}","pf");

    ci = TColor::GetColor("#ffcc00");
    entry->SetFillColor(ci);
    entry->SetFillStyle(1001);
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1.3);
    entry->SetTextFont(62);
    leg->Draw();
   
    leg = new TLegend(0.4711779,0.72,0.783208,0.8553846,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04307692);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    entry=leg->AddEntry("NULL","T_{pA} uncertainty","f");

    ci = TColor::GetColor("#8fbc8f");
    entry->SetFillColor(ci);
    entry->SetFillStyle(1001);
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1);
    entry->SetTextFont(62);
    entry=leg->AddEntry("NULL","T_{AA} uncertainty","f");

    ci = TColor::GetColor("#9999ff");
    entry->SetFillColor(ci);
    entry->SetFillStyle(1001);
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1);
    entry->SetTextFont(62);
    entry=leg->AddEntry("NULL","pp lumi. uncertainty","f");
    entry->SetFillColor(17);
    entry->SetFillStyle(1001);
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1);
    entry->SetTextFont(62);
    leg->Draw();
    tex = new TLatex(0.18,0.75,"|#eta_{CM}| < 1");
    tex->SetNDC();
    tex->SetTextAlign(13);
    tex->SetTextSize(0.047);
    tex->SetLineWidth(2);
    tex->Draw();
    Double_t xAxis4[39] = {0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24, 28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400, 500}; 
   
    TH1D *RAA_0_100584__3__3__3 = new TH1D("RAA_0_100584__3__3__3","",38, xAxis4);
    RAA_0_100584__3__3__3->SetBinContent(1,-1);
    RAA_0_100584__3__3__3->SetBinContent(2,-1);
    RAA_0_100584__3__3__3->SetBinContent(3,0.3605047);
    RAA_0_100584__3__3__3->SetBinContent(4,0.3850582);
    RAA_0_100584__3__3__3->SetBinContent(5,0.4188398);
    RAA_0_100584__3__3__3->SetBinContent(6,0.4411318);
    RAA_0_100584__3__3__3->SetBinContent(7,0.4564881);
    RAA_0_100584__3__3__3->SetBinContent(8,0.4798668);
    RAA_0_100584__3__3__3->SetBinContent(9,0.4985427);
    RAA_0_100584__3__3__3->SetBinContent(10,0.516094);
    RAA_0_100584__3__3__3->SetBinContent(11,0.5160761);
    RAA_0_100584__3__3__3->SetBinContent(12,0.5157593);
    RAA_0_100584__3__3__3->SetBinContent(13,0.5090351);
    RAA_0_100584__3__3__3->SetBinContent(14,0.4718224);
    RAA_0_100584__3__3__3->SetBinContent(15,0.381288);
    RAA_0_100584__3__3__3->SetBinContent(16,0.3016595);
    RAA_0_100584__3__3__3->SetBinContent(17,0.2568378);
    RAA_0_100584__3__3__3->SetBinContent(18,0.2362394);
    RAA_0_100584__3__3__3->SetBinContent(19,0.2313929);
    RAA_0_100584__3__3__3->SetBinContent(20,0.2397292);
    RAA_0_100584__3__3__3->SetBinContent(21,0.2795414);
    RAA_0_100584__3__3__3->SetBinContent(22,0.3109757);
    RAA_0_100584__3__3__3->SetBinContent(23,0.359047);
    RAA_0_100584__3__3__3->SetBinContent(24,0.4103232);
    RAA_0_100584__3__3__3->SetBinContent(25,0.4535739);
    RAA_0_100584__3__3__3->SetBinContent(26,0.4871136);
    RAA_0_100584__3__3__3->SetBinContent(27,0.5390386);
    RAA_0_100584__3__3__3->SetBinContent(28,0.6064001);
    RAA_0_100584__3__3__3->SetBinContent(29,0.6583816);
    RAA_0_100584__3__3__3->SetBinContent(30,0.6914895);
    RAA_0_100584__3__3__3->SetBinContent(31,0.7510186);
    RAA_0_100584__3__3__3->SetBinContent(32,0.8027461);
    RAA_0_100584__3__3__3->SetBinContent(33,0.8064291);
    RAA_0_100584__3__3__3->SetBinContent(34,0.8144861);
    RAA_0_100584__3__3__3->SetBinContent(35,0.8835603);
    RAA_0_100584__3__3__3->SetBinContent(36,0.9077175);
    RAA_0_100584__3__3__3->SetBinContent(37,0.8767255);
    RAA_0_100584__3__3__3->SetBinError(1,1.279903e-005);
    RAA_0_100584__3__3__3->SetBinError(2,1.387861e-005);
    RAA_0_100584__3__3__3->SetBinError(3,1.439234e-005);
    RAA_0_100584__3__3__3->SetBinError(4,1.673993e-005);
    RAA_0_100584__3__3__3->SetBinError(5,1.996438e-005);
    RAA_0_100584__3__3__3->SetBinError(6,2.296277e-005);
    RAA_0_100584__3__3__3->SetBinError(7,2.652961e-005);
    RAA_0_100584__3__3__3->SetBinError(8,2.303589e-005);
    RAA_0_100584__3__3__3->SetBinError(9,2.956636e-005);
    RAA_0_100584__3__3__3->SetBinError(10,3.708897e-005);
    RAA_0_100584__3__3__3->SetBinError(11,4.504344e-005);
    RAA_0_100584__3__3__3->SetBinError(12,5.385318e-005);
    RAA_0_100584__3__3__3->SetBinError(13,6.332594e-005);
    RAA_0_100584__3__3__3->SetBinError(14,4.278537e-005);
    RAA_0_100584__3__3__3->SetBinError(15,6.266848e-005);
    RAA_0_100584__3__3__3->SetBinError(16,8.299302e-005);
    RAA_0_100584__3__3__3->SetBinError(17,0.0001099248);
    RAA_0_100584__3__3__3->SetBinError(18,0.0001477243);
    RAA_0_100584__3__3__3->SetBinError(19,0.0002002684);
    RAA_0_100584__3__3__3->SetBinError(20,0.0001962782);
    RAA_0_100584__3__3__3->SetBinError(21,0.00043599);
    RAA_0_100584__3__3__3->SetBinError(22,0.0007955915);
    RAA_0_100584__3__3__3->SetBinError(23,0.001121086);
    RAA_0_100584__3__3__3->SetBinError(24,0.002349183);
    RAA_0_100584__3__3__3->SetBinError(25,0.003934258);
    RAA_0_100584__3__3__3->SetBinError(26,0.005431548);
    RAA_0_100584__3__3__3->SetBinError(27,0.007593713);
    RAA_0_100584__3__3__3->SetBinError(28,0.008532134);
    RAA_0_100584__3__3__3->SetBinError(29,0.006423681);
    RAA_0_100584__3__3__3->SetBinError(30,0.009668077);
    RAA_0_100584__3__3__3->SetBinError(31,0.009013466);
    RAA_0_100584__3__3__3->SetBinError(32,0.01230905);
    RAA_0_100584__3__3__3->SetBinError(33,0.01985833);
    RAA_0_100584__3__3__3->SetBinError(34,0.02497965);
    RAA_0_100584__3__3__3->SetBinError(35,0.04145952);
    RAA_0_100584__3__3__3->SetBinError(36,0.0459072);
    RAA_0_100584__3__3__3->SetBinError(37,0.152905);
    RAA_0_100584__3__3__3->SetMinimum(0);
    RAA_0_100584__3__3__3->SetMaximum(1.8);
    RAA_0_100584__3__3__3->SetEntries(11078.9);
    RAA_0_100584__3__3__3->SetDirectory(0);

    ci = TColor::GetColor("#ffcc00");
    RAA_0_100584__3__3__3->SetFillColor(ci);

    ci = TColor::GetColor("#000099");
    RAA_0_100584__3__3__3->SetLineColor(ci);
    RAA_0_100584__3__3__3->SetLineWidth(2);
    RAA_0_100584__3__3__3->SetMarkerStyle(20);
    RAA_0_100584__3__3__3->SetMarkerSize(1.3);
    RAA_0_100584__3__3__3->GetXaxis()->SetTitle("p_{T} (GeV)");
    RAA_0_100584__3__3__3->GetXaxis()->SetRange(1,37);
    RAA_0_100584__3__3__3->GetXaxis()->CenterTitle(true);
    RAA_0_100584__3__3__3->GetXaxis()->SetLabelFont(42);
    RAA_0_100584__3__3__3->GetXaxis()->SetLabelOffset(-0.005);
    RAA_0_100584__3__3__3->GetXaxis()->SetLabelSize(0.035);
    RAA_0_100584__3__3__3->GetXaxis()->SetTitleSize(0.05);
    RAA_0_100584__3__3__3->GetXaxis()->SetTitleFont(42);
    RAA_0_100584__3__3__3->GetYaxis()->SetTitle("R_{AA}, R_{pA}");
    RAA_0_100584__3__3__3->GetYaxis()->CenterTitle(true);
    RAA_0_100584__3__3__3->GetYaxis()->SetNdivisions(508);
    RAA_0_100584__3__3__3->GetYaxis()->SetLabelFont(42);
    RAA_0_100584__3__3__3->GetYaxis()->SetLabelSize(0.035);
    RAA_0_100584__3__3__3->GetYaxis()->SetTitleSize(0.035);
    RAA_0_100584__3__3__3->GetYaxis()->SetTitleOffset(1.6);
    RAA_0_100584__3__3__3->GetYaxis()->SetTitleFont(42);
    RAA_0_100584__3__3__3->GetZaxis()->SetLabelFont(42);
    RAA_0_100584__3__3__3->GetZaxis()->SetLabelSize(0.05);
    RAA_0_100584__3__3__3->GetZaxis()->SetTitleSize(0.09);
    RAA_0_100584__3__3__3->GetZaxis()->SetTitleFont(42);
    RAA_0_100584__3__3__3->Draw("sameaxis");
    RAA->Modified();
    RAA->cd();
    RAA->SetSelected(RAA);
  */
}
