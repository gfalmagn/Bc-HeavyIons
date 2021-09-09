#include <TGraphErrors.h>

namespace npjpsi2
{
  TGraphErrors* gre2TeV;
  TGraphErrors* gre2TeV2;
  TGraphErrors* gre2TeVstat;
  TGraphErrors* gre2TeV2stat;

  void expBeautyCMS_20170201_marker()
  {
    gre2TeVstat->Draw("p same");
    gre2TeV2stat->Draw("p same");
  }

  void expBeautyCMS_20170201()
  {
    //=========Macro generated from canvas: c2/c2
    //=========  (Wed Feb  1 19:16:45 2017) by ROOT version5.34/36

    Color_t fcolor = kViolet-8; // kGray+2
    Float_t falpha = 0.3; // 0.5
    Color_t mcolor = kViolet+4; // kGray+3
    Float_t logxerrfactor = 0.08;

    gre2TeV = new TGraphErrors(6);
    gre2TeV->SetName("Graph0");
    gre2TeV->SetTitle("Graph0");

    gre2TeV->SetFillColor(fcolor);
    gre2TeV->SetFillColorAlpha(fcolor, falpha);
    gre2TeV->SetMarkerColor(mcolor);
    gre2TeV->SetMarkerStyle(29);
    gre2TeV->SetMarkerSize(2);
    gre2TeV->SetLineColor(0);
    gre2TeV->SetPoint(0,      7.5,                 0.4664189);
    gre2TeV->SetPointError(0, 7.5*logxerrfactor,   0.05427119); // 0.5
    gre2TeV->SetPoint(1,      9,                   0.4856157);
    gre2TeV->SetPointError(1, 9*logxerrfactor,     0.04438349); // 0.5
    gre2TeV->SetPoint(2,      10.25,               0.4065521);
    gre2TeV->SetPointError(2, 10.25*logxerrfactor, 0.03578226); // 0.5
    gre2TeV->SetPoint(3,      12,                  0.4536073);
    gre2TeV->SetPointError(3, 12*logxerrfactor,    0.05860495); // 0.5
    gre2TeV->SetPoint(4,      14.5,                0.3694549);
    gre2TeV->SetPointError(4, 14.5*logxerrfactor,  0.0518269);  // 0.5
    gre2TeV->SetPoint(5,      23,                  0.3521646);
    gre2TeV->SetPointError(5, 23*logxerrfactor,    0.06848064); // 0.5
    gre2TeV->Draw("2 same");

    gre2TeVstat = new TGraphErrors(6);
    gre2TeVstat->SetName("Graph0stat");
    gre2TeVstat->SetTitle("Graph0stat");
    gre2TeVstat->SetMarkerColor(mcolor);
    gre2TeVstat->SetMarkerStyle(29);
    gre2TeVstat->SetMarkerSize(2);
    gre2TeVstat->SetLineWidth(1);
    gre2TeVstat->SetPoint(0,      7.5,   0.4664189);
    gre2TeVstat->SetPointError(0, 1,     0.02632702);
    gre2TeVstat->SetPoint(1,      9,     0.4856157);
    gre2TeVstat->SetPointError(1, 0.5,   0.0353226);
    gre2TeVstat->SetPoint(2,      10.25, 0.4065521);
    gre2TeVstat->SetPointError(2, 0.75,  0.0288783);
    gre2TeVstat->SetPoint(3,      12,    0.4536073);
    gre2TeVstat->SetPointError(3, 1,     0.03347202);
    gre2TeVstat->SetPoint(4,      14.5,  0.3694549);
    gre2TeVstat->SetPointError(4, 1.5,   0.03358354);
    gre2TeVstat->SetPoint(5,      23,    0.3521646);
    gre2TeVstat->SetPointError(5, 7,     0.0311077);
    // gre2TeVstat->Draw("p same");
   
    gre2TeV2 = new TGraphErrors(6);
    gre2TeV2->SetName("Graph1");
    gre2TeV2->SetTitle("Graph1");
    gre2TeV2->SetFillColor(fcolor);
    gre2TeV2->SetFillColorAlpha(fcolor, falpha);
    gre2TeV2->SetMarkerColor(mcolor);
    gre2TeV2->SetMarkerStyle(34);
    gre2TeV2->SetMarkerSize(1.7);
    gre2TeV2->SetLineColor(0);
    gre2TeV2->SetPoint(0,      3.75,               0.7535804);
    gre2TeV2->SetPointError(0, 3.75*logxerrfactor, 0.21432);    // 0.4
    gre2TeV2->SetPoint(1,      5,                  0.6891029);
    gre2TeV2->SetPointError(1, 5*logxerrfactor,    0.1218303);  // 0.4
    gre2TeV2->SetPoint(2,      6,                  0.5425812);
    gre2TeV2->SetPointError(2, 6*logxerrfactor,    0.08921681); // 0.4
    gre2TeV2->Draw("2 same");

    gre2TeV2stat = new TGraphErrors(3);
    gre2TeV2stat->SetName("Graph1stat");
    gre2TeV2stat->SetTitle("Graph1stat");
    gre2TeV2stat->SetMarkerColor(mcolor);
    gre2TeV2stat->SetMarkerStyle(34);
    gre2TeV2stat->SetMarkerSize(1.7);
    gre2TeV2stat->SetLineWidth(1);
    gre2TeV2stat->SetPoint(0,      3.75, 0.7535804);
    gre2TeV2stat->SetPointError(0, 0.75, 0.1209346);
    gre2TeV2stat->SetPoint(1,      5,    0.6891029);
    gre2TeV2stat->SetPointError(1, 0.5,  0.1059286);
    gre2TeV2stat->SetPoint(2,      6,    0.5425812);
    gre2TeV2stat->SetPointError(2, 0.5,  0.07363414);
    // gre2TeV2stat->Draw("p same");
  }
}
