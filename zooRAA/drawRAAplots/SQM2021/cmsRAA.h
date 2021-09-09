#ifndef _CMSRAA_H_
#define _CMSRAA_H_

#include <TLatex.h>
#include <TLine.h>

namespace cmsRAA
{
  template <class T>
  void sethemptystyle(T* hempty, Float_t xtitleoffset=-1, Float_t ytitleoffset=-1, Float_t xtitlesize=-1, Float_t ytitlesize=-1, Float_t xlabelsize=-1, Float_t ylabelsize=-1)
  {
    hempty->GetXaxis()->CenterTitle();
    hempty->GetYaxis()->CenterTitle();
    hempty->GetXaxis()->SetTitleFont(42);
    hempty->GetYaxis()->SetTitleFont(42);
    if(xtitleoffset>=0) hempty->GetXaxis()->SetTitleOffset(xtitleoffset);
    if(ytitleoffset>=0) hempty->GetYaxis()->SetTitleOffset(ytitleoffset);
    if(xtitlesize>=0) hempty->GetXaxis()->SetTitleSize(xtitlesize);
    if(ytitlesize>=0) hempty->GetYaxis()->SetTitleSize(ytitlesize);
    hempty->GetXaxis()->SetLabelFont(42);
    hempty->GetYaxis()->SetLabelFont(42);
    if(xlabelsize>=0) hempty->GetXaxis()->SetLabelSize(xlabelsize);
    if(ylabelsize>=0) hempty->GetYaxis()->SetLabelSize(ylabelsize);
    hempty->SetStats(0);
  }

  template <class T>
  void setthgrstyle(T* h, Color_t mcolor=-1, Style_t mstyle=-1, Size_t msize=-1, Color_t lcolor=-1, Style_t lstyle=-1, Width_t lwidth=-1, Color_t fcolor=-1, Float_t falpha=-1, Style_t fstyle=-1)
  {
    if(mcolor>=0) h->SetMarkerColor(mcolor);
    if(mstyle>=0) h->SetMarkerStyle(mstyle);
    if(msize>=0)  h->SetMarkerSize(msize);
    if(lcolor>=0) h->SetLineColor(lcolor);
    if(lstyle>=0) h->SetLineStyle(lstyle);
    if(lwidth>=0) h->SetLineWidth(lwidth);
    if(fcolor>=0) h->SetFillColor(fcolor);
    if(falpha>=0) h->SetFillColorAlpha(fcolor, falpha);
    if(fstyle>=0) h->SetFillStyle(fstyle);
  }

  void settex(TLatex* tex, Float_t tsize=0.04, Short_t align=11, Font_t tfont=42, Color_t tcolor=kBlack)
  {
    tex->SetNDC();
    if(tfont>0) tex->SetTextFont(tfont);
    if(align>0) tex->SetTextAlign(align);
    if(tsize>0) tex->SetTextSize(tsize);
    if(tcolor>0) tex->SetTextColor(tcolor);
  }

  void drawtex(Double_t x, Double_t y, const char* text, Float_t tsize=0.04, Short_t align=12, Style_t font=42, Color_t tcolor=kBlack)
  {
    TLatex* tex = new TLatex(x, y, text);
    settex(tex, tsize, align, font, tcolor);
    tex->Draw();
  }

  void setleg(TLegend* leg, Float_t size=0.04)
  {
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(size);
  }

  void drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor=kBlack, Style_t lstyle=1, Width_t lwidth=2)
  {
    TLine* l = new TLine(x1, y1, x2, y2);
    l->SetLineColor(lcolor);
    l->SetLineStyle(lstyle);
    l->SetLineWidth(lwidth);
    l->Draw();
  }

  void drawbox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t fcolor=kGray, Float_t falpha=0.4, Style_t fstyle=1001, Color_t lcolor=0, Style_t lstyle=1, Width_t lwidth=0)
  {
    TBox* b = new TBox(x1, y1, x2, y2);
    b->SetFillColor(fcolor);
    b->SetFillColorAlpha(fcolor, falpha);
    b->SetFillStyle(fstyle);
    b->SetLineColor(lcolor);
    b->SetLineStyle(lstyle);
    b->SetLineWidth(lwidth);
    b->Draw();
  }
}

#endif

#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TBox.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>

