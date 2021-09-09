#include "Tools.h"
#include "ChargedHad/RAA_0_10.C"
#include "ChargedHad/RAA_0_100.C"
#include "BmesonRaa/canvasRAAPbPb_0_100_ThmRAA.C"
#include "BtoDmesonRaa/canvasRAAPbPb_0_100_BtoDRAA.C"
#include "NonpromptJpsi/nonPrompt_276raa_20170201.h"
#include "NonpromptJpsi/nonPrompt_502raa_20170712.h"
#include "NonpromptJpsi/Prompt_502raa.h"
#include "BcRAA/BcRAA.C"
#include "Bs/BsRAA.C"
#include "systematics.h"
#include "nicePalette.h"

#include "cmsRAA.h"

void cmsRAAdraw(TString fileMB, TString file, Float_t centmin, Float_t centmax, 
                int isD, Int_t isHad, Int_t isB, Int_t isNjpsi, Int_t isND, Int_t isPjpsi, int isBc, int isBs, int isPsi2, int isUps,
                bool savepng=false)
{
  if(centmin!=0 || (centmax!=100 && (isB || isNjpsi || isND || isPjpsi))) { return; }

  gStyle->SetOptTitle(0);   
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);

  bool rewriteLumis = isBc && false;

  // prompt D0 -->  
  TFile* filePPMB = new TFile(fileMB.Data());  
  TGraphAsymmErrors* gNuclearModificationMB = (TGraphAsymmErrors*)filePPMB->Get("gNuclearModification");
  TH1D* hNuclearModificationMB = (TH1D*)filePPMB->Get("hNuclearModification");
  TFile* filePP = new TFile(file.Data());  
  TGraphAsymmErrors* gNuclearModification = (TGraphAsymmErrors*)filePP->Get("gNuclearModification");
  TH1D* hNuclearModification = (TH1D*)filePP->Get("hNuclearModification");
  cmsRAA::setthgrstyle(gNuclearModification, kGreen+3, 21, 1, kGreen+3, 1, 0, kGreen-9, 0.7, 1001);
  cmsRAA::setthgrstyle(gNuclearModificationMB, kGreen+3, 21, 1, kGreen+3, 1, 0, kGreen-9, 0.7, 1001);
  cmsRAA::setthgrstyle(hNuclearModification, kGreen+3, 21, 1., kGreen+3, 1, 3);
  cmsRAA::setthgrstyle(hNuclearModificationMB, kGreen+3, 21, 1., kGreen+3, 1, 3);
  // <-- prompt D0

  // Prompt J/psi + Upsi -->  
  TFile* fileJpsiUps = new TFile("RAA_vs_pt_BottomoniaCharmonia_5p02TeV_PbPb.root","read");  
  TGraphErrors* Jpsi_FullRap_stat = (isPjpsi && isPsi2)?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Jpsi_FullRap_Pt6p5_stat")):NULL;
  TGraphErrors* Jpsi_FullRap_syst = (isPjpsi && isPsi2)?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Jpsi_FullRap_Pt6p5_sys")):NULL;
  TGraphErrors* Jpsi_ForwRap_stat = (isPjpsi && isPsi2)?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Jpsi_ForwardRap_Pt3p0_stat")):NULL;
  TGraphErrors* Jpsi_ForwRap_syst = (isPjpsi && isPsi2)?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Jpsi_ForwardRap_Pt3p0_sys")):NULL;
  TGraphErrors* Psi2S_stat = isPsi2?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_psi2S_Rap1p6_Pt6p5_stat")):NULL;
  TGraphErrors* Psi2S_syst = isPsi2?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_psi2S_Rap1p6_Pt6p5_sys")):NULL;
  TGraphErrors* Ups1S_stat = isUps?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Ups1S_stat")):NULL;
  TGraphErrors* Ups1S_syst = isUps?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Ups1S_sys")):NULL;
  TGraphErrors* Ups2S_stat = isUps?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Ups2S_stat")):NULL;
  TGraphErrors* Ups2S_syst = isUps?((TGraphErrors*)fileJpsiUps->Get("graph_RAA_Ups2S_sys")):NULL;
  if(isPjpsi && isPsi2){
    cmsRAA::setthgrstyle(Jpsi_FullRap_stat, kviolet, 29, 1.5, kviolet, 1, 3, kvioletLight, 0.5);
    cmsRAA::setthgrstyle(Jpsi_FullRap_syst, kviolet, 29, 1.5, kviolet, 1, 0, kvioletLight, 0.5);
    cmsRAA::setthgrstyle(Jpsi_ForwRap_stat, kviolet, 30, 1.5, kviolet, 1, 3, kvioletLight, 0.5);
    cmsRAA::setthgrstyle(Jpsi_ForwRap_syst, kviolet, 30, 1.5, kviolet, 1, 0, kvioletLight, 0.5);
    cmsRAA::setthgrstyle(Psi2S_stat, korange, 47, 1.2, korange, 1, 3, korangeLight, 0.5);
    cmsRAA::setthgrstyle(Psi2S_syst, korange, 47, 1.2, korange, 1, 0, korangeLight, 0.5);
    cmsRAA::setthgrstyle(Ups1S_stat, kyellow, 21, 1, kyellow, 1, 3, kyellowLight, 0.5);
    cmsRAA::setthgrstyle(Ups1S_syst, kyellow, 21, 1, kyellow, 1, 0, kyellowLight, 0.5);
    cmsRAA::setthgrstyle(Ups2S_stat, kspring, 22, 1.2, kspring, 1, 3, kspringLight, 0.5);
    cmsRAA::setthgrstyle(Ups2S_syst, kspring, 22, 1.2, kspring, 1, 0, kspringLight, 0.5);
  }
  // <-- Prompt J/psi + Upsi

  TCanvas* canvasRAA = new TCanvas("canvasRAA", "canvasRAA", 600, 600);
  canvasRAA->cd();
  canvasRAA->SetFillColor(0);
  canvasRAA->SetBorderMode(0);
  canvasRAA->SetBorderSize(2);
  canvasRAA->SetLeftMargin(0.13);
  canvasRAA->SetRightMargin(0.025);
  canvasRAA->SetTopMargin(0.080);
  canvasRAA->SetBottomMargin(0.13);
  canvasRAA->SetFrameBorderMode(0);
  if(!(isBc && isUps && isPsi2)) canvasRAA->SetLogx();

  Float_t xaxismin = isBc?((isUps || isPsi2)?0.:1.):0.7; // 0.7, 1.0
  Float_t xaxismax = isBc?((isUps || isPsi2)?50:400):400; // 150, 400
  float yaxismax = isBc?((isUps || isPsi2)?3.1:3.1):1.75;
  TH2F* hemptyRAA = new TH2F("hemptyRAA", ";p_{T} [GeV];R_{AA}", 50, xaxismin, xaxismax, 10, 0, yaxismax);
  cmsRAA::sethemptystyle(hemptyRAA, 0.9, 1.0, 0.06, 0.06, 0.045, 0.045);
  hemptyRAA->GetXaxis()->SetLabelOffset(0.0);
  hemptyRAA->Draw();

  cmsRAA::drawline(xaxismin, 1, xaxismax, 1, kBlack, 2, 2);

  /* <-- syst --> */
  // charged particles
  if(isHad)
    { if(centmin==0 && centmax==100) { RAA_0_100(); }
      if(centmin==0 && centmax==10) { RAA_0_10(); } }
  // prompt D0
  if(isD)
    { gNuclearModification->Draw("2same");
      gNuclearModificationMB->Draw("2same"); }
  // B+
  if(isB==1 && centmin==0 && centmax==100)
    { bplus::canvasRAAPbPb_0_100_ThmRAA(); }
  // nonprompt jpsi 2.76 TeV
  if(isNjpsi==2 && centmin==0 && centmax==100)
    { npjpsi2::expBeautyCMS_20170201(); }
  // nonprompt jpsi 5.02 TeV
  if(isNjpsi==1 && centmin==0 && centmax==100)
    { npjpsi5::expBeautyCMS(); }
  // nonprompt D0
  if(isND==1 && centmin==0 && centmax==100)
    { npd::canvasRAAPbPb_0_100_BtoDRAA(); }
  // prompt jpsi 5.02 TeV
  if(isPjpsi && centmin==0 && centmax==100) { 
    if(isPsi2){
      Jpsi_FullRap_syst->Draw("2same");
      Jpsi_ForwRap_syst->Draw("2same");
    }
    else
      pjpsi5::expBeautyCMS(); 
  }
  //Psi(2S)
  if(isPsi2){
    Psi2S_syst->Draw("2same");
  }
  //Upsilon(nS)
  if(isUps){
    Ups1S_syst->Draw("2same");
    Ups2S_syst->Draw("2same");
  }
  // Bs
  if(isBs)
    { Bs::RAA(isBc); }
  // Bc
  if(isBc)
    { Bc::RAA(); }

  /* <-- markers --> */
  if(isHad==1)
    { if(centmin==0 && centmax==100) { RAA_0_100_marker(isBc); }
      if(centmin==0 && centmax==10) { RAA_0_10_marker(); } }
  if(isD)
    { hNuclearModification->Draw(isBc?"2EX0same":"2same");
      hNuclearModificationMB->Draw(isBc?"2EX0same":"2same"); }
  if(isB==1 && centmin==0 && centmax==100)
    { bplus::canvasRAAPbPb_0_100_ThmRAA_marker(isBc); }
  if(isNjpsi==2 && centmin==0 && centmax==100)
    { npjpsi2::expBeautyCMS_20170201_marker(); }
  if(isNjpsi==1 && centmin==0 && centmax==100)
    { npjpsi5::expBeautyCMS_marker(); }
  if(isND==1 && centmin==0 && centmax==100)
    { npd::canvasRAAPbPb_0_100_BtoDRAA_marker(); }
  if(isPjpsi && centmin==0 && centmax==100){ 
    if(isPsi2){
      Jpsi_FullRap_stat->Draw("psame");
      Jpsi_ForwRap_stat->Draw("psame");
    }
    else
      pjpsi5::expBeautyCMS_marker(); 
  }
  if(isPsi2){
    Psi2S_stat->Draw("psame");
  }
  if(isUps){
    Ups1S_stat->Draw("psame");
    Ups2S_stat->Draw("psame");
  }
  if(isBs)
    Bs::RAA_markers();

  Float_t systnormup = normalizationUncertaintyForRAA(centmin, centmax, true)*1.e-2;
  Float_t systnormlo = normalizationUncertaintyForRAA(centmin, centmax, false)*1.e-2;
  TBox* bSystnorm = new TBox(isUps?xaxismax:xaxismin, 1-systnormlo, isUps?(xaxismax-1.3):(xaxismin+0.2), 1+systnormup);
  bSystnorm->SetLineColor(16);
  bSystnorm->SetFillColor(16);
  bSystnorm->Draw();

  // legend preset -- >
  const std::string theader[] = {"Light", "Charm", "Beauty", "Open charm", "Beauty-Strange", "Beauty-Charm"};
  int nlinel = -1, nliner = -1, itheaderl = -1, itheaderll = -1, itheaderr = -1;
  // charged + charm + beauty
  if(isHad && isD && isB && isNjpsi && isND && !isPjpsi) { nlinel = 2; nliner = 5; itheaderl = 0; itheaderll = 1; itheaderr = 2; }
  // charm only
  if(!isHad && isD && !isB && !isNjpsi && !isND && !isPjpsi) { nlinel = 1; nliner = 0; itheaderl = 1; itheaderll = -1; itheaderr = -1; }
  // beauty only
  if(!isHad && !isD && isB && isNjpsi && isND && !isPjpsi) { nlinel = 0; nliner = 5; itheaderl = -1; itheaderll = -1; itheaderr = 2; }
  // charm + beauty
  if(!isHad && isD && isB && isNjpsi && isND && !isPjpsi) { nlinel = 1; nliner = 5; itheaderl = 1; itheaderll = -1; itheaderr = 2; }
  // charm + charged
  if(isHad && isD && !isB && !isNjpsi && !isND && !isPjpsi) { nlinel = 3; nliner = 0; itheaderl = 0; itheaderll = 1; itheaderr = -1; }
  // beauty + charged
  if(isHad && !isD && isB && isNjpsi && isND && !isPjpsi) { nlinel = 1; nliner = 5; itheaderl = 0; itheaderll = -1; itheaderr = 2; }
  // charm + nonprompt jpsi
  if(!isHad && isD && !isB && isNjpsi && !isND && !isPjpsi) { nlinel = 1; nliner = 3; itheaderl = 1; itheaderll = -1; itheaderr = 2; }
  // charm + nonprompt D
  if(!isHad && isD && !isB && !isNjpsi && isND && !isPjpsi) { nlinel = 1; nliner = 1; itheaderl = 1; itheaderll = -1; itheaderr = 2; }
  // charm + prompt jpsi
  if(!isHad && isD && !isB && !isNjpsi && !isND && isPjpsi) { nlinel = 0; nliner = 5; itheaderl = -1; itheaderll = -1; itheaderr = 3; }
  // 
  // charged + charm + B + NPD
  if(isHad && isD && isB && !isNjpsi && isND && !isPjpsi) { nlinel = 3; nliner = 2; itheaderl = 0; itheaderll = 1; itheaderr = 2; }
  // charged + charm + B + NPjpsi
  if(isHad && isD && isB && isNjpsi && !isND && !isPjpsi) { nlinel = 3; nliner = 4; itheaderl = 0; itheaderll = 1; itheaderr = 2; }
  // charged + charm + B
  if(isHad && isD && isB && !isNjpsi && !isND && !isPjpsi) { nlinel = 3; nliner = 1; itheaderl = 0; itheaderll = 1; itheaderr = 2; }
  // charged + charm + NPjpsi
  if(isHad && isD && !isB && isNjpsi && !isND && !isPjpsi) { nlinel = 3; nliner = 3; itheaderl = 0; itheaderll = 1; itheaderr = 2; }
  // charged + charm + NPD
  if(isHad && isD && !isB && !isNjpsi && isND && !isPjpsi) { nlinel = 3; nliner = 1; itheaderl = 0; itheaderll = 1; itheaderr = 2; }
  // Bc + charged + Bs + D + B + Jpsi
  if(isBc && isHad && isBs && isD && isB) { nlinel = 0; nliner = rewriteLumis?10:9; itheaderl = -1; itheaderll = -1; itheaderr = -1; }
  // Bc + pJpsi + Psi2S + Upsilon
  if(isBc && isPjpsi && isPsi2 && isUps) { nlinel = 0; nliner = rewriteLumis?10:9; itheaderl = -1; itheaderll = -1; itheaderr = -1; }
  //
  if(nlinel < 0 || nliner < 0) { std::cout << " \033[31;1mwarning: this combination is not predefiend: " 
                                           << isD << " " << isHad << " " << isB << " " << isNjpsi << " " << isND << " " << isPjpsi<< "\033[0m" << std::endl; return; }
  // <-- legend preset

  double lineHei = isBc?(rewriteLumis?0.0465:0.0525):0.0486;
  double legx1_r = isBc?(isUps?0.65:(rewriteLumis?0.6:0.625)):0.62, legx2_r = legx1_r + 0.460;
  double legy2_r = 0.860 - (lineHei/2)*(5-nliner)*0.8, legy1_r = 0.617 + (lineHei/2)*(5-nliner)*1.2;
  if(isBc){ legy2_r = rewriteLumis?0.89:0.925; legy1_r = legy2_r - lineHei*nliner; }
  double yheadderr = itheaderr>=0?legy2_r + 0.01:-1;
  double legx1_l = 0.22, legx2_l = legx1_l + isBc?0.43:0.45, legy2_l = isBc?0.66:(itheaderl>=0?0.78:0.82), legy1_l = legy2_l - 0.0468*nlinel;
  double yheadderl = itheaderl>=0?legy2_l + 0.01:-1, yheadderll = itheaderll>=0?yheadderl-lineHei*2:-1;

  //Lumi and centrality
  TString tlumi_ = Form("%s%s%s", ((isD||isND)?"0.53/":""), (isHad?"0.40/":""), (isB||isNjpsi||isPjpsi?"0.37/":""));
  TString tlumi(tlumi_, tlumi_.Length()-1);
  if(isBc) tlumi = "0.37-1.6";
  TString tlumipp_ = Form("%s%s", ((isD||isHad||isPjpsi||isB)?"27-":""), (isBc?"302-":""));
  TString tlumipp(tlumipp_, tlumipp_.Length()-1);
  cmsRAA::drawtex(0.96, 0.936, Form("5.02 TeV PbPb (%s nb^{-1}) + pp (%s pb^{-1})", tlumi.Data(), tlumipp.Data()), 0.038, 31);
  if(!isBc){
    cmsRAA::drawtex(0.165, 0.89, "#bf{CMS}#scale[0.6]{#it{ Supplementary}}", 0.062, 13);//Supplementary
    if(!isB && !isPjpsi && !isNjpsi) cmsRAA::drawtex(0.95, 0.27, "|y| < 1", 0.04, 32);
    cmsRAA::drawtex(0.955, 0.22, Form("Cent. %.0f-%.0f%s", centmin, centmax, "%"), 0.04, 32);
  }
  else {
    cmsRAA::drawtex(0.165, 0.9, "#bf{CMS}", 0.062, 13);
    cmsRAA::drawtex(0.155, 0.84, "#scale[0.6]{#it{ Supplementary}}", 0.062, 13);//Supplementary
    cmsRAA::drawtex(legx1_r-0.013, legy1_r+lineHei*(rewriteLumis?9:7.25)+0.046, Form("#it{2015, centrality %.0f-%.0f%s}", centmin, centmax, "%"), 0.03, 11);
    if(rewriteLumis) cmsRAA::drawtex(legx1_r-0.013, legy1_r+lineHei*9+0.013, "#scale[0.93]{PbPb ("+(TString)(isUps?"0.37":"0.37-0.53")+" nb^{-1}), pp (27 pb^{-1})}", 0.03, 11);
    cmsRAA::drawtex(legx1_r-0.012, legy1_r+lineHei*(rewriteLumis?3:2.22)+(isUps?-0.:0.051), "#it{2017-18, centrality 0-90%}",0.03,11);
    if(rewriteLumis) cmsRAA::drawtex(legx1_r-0.012, legy1_r+lineHei*3+(isUps?-0.033:0.018), "#scale[0.93]{PbPb (1.6 nb^{-1}), pp (302 pb^{-1})}",0.03,11);
  }

  //Legends
  TLegend* legendRAA_r = new TLegend(legx1_r, legy1_r, legx2_r, legy2_r, "");
  cmsRAA::setleg(legendRAA_r, 0.034);
  TLegend* legendRAA_l = new TLegend(legx1_l, legy1_l, legx2_l, legy2_l, "");
  cmsRAA::setleg(legendRAA_l, 0.034);
  
  if(isHad) 
    {
      TLegend* leg = isBc?legendRAA_r:legendRAA_l;
      if(isBc) { leg->AddEntry((TObject*)0, "", NULL); }
      leg->AddEntry(gTrackPt_leg, (TString)(isBc?"#scale[1.2]{#bf{h^{+}}}":"h^{+}")+", |#eta| < 1", "pf");//h^{#pm}
      if(isD && !isBc) { leg->AddEntry((TObject*)0, "", NULL); }
    }
  if(isD)
    {
      TLegend* leg = (isBc || isPjpsi)?legendRAA_r:legendRAA_l;
      leg->AddEntry(gNuclearModification, (TString)(isBc?"#scale[1.2]{#bf{D#scale[0.6]{#lower[-0.7]{0}}}}":"D#scale[0.6]{#lower[-0.7]{0}}")+", |y| < 1", "pf");// + #bar{D}#scale[0.6]{#lower[-0.7]{0}}
    }
  if(isB && centmin==0 && centmax==100)
    { legendRAA_r->AddEntry(bplus::grae, (TString)(isBc?"#scale[1.2]{#bf{B^{+}}}":"B^{+}")+", |y| < 2.4", "pf"); } //B^{#pm}

  if(isND)
    { legendRAA_r->AddEntry(npd::grae2, "(b #rightarrow) D^{0}", "pf"); }

  if(isNjpsi==2 && centmin==0 && centmax==100)
    {
      legendRAA_r->AddEntry((TObject*)0, "", NULL);
      legendRAA_r->AddEntry(npjpsi2::gre2TeV2, "1.6 < |y| < 2.4", "pf");
      legendRAA_r->AddEntry(npjpsi2::gre2TeV, "|y| < 2.4", "pf");
      cmsRAA::drawtex(legx1_r+0.02, legy1_r+lineHei*2+0.014, "(b #rightarrow) J/#psi (2.76 TeV)", 0.034, 11);
    }

  if(isNjpsi==1 && centmin==0 && centmax==100)
    {
      legendRAA_r->AddEntry((TObject*)0, "", NULL);
      legendRAA_r->AddEntry(npjpsi5::gre5TeV2, "1.8 < |y| < 2.4", "pf");
      legendRAA_r->AddEntry(npjpsi5::gre5TeV, "|y| < 2.4", "pf");
      cmsRAA::drawtex(legx1_r+0.02, legy1_r+lineHei*2+0.014, "(b #rightarrow) J/#psi", 0.034, 11);
    }

  if(isPjpsi && centmin==0 && centmax==100)
    {
      legendRAA_r->AddEntry((TObject*)0, "", NULL);
      //legendRAA_r->AddEntry((TObject*)0, "", NULL);
      if(isBc) {
	cmsRAA::drawtex(legx1_r+0.0, legy1_r+lineHei*(rewriteLumis?7.5:6.55)+0.012, "prompt #scale[1.2]{#bf{J/#psi}}", 0.034, 31);
      } else{
	cmsRAA::drawtex(legx1_r+0.02, legy1_r+lineHei*3, "Hidden charm", 0.034, 11, 72);
	cmsRAA::drawtex(legx1_r+0.02, legy1_r+lineHei*2+0.01, "prompt J/#psi", 0.034, 11);
      }
      legendRAA_r->AddEntry(isPsi2?Jpsi_ForwRap_syst:(pjpsi5::gre5TeV2), "1.8 < |y| < 2.4", "pf");
      legendRAA_r->AddEntry(isPsi2?Jpsi_FullRap_syst:(pjpsi5::gre5TeV),  "|y| < 2.4", "pf");
    }
  if(isPsi2){
    legendRAA_r->AddEntry(Psi2S_syst, "|y| < 1.6", "pf");
    cmsRAA::drawtex(legx1_r+0.0, legy1_r+lineHei*(rewriteLumis?6:5.05)+0.012, "#scale[1.2]{#bf{#psi(2S)}}", 0.034, 31);
  }
  if(isUps){
    legendRAA_r->AddEntry(Ups1S_syst, "|y| < 2.4", "pf");
    legendRAA_r->AddEntry(Ups2S_syst, "|y| < 2.4", "pf");
    cmsRAA::drawtex(legx1_r+0.0, legy1_r+lineHei*(rewriteLumis?5:4.05)+0.012, "#scale[1.2]{#bf{#Upsilon(1S)}}", 0.034, 31);
    cmsRAA::drawtex(legx1_r+0.0, legy1_r+lineHei*(rewriteLumis?4:3.05)+0.012, "#scale[1.2]{#bf{#Upsilon(2S)}}", 0.034, 31);
  }
  if(isBs) 
    legendRAA_r->AddEntry("Bs_graphSyst", "#bf{#scale[1.2]{B#scale[0.6]{#lower[0.7]{s}}#kern[-1.]{#scale[0.6]{#lower[-0.7]{0}}}}}, |y| < 2.4", "pf");// + #bar{B}#scale[0.6]{#lower[0.7]{s}}#kern[-1.]{#scale[0.6]{#lower[-0.7]{0}}}
  if(isBc) 
    {
      legendRAA_r->AddEntry((TObject*)0, "", NULL);
      if(rewriteLumis) legendRAA_r->AddEntry((TObject*)0, "", NULL);
      if(!isUps) legendRAA_r->AddEntry((TObject*)0, "", NULL);
      legendRAA_r->AddEntry("Graph1", "1.3 < |y| < 2.3", "pf");//|y| < 2.3 & (|y| > 1.3 or p_{T} > 11 GeV)
      legendRAA_r->AddEntry("Graph2", "|y| < 2.3 ", "pf");//|y| < 2.3 & (|y| > 1.3 or p_{T} > 11 GeV)
      cmsRAA::drawtex(legx1_r+(isUps?-0.04:0.115), legy1_r + (isUps?(lineHei*1.2):(lineHei*2+0.015)), "#scale[1.2]{#bf{B_{c}^{+}}}"+(TString)(isUps?"":" (visible kin.)"), 0.034, isUps?31:11);//B_{c}^{#pm}
      if(isUps) cmsRAA::drawtex(legx1_r, legy1_r + lineHei*0.34, "(visible kin.)", 0.034, isUps?31:11);//B_{c}^{#pm}
    }

  if(nliner) { legendRAA_r->Draw(); }
  if(nlinel) { legendRAA_l->Draw(); }

  if(yheadderr > 0) cmsRAA::drawtex(legx1_r+0.017, yheadderr, theader[itheaderr].c_str(), 0.034, -1, 72);
  if(yheadderl > 0) cmsRAA::drawtex(legx1_l, yheadderl, theader[itheaderl].c_str(), 0.034, -1, 72);//+0.017
  if(yheadderll > 0) cmsRAA::drawtex(legx1_l+0.017, yheadderll, theader[itheaderll].c_str(), 0.034, -1, 72);
  
  cmsRAA::drawtex(isUps?0.815:0.16, (isUps?0.114:0.184)+0.77/yaxismax, "T_{AA} and lumi.", 0.027, 12, 42, kGray+1);//0.61
  cmsRAA::drawtex(isUps?0.815:0.16, (isUps?0.084:0.154)+0.77/yaxismax, "uncert. (2015)", 0.027, 12, 42, kGray+1);//0.58

  canvasRAA->Update();
  canvasRAA->RedrawAxis();

  TString texHad = isHad?"_charged":"";
  TString texD = isD?"_D":"";
  TString texB = isB?"_B":"";
  TString texNjpsi = isNjpsi>0?(isNjpsi==2?"_Njpsi2TeV":"_Njpsi5TeV"):"";
  TString texND = isND?"_BtoD":"";    
  TString texPjpsi = isPjpsi?"_Pjpsi5TeV":"";
  TString texBc = isBc?"_Bc":"";
  TString texBs = isBs?"_Bs":"";
  TString texPsi2 = isPsi2?"_Psi2":"";
  TString texUps = isUps?"_UpsNS":"";

  TString filename = Form("canvasRAA%s%s%s%s%s%s%s%s%s%s_cent_%.0f_%.0f", texBc.Data(), texBs.Data(), texPsi2.Data(), texUps.Data(), texHad.Data(), texD.Data(), texB.Data(), texNjpsi.Data(), texND.Data(), texPjpsi.Data(), centmin, centmax);
  canvasRAA->SaveAs(Form("plotRAA/%s.pdf", filename.Data()));
  if(savepng) canvasRAA->SaveAs(Form("plotRAA/%s.png", filename.Data()));

  std::cout<<"\033[33;1m| [[%ATTACHURL%/"<<filename<<".pdf][<img alt=\""<<filename<<".png\" src=\"%ATTACHURLPATH%/"<<filename<<".png\" width=\"200\" />]]\033[0m"<<std::endl;
}

int main(int argc, char *argv[])
{
  if(argc==15)
    {
      cmsRAAdraw(argv[1], argv[2], atof(argv[3]), atof(argv[4]), 
                 atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atoi(argv[14]));
      return 0;
    }
  std::cout << "Wrong number of inputs" << std::endl;
  return 1;
}
