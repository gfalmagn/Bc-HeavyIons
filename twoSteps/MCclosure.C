#include "../helpers/hub.cpp"

void MCclosure(){

  //Create Hub gathering various data 
  Hub H = Hub(true, true);
  //H.SetxLW();
  H.SetMCclosure();

  vector<TH1F*> Ncor_true;
  vector<TH1F*> Ncor_step1;
  vector<TH1F*> Ncor_step2;
  vector<TH1F*> Ncor_step3;

  for(int t=0;t<_nMCclos;t++){
    //initialise hists
    Ncor_true.push_back( new TH1F("Ncor_true_toy"+(TString)to_string(t),"Ncor_true_toy"+(TString)to_string(t),3,0,3) );
    Ncor_step1.push_back( new TH1F("Ncor_step1_toy"+(TString)to_string(t),"Ncor_step1_toy"+(TString)to_string(t),3,0,3) );
    Ncor_step2.push_back( new TH1F("Ncor_step2_toy"+(TString)to_string(t),"Ncor_step2_toy"+(TString)to_string(t),3,0,3) );
    Ncor_step3.push_back( new TH1F("Ncor_step3_toy"+(TString)to_string(t),"Ncor_step3_toy"+(TString)to_string(t),3,0,3) );

    //fill hists
    for(int b=0;b<=_NanaBins;b++){
      Ncor_true[t]->SetBinContent(b+1,H.Ycorr_MCclos[t].Val[b]);
      Ncor_step1[t]->SetBinContent(b+1,H.Ycorr_MCclos_step1[t].Val[b]);
      Ncor_step2[t]->SetBinContent(b+1,H.Ycorr_MCclos_step2[t].Val[b]);
      Ncor_step3[t]->SetBinContent(b+1,H.Ycorr_MCclos_step3[t].Val[b]);
    }

    //style of hists
    Ncor_true[t]->SetLineColor((t==0)?kRed:kBlue);
    Ncor_step1[t]->SetLineColor((t==0)?kRed:kBlue);
    Ncor_step2[t]->SetLineColor((t==0)?kRed:kBlue);
    Ncor_step3[t]->SetLineColor((t==0)?kRed:kBlue);
    Ncor_true[t]->SetLineStyle(1);
    Ncor_step1[t]->SetLineStyle(2);
    Ncor_step2[t]->SetLineStyle(9);
    Ncor_step3[t]->SetLineStyle(10);
    Ncor_true[t]->SetLineWidth(3);
    Ncor_step1[t]->SetLineWidth(3);
    Ncor_step2[t]->SetLineWidth(3);
    Ncor_step3[t]->SetLineWidth(3);
    Ncor_true[t]->SetMarkerSize(0);
    Ncor_step1[t]->SetMarkerSize(0);
    Ncor_step2[t]->SetMarkerSize(0);
    Ncor_step3[t]->SetMarkerSize(0);
    Ncor_true[t]->GetXaxis()->SetBinLabel(1,"PbPb integrated");
    Ncor_true[t]->GetXaxis()->SetBinLabel(2,"PbPb p_{T} bin1");
    Ncor_true[t]->GetXaxis()->SetBinLabel(3,"PbPb p_{T} bin2");
  }

  TLegend* leg = new TLegend(0.49,0.6,0.95,0.95);
  for(int t=0;t<_nMCclos;t++){
    leg->AddEntry(Ncor_true[t],"true, toy"+(TString)to_string(t+1),"l");
    leg->AddEntry(Ncor_step1[t],"original MC, toy"+(TString)to_string(t+1),"l");
    leg->AddEntry(Ncor_step2[t],"MC 1st-step p_{T} correction, toy"+(TString)to_string(t+1),"l");
    leg->AddEntry(Ncor_step3[t],"MC 2nd-step p_{T} correction, toy"+(TString)to_string(t+1),"l");
  }
  leg->SetTextSize(0.029);
  leg->SetBorderSize(0);

  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("c1","c1",1500,1500);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.15);
  Ncor_true[1]->GetYaxis()->SetRangeUser(0,1.5*Ncor_true[1]->GetMaximum());
  Ncor_true[1]->GetYaxis()->SetTitle("corrected yield");
  Ncor_true[1]->SetTitle("");

  for(int t=_nMCclos-1;t>=0;t--){
    Ncor_true[t]->Draw((t==_nMCclos)?"E":"Esame");
    Ncor_step1[t]->Draw("Esame");
    Ncor_step2[t]->Draw("Esame");
    Ncor_step3[t]->Draw("Esame");
  }
  leg->Draw("same");

  c1->SaveAs("figs/MCclosureTest_forTwoSteps.pdf");
  
}
