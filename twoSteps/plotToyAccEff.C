#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH2Poly.h"
#include "TLine.h"
#include "TClonesArray.h"
#include "TStyle.h"
#include "TLegend.h"
#include <TVirtualFitter.h>
#include "../PbPb18/Utilities/EVENTUTILS.h"
#include "../helpers/Cuts_BDT.h"
#include "../helpers/Cuts.h"
#include "../acceptance/SgMuonAcceptanceCuts.h"
#include "../helpers/AccEff2DBinning.h"

void plotToyAccEff(bool ispp=true){

  //Grab acc and eff for all toys
  vector<vector<vector<float> > > *acc_biased; //(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) ))
  TFile *facc = new TFile("../acceptance/acceptanceMap.root","READ");
  facc->GetObject("acceptance_oneBinned_biased", acc_biased); 

  //eff
  vector<vector<float> > *eff_biased;
  TFile *feff = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
  feff->GetObject("efficiency_oneBinned_biased"+(TString)(ispp?"_pp":"_PbPb"), eff_biased); 
  //names of "fit" methods of pT distributions
  

  //Grab nominal acc and eff, without second-step pT biasing of MC
  vector<float> *acc_oneBinned;
  facc->GetObject("acceptance_oneBinned", acc_oneBinned);
  vector<vector<float> > *eff_oneBinned;
  feff->GetObject("efficiency_oneBinned"+(TString)(ispp?"_pp":"_PbPb"), eff_oneBinned);

  //histograms of acc and eff filled with toys
  vector<vector<TH1F*> > acc_toys(_biasNmeth);
  vector<vector<TH1F*> > eff_toys(_biasNmeth);
  vector<vector<TH1F*> > acceff_toys(_biasNmeth);

  int col = ispp?0:1;
  for(int b=1;b<=_NanaBins;b++){
    
    //histograms definition
    float typacc = (*acc_oneBinned)[b];
    float typeff = (*eff_oneBinned)[b][0];
    float typacceff = typacc*typeff;
    for(int m=0;m<_biasNmeth;m++){

      acc_toys[m].push_back(new TH1F("acc_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "acceptance toys "+(TString)(ispp?"pp":"PbPb")+" bin"+(TString)to_string(b)+";acceptance;n_{toys}", 22, (ispp?0.98:0.75)*typacc, (ispp?1.12:1.2)*typacc));
      eff_toys[m].push_back(new TH1F("eff_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "efficiency toys "+(TString)(ispp?"pp":"PbPb")+" bin"+(TString)to_string(b)+";efficiency;n_{toys}", 22, (ispp?0.99:0.8)*typeff, (ispp?1.04:1.2)*typeff));
      acceff_toys[m].push_back(new TH1F("acceff_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "acc #times eff toys "+(TString)(ispp?"pp":"PbPb")+" bin"+(TString)to_string(b)+";acc #times eff ;n_{toys}", 22, (ispp?0.98:0.65)*typacceff, (ispp?1.15:1.2)*typacceff));
      
      //Run over all toys
      for(int t=0;t<_biasNtoys;t++){
	int ntoy = _biasNmeth + m*_biasNtoys + t;
	//cout<<"b toyNumber acc eff = "<<b<<" "<<t<<" "<< (*acc_biased)[col][b][ntoy] <<" "<<(*eff_biased)[b][ntoy]<<endl;
	acc_toys[m][b-1]->Fill((*acc_biased)[col][b][ntoy]);
	eff_toys[m][b-1]->Fill((*eff_biased)[b][ntoy]);
	acceff_toys[m][b-1]->Fill((*acc_biased)[col][b][ntoy] * (*eff_biased)[b][ntoy]);
      }
    }

    //lines for nominal values
    TLine *unbiased = new TLine();
    unbiased->SetLineColor(kBlack);
    unbiased->SetLineWidth(3);
    unbiased->SetLineStyle(7);
    TH1F *dum = new TH1F("line_unbiased_bin"+(TString)to_string(b),"line_unbiased",10,0,1);
    dum->SetLineColor(kBlack);
    dum->SetLineWidth(3);
    dum->SetLineStyle(7);

    vector<Color_t> cols = {kRed,kBlue,kGreen+2,kMagenta,kOrange};
    vector<TLine> biasedNom(_biasNmeth, TLine());
    for(int m=0;m<_biasNmeth;m++){
      biasedNom[m].SetLineColor(cols[m]);
      biasedNom[m].SetLineWidth(3);
      biasedNom[m].SetLineStyle(7);
    }

    //plot
    gStyle->SetOptStat(0);
    TLegend *leg = new TLegend(0.62,0.67,0.9,0.89);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);

    TCanvas *c1 = new TCanvas("c1_"+(TString)to_string(b),"c1",3200,1200);
    c1->Divide(3,1);

    c1->cd(1)->SetLeftMargin(0.09);
    c1->cd(1)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      acc_toys[m][b-1]->SetLineColor(cols[m]);
      acc_toys[m][b-1]->SetLineWidth(2);
      acc_toys[m][b-1]->GetYaxis()->SetRangeUser(0,1.1*acc_toys[m][b-1]->GetMaximum());
      acc_toys[m][b-1]->GetYaxis()->SetTitleOffset(1.15);
      acc_toys[m][b-1]->Draw((TString)((m==0)?"hist":"histsame"));
      leg->AddEntry(acc_toys[m][b-1], _biasMethName[m]);
      biasedNom[m].DrawLine((*acc_biased)[col][b][m],0, (*acc_biased)[col][b][m], 0.4*acc_toys[0][b-1]->GetMaximum()); //positioned after SetRangeUser
    }
    unbiased->DrawLine(typacc,0, typacc, 0.4*acc_toys[0][b-1]->GetMaximum());
    leg->AddEntry(dum, "unbiased MC");
    leg->Draw("same");

    c1->cd(2)->SetLeftMargin(0.09);
    c1->cd(2)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      eff_toys[m][b-1]->SetLineColor(cols[m]);
      eff_toys[m][b-1]->SetLineWidth(2);
      eff_toys[m][b-1]->GetYaxis()->SetRangeUser(0,1.1*eff_toys[m][b-1]->GetMaximum());
      eff_toys[m][b-1]->GetYaxis()->SetTitleOffset(1.15);
      eff_toys[m][b-1]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine((*eff_biased)[b][m],0, (*eff_biased)[b][m], 0.4*eff_toys[0][b-1]->GetMaximum());
    }
    unbiased->DrawLine(typeff,0, typeff, 0.4*eff_toys[0][b-1]->GetMaximum());

    c1->cd(3)->SetLeftMargin(0.09);
    c1->cd(3)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      acceff_toys[m][b-1]->SetLineColor(cols[m]);
      acceff_toys[m][b-1]->SetLineWidth(2);
      acceff_toys[m][b-1]->GetYaxis()->SetRangeUser(0,1.1*acceff_toys[m][b-1]->GetMaximum());
      acceff_toys[m][b-1]->GetYaxis()->SetTitleOffset(1.15);
      acceff_toys[m][b-1]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine((*acc_biased)[col][b][m] * (*eff_biased)[b][m],0, (*acc_biased)[col][b][m] * (*eff_biased)[b][m], 0.4*acceff_toys[0][b-1]->GetMaximum());
    }
    unbiased->DrawLine(typacceff,0, typacceff, 0.4*acceff_toys[0][b-1]->GetMaximum());
    c1->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)(ispp?"pp":"PbPb")+"_bin"+(TString)to_string(b)+".pdf");

    cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" unbiased one-binned AccEff = "<<typacceff<<endl;
    for(int m=0;m<_biasNmeth;m++)
      cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" mean, rms, and rms/mean of AccEff for toys of biasing method "<<_biasMethName[m]<<" = "<<acceff_toys[m][b-1]->GetMean()<<" "<<acceff_toys[m][b-1]->GetRMS()<<" "<<acceff_toys[m][b-1]->GetRMS() / acceff_toys[m][b-1]->GetMean()<<endl;

  }


}
