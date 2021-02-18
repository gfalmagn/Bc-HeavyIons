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
#include "../helpers/Tools.h"

void plotToyAccEff(bool ispp=true, bool secondStep=true){

  //Grab acc and eff for all toys
  vector<vector<vector<float> > > *acc_biased; //(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) ))
  TFile *facc = new TFile("../acceptance/acceptanceMap.root","READ");
  facc->GetObject("acceptance_oneBinned_biased"+(TString)(secondStep?"_2ndStep":""), acc_biased); 

  //eff
  vector<vector<float> > *eff_biased;
  TFile *feff = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
  feff->GetObject("efficiency_oneBinned_biased"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?"_2ndStep":""), eff_biased); 
  //names of "fit" methods of pT distributions
  

  //Grab nominal acc and eff, without second-step pT biasing of MC
  vector<vector<float> > *acc_oneBinned;
  facc->GetObject("acceptance_oneBinned", acc_oneBinned);
  vector<vector<float> > *eff_oneBinned;
  feff->GetObject("efficiency_oneBinned"+(TString)(ispp?"_pp":"_PbPb"), eff_oneBinned);
  vector<vector<float> > *acc_oneBinned_2ndStep;
  vector<vector<float> > *eff_oneBinned_2ndStep;
  if(secondStep){
    facc->GetObject("acceptance_oneBinned"+(TString)(secondStep?"_2ndStep":""), acc_oneBinned_2ndStep);
    feff->GetObject("efficiency_oneBinned"+(TString)(ispp?"_pp":"_PbPb")+(TString)(secondStep?"_2ndStep":""), eff_oneBinned_2ndStep);
  }

  //histograms of acc and eff filled with toys
  vector<vector<TH1F*> > accI_toys(_biasNmeth);
  vector<vector<TH1F*> > effI_toys(_biasNmeth);
  vector<vector<TH1F*> > acceffI_toys(_biasNmeth);
  
  //output
  vector<vector<float> > acceffI_biased(_NanaBins);
  vector<float> DevToMean(_NanaBins);
  vector<vector<double> > nomiAndErr(_NanaBins);

  int col = ispp?0:1;
  for(int b=1;b<=_NanaBins;b++){
    
    //histograms definition
    float accI_unbiased = 1/(*acc_oneBinned)[1-(int)ispp][b];
    float effI_unbiased = 1/(*eff_oneBinned)[b][0];
    float acceffI_unbiased = accI_unbiased*effI_unbiased;
    float accI_1stStep = secondStep?(1/(*acc_oneBinned_2ndStep)[1-(int)ispp][b]):1;
    float effI_1stStep = secondStep?(1/(*eff_oneBinned_2ndStep)[b][0]):1;
    float acceffI_1stStep = accI_1stStep*effI_1stStep;
    float maxacc = MaxVec((*acc_biased)[col][b], 1/accI_unbiased, ispp?1e5:(_biasNmeth+2*_biasNtoys));//exclude outliers (especially 4) from third method
    float minacc = MinVec((*acc_biased)[col][b], 1/accI_unbiased, ispp?(_biasNmeth+1*_biasNtoys):1e5);
    float maxeff = MaxVec((*eff_biased)[b], 1/effI_unbiased, ispp?1e5:(_biasNmeth+2*_biasNtoys));
    float mineff = MinVec((*eff_biased)[b], 1/effI_unbiased, ispp?(_biasNmeth+1*_biasNtoys):1e5);
    float maxacceff = maxacc*maxeff;
    float minacceff = minacc*mineff;
    for(int m=0;m<_biasNmeth;m++) acceffI_biased[b-1].push_back(1/((*acc_biased)[col][b][m] * (*eff_biased)[b][m]));

    for(int m=0;m<_biasNmeth;m++){

      accI_toys[m].push_back(new TH1F("acc_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "acceptance toys "+(TString)(ispp?"pp":"PbPb")+" bin"+(TString)to_string(b)+";1/acceptance;n_{toys}", 23, 1/maxacc-0.05*(1/minacc - 1/maxacc) , 1/minacc+0.1*(1/minacc - 1/maxacc)));
      effI_toys[m].push_back(new TH1F("eff_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "efficiency toys "+(TString)(ispp?"pp":"PbPb")+" bin"+(TString)to_string(b)+";1/efficiency;n_{toys}", 23, 1/maxeff-0.05*(1/mineff - 1/maxeff) , 1/mineff+0.1*(1/mineff - 1/maxeff)));
      acceffI_toys[m].push_back(new TH1F("acceff_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "acc #times eff toys "+(TString)(ispp?"pp":"PbPb")+" bin"+(TString)to_string(b)+";1/(acc #times eff) ;n_{toys}", 23, 1/maxacceff-0.05*(1/minacceff - 1/maxacceff) , 1/minacceff+0.1*(1/minacceff - 1/maxacceff)));

      //Run over all toys
      for(int t=0;t<_biasNtoys;t++){
	int ntoy = _biasNmeth + m*_biasNtoys + t;
	//cout<<"b toyNumber acc eff = "<<b<<" "<<t<<" "<< (*acc_biased)[col][b][ntoy] <<" "<<(*eff_biased)[b][ntoy]<<endl;
	accI_toys[m][b-1]->Fill(1/((*acc_biased)[col][b][ntoy]));
	effI_toys[m][b-1]->Fill(1/((*eff_biased)[b][ntoy]));
	acceffI_biased[b-1].push_back(1/((*acc_biased)[col][b][ntoy] * (*eff_biased)[b][ntoy]));
	acceffI_toys[m][b-1]->Fill(acceffI_biased[b-1][ntoy]);
      }
    }

    //lines for nominal values
    TLine *unbiased = new TLine();
    unbiased->SetLineColor(kBlack);
    unbiased->SetLineWidth(3);
    unbiased->SetLineStyle(7);
    TLine *corr1stStep = new TLine();
    corr1stStep->SetLineColor(kMagenta+3);
    corr1stStep->SetLineWidth(3);
    corr1stStep->SetLineStyle(3);
    TH1F *dum = new TH1F("line_unbiased_bin"+(TString)to_string(b),"line_unbiased",10,0,1);
    dum->SetLineColor(kBlack);
    dum->SetLineWidth(3);
    dum->SetLineStyle(7);
    TH1F *dum2 = new TH1F("line_1stStep_bin"+(TString)to_string(b),"line_1stStep",10,0,1);
    dum2->SetLineColor(kMagenta+3);
    dum2->SetLineWidth(3);
    dum2->SetLineStyle(3);

    vector<Color_t> cols = {kRed,kBlue,kGreen+2,kMagenta,kOrange};
    vector<TLine> biasedNom(_biasNmeth, TLine());
    for(int m=0;m<_biasNmeth;m++){
      biasedNom[m].SetLineColor(cols[m]);
      biasedNom[m].SetLineWidth(3);
      biasedNom[m].SetLineStyle(7);
    }

    //plot
    gStyle->SetOptStat(0);
    TLegend *leg = new TLegend(0.12,0.62,0.45,0.89);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);

    TCanvas *c1 = new TCanvas("c1_"+(TString)to_string(b),"c1",3200,1200);
    c1->Divide(3,1);

    c1->cd(1)->SetLeftMargin(0.09);
    c1->cd(1)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      accI_toys[m][b-1]->SetLineColor(cols[m]);
      accI_toys[m][b-1]->SetLineWidth(2);
      accI_toys[m][b-1]->GetYaxis()->SetRangeUser(0,1.15*accI_toys[m][b-1]->GetMaximum());
      accI_toys[m][b-1]->GetYaxis()->SetTitleOffset(1.15);
      accI_toys[m][b-1]->Draw((TString)((m==0)?"hist":"histsame"));
      leg->AddEntry(accI_toys[m][b-1], _biasMethName[m]);
      biasedNom[m].DrawLine(1/(*acc_biased)[col][b][m],0, 1/(*acc_biased)[col][b][m], 0.4*accI_toys[0][b-1]->GetMaximum()); //positioned after SetRangeUser
    }
    unbiased->DrawLine(accI_unbiased,0, accI_unbiased, 0.4*accI_toys[0][b-1]->GetMaximum());
    corr1stStep->DrawLine(accI_1stStep,0, accI_1stStep, 0.4*accI_toys[0][b-1]->GetMaximum());
    leg->AddEntry(dum, "unbiased MC");
    leg->AddEntry(dum2, "1^{st}-step corrected");

    c1->cd(2)->SetLeftMargin(0.09);
    c1->cd(2)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      effI_toys[m][b-1]->SetLineColor(cols[m]);
      effI_toys[m][b-1]->SetLineWidth(2);
      effI_toys[m][b-1]->GetYaxis()->SetRangeUser(0,1.15*effI_toys[m][b-1]->GetMaximum());
      effI_toys[m][b-1]->GetYaxis()->SetTitleOffset(1.15);
      effI_toys[m][b-1]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(1/(*eff_biased)[b][m],0, 1/(*eff_biased)[b][m], 0.4*effI_toys[0][b-1]->GetMaximum());
    }
    unbiased->DrawLine(effI_unbiased,0, effI_unbiased, 0.4*effI_toys[0][b-1]->GetMaximum());
    corr1stStep->DrawLine(effI_1stStep,0, effI_1stStep, 0.4*effI_toys[0][b-1]->GetMaximum());

    c1->cd(3)->SetLeftMargin(0.09);
    c1->cd(3)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      acceffI_toys[m][b-1]->SetLineColor(cols[m]);
      acceffI_toys[m][b-1]->SetLineWidth(2);
      acceffI_toys[m][b-1]->GetYaxis()->SetRangeUser(0,1.15*acceffI_toys[m][b-1]->GetMaximum());
      acceffI_toys[m][b-1]->GetYaxis()->SetTitleOffset(1.15);
      acceffI_toys[m][b-1]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(acceffI_biased[b-1][m],0, acceffI_biased[b-1][m], 0.4*acceffI_toys[0][b-1]->GetMaximum());
    }
    unbiased->DrawLine(acceffI_unbiased,0, acceffI_unbiased, 0.4*acceffI_toys[0][b-1]->GetMaximum());
    corr1stStep->DrawLine(acceffI_1stStep,0, acceffI_1stStep, 0.4*acceffI_toys[0][b-1]->GetMaximum());
    leg->Draw("same");

    c1->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+"_bin"+(TString)to_string(b)+".pdf");

    TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
    c3->SetLeftMargin(0.09);
    c3->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      acceffI_toys[m][b-1]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(acceffI_biased[b-1][m],0, acceffI_biased[b-1][m], 0.4*acceffI_toys[0][b-1]->GetMaximum());
    }
    unbiased->DrawLine(acceffI_unbiased,0, acceffI_unbiased, 0.4*acceffI_toys[0][b-1]->GetMaximum());
    corr1stStep->DrawLine(acceffI_1stStep,0, acceffI_1stStep, 0.4*acceffI_toys[0][b-1]->GetMaximum());
    leg->Draw("same");

    c3->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+"_bin"+(TString)to_string(b)+"_onlyAccEffProduct.pdf");

    //Get mean of the three possible nominal values, and low and high RMS from this mean
    //mean
    double mean = SimpleMean(acceffI_biased[b-1],0, _biasNmeth-1);
    //Deviation of the three nominals to the mean
    double DevToMeanLo = 0, DevToMeanHi = 0;
    for(int m=0;m<_biasNmeth;m++){
      if(acceffI_biased[b-1][m] < mean - DevToMeanLo) DevToMeanLo = fabs(acceffI_biased[b-1][m]-mean);
      if(acceffI_biased[b-1][m] > mean + DevToMeanHi) DevToMeanHi = fabs(acceffI_biased[b-1][m]-mean);
    }
    //The RMS's for the three methods, from their respective nominal
    vector<vector<double> > meanAndRMS; 
    for(int m=0;m<_biasNmeth;m++) meanAndRMS.push_back( DoubleSidedRMS(acceffI_biased[b-1], 
								       _biasNmeth+m*_biasNtoys, _biasNmeth+(m+1)*_biasNtoys-1,
								       acceffI_biased[b-1][m]) );
    //Find max RMS of the three methods
    meanAndRMS.push_back(vector<double>{mean,0,0,0});
    for(int m=0;m<_biasNmeth;m++){
      if(meanAndRMS[m][1] > meanAndRMS[_biasNmeth][1]) meanAndRMS[_biasNmeth][1] = meanAndRMS[m][1];//total
      if(meanAndRMS[m][2] > meanAndRMS[_biasNmeth][2]) meanAndRMS[_biasNmeth][2] = meanAndRMS[m][2];//lo
      if(meanAndRMS[m][3] > meanAndRMS[_biasNmeth][3]) meanAndRMS[_biasNmeth][3] = meanAndRMS[m][3];//hi
    }

    //Write output
    cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" unbiased one-binned AccEffInverse = "<<acceffI_unbiased<<endl;
    cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" 1st-step-corrected one-binned AccEffInverse = "<<acceffI_1stStep<<endl;
    for(int m=0;m<_biasNmeth;m++)
      cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" mean, lo rms/mean, up rms/mean, relative_difference_with_unbiased of AccEffInverse for toys of biasing method "<<_biasMethName[m]<<" = "<<meanAndRMS[m][0]<<" "<<meanAndRMS[m][2]/meanAndRMS[m][0]<<" "<<meanAndRMS[m][3]/meanAndRMS[m][0]<<" "<<(meanAndRMS[m][0] - acceffI_unbiased)/acceffI_unbiased <<endl;
    cout<<"Mean of the three nominals (final nominal), max rel Up deviation, max rel Lo deviation = "<<meanAndRMS[_biasNmeth][0]<<" "<<DevToMeanHi/meanAndRMS[_biasNmeth][0]<<" "<<DevToMeanLo/meanAndRMS[_biasNmeth][0]<<endl;
    cout<<"Maximal (among three methods) rel RMS total, lo, hi = "<<meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0]<<endl;

    //Add the deviation of the three nominal to their means into the final syst error
    DevToMean[b-1] = (DevToMeanLo+DevToMeanHi)/2;
    meanAndRMS[_biasNmeth][1] = sqrt(pow(meanAndRMS[_biasNmeth][1],2) + pow(DevToMean[b-1],2));
    meanAndRMS[_biasNmeth][2] = sqrt(pow(meanAndRMS[_biasNmeth][2],2) + pow(DevToMeanLo,2));
    meanAndRMS[_biasNmeth][3] = sqrt(pow(meanAndRMS[_biasNmeth][3],2) + pow(DevToMeanHi,2));
    cout<<"Final acceffInverse rel uncertainty total, lo, hi = "<<meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0]<<endl;

    nomiAndErr[b-1] = meanAndRMS[_biasNmeth];
  }

  //Calculate ana bins correlation factor (result>0.999)
  //Will assume 0.3 correlation for subdominant DevToMean 
  vector<float> corrFactors;
  for(int m=0;m<_biasNmeth;m++){
    for(int b=0;b<_NanaBins;b++){
      for(int b2=b+1;b2<_NanaBins;b2++){
	float sumcor = 0, sum1=0, sum2=0;
	for(int i=_biasNmeth+m*_biasNtoys;i<_biasNmeth+(m+1)*_biasNtoys;i++){
	  sumcor += (acceffI_biased[b][i]-acceffI_biased[b][m]) * (acceffI_biased[b2][i]-acceffI_biased[b2][m]);
	  sum1 += pow((acceffI_biased[b][i]-acceffI_biased[b][m]),2);
	  sum2 += pow((acceffI_biased[b2][i]-acceffI_biased[b2][m]),2);
	}
	corrFactors.push_back(sumcor/sqrt(sum1*sum2));
      }
    }
  }

  float corrDevToMean = 0.3; //arbitrary but subdominant
  float corrfact = SimpleMean(corrFactors,0,corrFactors.size()-1);
  cout<<"correlation factors from the three methods (rms only), and mean = "<<corrFactors[0]<<" "<<corrFactors[1]<<" "<<corrFactors[2]<<" "<<corrfact<<endl;
  vector<float> corrfactFinal;
  corrfactFinal.push_back((corrDevToMean*DevToMean[0]*DevToMean[1] 
			   + corrfact * sqrt(pow(nomiAndErr[0][1],2) - pow(DevToMean[0],2)) * sqrt(pow(nomiAndErr[1][1],2) - pow(DevToMean[1],2))
			   )/ (nomiAndErr[0][1]*nomiAndErr[1][1]) );
  cout<<"final correlation factor = "<<corrfactFinal[0]<<endl;

  //Record the final correction 
  TFile *fout = new TFile("AccEffFrom2ndStepToys.root","UPDATE");  
  fout->WriteObject(&nomiAndErr,"InvAccEffFromCorrMC_withSystErr_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":""));
  fout->WriteObject(&corrfactFinal,"InvAccEffFromCorrMC_LinearisedCorrelationFactor_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":""));

}
