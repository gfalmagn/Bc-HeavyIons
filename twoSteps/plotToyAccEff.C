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
#include "../helpers/hub.cpp"

void ToyAccEff(Hub H, bool ispp=true, bool secondStep=true, bool inCentBins=false, bool keepfile=true){
  
  if(ispp) inCentBins=false;

  //for error on yields, without acceff
  Result ResAll;
  ResAll = ispp?SumResult({H.fit.pp_pt,H.metafit.pp_pt,H.TnP.pp_pt} , H.Ycorr.pp_pt.Val,true):(
     inCentBins?SumResult({H.fit.PbPb_cent,H.metafit.PbPb_cent,H.TnP.PbPb_cent} , H.Ycorr.PbPb_cent.Val,true):(
		SumResult({H.fit.PbPb_pt,H.metafit.PbPb_pt,H.TnP.PbPb_pt} , H.Ycorr.PbPb_pt.Val,true)));

  //eff
  vector<vector<vector<float> > > eff_biased = inCentBins?H.eff_biased_cent:H.eff_biased;
  vector<vector<float> > eff_oneBinned = inCentBins?H.eff_oneBinned_cent:H.eff_oneBinned;
  vector<vector<float> > eff_oneBinned_2ndStep = inCentBins?H.eff_oneBinned_cent_2ndStep:H.eff_oneBinned_2ndStep;

  //nominal corrected yields
  vector<float> ycorr = ispp?H.Ycorr.pp_pt.Val:(inCentBins?H.Ycorr.PbPb_cent.Val:H.Ycorr.PbPb_pt.Val); //if we take the simple nominal method and not an average of methods (cf Draw_metafit), this is supposed to = yield_postfit / (acc_oneBinned*eff_oneBinned), with _2ndStep in all component if(secondStep)

  //histograms of acc and eff filled with toys
  vector<vector<TH1F*> > accI_toys(_biasNmeth);
  vector<vector<TH1F*> > effI_toys(_biasNmeth);
  vector<vector<TH1F*> > acceffI_toys(_biasNmeth);
  vector<vector<TH1F*> > ycorr_toys(_biasNmeth);
  
  //output
  vector<vector<float> > acceffI_biased(_NanaBins+1);
  vector<float> DevToMean(_NanaBins+1);
  vector<float> DevToMean_y(_NanaBins+1);
  vector<vector<double> > nomiAndErr(_NanaBins+1);
  vector<vector<double> > nomiAndErr_y(_NanaBins+1);
  vector<vector<float> > ycorr_biased(_NanaBins+1);

  int col = ispp?0:1;
  for(int b=0;b<=_NanaBins;b++){

    cout<<"\n***** Bin "<<(TString)((b==0)?"integrated":to_string(b))<<endl<<endl;

    //AccAndEff values
    int bacc = inCentBins?0:b; //bin number for acceptance (integrated when doing centrality dependence)
    float accI_unbiased = 1/H.acc_oneBinned[1-(int)ispp][bacc];
    float effI_unbiased = 1/eff_oneBinned[1-(int)ispp][b];
    float acceffI_unbiased = accI_unbiased*effI_unbiased;
    float accI_1stStep = secondStep?(1/H.acc_oneBinned_2ndStep[1-(int)ispp][bacc]):1;
    float effI_1stStep = secondStep?(1/H.eff_oneBinned_2ndStep[1-(int)ispp][b]):1;
    float acceffI_1stStep = accI_1stStep*effI_1stStep;
    float maxacc = MaxVec(H.acc_biased[col][bacc], 1/accI_unbiased, ispp?1e5:(_biasNmeth+2*_biasNtoys));//exclude outliers (especially 4) from third method
    float minacc = MinVec(H.acc_biased[col][bacc], 1/accI_unbiased, ispp?(_biasNmeth+1*_biasNtoys):1e5);
    float maxeff = MaxVec(eff_biased[col][b], 1/effI_unbiased, ispp?1e5:(_biasNmeth+2*_biasNtoys));
    float mineff = MinVec(eff_biased[col][b], 1/effI_unbiased, ispp?(_biasNmeth+1*_biasNtoys):1e5);
    float maxacceff = maxacc*maxeff;
    float minacceff = minacc*mineff;
    float acceffI_default = secondStep?acceffI_1stStep:acceffI_unbiased ;

    //corrected yields values
    for(int t=0;t<_biasNmeth*(_biasNtoys+1);t++){
      //yield-ratio-variation * nominal-yield * cancel-nominal-AccEff / varied-AccEff
      ycorr_biased[b].push_back( H.corrYr_varied[col][b][t] * ycorr[b] / acceffI_default / (H.acc_biased[col][bacc][t] * eff_biased[col][b][t]) );
    }
    float maxycorr = MaxVec(ycorr_biased[b], ycorr[b], 1e5);
    float minycorr = MinVec(ycorr_biased[b], ycorr[b], 1e5);
    float ycorr_unbiased = ycorr[b] * acceffI_unbiased / acceffI_default;
    float ycorr_1stStep = ycorr[b] * acceffI_1stStep / acceffI_default;
    
    for(int m=0;m<_biasNmeth;m++) acceffI_biased[b].push_back(1/(H.acc_biased[col][bacc][m] * eff_biased[col][b][m]));

    //histograms definition
    for(int m=0;m<_biasNmeth;m++){
      accI_toys[m].push_back(new TH1F("acc_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "acceptance toys "+(TString)(ispp?"pp ":"PbPb ")+(TString)((b==0 || inCentBins)?"integrated":("p_{T} bin"+(TString)to_string(b)))+";1/acceptance;n_{toys}", 23, 1/maxacc-0.05*(1/minacc - 1/maxacc) , 1/minacc+0.1*(1/minacc - 1/maxacc)));
      effI_toys[m].push_back(new TH1F("eff_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "efficiency toys "+(TString)(ispp?"pp ":"PbPb ")+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";1/efficiency;n_{toys}", 23, 1/maxeff-0.05*(1/mineff - 1/maxeff) , 1/mineff+0.1*(1/mineff - 1/maxeff)));
      acceffI_toys[m].push_back(new TH1F("acceff_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "acc #times eff toys "+(TString)(ispp?"pp ":"PbPb ")+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";1/(acc #times eff) ;n_{toys}", 23, 1/maxacceff-0.05*(1/minacceff - 1/maxacceff) , 1/minacceff+0.1*(1/minacceff - 1/maxacceff)));
      ycorr_toys[m].push_back(new TH1F("ycorr_toys"+(TString)(ispp?"_pp":"_PbPb")+"_bin"+(TString)to_string(b)+"_meth"+(TString)to_string(m), "acc #times eff toys "+(TString)(ispp?"pp ":"PbPb ")+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";N_{corrected} ;n_{toys}", 23, minycorr-0.05*(maxycorr - minycorr) , maxycorr+0.1*(maxycorr - minycorr)));

      //Run over all toys
      for(int t=0;t<_biasNtoys;t++){
	int ntoy = _biasNmeth + m*_biasNtoys + t;
	//cout<<"b toyNumber acc eff = "<<b<<" "<<t<<" "<< H.acc_biased[col][bacc][ntoy] <<" "<<eff_biased[col][b][ntoy]<<endl;
	accI_toys[m][b]->Fill(1/(H.acc_biased[col][bacc][ntoy]));
	effI_toys[m][b]->Fill(1/(eff_biased[col][b][ntoy]));
	acceffI_biased[b].push_back(1/(H.acc_biased[col][bacc][ntoy] * eff_biased[col][b][ntoy]));
	acceffI_toys[m][b]->Fill(acceffI_biased[b][ntoy]);
	ycorr_toys[m][b]->Fill(ycorr_biased[b][ntoy]);
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

    TCanvas *c1 = new TCanvas("c1_"+(TString)to_string(b),"c1",3500,1200);
    c1->Divide(4,1);

    c1->cd(1)->SetLeftMargin(0.09);
    c1->cd(1)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      accI_toys[m][b]->SetLineColor(cols[m]);
      accI_toys[m][b]->SetLineWidth(2);
      accI_toys[m][b]->GetYaxis()->SetRangeUser(0,1.15*accI_toys[m][b]->GetMaximum());
      accI_toys[m][b]->GetYaxis()->SetTitleOffset(1.15);
      accI_toys[m][b]->Draw((TString)((m==0)?"hist":"histsame"));
      leg->AddEntry(accI_toys[m][b], _biasMethName[m]);
      biasedNom[m].DrawLine(1/H.acc_biased[col][bacc][m],0, 1/H.acc_biased[col][bacc][m], 0.4*accI_toys[0][b]->GetMaximum()); //positioned after SetRangeUser
    }
    unbiased->DrawLine(accI_unbiased,0, accI_unbiased, 0.4*accI_toys[0][b]->GetMaximum());
    corr1stStep->DrawLine(accI_1stStep,0, accI_1stStep, 0.4*accI_toys[0][b]->GetMaximum());
    leg->AddEntry(dum, "unbiased MC");
    leg->AddEntry(dum2, "1^{st}-step corrected");

    c1->cd(2)->SetLeftMargin(0.09);
    c1->cd(2)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      effI_toys[m][b]->SetLineColor(cols[m]);
      effI_toys[m][b]->SetLineWidth(2);
      effI_toys[m][b]->GetYaxis()->SetRangeUser(0,1.15*effI_toys[m][b]->GetMaximum());
      effI_toys[m][b]->GetYaxis()->SetTitleOffset(1.15);
      effI_toys[m][b]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(1/eff_biased[col][b][m],0, 1/eff_biased[col][b][m], 0.4*effI_toys[0][b]->GetMaximum());
    }
    unbiased->DrawLine(effI_unbiased,0, effI_unbiased, 0.4*effI_toys[0][b]->GetMaximum());
    corr1stStep->DrawLine(effI_1stStep,0, effI_1stStep, 0.4*effI_toys[0][b]->GetMaximum());

    c1->cd(3)->SetLeftMargin(0.09);
    c1->cd(3)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      acceffI_toys[m][b]->SetLineColor(cols[m]);
      acceffI_toys[m][b]->SetLineWidth(2);
      acceffI_toys[m][b]->GetYaxis()->SetRangeUser(0,1.15*acceffI_toys[m][b]->GetMaximum());
      acceffI_toys[m][b]->GetYaxis()->SetTitleOffset(1.15);
      acceffI_toys[m][b]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(acceffI_biased[b][m],0, acceffI_biased[b][m], 0.4*acceffI_toys[0][b]->GetMaximum());
    }
    unbiased->DrawLine(acceffI_unbiased,0, acceffI_unbiased, 0.4*acceffI_toys[0][b]->GetMaximum());
    corr1stStep->DrawLine(acceffI_1stStep,0, acceffI_1stStep, 0.4*acceffI_toys[0][b]->GetMaximum());
    leg->Draw("same");

    c1->cd(4)->SetLeftMargin(0.09);
    c1->cd(4)->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      ycorr_toys[m][b]->SetLineColor(cols[m]);
      ycorr_toys[m][b]->SetLineWidth(2);
      ycorr_toys[m][b]->GetYaxis()->SetRangeUser(0,1.15*ycorr_toys[m][b]->GetMaximum());
      ycorr_toys[m][b]->GetYaxis()->SetTitleOffset(1.15);
      ycorr_toys[m][b]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(ycorr_biased[b][m],0, ycorr_biased[b][m], 0.4*ycorr_toys[0][b]->GetMaximum());
    }
    unbiased->DrawLine(ycorr_unbiased, 0, ycorr_unbiased, 0.4*ycorr_toys[0][b]->GetMaximum());
    corr1stStep->DrawLine(ycorr_1stStep,0, ycorr_1stStep, 0.4*ycorr_toys[0][b]->GetMaximum());
    leg->Draw("same");

    c1->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+".pdf");

    TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
    c3->SetLeftMargin(0.09);
    c3->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      acceffI_toys[m][b]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(acceffI_biased[b][m],0, acceffI_biased[b][m], 0.4*acceffI_toys[0][b]->GetMaximum());
    }
    unbiased->DrawLine(acceffI_unbiased,0, acceffI_unbiased, 0.4*acceffI_toys[0][b]->GetMaximum());
    corr1stStep->DrawLine(acceffI_1stStep,0, acceffI_1stStep, 0.4*acceffI_toys[0][b]->GetMaximum());
    leg->Draw("same");

    c3->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+"_onlyAccEffProduct.pdf");

    TCanvas *c4 = new TCanvas("c4","c4",1500,1500);
    c4->SetLeftMargin(0.09);
    c4->SetRightMargin(0.02);
    for(int m=0;m<_biasNmeth;m++){
      ycorr_toys[m][b]->Draw((TString)((m==0)?"hist":"histsame"));
      biasedNom[m].DrawLine(ycorr_biased[b][m],0, ycorr_biased[b][m], 0.4*ycorr_toys[0][b]->GetMaximum());
    }
    unbiased->DrawLine(ycorr_unbiased, 0, ycorr_unbiased, 0.4*ycorr_toys[0][b]->GetMaximum());
    corr1stStep->DrawLine(ycorr_1stStep,0, ycorr_1stStep, 0.4*ycorr_toys[0][b]->GetMaximum());
    leg->Draw("same");

    c4->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)(ispp?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+"_onlyCorrYields.pdf");

    //Get mean of the three possible nominal values, and low and high RMS from this mean
    //mean
    double mean = SimpleMean(acceffI_biased[b],0, _biasNmeth-1);
    double mean_y = SimpleMean(ycorr_biased[b],0, _biasNmeth-1);//ycorr[b] * mean / acceffI_default;
    //Deviation of the three nominals to the mean
    double DevToMeanLo = 0, DevToMeanHi = 0;
    for(int m=0;m<_biasNmeth;m++){
      if(acceffI_biased[b][m] < mean - DevToMeanLo) DevToMeanLo = fabs(acceffI_biased[b][m]-mean);
      if(acceffI_biased[b][m] > mean + DevToMeanHi) DevToMeanHi = fabs(acceffI_biased[b][m]-mean);
    }
    float DevToMeanLo_y = DevToMeanLo * ycorr[b] / acceffI_default;
    float DevToMeanHi_y = DevToMeanHi * ycorr[b] / acceffI_default;
    //The RMS's for the three methods, from their respective nominal
    vector<vector<double> > meanAndRMS,meanAndRMS_y; 
    for(int m=0;m<_biasNmeth;m++) {
      meanAndRMS.push_back( DoubleSidedRMS(acceffI_biased[b], 
					   _biasNmeth+m*_biasNtoys, _biasNmeth+(m+1)*_biasNtoys-1,
					   acceffI_biased[b][m]) );
      meanAndRMS_y.push_back( DoubleSidedRMS(ycorr_biased[b], 
					     _biasNmeth+m*_biasNtoys, _biasNmeth+(m+1)*_biasNtoys-1,
					     ycorr_biased[b][m]) );
    }

    //Find max RMS of the three methods
    meanAndRMS.push_back(vector<double>{mean,0,0,0});
    meanAndRMS_y.push_back(vector<double>{mean_y,0,0,0});
    for(int m=0;m<_biasNmeth;m++){
      if(meanAndRMS[m][1] > meanAndRMS[_biasNmeth][1]) meanAndRMS[_biasNmeth][1] = meanAndRMS[m][1];//total
      if(meanAndRMS[m][2] > meanAndRMS[_biasNmeth][2]) meanAndRMS[_biasNmeth][2] = meanAndRMS[m][2];//lo
      if(meanAndRMS[m][3] > meanAndRMS[_biasNmeth][3]) meanAndRMS[_biasNmeth][3] = meanAndRMS[m][3];//hi
      if(meanAndRMS_y[m][1] > meanAndRMS_y[_biasNmeth][1]) meanAndRMS_y[_biasNmeth][1] = meanAndRMS_y[m][1];//total
      if(meanAndRMS_y[m][2] > meanAndRMS_y[_biasNmeth][2]) meanAndRMS_y[_biasNmeth][2] = meanAndRMS_y[m][2];//lo
      if(meanAndRMS_y[m][3] > meanAndRMS_y[_biasNmeth][3]) meanAndRMS_y[_biasNmeth][3] = meanAndRMS_y[m][3];//hi
    }

    //Write output
    cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" unbiased one-binned AccEffInverse = "<<acceffI_unbiased<<endl;
    cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" 1st-step-corrected one-binned AccEffInverse = "<<acceffI_1stStep<<endl;
    for(int m=0;m<_biasNmeth;m++)
      cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" mean, lo rms/mean, up rms/mean, relative_difference_with_unbiased of AccEffInverse for toys of biasing method "<<_biasMethName[m]<<" = "<<meanAndRMS[m][0]<<" "<<meanAndRMS[m][2]/meanAndRMS[m][0]<<" "<<meanAndRMS[m][3]/meanAndRMS[m][0]<<" "<<(meanAndRMS[m][0] - acceffI_unbiased)/acceffI_unbiased <<endl;
    cout<<"Mean of the three nominals (final nominal), max rel Up deviation, max rel Lo deviation = "<<meanAndRMS[_biasNmeth][0]<<" "<<DevToMeanHi/meanAndRMS[_biasNmeth][0]<<" "<<DevToMeanLo/meanAndRMS[_biasNmeth][0]<<endl;
    cout<<"Maximal (among three methods) rel RMS total, lo, hi = "<<meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0]<<endl;

    cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" ycorr from unbiased one-binned AccEff = "<<ycorr_unbiased<<endl;
    cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" ycorr from 1st-step-corrected one-binned AccEff = "<<ycorr_1stStep<<endl;
    for(int m=0;m<_biasNmeth;m++)
      cout<<(TString)(ispp?"pp ":"PbPb ")<<"bin "<<b<<" mean, lo rms/mean, up rms/mean, relative_difference_with_unbiased of ycorr for toys of biasing method "<<_biasMethName[m]<<" = "<<meanAndRMS_y[m][0]<<" "<<meanAndRMS_y[m][2]/meanAndRMS_y[m][0]<<" "<<meanAndRMS_y[m][3]/meanAndRMS_y[m][0]<<" "<<(meanAndRMS_y[m][0] - ycorr_unbiased)/ycorr_unbiased <<endl;
    cout<<"Mean of the three nominals (final nominal), max rel Up deviation, max rel Lo deviation = "<<meanAndRMS_y[_biasNmeth][0]<<" "<<DevToMeanHi_y/meanAndRMS_y[_biasNmeth][0]<<" "<<DevToMeanLo_y/meanAndRMS_y[_biasNmeth][0]<<endl;
    cout<<"Maximal (among three methods) rel RMS total, lo, hi = "<<meanAndRMS_y[_biasNmeth][1]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][2]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][3]/meanAndRMS_y[_biasNmeth][0]<<endl;

    //Add the deviation of the three nominal to their means into the final syst error
    DevToMean[b] = (DevToMeanLo+DevToMeanHi)/2;
    meanAndRMS[_biasNmeth][1] = sqrt(pow(meanAndRMS[_biasNmeth][1],2) + pow(DevToMean[b],2));
    meanAndRMS[_biasNmeth][2] = sqrt(pow(meanAndRMS[_biasNmeth][2],2) + pow(DevToMeanLo,2));
    meanAndRMS[_biasNmeth][3] = sqrt(pow(meanAndRMS[_biasNmeth][3],2) + pow(DevToMeanHi,2));
    cout<<"Final acceffInverse rel uncertainty total, lo, hi = "<<meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0]<<endl;
    // cout<<"Fit rel uncertainty total, lo, hi = "<<(ispp?H.fit.pp_pt:(inCentBins?H.fit.PbPb_cent:H.fit.PbPb_pt)).RelErr[b]<<" "<<(ispp?H.fit.pp_pt:(inCentBins?H.fit.PbPb_cent:H.fit.PbPb_pt)).RelErrLo[b]<<" "<<(ispp?H.fit.pp_pt:(inCentBins?H.fit.PbPb_cent:H.fit.PbPb_pt)).RelErrHi[b]<<endl;
    // cout<<"Metafit rel uncertainty total, lo, hi = "<<(ispp?H.metafit.pp_pt:(inCentBins?H.metafit.PbPb_cent:H.metafit.PbPb_pt)).RelErr[b]<<" "<<(ispp?H.metafit.pp_pt:(inCentBins?H.metafit.PbPb_cent:H.metafit.PbPb_pt)).RelErrLo[b]<<" "<<(ispp?H.metafit.pp_pt:(inCentBins?H.metafit.PbPb_cent:H.metafit.PbPb_pt)).RelErrHi[b]<<endl;
    cout<<"Full rel uncertainty without AccEff total, lo, hi = "<<ResAll.RelErr[b]<<" "<<ResAll.RelErrLo[b]<<" "<<ResAll.RelErrHi[b]<<endl;
    cout<<"Simple quadSum of AccEff and others rel uncertainties, total, lo, hi = "<<quadSum(ResAll.RelErr[b],meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0])<<" "<<quadSum(ResAll.RelErrLo[b],meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0])<<" "<<quadSum(ResAll.RelErrHi[b],meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0])<<endl;
    DevToMean_y[b] = (DevToMeanLo_y+DevToMeanHi_y)/2;
    meanAndRMS_y[_biasNmeth][1] = sqrt(pow(meanAndRMS_y[_biasNmeth][1],2) + pow(DevToMean_y[b],2));
    meanAndRMS_y[_biasNmeth][2] = sqrt(pow(meanAndRMS_y[_biasNmeth][2],2) + pow(DevToMeanLo_y,2));
    meanAndRMS_y[_biasNmeth][3] = sqrt(pow(meanAndRMS_y[_biasNmeth][3],2) + pow(DevToMeanHi_y,2));
    cout<<"Final ycorr rel uncertainty total, lo, hi = "<<meanAndRMS_y[_biasNmeth][1]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][2]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][3]/meanAndRMS_y[_biasNmeth][0]<<endl;

    nomiAndErr[b] = meanAndRMS[_biasNmeth];
    nomiAndErr_y[b] = meanAndRMS_y[_biasNmeth];
  }

  //Calculate ana bins correlation factor (result>0.999)
  //Will assume 0.3 correlation for subdominant DevToMean 
  vector<float> corrFactors,corrFactors_y;
  for(int m=0;m<_biasNmeth;m++){
    for(int b=1;b<=_NanaBins;b++){
      for(int b2=b+1;b2<=_NanaBins;b2++){
	float sumcor = 0, sum1=0, sum2=0;
	float sumcor_y = 0, sum1_y=0, sum2_y=0;
	for(int i=_biasNmeth+m*_biasNtoys;i<_biasNmeth+(m+1)*_biasNtoys;i++){
	  sumcor += (acceffI_biased[b][i]-acceffI_biased[b][m]) * (acceffI_biased[b2][i]-acceffI_biased[b2][m]);
	  sum1 += pow((acceffI_biased[b][i]-acceffI_biased[b][m]),2);
	  sum2 += pow((acceffI_biased[b2][i]-acceffI_biased[b2][m]),2);
	  sumcor_y += (ycorr_biased[b][i]-ycorr_biased[b][m]) * (ycorr_biased[b2][i]-ycorr_biased[b2][m]);
	  sum1_y += pow((ycorr_biased[b][i]-ycorr_biased[b][m]),2);
	  sum2_y += pow((ycorr_biased[b2][i]-ycorr_biased[b2][m]),2);
	}
	corrFactors.push_back(sumcor/sqrt(sum1*sum2));
	corrFactors_y.push_back(sumcor_y/sqrt(sum1_y*sum2_y));
      }
    }
  }

  float corrDevToMean = 0.3; //arbitrary but subdominant
  float corrfact = SimpleMean(corrFactors,0,corrFactors.size()-1);
  cout<<"correlation factors from the three methods (rms only), and mean = "<<corrFactors[0]<<" "<<corrFactors[1]<<" "<<corrFactors[2]<<" "<<corrfact<<endl;
  float corrfact_y = SimpleMean(corrFactors_y,0,corrFactors_y.size()-1);
  cout<<"correlation factors (for yields) from the three methods (rms only), and mean = "<<corrFactors_y[0]<<" "<<corrFactors_y[1]<<" "<<corrFactors_y[2]<<" "<<corrfact_y<<endl;
  vector<float> corrfactFinal,corrfactFinal_y;
  corrfactFinal.push_back((corrDevToMean*DevToMean[1]*DevToMean[2] 
			   + corrfact * sqrt(pow(nomiAndErr[1][1],2) - pow(DevToMean[1],2)) * sqrt(pow(nomiAndErr[2][1],2) - pow(DevToMean[2],2))
			   )/ (nomiAndErr[1][1]*nomiAndErr[2][1]) );
  cout<<"final correlation factor = "<<corrfactFinal[0]<<endl;
  corrfactFinal_y.push_back((corrDevToMean*DevToMean_y[1]*DevToMean_y[2] 
			     + corrfact_y * sqrt(pow(nomiAndErr_y[1][1],2) - pow(DevToMean_y[1],2)) * sqrt(pow(nomiAndErr_y[2][1],2) - pow(DevToMean_y[2],2))
			   )/ (nomiAndErr_y[1][1]*nomiAndErr_y[2][1]) );
  cout<<"final correlation factor (for yields) = "<<corrfactFinal_y[0]<<endl;

  //Record the final correction 
  TFile *fout = new TFile("AccEffFrom2ndStepToys.root",keepfile?"UPDATE":"RECREATE");  
  fout->WriteObject(&nomiAndErr,"InvAccEffFromCorrMC_withSystErr_"+(TString)(ispp?"pp":"PbPb")+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
  fout->WriteObject(&corrfactFinal,"InvAccEffFromCorrMC_LinearisedCorrelationFactor_"+(TString)(ispp?"pp":"PbPb")+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
  fout->WriteObject(&nomiAndErr_y,"CorrYieldsFromCorrMC_withSystErr_"+(TString)(ispp?"pp":"PbPb")+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
  fout->WriteObject(&corrfactFinal_y,"CorrYieldsFromCorrMC_LinearisedCorrelationFactor_"+(TString)(ispp?"pp":"PbPb")+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
  fout->Close();
}

void plotToyAccEff(bool secondStep=true){

  //Create Hub gathering various data 
  Hub H = Hub(secondStep,false);
  H.SetAccEff();
  H.SetFit(false);
  H.SetMetafit();
  H.SetTnP();

  cout<<"\n**** pp pT dependence ****\n"<<endl;
  ToyAccEff(H,true,secondStep,false,false);
  cout<<"\n**** PbPb pT dependence ****\n"<<endl;
  ToyAccEff(H,false,secondStep,false,true);
  cout<<"\n**** PbPb centrality dependence ****\n"<<endl;
  ToyAccEff(H,false,secondStep,true,true); //centrality dep
}
