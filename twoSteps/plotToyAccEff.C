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

void ToyAccEff(Hub H, bool secondStep=true, bool inCentBins=false, bool keepfile=true, bool maxRMS=false){
  
  //for error on yields, without acceff
  vector<Result> ResAll = vector<Result>{ SumResult({H.fit.pp_pt,H.metafit.pp_pt,H.TnP.pp_pt} , H.Ycorr.pp_pt.Val,true) , 
					  inCentBins?SumResult({H.fit.PbPb_cent,H.metafit.PbPb_cent,H.TnP.PbPb_cent} , H.Ycorr.PbPb_cent.Val,true):
					  (SumResult({H.fit.PbPb_pt,H.metafit.PbPb_pt,H.TnP.PbPb_pt} , H.Ycorr.PbPb_pt.Val,true)) ,
					  inCentBins?SumResult({H.fit.RAA_cent,H.metafit.RAA_cent,H.TnP.RAA_cent} , H.Ycorr.RAA_cent.Val,true):
					  (SumResult({H.fit.RAA_pt,H.metafit.RAA_pt,H.TnP.RAA_pt} , H.Ycorr.RAA_pt.Val,true)) };

  //eff
  vector<vector<vector<float> > > eff_biased = inCentBins?H.eff_biased_cent:H.eff_biased;
  vector<vector<float> > eff_oneBinned = inCentBins?H.eff_oneBinned_cent:H.eff_oneBinned;
  vector<vector<float> > eff_oneBinned_2ndStep = inCentBins?H.eff_oneBinned_cent_2ndStep:H.eff_oneBinned_2ndStep;
  if(inCentBins){ //take non-centrality for integrated pp
    eff_biased[0].push_back( H.eff_biased[0][0] );
    eff_oneBinned[0].push_back( H.eff_oneBinned[0][0] );
    eff_oneBinned_2ndStep[0].push_back( H.eff_oneBinned_2ndStep[0][0] );
  }

  //nominal corrected yields
  vector<vector<float> > ycorr = vector<vector<float> >{ H.Ycorr.pp_pt.Val, inCentBins?H.Ycorr.PbPb_cent.Val:H.Ycorr.PbPb_pt.Val, inCentBins?H.Ycorr.RAA_cent.Val:H.Ycorr.RAA_pt.Val }; //if we take the simple nominal method and not an average of methods (cf Draw_metafit), this is supposed to = yield_postfit / (acc_oneBinned*eff_oneBinned), with _2ndStep in all component if(secondStep)

  //histograms of acc and eff filled with toys
  vector<vector<vector<TH1F*> > > accI_toys(3, vector<vector<TH1F*> >(_biasNmeth));
  vector<vector<vector<TH1F*> > > effI_toys(3, vector<vector<TH1F*> >(_biasNmeth));
  vector<vector<vector<TH1F*> > > acceffI_toys(3, vector<vector<TH1F*> >(_biasNmeth));
  vector<vector<vector<TH1F*> > > ycorr_toys(3, vector<vector<TH1F*> >(_biasNmeth));

  //mins and maxs
  vector<vector<float> > minacc(3);
  vector<vector<float> > maxacc(3);
  vector<vector<float> > mineff(3);
  vector<vector<float> > maxeff(3);
  vector<vector<float> > minacceff(3);
  vector<vector<float> > maxacceff(3);
  vector<vector<float> > minycorr(3);
  vector<vector<float> > maxycorr(3);

  //output
  vector<vector<vector<float> > > acceffI_biased(3,vector<vector<float> >(_NanaBins+1));
  vector<vector<float> > DevToMean(3,vector<float> (_NanaBins+1));
  vector<vector<float> > DevToMean_y(3,vector<float> (_NanaBins+1));
  vector<vector<vector<double> > > nomiAndErr(3,vector<vector<double> >(_NanaBins+1));
  vector<vector<vector<double> > > nomiAndErr_y(3,vector<vector<double> >(_NanaBins+1));
  vector<vector<vector<float> > > ycorr_biased(3,vector<vector<float> >(_NanaBins+1));

  for(int col=0;col<3;col++){
    for(int b=0;b<=((inCentBins && col==0)?0:_NanaBins);b++){ //only integrated for pp when inCentBins=true

      cout<<"\n\n***** "<<(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))<<endl<<endl;
      cout<<  "\n***** Bin "<<(TString)((b==0)?"integrated":to_string(b))<<endl<<endl;

      //AccAndEff values
      int bacc = inCentBins?0:b; //bin number for acceptance (integrated when doing centrality dependence)
      float accI_unbiased = (col<2)?( 
				     1/H.acc_oneBinned[col][bacc] //main value
				     ):( H.acc_oneBinned[0][bacc]/H.acc_oneBinned[1][bacc] ); //for RAA: (1/acc(PbPb)) / (1/acc(pp))
      float effI_unbiased = (col<2)?(
				     1/eff_oneBinned[col][b] //main value
				     ):( eff_oneBinned[0][bacc]/eff_oneBinned[1][b] );
      float acceffI_unbiased = accI_unbiased*effI_unbiased;

      float accI_1stStep = secondStep?( (col<2)?( 
						 1/H.acc_oneBinned_2ndStep[col][bacc] //main value
						 ):( H.acc_oneBinned_2ndStep[0][bacc]/H.acc_oneBinned_2ndStep[1][bacc] ) ):1;
      float effI_1stStep = secondStep?( (col<2)?(
						 1/eff_oneBinned_2ndStep[col][b] //main value
						 ):( eff_oneBinned_2ndStep[0][bacc]/eff_oneBinned_2ndStep[1][b] ) ):1;
      float acceffI_1stStep = accI_1stStep*effI_1stStep;

      maxacc[col].push_back( (col<2)?(
				      MaxVec(H.acc_biased[col][bacc], 1/accI_unbiased, (col==0)?1e5:(_biasNmeth+1*_biasNtoys)) //exclude outliers (especially 4) from third method
				  ):( MaxVec(H.acc_biased[1][bacc], H.acc_oneBinned[1][bacc]/H.acc_oneBinned[0][bacc], (col==0)?1e5:(_biasNmeth+1*_biasNtoys), H.acc_biased[0][bacc]) ) ); //max of the vector of PbPb/pp acceptances
      minacc[col].push_back( (col<2)?(
				      MinVec(H.acc_biased[col][bacc], 1/accI_unbiased, (col!=1)?(_biasNmeth+1*_biasNtoys):1e5)
				  ):( MinVec(H.acc_biased[1][bacc], H.acc_oneBinned[1][bacc]/H.acc_oneBinned[0][bacc], (col!=1)?(_biasNmeth+1*_biasNtoys):1e5, H.acc_biased[0][bacc]) ) );
      maxeff[col].push_back( (col<2)?(
				      MaxVec(eff_biased[col][b], 1/effI_unbiased, (col==0)?1e5:(_biasNmeth+1*_biasNtoys))
				  ):( MaxVec(eff_biased[1][b], eff_oneBinned[1][b]/eff_oneBinned[0][bacc], (col==0)?1e5:(_biasNmeth+1*_biasNtoys), eff_biased[0][bacc]) ) );
      mineff[col].push_back( (col<2)?(
				      MinVec(eff_biased[col][b], 1/effI_unbiased, (col!=1)?(_biasNmeth+1*_biasNtoys):1e5)
				  ):( MinVec(eff_biased[1][b], eff_oneBinned[1][b]/eff_oneBinned[0][bacc], (col!=1)?(_biasNmeth+1*_biasNtoys):1e5, eff_biased[0][bacc]) ) );
      maxacceff[col].push_back( maxacc[col][b]*maxeff[col][b] );
      minacceff[col].push_back( minacc[col][b]*mineff[col][b] );
      float acceffI_default = secondStep?acceffI_1stStep:acceffI_unbiased ;

      //corrected yields values
      for(int t=0;t<_biasNmeth*(_biasNtoys+1);t++){
	//yield-ratio-variation * nominal-yield * cancel-nominal-AccEff / varied-AccEff
	if(col<2)
	  ycorr_biased[col][b].push_back( ycorr[col][b] * H.corrYr_varied[col][bacc][t] / acceffI_default / (H.acc_biased[col][bacc][t] * eff_biased[col][b][t]) ); //take pt-integrated bin for variations of yields in centrality bins
	else
	  ycorr_biased[col][b].push_back( ycorr[col][b] / acceffI_default
					                * (H.corrYr_varied[1][bacc][t] / (H.acc_biased[1][bacc][t] * eff_biased[1][b][t]))
					                / (H.corrYr_varied[0][bacc][t] / (H.acc_biased[0][bacc][t] * eff_biased[0][bacc][t])) );
      }
      maxycorr[col].push_back( MaxVec(ycorr_biased[col][b], ycorr[col][b], 1e5) );
      minycorr[col].push_back( MinVec(ycorr_biased[col][b], ycorr[col][b], 1e5) );
      float ycorr_unbiased = ycorr[col][b] * acceffI_unbiased / acceffI_default;
      float ycorr_1stStep = ycorr[col][b] * acceffI_1stStep / acceffI_default;
    
      for(int m=0;m<_biasNmeth;m++){
	//nominal
	acceffI_biased[col][b].push_back( (col<2)?( 1/(H.acc_biased[col][bacc][m] * eff_biased[col][b][m])
						):( (H.acc_biased[0][bacc][m] * eff_biased[0][bacc][m]) / (H.acc_biased[1][bacc][m] * eff_biased[1][b][m]) ) );
      }

      for(int m=0;m<_biasNmeth;m++){
	//histograms definition
	accI_toys[col][m].push_back(new TH1F("acc_toys"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":"")+"_meth"+(TString)to_string(m), "acceptance toys "+(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))+(TString)((b==0 || inCentBins)?"integrated":("p_{T} bin"+(TString)to_string(b)))+";"+(TString)((col<2)?"1/acceptance":"acceptance(pp) / acceptance(PbPb)")+";n_{toys}", 23, max(0., 1/maxacc[col][b]-0.05*(1/minacc[col][b] - 1/maxacc[col][b])) , 1/minacc[col][b]+0.1*(1/minacc[col][b] - 1/maxacc[col][b])));
	effI_toys[col][m].push_back(new TH1F("eff_toys"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":"")+"_meth"+(TString)to_string(m), "efficiency toys "+(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";"+(TString)((col<2)?"1/efficiency":"efficiency(pp) / efficiency(PbPb)")+";n_{toys}", 23, max(0.,1/maxeff[col][b]-0.05*(1/mineff[col][b] - 1/maxeff[col][b])) , 1/mineff[col][b]+0.1*(1/mineff[col][b] - 1/maxeff[col][b])));
	acceffI_toys[col][m].push_back(new TH1F("acceff_toys"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":"")+"_meth"+(TString)to_string(m), "acc #times eff toys "+(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";"+(TString)((col<2)?"1/(acc #times eff)":"(acc #times eff)(pp) / (acc #times eff) (PbPb)")+";n_{toys}", 23, max(0.,1/maxacceff[col][b]-0.05*(1/minacceff[col][b] - 1/maxacceff[col][b])) , 1/minacceff[col][b]+0.1*(1/minacceff[col][b] - 1/maxacceff[col][b])));
	ycorr_toys[col][m].push_back(new TH1F("ycorr_toys"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":"")+"_meth"+(TString)to_string(m), (TString)((col<2)?"yields":"R_{PbPb}")+" for acc #times eff toys "+(TString)((col==0)?"pp ":((col==1)?"PbPb ":""))+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";"+(TString)((col<2)?"N_{corrected}":"R_{PbPb}")+";n_{toys}", 23, max(0.,minycorr[col][b]-0.05*(maxycorr[col][b] - minycorr[col][b])) , maxycorr[col][b]+0.1*(maxycorr[col][b] - minycorr[col][b])));

	//Run over all toys
	for(int t=0;t<_biasNtoys;t++){
	  int ntoy = _biasNmeth + m*_biasNtoys + t;
	  //cout<<"b toyNumber acc eff = "<<b<<" "<<t<<" "<< H.acc_biased[col][bacc][ntoy] <<" "<<eff_biased[col][b][ntoy]<<endl;
	  accI_toys[col][m][b]->Fill( (col<2)?( 1/(H.acc_biased[col][bacc][ntoy])
					    ):(    H.acc_biased[0][bacc][ntoy] / H.acc_biased[1][bacc][ntoy] ) );
	  effI_toys[col][m][b]->Fill( (col<2)?( 1/(eff_biased[col][b][ntoy])
					    ):(    eff_biased[0][bacc][ntoy] / eff_biased[1][b][ntoy] ) );
	  acceffI_biased[col][b].push_back( (col<2)?( 1/(H.acc_biased[col][bacc][ntoy] * eff_biased[col][b][ntoy])
						  ):(   (H.acc_biased[0][bacc][ntoy] * eff_biased[0][bacc][ntoy]) / (H.acc_biased[1][bacc][ntoy] * eff_biased[1][b][ntoy]) ) );
	  acceffI_toys[col][m][b]->Fill(acceffI_biased[col][b][ntoy]);
	  ycorr_toys[col][m][b]->Fill(ycorr_biased[col][b][ntoy]);
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
      TH1F *dum = new TH1F("line_unbiased"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":""),"line_unbiased",10,0,1);
      dum->SetLineColor(kBlack);
      dum->SetLineWidth(3);
      dum->SetLineStyle(7);
      TH1F *dum2 = new TH1F("line_1stStep"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":""),"line_1stStep",10,0,1);
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
      if(b>0 || !inCentBins){
	gStyle->SetOptStat(0);
	TLegend *leg = new TLegend(0.12,0.54,0.45,0.89);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);

	TCanvas *c1 = new TCanvas("c1_"+(TString)to_string(b),"c1",3500,1200);
	c1->Divide(4,1);

	c1->cd(1)->SetLeftMargin(0.09);
	c1->cd(1)->SetRightMargin(0.02);
	for(int m=0;m<_biasNmeth;m++){
	  accI_toys[col][m][b]->SetLineColor(cols[m]);
	  accI_toys[col][m][b]->SetLineWidth(2);
	  accI_toys[col][m][b]->GetYaxis()->SetRangeUser(0,1.15*accI_toys[col][m][b]->GetMaximum());
	  accI_toys[col][m][b]->GetYaxis()->SetTitleOffset(1.15);
	  accI_toys[col][m][b]->Draw((TString)((m==0)?"hist":"histsame"));
	  leg->AddEntry(accI_toys[col][m][b], _biasMethName[m]);
	  double accI_biased = (col<2)?( 1/H.acc_biased[col][bacc][m] ):( H.acc_biased[0][bacc][m] / H.acc_biased[1][bacc][m] );
	  biasedNom[m].DrawLine(accI_biased,0, accI_biased, 0.4*accI_toys[col][0][b]->GetMaximum()); //positioned after SetRangeUser
	}
	unbiased->DrawLine(accI_unbiased,0, accI_unbiased, 0.4*accI_toys[col][0][b]->GetMaximum());
	corr1stStep->DrawLine(accI_1stStep,0, accI_1stStep, 0.4*accI_toys[col][0][b]->GetMaximum());
	leg->AddEntry(dum, "unbiased MC");
	leg->AddEntry(dum2, "1^{st}-step corrected");

	c1->cd(2)->SetLeftMargin(0.09);
	c1->cd(2)->SetRightMargin(0.02);
	for(int m=0;m<_biasNmeth;m++){
	  effI_toys[col][m][b]->SetLineColor(cols[m]);
	  effI_toys[col][m][b]->SetLineWidth(2);
	  effI_toys[col][m][b]->GetYaxis()->SetRangeUser(0,1.15*effI_toys[col][m][b]->GetMaximum());
	  effI_toys[col][m][b]->GetYaxis()->SetTitleOffset(1.15);
	  effI_toys[col][m][b]->Draw((TString)((m==0)?"hist":"histsame"));
	  double effI_biased = (col<2)?( 1/eff_biased[col][b][m] ):( eff_biased[0][bacc][m] / eff_biased[1][b][m] );
	  biasedNom[m].DrawLine(effI_biased,0, effI_biased, 0.4*effI_toys[col][0][b]->GetMaximum());
	}
	unbiased->DrawLine(effI_unbiased,0, effI_unbiased, 0.4*effI_toys[col][0][b]->GetMaximum());
	corr1stStep->DrawLine(effI_1stStep,0, effI_1stStep, 0.4*effI_toys[col][0][b]->GetMaximum());

	c1->cd(3)->SetLeftMargin(0.09);
	c1->cd(3)->SetRightMargin(0.02);
	for(int m=0;m<_biasNmeth;m++){
	  acceffI_toys[col][m][b]->SetLineColor(cols[m]);
	  acceffI_toys[col][m][b]->SetLineWidth(2);
	  acceffI_toys[col][m][b]->GetYaxis()->SetRangeUser(0,1.15*acceffI_toys[col][m][b]->GetMaximum());
	  acceffI_toys[col][m][b]->GetYaxis()->SetTitleOffset(1.15);
	  acceffI_toys[col][m][b]->Draw((TString)((m==0)?"hist":"histsame"));
	  biasedNom[m].DrawLine(acceffI_biased[col][b][m],0, acceffI_biased[col][b][m], 0.4*acceffI_toys[col][0][b]->GetMaximum());
	}
	unbiased->DrawLine(acceffI_unbiased,0, acceffI_unbiased, 0.4*acceffI_toys[col][0][b]->GetMaximum());
	corr1stStep->DrawLine(acceffI_1stStep,0, acceffI_1stStep, 0.4*acceffI_toys[col][0][b]->GetMaximum());
	leg->Draw("same");

	c1->cd(4)->SetLeftMargin(0.09);
	c1->cd(4)->SetRightMargin(0.02);
	for(int m=0;m<_biasNmeth;m++){
	  ycorr_toys[col][m][b]->SetLineColor(cols[m]);
	  ycorr_toys[col][m][b]->SetLineWidth(2);
	  ycorr_toys[col][m][b]->GetYaxis()->SetRangeUser(0,1.15*ycorr_toys[col][m][b]->GetMaximum());
	  ycorr_toys[col][m][b]->GetYaxis()->SetTitleOffset(1.15);
	  ycorr_toys[col][m][b]->Draw((TString)((m==0)?"hist":"histsame"));
	  biasedNom[m].DrawLine(ycorr_biased[col][b][m],0, ycorr_biased[col][b][m], 0.4*ycorr_toys[col][0][b]->GetMaximum());
	}
	unbiased->DrawLine(ycorr_unbiased, 0, ycorr_unbiased, 0.4*ycorr_toys[col][0][b]->GetMaximum());
	corr1stStep->DrawLine(ycorr_1stStep,0, ycorr_1stStep, 0.4*ycorr_toys[col][0][b]->GetMaximum());
	
	c1->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+".pdf");

	TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
	c3->SetLeftMargin(0.09);
	c3->SetRightMargin(0.02);
	for(int m=0;m<_biasNmeth;m++){
	  acceffI_toys[col][m][b]->Draw((TString)((m==0)?"hist":"histsame"));
	  biasedNom[m].DrawLine(acceffI_biased[col][b][m],0, acceffI_biased[col][b][m], 0.4*acceffI_toys[col][0][b]->GetMaximum());
	}
	unbiased->DrawLine(acceffI_unbiased,0, acceffI_unbiased, 0.4*acceffI_toys[col][0][b]->GetMaximum());
	corr1stStep->DrawLine(acceffI_1stStep,0, acceffI_1stStep, 0.4*acceffI_toys[col][0][b]->GetMaximum());
	leg->Draw("same");

	c3->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+"_onlyAccEffProduct.pdf");

	TCanvas *c4 = new TCanvas("c4","c4",1500,1500);
	c4->SetLeftMargin(0.09);
	c4->SetRightMargin(0.02);
	for(int m=0;m<_biasNmeth;m++){
	  ycorr_toys[col][m][b]->Draw((TString)((m==0)?"hist":"histsame"));
	  biasedNom[m].DrawLine(ycorr_biased[col][b][m],0, ycorr_biased[col][b][m], 0.4*ycorr_toys[col][0][b]->GetMaximum());
	}
	unbiased->DrawLine(ycorr_unbiased, 0, ycorr_unbiased, 0.4*ycorr_toys[col][0][b]->GetMaximum());
	corr1stStep->DrawLine(ycorr_1stStep,0, ycorr_1stStep, 0.4*ycorr_toys[col][0][b]->GetMaximum());
	leg->Draw("same");
	gPad->Update();
	leg->SetX1NDC(0.6);
	leg->SetX2NDC(0.98);
	gPad->Modified();

	c4->SaveAs("figs/AcceptanceEfficiencyToys_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+"_onlyCorrYields.pdf");
      }

      //Get mean of the three possible nominal values, and low and high RMS from this mean
      //mean
      double mean = SimpleMean(acceffI_biased[col][b],0, _biasNmeth-1);
      double mean_y = SimpleMean(ycorr_biased[col][b],0, _biasNmeth-1);//ycorr[col][b] * mean / acceffI_default;
      //Deviation of the three nominals to the mean
      double DevToMeanLo = 0, DevToMeanHi = 0;
      for(int m=0;m<_biasNmeth;m++){
	if(acceffI_biased[col][b][m] < mean - DevToMeanLo) DevToMeanLo = fabs(acceffI_biased[col][b][m]-mean);
	if(acceffI_biased[col][b][m] > mean + DevToMeanHi) DevToMeanHi = fabs(acceffI_biased[col][b][m]-mean);
      }
      float DevToMeanLo_y = DevToMeanLo * ycorr[col][b] / acceffI_default;
      float DevToMeanHi_y = DevToMeanHi * ycorr[col][b] / acceffI_default;
      //The RMS's for the three methods, from their respective nominal
      vector<vector<double> > meanAndRMS,meanAndRMS_y; 
      for(int m=0;m<_biasNmeth;m++) {
	int start = m, end = m+1;
	if(!maxRMS){ start = 0; end = _biasNmeth;}
	meanAndRMS.push_back( DoubleSidedRMS(acceffI_biased[col][b],  //dummy if m<__biasNmeth-1 and maxRMS=false
					     _biasNmeth + start*_biasNtoys, _biasNmeth + end*_biasNtoys-1,
					     maxRMS ? acceffI_biased[col][b][m] : mean) );
	meanAndRMS_y.push_back( DoubleSidedRMS(ycorr_biased[col][b], 
					       _biasNmeth + start*_biasNtoys, _biasNmeth + end*_biasNtoys-1,
					       maxRMS ? ycorr_biased[col][b][m] : mean_y) );
      }

      //Find max RMS of the three methods
      meanAndRMS.push_back(vector<double>{mean,0,0,0});
      meanAndRMS_y.push_back(vector<double>{mean_y,0,0,0});
      if(maxRMS){
	for(int m=0;m<_biasNmeth;m++){
	  if(meanAndRMS[m][1] > meanAndRMS[_biasNmeth][1]) meanAndRMS[_biasNmeth][1] = meanAndRMS[m][1];//total
	  if(meanAndRMS[m][2] > meanAndRMS[_biasNmeth][2]) meanAndRMS[_biasNmeth][2] = meanAndRMS[m][2];//lo
	  if(meanAndRMS[m][3] > meanAndRMS[_biasNmeth][3]) meanAndRMS[_biasNmeth][3] = meanAndRMS[m][3];//hi
	  if(meanAndRMS_y[m][1] > meanAndRMS_y[_biasNmeth][1]) meanAndRMS_y[_biasNmeth][1] = meanAndRMS_y[m][1];//total
	  if(meanAndRMS_y[m][2] > meanAndRMS_y[_biasNmeth][2]) meanAndRMS_y[_biasNmeth][2] = meanAndRMS_y[m][2];//lo
	  if(meanAndRMS_y[m][3] > meanAndRMS_y[_biasNmeth][3]) meanAndRMS_y[_biasNmeth][3] = meanAndRMS_y[m][3];//hi
	}
      }
      else{
	for(int e=1;e<=3;e++){
	  meanAndRMS[_biasNmeth][e] = meanAndRMS[_biasNmeth-1][e];
	  meanAndRMS_y[_biasNmeth][e] = meanAndRMS_y[_biasNmeth-1][e];
	}
      }

      //Write output
      cout<<(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))<<"bin "<<b<<" unbiased one-binned AccEffInverse = "<<acceffI_unbiased<<endl;
      cout<<(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))<<"bin "<<b<<" 1st-step-corrected one-binned AccEffInverse = "<<acceffI_1stStep<<endl;
      for(int m=0;m<(maxRMS?_biasNmeth:1);m++)
	cout<<(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))<<"bin "<<b<<" mean, lo rms/mean, up rms/mean, relative_difference_with_unbiased of AccEffInverse for "<<(TString)(maxRMS?("toys of biasing method "+_biasMethName[m]):"all toys")<<" = "<<meanAndRMS[m][0]<<" "<<meanAndRMS[m][2]/meanAndRMS[m][0]<<" "<<meanAndRMS[m][3]/meanAndRMS[m][0]<<" "<<(meanAndRMS[m][0] - acceffI_unbiased)/acceffI_unbiased <<endl;
      cout<<"Mean of the three nominals (final nominal), max rel Up deviation, max rel Lo deviation = "<<meanAndRMS[_biasNmeth][0]<<" "<<DevToMeanHi/meanAndRMS[_biasNmeth][0]<<" "<<DevToMeanLo/meanAndRMS[_biasNmeth][0]<<endl;
      cout<<(TString)(maxRMS?"Maximal (among three methods) ":"")<<"rel RMS total, lo, hi = "<<meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0]<<endl;

      cout<<(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))<<"bin "<<b<<" ycorr from unbiased one-binned AccEff = "<<ycorr_unbiased<<endl;
      cout<<(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))<<"bin "<<b<<" ycorr from 1st-step-corrected one-binned AccEff = "<<ycorr_1stStep<<endl;
      for(int m=0;m<(maxRMS?_biasNmeth:1);m++)
	cout<<(TString)((col==0)?"pp ":((col==1)?"PbPb ":"RAA "))<<"bin "<<b<<" mean, lo rms/mean, up rms/mean, relative_difference_with_unbiased of ycorr for "<<(TString)(maxRMS?("toys of biasing method "+_biasMethName[m]):"all \
toys")<<" = "<<meanAndRMS_y[m][0]<<" "<<meanAndRMS_y[m][2]/meanAndRMS_y[m][0]<<" "<<meanAndRMS_y[m][3]/meanAndRMS_y[m][0]<<" "<<(meanAndRMS_y[m][0] - ycorr_unbiased)/ycorr_unbiased <<endl;
      cout<<"Mean of the three nominals (final nominal), max rel Up deviation, max rel Lo deviation = "<<meanAndRMS_y[_biasNmeth][0]<<" "<<DevToMeanHi_y/meanAndRMS_y[_biasNmeth][0]<<" "<<DevToMeanLo_y/meanAndRMS_y[_biasNmeth][0]<<endl;
      cout<<(TString)(maxRMS?"Maximal (among three methods) ":"")<<"rel RMS total, lo, hi = "<<meanAndRMS_y[_biasNmeth][1]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][2]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][3]/meanAndRMS_y[_biasNmeth][0]<<endl;

      //Add the deviation of the three nominal to their means into the final syst error
      DevToMean[col][b] = (DevToMeanLo+DevToMeanHi)/2;
      meanAndRMS[_biasNmeth][1] = sqrt(pow(meanAndRMS[_biasNmeth][1],2) + pow(DevToMean[col][b],2));
      meanAndRMS[_biasNmeth][2] = sqrt(pow(meanAndRMS[_biasNmeth][2],2) + pow(DevToMeanLo,2));
      meanAndRMS[_biasNmeth][3] = sqrt(pow(meanAndRMS[_biasNmeth][3],2) + pow(DevToMeanHi,2));
      cout<<"Final acceffInverse rel uncertainty total, lo, hi = "<<meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0]<<" "<<meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0]<<endl;
      // cout<<"Fit rel uncertainty total, lo, hi = "<<(ispp?H.fit.pp_pt:(inCentBins?H.fit.PbPb_cent:H.fit.PbPb_pt)).RelErr[b]<<" "<<(ispp?H.fit.pp_pt:(inCentBins?H.fit.PbPb_cent:H.fit.PbPb_pt)).RelErrLo[b]<<" "<<(ispp?H.fit.pp_pt:(inCentBins?H.fit.PbPb_cent:H.fit.PbPb_pt)).RelErrHi[b]<<endl;
      // cout<<"Metafit rel uncertainty total, lo, hi = "<<(ispp?H.metafit.pp_pt:(inCentBins?H.metafit.PbPb_cent:H.metafit.PbPb_pt)).RelErr[b]<<" "<<(ispp?H.metafit.pp_pt:(inCentBins?H.metafit.PbPb_cent:H.metafit.PbPb_pt)).RelErrLo[b]<<" "<<(ispp?H.metafit.pp_pt:(inCentBins?H.metafit.PbPb_cent:H.metafit.PbPb_pt)).RelErrHi[b]<<endl;
      cout<<"Full rel uncertainty without AccEff total, lo, hi = "<<ResAll[col].RelErr[b]<<" "<<ResAll[col].RelErrLo[b]<<" "<<ResAll[col].RelErrHi[b]<<endl;
      cout<<"Simple quadSum of AccEff and others rel uncertainties, total, lo, hi = "<<quadSum(ResAll[col].RelErr[b],meanAndRMS[_biasNmeth][1]/meanAndRMS[_biasNmeth][0])<<" "<<quadSum(ResAll[col].RelErrLo[b],meanAndRMS[_biasNmeth][2]/meanAndRMS[_biasNmeth][0])<<" "<<quadSum(ResAll[col].RelErrHi[b],meanAndRMS[_biasNmeth][3]/meanAndRMS[_biasNmeth][0])<<endl;
      DevToMean_y[col][b] = (DevToMeanLo_y+DevToMeanHi_y)/2;
      meanAndRMS_y[_biasNmeth][1] = sqrt(pow(meanAndRMS_y[_biasNmeth][1],2) + pow(DevToMean_y[col][b],2));
      meanAndRMS_y[_biasNmeth][2] = sqrt(pow(meanAndRMS_y[_biasNmeth][2],2) + pow(DevToMeanLo_y,2));
      meanAndRMS_y[_biasNmeth][3] = sqrt(pow(meanAndRMS_y[_biasNmeth][3],2) + pow(DevToMeanHi_y,2));
      cout<<"Final ycorr rel uncertainty total, lo, hi = "<<meanAndRMS_y[_biasNmeth][1]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][2]/meanAndRMS_y[_biasNmeth][0]<<" "<<meanAndRMS_y[_biasNmeth][3]/meanAndRMS_y[_biasNmeth][0]<<endl;

      nomiAndErr[col][b] = meanAndRMS[_biasNmeth];
      nomiAndErr_y[col][b] = meanAndRMS_y[_biasNmeth];
      
      //
      // 2D+3D plot AccEff - observed yield
      //
      vector<float> obsyield;
      float obsYmin=1e20, obsYmax=0;

      for(int t=0;t<_biasNmeth*(1+_biasNtoys);t++){
	if(col<2) obsyield.push_back( ((col==0)?H.YieldsPostfit_pp:(inCentBins?H.YieldsPostfit_PbPbcent:H.YieldsPostfit_PbPb))[b] * H.corrYr_varied[col][bacc][t] ); //position 2-3 for bins 1-2
	else obsyield.push_back( ((inCentBins?H.YieldsPostfit_PbPbcent:H.YieldsPostfit_PbPb)[b] / H.YieldsPostfit_pp[b]) * H.corrYr_varied[1][bacc][t]/H.corrYr_varied[0][bacc][t] );
	
	if(obsyield[t]<obsYmin) obsYmin = obsyield[t];
	if(obsyield[t]>obsYmax) obsYmax = obsyield[t];
      }

      TH2F* AEvsYield = new TH2F("AccEffVsYields_toys"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":""), " 1/(acc #times eff) VS obs_yields "+(TString)((col==0)?"pp ":((col==1)?"PbPb ":""))+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";"+(TString)((col<2)?"observed yield":"ratio of observed yields")+(TString)((col<2)?";1/(acc #times eff)":";(acc #times eff)(pp)/(acc #times eff)(PbPb)"), 40, max(0., obsYmin-0.05*(obsYmax-obsYmin)), obsYmax+0.05*(obsYmax-obsYmin), 40, max(0.,1/maxacceff[col][b]-0.05*(1/minacceff[col][b] - 1/maxacceff[col][b])) , 1/minacceff[col][b]+0.1*(1/minacceff[col][b] - 1/maxacceff[col][b]));
      TH2F* AEvsYield3D = new TH2F("AccEffVsYields3D_toys"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":""), "(AccEff_correction) #times  obs_yields for "+(TString)((col==0)?"pp ":((col==1)?"PbPb ":""))+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";"+(TString)((col<2)?"observed yield":"ratio of observed yields")+(TString)((col<2)?";1/(acc #times eff);corrected yield":";(acc #times eff)(pp)/(acc #times eff)(PbPb);R_{PbPb}"), 180, max(0., (double)obsYmin), obsYmax, 180, max(0.,(double)1/maxacceff[col][b]) , 1/minacceff[col][b]);
      TH2F* CorrVsObsYield3D = new TH2F("CorrVsObsYield3D_toys"+(TString)((col==0)?"_pp":((col==1)?"_PbPb":"_RAA"))+"_bin"+(TString)to_string(b)+(TString)(inCentBins?"cent":""), "Corrected yield VS obs_yields for "+(TString)((col==0)?"pp ":((col==1)?"PbPb ":""))+(TString)((b==0)?"integrated":((TString)(inCentBins?"cent ":"p_{T} ")+"bin"+(TString)to_string(b)))+";"+(TString)((col<2)?"observed yield":"ratio of observed yields")+(TString)((col<2)?";corrected yield;1/(acc #times eff)":"R_{PbPb};(acc #times eff)(pp)/(acc #times eff)(PbPb)"), 180, max(0., obsYmin-0.05*(obsYmax-obsYmin)), obsYmax+0.05*(obsYmax-obsYmin), 180, max(0.,(double)minycorr[col][b]) , maxycorr[col][b]);
      //for RAA only
      //AccEff when inCentBins
      TH2F* AEppVsPbPb3D = (col<2 || !inCentBins)?(new TH2F()):( new TH2F("AEppVsPbPb3D_toys_RAA_bin"+(TString)to_string(b)+"cent", "AccEff pp integrated Vs PbPb cent bin"+(TString)to_string(b)+";1/acc #times eff (pp);1/acc #times eff (PbPb);acc #times eff pp/PbPb", 180, max(0.,(double)1/maxacceff[0][0]) , 1/minacceff[0][0], 180, max(0.,(double)1/maxacceff[1][b]) , 1/minacceff[1][b]) );
      TH2F* CorrYppVsPbPb3D = (col<2 || inCentBins)?(new TH2F()):( new TH2F("CorrYppVsPbPb3D_toys_RAA_bin"+(TString)to_string(b), "Corrected yield PbPb vs pp pT bin"+(TString)to_string(b)+";corrected yield (pp);corrected yield (PbPb);R_{PbPb}", 180, max(0.,(double)minycorr[0][b]) , maxycorr[0][b], 180, max(0.,(double)minycorr[1][b]) , maxycorr[1][b] ) );
      TH2F* AEppVsPbPb = (col<2 || !inCentBins)?(new TH2F()):( new TH2F("AEppVsPbPb_toys_RAA_bin"+(TString)to_string(b)+"cent", "AccEff pp integrated Vs PbPb cent bin"+(TString)to_string(b)+";1/acc #times eff (pp);1/acc #times eff (PbPb);acc #times eff pp/PbPb", 40, 0.99*max(0.,(double)1/maxacceff[0][0]) , 1.01/minacceff[0][0], 40, 0.99*max(0.,(double)1/maxacceff[1][b]) , 1.01/minacceff[1][b]) );
      TH2F* CorrYppVsPbPb = (col<2 || inCentBins)?(new TH2F()):( new TH2F("CorrYppVsPbPb_toys_RAA_bin"+(TString)to_string(b), "Corrected yield PbPb vs pp pT bin"+(TString)to_string(b)+";corrected yield (pp);corrected yield (PbPb);R_{PbPb}", 40, max(0.,(double)0.99*minycorr[0][b]) , 1.01*maxycorr[0][b], 40, 0.99*max(0.,(double)minycorr[1][b]) , 1.01*maxycorr[1][b] ) );

      for(int t=_biasNmeth;t<_biasNmeth*(1+_biasNtoys);t++){
	AEvsYield->Fill(obsyield[t], acceffI_biased[col][b][t]);
	if(AEvsYield3D->GetBinContent(AEvsYield3D->FindBin(obsyield[t], acceffI_biased[col][b][t])) == 0) //goal here is to fill only once each non-empty bin
	  AEvsYield3D->Fill( obsyield[t], acceffI_biased[col][b][t], ycorr_biased[col][b][t]);
	if(CorrVsObsYield3D->GetBinContent(CorrVsObsYield3D->FindBin(obsyield[t], ycorr_biased[col][b][t])) == 0) //goal here is to fill only once each non-empty bin
	  CorrVsObsYield3D->Fill( obsyield[t], ycorr_biased[col][b][t], acceffI_biased[col][b][t]);

	if(col==2 && inCentBins){
	  if(AEppVsPbPb3D->GetBinContent(AEppVsPbPb3D->FindBin(acceffI_biased[0][0][t], acceffI_biased[1][b][t])) == 0) //goal here is to fill only once each non-empty bin
	    AEppVsPbPb3D->Fill( acceffI_biased[0][0][t], acceffI_biased[1][b][t], acceffI_biased[1][b][t] / acceffI_biased[0][0][t] );
	  AEppVsPbPb->Fill( acceffI_biased[0][0][t], acceffI_biased[1][b][t]);
	}     
	if(col==2 && !inCentBins){
	  if(CorrYppVsPbPb3D->GetBinContent(CorrYppVsPbPb3D->FindBin(ycorr_biased[0][b][t], ycorr_biased[1][b][t])) == 0) //goal here is to fill only once each non-empty bin
	    CorrYppVsPbPb3D->Fill( ycorr_biased[0][b][t], ycorr_biased[1][b][t], ycorr_biased[2][b][t] );
	  CorrYppVsPbPb->Fill( ycorr_biased[0][b][t], ycorr_biased[1][b][t]);
	}
      } 

      //graph with nominal
      double AEInom[] = {nomiAndErr[col][b][0]};
      double obsynom[] = {obsyield[0]};
      double corrynom[] = {ycorr_biased[col][b][0]};
      TGraph* g_nom = new TGraph(1,obsynom, AEInom);//"nominal AccEff vs observed yield"
      g_nom->SetMarkerSize(5);
      g_nom->SetMarkerStyle(34);
      TGraph* g_nom2 = new TGraph(1,obsynom, corrynom);
      g_nom2->SetMarkerSize(5);
      g_nom2->SetMarkerStyle(34);

      gStyle->SetPalette(kThermometer);

      if(b==0 && inCentBins) continue;

      TCanvas *c5 = new TCanvas("c5","c5",1500,1500);
      AEvsYield->Draw("COLZ");
      g_nom->Draw("P");
      c5->SaveAs("figs/AcceptanceEfficiencyVsObsYields_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+".pdf");

      TCanvas *c6 = new TCanvas("c6","c6",1500,1500);
      c6->SetRightMargin(0.19);
      c6->SetLeftMargin(0.15);
      CorrVsObsYield3D->GetZaxis()->SetRangeUser(max(0.,(double)(1/maxacceff[col][b])) , 1/minacceff[col][b]);
      CorrVsObsYield3D->SetContour(100);
      CorrVsObsYield3D->GetZaxis()->SetTitleOffset(1.5);
      CorrVsObsYield3D->Draw("COLZ");
      g_nom2->Draw("P");
      c6->SaveAs("figs/CorrectedVsObsYields3D_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+".pdf");

      TCanvas *c7 = new TCanvas("c7","c7",1500,1500);
      c7->SetRightMargin(0.19);
      c7->SetLeftMargin(0.12);
      AEvsYield3D->GetZaxis()->SetRangeUser(max(0.,(double)minycorr[col][b]) , maxycorr[col][b]);
      AEvsYield3D->SetContour(100);
      AEvsYield3D->GetZaxis()->SetTitleOffset(1.9);
      AEvsYield3D->Draw("COLZ");
      g_nom->Draw("P");

      float corrc7 = AEvsYield->GetCorrelationFactor();
      TLatex corrT;
      corrT.SetNDC();
      corrT.SetTextFont(42);
      corrT.SetTextSize(0.043);
      corrT.DrawLatex(0.2,0.14,Form("#rho = %.2f",corrc7));
      cout<<"correlation between AccEff correction and other uncertainties on the yields = "<<corrc7<<endl<<endl;

      c7->SaveAs("figs/AcceptanceEfficiencyTimesObsYields3D_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(TString)(secondStep?"_2ndStep":"")+"_kinBin"+(TString)to_string(bacc)+(inCentBins?("_centBin"+(TString)to_string(b)):"")+".pdf");

      if(col==2 && inCentBins){
	double AEInompp[] = {nomiAndErr[0][0][0]};
	double AEInomPbPb[] = {nomiAndErr[1][b][0]};

	TGraph* g_nom3 = new TGraph(1,AEInompp,AEInomPbPb);
	g_nom3->SetMarkerSize(5);
	g_nom3->SetMarkerStyle(34);

	TCanvas *c8 = new TCanvas("c8","c8",1500,1500);
	c8->SetRightMargin(0.19);
	c8->SetLeftMargin(0.15);
	AEppVsPbPb3D->GetZaxis()->SetRangeUser(max(0.,(double)1/maxacceff[2][b]) , 1/minacceff[2][b]);
	AEppVsPbPb3D->SetContour(100);
	AEppVsPbPb3D->GetZaxis()->SetTitleOffset(1.5);
	AEppVsPbPb3D->Draw("COLZ");
	g_nom3->Draw("P");

	float corrc8 = AEppVsPbPb->GetCorrelationFactor();
	TLatex corrT8;
	corrT8.SetNDC();
	corrT8.SetTextFont(42);
	corrT8.SetTextSize(0.043);
	corrT8.DrawLatex(0.2,0.14,Form("#rho = %.2f",corrc8));
	cout<<"correlation between AccEff correction for pp integrated and for PbPb cent bin "<<b<<" = "<<corrc8<<endl;

	c8->SaveAs("figs/AcceptanceEfficiency3D_ppintegrated_Vs_PbPb_centBin"+(TString)to_string(b)+(TString)(secondStep?"_2ndStep":"")+".pdf");
      }

      if(col==2 && !inCentBins){
	double corYnompp[] = {nomiAndErr_y[0][b][0]};
	double corYnomPbPb[] = {nomiAndErr_y[1][b][0]};

	TGraph* g_nom4 = new TGraph(1,corYnompp,corYnomPbPb);
	g_nom4->SetMarkerSize(5);
	g_nom4->SetMarkerStyle(34);

	TCanvas *c9 = new TCanvas("c9","c9",1500,1500);
	c9->SetRightMargin(0.19);
	c9->SetLeftMargin(0.15);
	CorrYppVsPbPb3D->GetZaxis()->SetRangeUser(max(0.,(double)minycorr[2][b]) , maxycorr[2][b]);
	CorrYppVsPbPb3D->SetContour(100);
	CorrYppVsPbPb3D->GetZaxis()->SetTitleOffset(1.5);
	CorrYppVsPbPb3D->Draw("COLZ");
	g_nom4->Draw("P");

	float corrc9 = CorrYppVsPbPb->GetCorrelationFactor();
	TLatex corrT9;
	corrT9.SetNDC();
	corrT9.SetTextFont(42);
	corrT9.SetTextSize(0.043);
	corrT9.DrawLatex(0.2,0.14,Form("#rho = %.2f",corrc9));
	cout<<"correlation between corrected yields for pp and for PbPb kinBin "<<b<<" = "<<corrc9<<endl;

	c9->SaveAs("figs/CorrectedYields3D_ppVsPbPb_kinBin"+(TString)to_string(b)+(TString)(secondStep?"_2ndStep":"")+".pdf");
      }

    }

    if(inCentBins && col==0) continue; //integrated pp was already done

    //Calculate ana bins correlation factor (result>0.999)
    //Will assume 0.3 correlation for subdominant DevToMean 
    vector<float> corrFactors,corrFactors_y;
    for(int m=0;m<_biasNmeth;m++){
      for(int b=1;b<=_NanaBins;b++){
	for(int b2=b+1;b2<=_NanaBins;b2++){
	  float sumcor = 0, sum1=0, sum2=0;
	  float sumcor_y = 0, sum1_y=0, sum2_y=0;
	  for(int i=_biasNmeth+m*_biasNtoys;i<_biasNmeth+(m+1)*_biasNtoys;i++){
	    sumcor += (acceffI_biased[col][b][i]-acceffI_biased[col][b][m]) * (acceffI_biased[col][b2][i]-acceffI_biased[col][b2][m]);
	    sum1 += pow((acceffI_biased[col][b][i]-acceffI_biased[col][b][m]),2);
	    sum2 += pow((acceffI_biased[col][b2][i]-acceffI_biased[col][b2][m]),2);
	    sumcor_y += (ycorr_biased[col][b][i]-ycorr_biased[col][b][m]) * (ycorr_biased[col][b2][i]-ycorr_biased[col][b2][m]);
	    sum1_y += pow((ycorr_biased[col][b][i]-ycorr_biased[col][b][m]),2);
	    sum2_y += pow((ycorr_biased[col][b2][i]-ycorr_biased[col][b2][m]),2);
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
    corrfactFinal.push_back((corrDevToMean*DevToMean[col][1]*DevToMean[col][2] 
			     + corrfact * sqrt(pow(nomiAndErr[col][1][1],2) - pow(DevToMean[col][1],2)) * sqrt(pow(nomiAndErr[col][2][1],2) - pow(DevToMean[col][2],2))
			     )/ (nomiAndErr[col][1][1]*nomiAndErr[col][2][1]) );
    cout<<"final correlation factor = "<<corrfactFinal[0]<<endl;
    corrfactFinal_y.push_back((corrDevToMean*DevToMean_y[col][1]*DevToMean_y[col][2] 
			       + corrfact_y * sqrt(pow(nomiAndErr_y[col][1][1],2) - pow(DevToMean_y[col][1],2)) * sqrt(pow(nomiAndErr_y[col][2][1],2) - pow(DevToMean_y[col][2],2))
			       )/ (nomiAndErr_y[col][1][1]*nomiAndErr_y[col][2][1]) );
    cout<<"final correlation factor (for yields) = "<<corrfactFinal_y[0]<<endl;

    //Record the final correction 
    TFile *fout = new TFile("AccEffFrom2ndStepToys.root",keepfile?"UPDATE":"RECREATE");  
    fout->WriteObject(&nomiAndErr[col],"InvAccEffFromCorrMC_withSystErr_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
    fout->WriteObject(&corrfactFinal,"InvAccEffFromCorrMC_LinearisedCorrelationFactor_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
    fout->WriteObject(&nomiAndErr_y[col],"CorrYieldsFromCorrMC_withSystErr_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
    fout->WriteObject(&corrfactFinal_y,"CorrYieldsFromCorrMC_LinearisedCorrelationFactor_"+(TString)((col==0)?"pp":((col==1)?"PbPb":"RAA"))+(inCentBins?"_inCentBins":"")+(TString)(secondStep?"_2ndStep":""));
    fout->Close();

  }

}

void plotToyAccEff(bool secondStep=true){

  //Create Hub gathering various data 
  Hub H = Hub(secondStep,false);
  H.SetAccEff();
  H.SetFit(false);
  H.SetMetafit();
  H.SetTnP();
  H.ScaleByLumi(true); //don't normalise pp and PbPb corrected yields to the cross sections

  cout<<"\n**** pT dependence ****\n"<<endl;
  ToyAccEff(H,secondStep,false,true);
  cout<<"\n**** centrality dependence ****\n"<<endl;
  ToyAccEff(H,secondStep,true,true); //centrality dep
}
