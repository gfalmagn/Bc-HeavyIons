#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"
#include "../helpers/SystematicsList.h"

void Draw_metafitSyst(bool ispp=true){
  
  vector<vector<vector<float> > > *Yields_postfit[2][_nMetafitSyst];
  vector<float> *rsig_relerr[3][_nMetafitSyst];
  vector<float> *eff_oneBinned[3];
  vector<float> *acc_oneBinned[3];

  vector<float> y_nom[3];
  vector<float> y_systErr[3];
  vector<float> y_systErr_maxDev[3];
  vector<float> y_systCorr[3];
  //vector<float> y_fitErr(nbins);

  TFile *infile = new TFile("corrected_yields.root","READ");

  vector<TGraphErrors*> gSyst[3], gSystRatio[3];
  const int nbins = _NanaBins;
  double xErr[nbins], yErr[3][_nMetafitSyst][nbins], yErr_ratio[3][_nMetafitSyst][nbins], zero[nbins], x[nbins];
  double y_oneBinned[3][_nMetafitSyst][nbins], y_ratio[3][_nMetafitSyst][nbins];

  for(int ispp=0;ispp<3;ispp++){
    cout<<"\n    ***********    "<<((ispp==1)?"pp":((ispp==0)?"PbPb":"RAA"))<<"    *************"<<endl;

    TFile *infile3, *infile2;
    if(ispp<2){
      infile3 = new TFile("../acceptance/acceptanceMap.root","READ");
      infile3->GetObject("acceptance_oneBinned", acc_oneBinned[ispp]);

      infile2 = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
      infile2->GetObject("efficiency_oneBinned"+(TString)((ispp==1)?"_pp":"_PbPb"), eff_oneBinned[ispp]);
    }
  
    for(int isys=0;isys<_nMetafitSyst;isys++){
      cout<<"\n  *** Building graph for systematic: "<<systPrettyName[isys]<<endl;
      if(ispp<2){
	infile->GetObject("Yields_postfit"+(TString)((ispp==1)?"_pp":"_PbPb")+systName[isys]+systExtName[isys], Yields_postfit[ispp][isys]);
	infile->GetObject("rsig_relerr"+(TString)((ispp==1)?"_pp":"_PbPb")+systName[isys]+systExtName[isys], rsig_relerr[ispp][isys]);

	for(int b=1;b<=nbins;b++){
	  cout<<"analysis bin #"<<b<<endl;
	  //      cout<<"one-binned acceptance, efficiency = "<<(*acc_oneBinned)[b]<<" "<<(*eff_oneBinned)[b]<<endl;

	  x[b-1] = (_BcPtmin[b]+_BcPtmax[b])/2;
	  xErr[b-1] = (-_BcPtmin[b]+_BcPtmax[b])/2;
	  zero[b-1] = 0;
	  y_oneBinned[ispp][isys][b-1] = (*Yields_postfit[ispp][isys])[0][b][0] / ((*acc_oneBinned[ispp])[b] * (*eff_oneBinned[ispp])[b]);
	  y_ratio[ispp][isys][b-1] = y_oneBinned[ispp][isys][b-1]/y_oneBinned[ispp][0][b-1];
	  yErr[ispp][isys][b-1] = (*rsig_relerr[ispp][isys])[b] * y_oneBinned[ispp][isys][b-1];
	  yErr_ratio[ispp][isys][b-1] = (*rsig_relerr[ispp][isys])[b];
	  cout<<"b, y, yErr = "<<b<<" "<<y_oneBinned[ispp][isys][b-1]<<" "<<yErr[ispp][isys][b-1]<<endl;
	}
      }
      else{
	for(int b=1;b<=nbins;b++){
	  y_oneBinned[ispp][isys][b-1] = y_oneBinned[0][isys][b-1] / y_oneBinned[1][isys][b-1];
	  y_ratio[ispp][isys][b-1] = y_ratio[0][isys][b-1] / y_ratio[1][isys][b-1];
	  yErr_ratio[ispp][isys][b-1] = sqrt(pow(yErr[0][isys][b-1]/y_oneBinned[0][isys][b-1],2) + pow(yErr[1][isys][b-1]/y_oneBinned[1][isys][b-1],2));
	  yErr[ispp][isys][b-1] = y_oneBinned[ispp][isys][b-1] * yErr_ratio[ispp][isys][b-1];
	}
      }
    
      gSyst[ispp].push_back(new TGraphErrors(nbins, x,y_oneBinned[ispp][isys], xErr,(isys==0)?yErr[ispp][0]:zero));
      gSystRatio[ispp].push_back(new TGraphErrors(nbins, x,y_ratio[ispp][isys], xErr,(isys==0)?yErr_ratio[ispp][0]:zero));

      gSyst[ispp][isys]->SetMarkerSize(3);
      gSyst[ispp][isys]->SetMarkerStyle(systStyle[isys]);
      gSyst[ispp][isys]->SetMarkerColor(systCol[isys]);
      gSyst[ispp][isys]->GetYaxis()->SetRangeUser(0.7*y_oneBinned[ispp][0][1],1.3*y_oneBinned[ispp][0][0]);
      gSyst[ispp][isys]->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);//+3*yErr[0]);
      gSyst[ispp][isys]->SetTitle("Metafit systematic variations of "+(TString)((ispp==2)?"R_{PbPb}":"corrected yields"));
      gSyst[ispp][isys]->GetYaxis()->SetTitle((ispp==2)?"R_{PbPb}":"corrected yield");
      gSyst[ispp][isys]->GetXaxis()->SetTitle("p_{T}(#mu#mu#mu) [GeV]");

      gSystRatio[ispp][isys]->SetMarkerSize(3);
      gSystRatio[ispp][isys]->SetMarkerStyle(systStyle[isys]);
      gSystRatio[ispp][isys]->SetMarkerColor(systCol[isys]);
      gSystRatio[ispp][isys]->GetYaxis()->SetRangeUser((ispp==1)?0.89:0.72,(ispp==1)?1.16:1.4);
      gSystRatio[ispp][isys]->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);//+3*yErr[0]);
      gSystRatio[ispp][isys]->SetTitle("Ratio of metafit systematics to nominal");
      gSystRatio[ispp][isys]->GetYaxis()->SetTitle((ispp==2)?"R_{PbPb}^{systematic}/R_{PbPb}^{nominal}":"yield^{systematic}/yield^{nominal}");
      gSystRatio[ispp][isys]->GetXaxis()->SetTitle("p_{T}(#mu#mu#mu) [GeV]");
    }

    //DRAW CMS Preliminary 
    TLatex CMStag;
    CMStag.SetNDC();
    CMStag.SetTextFont(42);
    CMStag.SetTextSize(0.035);
    
    if(ispp<2){
      TCanvas *c1 = new TCanvas("c1"+(TString)ispp,"c1",2000,2000);

      TLegend *leg = new TLegend(0.4,0.6,0.9,0.9);
      for(int isys=0;isys<_nMetafitSyst;isys++){
	leg->AddEntry(gSyst[ispp][isys],(isys==0)?"nominal, fit error only":systPrettyName[isys],(isys==0)?"lpe":"lp");
      }
      leg->SetTextSize(0.023);
      leg->SetBorderSize(0);

      gSyst[ispp][0]->Draw("AP");
      for(int isys=0;isys<_nMetafitSyst;isys++){
	gSyst[ispp][isys]->Draw("Psame");
      }
      gSyst[ispp][0]->Draw("Psame");
      leg->Draw("same");
      //gPad->SetLogy();

      CMStag.DrawLatex(0.19,0.85,"#font[61]{CMS "+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":""))+"}");
      CMStag.DrawLatex(0.19,0.80,"#font[52]{Work in progress}");

      c1->SetLeftMargin(0.15);
      c1->SaveAs("correctedYields_metafitSyst"+(TString)((ispp==1)?"_pp":((ispp==0)?"_PbPb":""))+".pdf");
    }

    TCanvas *c2 = new TCanvas("c2"+(TString)ispp,"c2",2000,2000);

    TLegend *leg2 = new TLegend(0.4,0.6,0.9,0.9);
    for(int isys=0;isys<_nMetafitSyst;isys++){
      leg2->AddEntry(gSystRatio[ispp][isys],(isys==0)?"nominal, fit error only":systPrettyName[isys],(isys==0)?"lpe":"lp");
    }
    leg2->SetTextSize(0.023);
    leg2->SetBorderSize(0);

    gSystRatio[ispp][0]->Draw("AP");
    for(int isys=0;isys<_nMetafitSyst;isys++){
      gSystRatio[ispp][isys]->Draw("Psame");
    }
    leg2->Draw("same");

    CMStag.DrawLatex(0.19,0.85,"#font[61]{CMS "+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":""))+"}");
    CMStag.DrawLatex(0.19,0.80,"#font[52]{Work in progress}");

    TLine *one = new TLine(_BcPtmin[0],1,_BcPtmax[0],1);
    one->SetLineWidth(2);
    one->SetLineStyle(2);
    one->Draw("same");

    c2->SetLeftMargin(0.15);
    c2->SaveAs((TString)((ispp<2)?"correctedYields":"RAA")+"_metafitSyst_ratioToNominal"+(TString)((ispp==1)?"_pp":((ispp==0)?"_PbPb":""))+".pdf");

    //******************************************************
    //Record the (nominal result +) systematic error
    //******************************************************

    for(int b=0;b<nbins;b++){
      y_nom[ispp].push_back( y_ratio[ispp][0][b] );
      y_systErr[ispp].push_back( 0 ); //rms of the deviations to nominal y
      y_systErr_maxDev[ispp].push_back( 0 );//max deviation from nominal
      for(int isys=1;isys<_nMetafitSyst;isys++){
	y_systErr_maxDev[ispp][b] = max(y_systErr[ispp][b] , (float)fabs(y_ratio[ispp][isys][b] - y_nom[ispp][b]) );
	y_systErr[ispp][b] += pow(y_ratio[ispp][isys][b] - y_nom[ispp][b] ,2); 
	cout<<"b, isys, y = "<<b<<" "<<isys<<" "<<y_ratio[ispp][isys][b]<<endl;
      }
      y_systErr[ispp][b] /= (_nMetafitSyst-1);
      y_systErr[ispp][b] = sqrt( y_systErr[ispp][b] );
      cout<<"******   b, y systematic error ="<<b<<" "<<y_systErr[ispp][b]<<endl;
      //y_fitErr[b] = (*rsig_relerr)[b+1] * y_nom[ispp][b];
    }

    int idx = 0;
    for(int b=0;b<nbins;b++){
      for(int b2=b+1;b2<nbins;b2++){
	y_systCorr[ispp].push_back(0);
	for(int isys=1;isys<_nMetafitSyst;isys++){
	  y_systCorr[ispp][idx] += (y_ratio[ispp][isys][b]-y_nom[ispp][b]) * (y_ratio[ispp][isys][b2]-y_nom[ispp][b2]);
	}
	y_systCorr[ispp][idx] /= (_nMetafitSyst-1);
	y_systCorr[ispp][idx] /= y_systErr[ispp][b];
	y_systCorr[ispp][idx] /= y_systErr[ispp][b2];
	cout<<"coorelation factor of metafit syst between bins "<<b<<" and "<<b2<<" = "<<y_systCorr[ispp][idx]<<endl;
	idx += 1;
      }
    }

    TFile * outf = new TFile("corrected_yields.root","UPDATE");
    //  outf->WriteObject(&y_nom[ispp],"FinalCorrectedYield"+(TString)((ispp==1)?"_pp":"_PbPb"));
    outf->WriteObject(&y_systErr[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitRelSystError"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    outf->WriteObject(&y_systErr_maxDev[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitRelSystErrorMaxDeviation"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    outf->WriteObject(&y_systCorr[ispp],(TString)((ispp<2)?"CorrectedYields":"RAA")+"_MetafitSyst_LinearizedCorrelationMatrix"+(TString)((ispp==1)?"pp":((ispp==0)?"PbPb":"")));
    //  outf->WriteObject(&y_fitErr,"FinalCorrectedYield_fitError"+(TString)((ispp==1)?"_pp":"_PbPb"));
    outf->Close();
  }

}
