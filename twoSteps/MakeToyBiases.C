#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TRandom.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "../helpers/Definitions.h"
#include "../helpers/Cuts.h"

TF1* f1, f2;

//to find the x_LW_i value satisfying f(x_LW_i) = integralOnBin_i(f)
//h2 is the histogram containing the bin integrals
//Assumes that h is decreasing and monotonous
vector<double> GetLWabciss(TH1F* h, TH1F* h2){
  int nB = h2->GetNbinsX(); // =2
  int iB = 1;
  vector<double> res;
  for(int bin=1;bin<=h->GetNbinsX();bin++){
    if(iB>nB) break;
    if(h->GetBinContent(bin)<h2->GetBinContent(iB)){ //found the crossing point with the integral
      res.push_back(h->GetBinLowEdge(bin));
      iB +=1;
    }
  }
  if(res.size()!=nB) cout<<"ERROR in GetLWabciss: did not find all the crossings"<<endl;
  return res;
}

//version where we calculate the integral
vector<double> GetLWabciss(TF1* f, vector<double> xbins){
  int nB = xbins.size()-1; // =2
  vector<double> integ;
  vector<double> res;
  int iB = 1;
  int nEval = 500;
  
  for(int i=1;i<=nB;i++){
    integ.push_back(f->Integral(xbins[i-1], xbins[i]));
    integ[i-1] /= xbins[i]-xbins[i-1];
  }

  for(int bin=0;bin<nEval;bin++){
    if(iB>nB) break;
    double xeval = xbins[0] + ((float)bin/nEval) * (xbins[nB]-xbins[0]);
    double ev = f->Eval( xeval );
    if(ev<integ[iB-1]){ //found the crossing point with the integral
      res.push_back(xeval);
      iB +=1;
    }
  }
  if(res.size()!=nB) cout<<"ERROR in GetLWabciss: did not find all the crossings"<<endl;
  return res;
}

TGraphAsymmErrors* CloneAndDivide(TGraphAsymmErrors* gr, TH1F* h, TString name){
  
  TGraphAsymmErrors* res = (TGraphAsymmErrors*)gr->Clone(name);
  for(int b=0;b<_NanaBins;b++){
    res->SetPoint(b, gr->GetX()[b] , 
		  (h->GetBinContent(b+1)!=0)?( gr->GetY()[b] / h->GetBinContent(b+1) ):1 );
    res->SetPointError(b, gr->GetEXlow()[b] , gr->GetEXhigh()[b] , 
		       (h->GetBinContent(b+1)!=0)?( gr->GetEYlow()[b] / h->GetBinContent(b+1) ):1 ,
		       (h->GetBinContent(b+1)!=0)?( gr->GetEYhigh()[b] / h->GetBinContent(b+1) ):1 );
  }
  
  res->GetHistogram()->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);
  res->GetHistogram()->GetYaxis()->SetTitle("ratio to MC fit");
  res->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  res->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  res->GetHistogram()->GetXaxis()->SetMoreLogLabels();
  res->GetHistogram()->GetYaxis()->SetMoreLogLabels();
  return res;
}

TGraph* CloneAndDivide(TGraph* gr, TH1F* h, TString name){
  
  TGraph* res = (TGraph*)gr->Clone(name);
  for(int b=0;b<_NanaBins;b++){
    res->SetPoint(b, gr->GetX()[b] , 
		  (h->GetBinContent(b+1)!=0)?( gr->GetY()[b] / h->GetBinContent(b+1) ):1 );
  }
  
  return res;
}

TH1F* CloneAndDivide(TF1* fct, TH1F* h, TString name){
  
  TH1F* res = (TH1F*)h->Clone(name);
  //  res->SetDirectory(0);
  for(int b=1;b<=h->GetNbinsX();b++){
    res->SetBinContent(b, 
		       (h->GetBinContent(b)!=0)?( fct->Eval(h->GetBinCenter(b)) / h->GetBinContent(b) ):1 );
    res->SetBinError(b, 0);
  }

  res->SetLineWidth(1);
  res->GetListOfFunctions()->Remove(res->GetFunction("bias_var"+(TString)(to_string(_biasNmeth-1)))); //remove fit function so it's not displayed
  return res;
}

TH1F* CloneAndDivide(TF1* fct, TF1* denom, TString name, int npts=700){
  
  TH1F *res = new TH1F(name,name+";p_{T} [GeV];MC ratio biased/default",npts,_BcPtmin[0],_BcPtmax[0]);
  for(int b=1;b<=npts;b++){
    res->SetBinContent(b,fct->Eval(res->GetBinCenter(b))/denom->Eval(res->GetBinCenter(b)));
  }
  return res;
}

int FindParamsPower(double *pars, double *y, vector<double> xbins, double *parLoLim, double *parHiLim, int verbose=1){
  ROOT::Math::MultiRootFinder r(0);
  r.SetPrintLevel(verbose);

  for(int i=0;i<xbins.size()-1;i++){
    if(parLoLim[i]==0 || parHiLim[i]==0){
      cout<<"ERROR: forbidden to set parameter limits = 0 in FindParamsPower!! Return now."<<endl;
      return 1;}
  }

  //Parameters [0]-[4]: 3 bin limits , then 2 y values
  //y=n and x=K(constant) are the function variables
  //Max functions are to implement parameter limits in the search. 1e3 penalty is a bit arbitrary
  TF2 * f1 = new TF2("f1",TString::Format("([0]**(y+1) - [1]**(y+1) - [2]*(y+1)/x) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) )",
					  parHiLim[0], parHiLim[0], parLoLim[0], parLoLim[0], parHiLim[1], parHiLim[1], parLoLim[1], parLoLim[1]));
  TF2 * f2 = new TF2("f2",TString::Format("([0]**(y+1) - [1]**(y+1) - [2]*(y+1)/x) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) )",
					  parHiLim[0], parHiLim[0], parLoLim[0], parLoLim[0], parHiLim[1], parHiLim[1], parLoLim[1], parLoLim[1]));
  f1->SetParameters(xbins[1], xbins[0], y[0]);
  f2->SetParameters(xbins[2], xbins[1], y[1]);

  // wrap the functions
  ROOT::Math::WrappedMultiTF1 g1(*f1,2);
  ROOT::Math::WrappedMultiTF1 g2(*f2,2);
  r.AddFunction(g1);
  r.AddFunction(g2);

  r.Solve(pars,1000,1e-4,1e-6);
  f1->Delete();
  f2->Delete();

  const double *x = r.X();
  for(int i=0;i<xbins.size()-1;i++) pars[i] = x[i];
  return r.Status();
}

double PowerToLogIntegral(double *x, double *p) {
  TF1* integrand = new TF1("integrand","[0]*(x**( [1]+[2]*TMath::Log(x) ))",_BcPtmin[0],_BcPtmax[0]);
  int fix = (int)(p[4]+0.01);
  integrand->SetParameter(0,(fix==0)?p[0]:x[0]);
  integrand->SetParameter(1,(fix==1)?p[0]:x[(fix==2)?1:0]);
  integrand->SetParameter(2,(fix==2)?p[0]:x[1]);

  double res = integrand->Integral(p[1],p[2]) - p[3];
  integrand->Delete();
  double parlimsPen = TMath::Max(1., (x[0]-p[5]) * 1e3/fabs(p[5])) //implementing limit parHiLim[0]
                    * TMath::Max(1., (-x[0]+p[6]) * 1e3/fabs(p[6])) //implementing limit parLoLim[0]
                    * TMath::Max(1., (x[1]-p[7]) * 1e3/fabs(p[7])) //implementing limit parHiLim[1]
                    * TMath::Max(1., (-x[1]+p[8]) * 1e3/fabs(p[8])); //implementing limit parLoLim[1]
  // if(parlimsPen-1 > 1e-6) {cout<<"Give penalty to limit parameters = "<<parlimsPen<<". PowerToLogIntegral evaluated at x0,x1 = "<<x[0]<<" "<<x[1]<<endl; 
  //   cout<<"separate penalties for parHiLim[0],parLoLim[0],parHiLim[1],parLoLim[1] = "<<TMath::Max(1., (x[0]-p[5]) * 1e3/fabs(p[5]))<<" "<<TMath::Max(1., (-x[0]+p[6]) * 1e3/fabs(p[6]))<<" "<<TMath::Max(1., (x[1]-p[7]) * 1e3/fabs(p[7]))<<" "<<TMath::Max(1., (-x[1]+p[8]) * 1e3/fabs(p[8]))<<endl;
  //   cout<<"fixed par = "<<p[0]<<endl;
  //   cout<<"parHiLim[0],parLoLim[0],parHiLim[1],parLoLim[1] = "<<p[5]<<" "<<p[6]<<" "<<p[7]<<" "<<p[8]<<endl;
  // }
  return res * parlimsPen;

}

int FindParamsPowerToLog(double *pars, double *y, vector<double> xbins, double fix, double fixedPar, double *parLoLim, double *parHiLim, int verbose=1){
  ROOT::Math::MultiRootFinder r(0);
  r.SetPrintLevel(verbose);

  //Parameters [0]-[4]: 3 bin limits , then 2 y values
  TF1 * f1 = new TF1("f1",PowerToLogIntegral,0,100,9); //0,100 are dummy borns
  TF1 * f2 = new TF1("f2",PowerToLogIntegral,0,100,9); //0,100 are dummy borns
  f1->SetParameters(fixedPar, xbins[0], xbins[1], y[0], fix, parHiLim[0], parLoLim[0], parHiLim[1], parLoLim[1]);
  f2->SetParameters(fixedPar, xbins[1], xbins[2], y[1], fix, parHiLim[0], parLoLim[0], parHiLim[1], parLoLim[1]);

  // wrap the functions
  ROOT::Math::WrappedMultiTF1 g1(*f1,2);
  ROOT::Math::WrappedMultiTF1 g2(*f2,2);
  r.AddFunction(g1);
  r.AddFunction(g2);

  // std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
  // if(verbose==0){
  //   std::ofstream   fout("/dev/null");
  //   std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
  // }
  r.Solve(pars,1000,1e-4,1e-6);
  //  if(verbose==0) std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
  f1->Delete();
  f2->Delete();

  const double *x = r.X();
  for(int i=0;i<xbins.size()-1;i++) pars[i] = x[i];
  return r.Status();
}

void MakeToyBiases(bool firstStep=true, bool plotOnly=false){

  gStyle->SetOptStat(0);
  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);
  TVirtualFitter::SetMaxIterations( 10000 );//default is 5000
  gRandom = new TRandom3(234); //some seed give one or two doubly-failed fits, this one does not

  if(plotOnly) _biasNtoys = 40;
  bool preAEmeanPt = true; // use the mean pT in acceptance uncut sample rather than of observed events
  bool useLWabciss = true; //Use Lafferty and Wyatt prescription for position of the bin centers

  //pT distribution of gen signal MC (acceptance sample), before any cuts = original pT spectrum
  cout<<"Getting pT distribution from original MC..."<<endl;
  TFile *fileMC = TFile::Open("/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/acceptance/BcToJpsiMuNu_BCVEGPY_PYTHIA8_pp5TeV_RunIIpp5Spring18DR-00093_acceptance_14092020_ONIATREE_3640k.root");
  TTree* T = (TTree*)fileMC->Get("hionia/myTree");
  int nentries = T->GetEntries();

  Short_t Gen_Bc_size;
  T->SetBranchAddress("Gen_Bc_size", &Gen_Bc_size);
  TClonesArray *Gen_3mu_4mom = new TClonesArray(); TBranch *b_Gen_3mu_4mom;
  T->SetBranchAddress("Gen_3mu_4mom", &Gen_3mu_4mom, &b_Gen_3mu_4mom);

  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("Gen_Bc_size",1);
  T->SetBranchStatus("Gen_3mu_4mom",1);

  int nbins = 150;
  auto h_pTMC_cont = new TH1F("h_pTMC_cont","h_pTMC_cont;p_{T} [GeV];#frac{dN_{MC}}{dp_{T} dy}",nbins,_BcPtmin[0],_BcPtmax[0]);
  float bwpt_cont = (_BcPtmax[0]-_BcPtmin[0])/nbins;
  float anabins[] = {_BcPtmin[1], _BcPtmin[2],_BcPtmax[2]};
  auto h_pTMC = new TH1F("h_pTMC","h_pTMC",_NanaBins,anabins);

  float ntot=0;
  vector<float> nbin = vector<float>(_NanaBins, 0);
  vector<float> ptSum = vector<float>(_NanaBins, 0);

  for(int j=0; j<nentries; j++){
    Gen_3mu_4mom->Clear();
    T->GetEntry(j);
    for(int iBc=0;iBc<Gen_Bc_size;iBc++){
      TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom->At(iBc);
      float pt = gen3mu->Pt();
      float rap = gen3mu->Rapidity();
      float mass = gen3mu->M();
      if(mass<_mBcMin || mass>_mBcMax) continue; //exclude the events excluded in the template fit
      
      //correct by pt and Y bin width
      float bwpt_2, bwY;
      bool infid=false;
      ntot+=1;
      for(int b=0;b<_NanaBins;b++){
	if(inFidCuts(b+1, pt, rap)){
	  bwpt_2 = _BcPtmax[b+1]-_BcPtmin[b+1];
	  bwY = _BcYmax[b+1]-_BcYmin[b+1];
	  ptSum[b] += pt;
	  nbin[b] += 1;
	  infid=true;
	}
      }
      if(!infid) continue;

      h_pTMC_cont->Fill(pt, 1/(bwpt_cont*bwY));
      h_pTMC->Fill(pt, 1/(bwpt_2*bwY));
    }
  }

  float avPt_preAE[] = {ptSum[0]/nbin[0], ptSum[1]/nbin[1]};
  vector<double> x_LW = GetLWabciss(h_pTMC_cont,h_pTMC);

  cout<<"in acc sample, ntot, nbin1, nbin2, fraction passing bin1 or bin2 = "<<ntot<<" "<<nbin[0]<<" "<<nbin[1]<<" "<<nbin[0]/ntot<<" "<<nbin[1]/ntot<<" "<<endl;
  cout<<"average pT in acceptance MC sample, bin1 and bin2 = "<<avPt_preAE[0]<<" "<<avPt_preAE[1]<<endl;
  cout<<"LW abcisses from MC spectrum, common to pp and PbPb, bin1 and 2 = "<<x_LW[0]<<" "<<x_LW[1]<<endl;

  TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
  c2->SetLogy();
  c2->SetLogx();
  h_pTMC_cont->SetLineWidth(2);
  h_pTMC->SetLineWidth(2);
  h_pTMC->SetLineColor(kBlack);
  h_pTMC_cont->GetXaxis()->SetMoreLogLabels();
  h_pTMC_cont->Draw("hist");
  h_pTMC->Draw("histsame");
  c2->SaveAs("figs/pT_originalMC.pdf");  

  //scaling to expected corrected yield in data (prefit)
  float tmp_scale = 1.17; //norm of MC seems a bit too low, temporary fix
  TH1F *h_pTMC_cont_pp = (TH1F*)h_pTMC_cont->Clone("h_pTMC_cont_pp");
  TH1F *h_pTMC_pp = (TH1F*)h_pTMC->Clone("h_pTMC_pp");
  h_pTMC_cont_pp->Scale(tmp_scale*_scaleMCsig[true] * 3000000 / nentries);
  h_pTMC_pp->Scale(tmp_scale*_scaleMCsig[true] * 3000000 / nentries);
  TH1F *h_pTMC_cont_PbPb = (TH1F*)h_pTMC_cont->Clone("h_pTMC_cont_PbPb");
  TH1F *h_pTMC_PbPb = (TH1F*)h_pTMC->Clone("h_pTMC_PbPb");
  h_pTMC_cont_PbPb->Scale(tmp_scale*_scaleMCsig[false] * Ncoll_MB * 4200000 / nentries);
  h_pTMC_PbPb->Scale(tmp_scale*_scaleMCsig[false] * Ncoll_MB * 4200000 / nentries);

  //Grab nominal acc and eff, without second-step pT biasing of MC, and postfit uncorrected yields
  cout<<"Grabbing info of all acc, eff, yields, ..."<<endl;
  //acc
  TFile *facc = new TFile("../acceptance/acceptanceMap.root","READ");
  vector<float> *acc_oneBinned;
  facc->GetObject("acceptance_oneBinned", acc_oneBinned);
  TFile *feff = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
  //eff
  vector<vector<float> > *eff_oneBinned_pp;
  feff->GetObject("efficiency_oneBinned_pp", eff_oneBinned_pp);
  vector<vector<float> > *eff_oneBinned_PbPb;
  feff->GetObject("efficiency_oneBinned_PbPb", eff_oneBinned_PbPb);

  //postfit (and prefit) yields
  TFile *infile = new TFile("../AccEffCorr/corrected_yields.root","READ");
  vector<vector<vector<float> > > *Yields_postfit_pp;
  infile->GetObject("Yields_postfit_pp", Yields_postfit_pp);
  vector<vector<vector<float> > > *Yields_postfit_PbPb;
  infile->GetObject("Yields_postfit_PbPb", Yields_postfit_PbPb);
  vector<vector<vector<float> > > *Yields_prefit_pp;
  infile->GetObject("Yields_prefit_pp", Yields_prefit_pp);
  vector<vector<vector<float> > > *Yields_prefit_PbPb;
  infile->GetObject("Yields_prefit_PbPb", Yields_prefit_PbPb);
  //nominal corrected yield
  vector<float> *y_nom_pp;
  infile->GetObject("FinalCorrectedYield"+(TString)(firstStep?"_1stStep":"_2ndStep")+"_pp", y_nom_pp);
  vector<float> *y_nom_PbPb;
  infile->GetObject("FinalCorrectedYield"+(TString)(firstStep?"_1stStep":"_2ndStep")+"_PbPb", y_nom_PbPb);

  //pp and PbPb corrected yields
  //  vector<vector<double> > yield = vector<vector<double> >(2, vector<double>(_NanaBins,0));
  vector<vector<double> > ycorr = vector<vector<double> >(2, vector<double>(_NanaBins,0));
  vector<double> ptwidth, Ywidth;
  for(int b=0;b<_NanaBins;b++){
    //yield[0][b] = (*Yields_postfit_pp)[0][b+1][0];
    //yield[1][b] = (*Yields_postfit_PbPb)[0][b+1][0];

    ptwidth.push_back( _BcPtmax[b+1]-_BcPtmin[b+1] );
    Ywidth.push_back( _BcYmax[b+1]-_BcYmin[b+1] );
    ycorr[0][b] = (*y_nom_pp)[b] / ( ptwidth[b] * Ywidth[b]); //yield[0][b] / ((*acc_oneBinned)[b+1] * (*eff_oneBinned_pp)[b+1][0]
    ycorr[1][b] = (*y_nom_PbPb)[b] / ( ptwidth[b] * Ywidth[b]);//yield[1][b] / ((*acc_oneBinned)[b+1] * (*eff_oneBinned_PbPb)[b+1][0]
  }

  //Grab signal normalisations
  vector<float> *rsig_pp;//[pt bin+1]
  infile->GetObject("rsig_pp",rsig_pp);
  vector<vector<float> > *rrelerr_pp; //[pt bin+1][sym err, lo err, hi err]
  infile->GetObject("rsig_relerr_pp",rrelerr_pp);
  cout<<"rrelerr_pp low bin1 bin2 = "<< (*rrelerr_pp)[1][1]<<" "<<(*rrelerr_pp)[2][1]<<endl;
  vector<float> *metafitErrCorr_pp;
  infile->GetObject("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixpp", metafitErrCorr_pp);
  vector<float> *r1r2Corr_pp;
  infile->GetObject("r1r2Correlation_pp", r1r2Corr_pp);

  vector<float> * rsig_PbPb;
  infile->GetObject("rsig_PbPb",rsig_PbPb);
  vector<vector<float> > * rrelerr_PbPb; //[pt bin+1][sym err, lo err, hi err]
  infile->GetObject("rsig_relerr_PbPb",rrelerr_PbPb);
  cout<<"rrelerr_PbPb low bin1 bin2 = "<< (*rrelerr_PbPb)[1][1]<<" "<<(*rrelerr_PbPb)[2][1]<<endl;
  vector<float> *metafitErrCorr_PbPb;
  infile->GetObject("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixPbPb", metafitErrCorr_PbPb);
  vector<float> *r1r2Corr_PbPb;
  infile->GetObject("r1r2Correlation_PbPb", r1r2Corr_PbPb);

  vector<float> *y_metafitRelErrLo_pp;
  vector<float> *y_metafitRelErrHi_pp;
  infile->GetObject("CorrectedYields_MetafitRelSystErrorLo_pp", y_metafitRelErrLo_pp);
  infile->GetObject("CorrectedYields_MetafitRelSystErrorHi_pp", y_metafitRelErrHi_pp); 
  vector<float> *y_metafitRelErrLo_PbPb;
  vector<float> *y_metafitRelErrHi_PbPb;
  infile->GetObject("CorrectedYields_MetafitRelSystErrorLo_PbPb", y_metafitRelErrLo_PbPb);
  infile->GetObject("CorrectedYields_MetafitRelSystErrorHi_PbPb", y_metafitRelErrHi_PbPb); 

  cout<<"Metafit relerr_pp low bin1 bin2 = "<< (*y_metafitRelErrLo_pp)[0]<<" "<<(*y_metafitRelErrLo_pp)[1]<<endl;
  cout<<"Metafit relerr_PbPb low bin1 bin2 = "<< (*y_metafitRelErrLo_PbPb)[0]<<" "<<(*y_metafitRelErrLo_PbPb)[1]<<endl;
 
  vector<vector<float> > rSig(2, vector<float>(_NanaBins,0) );//[pp or PbPb][pt bins from 0]
  vector<vector<float> > rErrHi(2, vector<float>(_NanaBins,0));//[pp or PbPb][pt bins from 0]
  vector<vector<float> > rErrLo(2, vector<float>(_NanaBins,0));
  vector<vector<float> > corr(2, vector<float>((int)(_NanaBins*(_NanaBins-1)/2), 0));
  //errors on corrected yield
  vector<vector<double> > yErrHi(2, vector<double>(_NanaBins,0));//[pp or PbPb][pt bins from 0]
  vector<vector<double> > yErrLo(2, vector<double>(_NanaBins,0));
  
  //vector<vector<float> > ptAverage
  //ptAverage[col][trueb]
  for(int b=0;b<_NanaBins;b++){
    cout<<"bin = "<<b+1<<endl;
    rSig[0][b] = (*rsig_pp)[b+1];
    cout<<"rSig[0][b] = "<<rSig[0][b]<<endl;
    rErrLo[0][b] = rSig[0][b] * sqrt(pow((*rrelerr_pp)[b+1][1],2) + pow((*y_metafitRelErrLo_pp)[b],2));
    cout<<"rErrLo[0][b] = "<<rErrLo[0][b]<<endl;
    rErrHi[0][b] = rSig[0][b] * sqrt(pow((*rrelerr_pp)[b+1][2],2) + pow((*y_metafitRelErrHi_pp)[b],2));
    yErrLo[0][b] = ycorr[0][b] * sqrt(pow((*rrelerr_pp)[b+1][1],2) + pow((*y_metafitRelErrLo_pp)[b],2));
    cout<<"yErrLo[0][b] = "<<yErrLo[0][b]<<endl;
    cout<<"relerrlo = "<<sqrt(pow((*rrelerr_pp)[b+1][1],2) + pow((*y_metafitRelErrLo_pp)[b],2))<<endl;
    yErrHi[0][b] = ycorr[0][b] * sqrt(pow((*rrelerr_pp)[b+1][2],2) + pow((*y_metafitRelErrHi_pp)[b],2));
    cout<<"corrected yield pp   bin "<<b+1<<" = "<<ycorr[0][b]<<endl;
    cout<<"corrected yield pp (MC)   bin "<<b+1<<" = "<<h_pTMC_pp->GetBinContent(b+1)<<endl;
    
    rSig[1][b] = (*rsig_PbPb)[b+1];
    cout<<"rSig[1][b] = "<<rSig[1][b]<<endl;
    rErrLo[1][b] = rSig[1][b] * sqrt(pow((*rrelerr_PbPb)[b+1][1],2) + pow((*y_metafitRelErrLo_PbPb)[b],2));
    cout<<"rErrLo[1][b] = "<<rErrLo[1][b]<<endl;
    rErrHi[1][b] = rSig[1][b] * sqrt(pow((*rrelerr_PbPb)[b+1][2],2) + pow((*y_metafitRelErrHi_PbPb)[b],2));
    yErrLo[1][b] = ycorr[1][b] * sqrt(pow((*rrelerr_PbPb)[b+1][1],2) + pow((*y_metafitRelErrLo_PbPb)[b],2));
    cout<<"yErrLo[1][b] = "<<yErrLo[1][b]<<endl;
    cout<<"relerrlo = "<<sqrt(pow((*rrelerr_PbPb)[b+1][1],2) + pow((*y_metafitRelErrLo_PbPb)[b],2))<<endl;
    yErrHi[1][b] = ycorr[1][b] * sqrt(pow((*rrelerr_PbPb)[b+1][2],2) + pow((*y_metafitRelErrHi_PbPb)[b],2));
    cout<<"corrected yield PbPb bin "<<b+1<<" = "<<ycorr[1][b]<<endl;
    cout<<"corrected yield PbPb (MC) bin "<<b+1<<" = "<<h_pTMC_PbPb->GetBinContent(b+1)<<endl;
  }

  corr[0][0] = ( (*r1r2Corr_pp)[0] * (*rrelerr_pp)[1][0] * (*rrelerr_pp)[2][0] //symmetric errors for correlations
		 +(*metafitErrCorr_pp)[0] * (((*y_metafitRelErrLo_pp)[0]+(*y_metafitRelErrHi_pp)[0])/2) * (((*y_metafitRelErrLo_pp)[1]+(*y_metafitRelErrHi_pp)[1])/2) 
		 )/ sqrt( (pow((*rrelerr_pp)[1][0],2) + pow(((*y_metafitRelErrLo_pp)[0]+(*y_metafitRelErrHi_pp)[0])/2,2)) * (pow((*rrelerr_pp)[2][0],2) + pow(((*y_metafitRelErrLo_pp)[1]+(*y_metafitRelErrHi_pp)[1])/2,2)) ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  corr[1][0] = ( (*r1r2Corr_PbPb)[0] * (*rrelerr_PbPb)[1][0] * (*rrelerr_PbPb)[2][0] //symmetric errors for correlations
		 +(*metafitErrCorr_PbPb)[0] * (((*y_metafitRelErrLo_PbPb)[0]+(*y_metafitRelErrHi_PbPb)[0])/2) * (((*y_metafitRelErrLo_PbPb)[1]+(*y_metafitRelErrHi_PbPb)[1])/2)
		 )/ sqrt( (pow((*rrelerr_PbPb)[1][0],2) + pow(((*y_metafitRelErrLo_PbPb)[0]+(*y_metafitRelErrHi_PbPb)[0])/2,2)) * (pow((*rrelerr_PbPb)[2][0],2) + pow(((*y_metafitRelErrLo_PbPb)[1]+(*y_metafitRelErrHi_PbPb)[1])/2,2)) ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))

  vector<vector<vector<double> > > dipoints = vector<vector<vector<double> > >(2, vector<vector<double> >(_biasNmeth*(_biasNtoys+1), vector<double>(4,0))); //[{x1,x2,r1,r2}] 
  vector<vector<vector<TF1*> > > mcfit = vector<vector<vector<TF1*> > >(2,vector<vector<TF1*> >(_biasNmeth, vector<TF1*>(2)));
  vector<vector<TF1*> > bias = vector<vector<TF1*> >(2,vector<TF1*>(_biasNmeth*(_biasNtoys+1)));
  vector<vector<TH1F*> > bias_r = vector<vector<TH1F*> >(2,vector<TH1F*>(_biasNmeth*(_biasNtoys+1)));
  vector<vector<TGraph*> > toysDipts = vector<vector<TGraph*> >(2 , vector<TGraph*>(_biasNmeth*(_biasNtoys+1)));
  vector<vector<TGraph*> > toysDipts_r = vector<vector<TGraph*> >(2 , vector<TGraph*>(_biasNmeth*(_biasNtoys+1)));
  vector<Color_t> color = vector<Color_t>(_biasNmeth*(_biasNtoys+1));
  vector<double> x_LW_MC;
  vector<vector<vector<double> > > x_LW_2ndstep = vector<vector<vector<double> > >(2, vector<vector<double> >(_biasNmeth));

  vector<double> xbins; xbins.push_back(_BcPtmin[1]);
  for(int b=1;b<=_NanaBins;b++) xbins.push_back(_BcPtmax[b]);

  cout<<"Running "<<_biasNtoys<<" toy biases for each of the "<<_biasNmeth<<" methods"<<endl;
  for(int col=0;col<2;col++){//pp or PbPb

    double yMC[] = {((col==0)?h_pTMC_pp:h_pTMC_PbPb)->GetBinContent(1) * ptwidth[0],
		    ((col==0)?h_pTMC_pp:h_pTMC_PbPb)->GetBinContent(2) * ptwidth[1]};

    for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
      //cout<<"variation v="<<v<<endl;
    failedVar:

      //make toys for dipoints
      double z1 = gRandom->Gaus(),z2 = gRandom->Gaus();
      if(v<3) {z1=0; z2=0;}
      double z22 = ( corr[col][0]*z1+sqrt(1-pow(corr[col][0],2))*z2 );
      double r1 = ycorr[col][0] + z1 * ((z1>0)?yErrHi:yErrLo)[col][0];
      double r2 = ycorr[col][1] + z22 * ((z22>0)?yErrHi:yErrLo)[col][1];
      if(r1<0.1*ycorr[col][0]) {r1=0.1*ycorr[col][0]; cout<<"using minimal r1"<<endl;}
      if(r2<0.1*ycorr[col][1]) {r2=0.1*ycorr[col][1]; cout<<"using minimal r2"<<endl;}
      if(v<_biasNmeth){
	x_LW_MC.push_back( useLWabciss?x_LW[0]:(preAEmeanPt?avPt_preAE[0]:ptAverage[col][0]) );
	x_LW_MC.push_back( useLWabciss?x_LW[1]:(preAEmeanPt?avPt_preAE[1]:ptAverage[col][1]) );
      }

      double xt[] = {x_LW_MC[0], x_LW_MC[1]};
      dipoints[col][v] = vector<double>{xt[0],xt[1],r1,r2};

      double yt[] = {r1, r2};
      toysDipts[col][v] = new TGraph(2,xt,yt);
      toysDipts[col][v]->SetMarkerSize(0.9);
      toysDipts[col][v]->SetMarkerStyle(20);

      //"fit" the dipoints
      for(int m=0;m<_biasNmeth;m++){
	if(v==m || (v>=_biasNmeth+m*_biasNtoys && v<_biasNmeth+(m+1)*_biasNtoys) ){//"_biasNtoys" number of toys per method
	  // if(m==0){//linear, limiting to 3sigma around a given point 
	  //   //from lowest pt edge to highest pt edge
	  //   // x1 = _BcPtmin[0];
	  //   // x2 = _BcPtmax[0];
	  //   bias[col][v] = new TF1("bias_var"+(TString)(to_string(v)),"max(min([0]*x + [1] , [2]) , [3])",_BcPtmin[0],_BcPtmax[0]); //[2] is the maximum value, [3] is the minimum
	  //   bias[col][v]->SetParameter(0, (r2-r1) / (x2-x1) );
	  //   bias[col][v]->SetParameter(1, r2 - bias[col][v]->GetParameter(0) * x2 );
	  //   bias[col][v]->SetParameter(2, max(r1*(1+2.5*yErrHi[col][0]/ycorr[col][0]) , r2*(1+2.5*yErrHi[col][1]/ycorr[col][1]))); //2.5 times the relative error is allowed
	  //   bias[col][v]->SetParameter(3, min(r1*(1-2.5*yErrLo[col][0]/ycorr[col][0]) , r2*(1-2.5*yErrLo[col][1]/ycorr[col][1])) );
	  //   bias[col][v]->SetLineColor(kRed);
	  //   toysDipts[col][v]->SetMarkerColor(kRed);
	  // }
	  // else if(m==1){ //quadratic, passing through average pt x-points and highest pt edge
	  //   double x3 = _BcPtmax[0];
	  //   bias[col][v] = new TF1("bias_var"+(TString)(to_string(v)),"[0]*x*x+[1]*x+[2]",_BcPtmin[0],_BcPtmax[0]);
	  //   bias[col][v]->SetParameter(0, (r1-r2) / ((x1-x2)*(x1-x3)) );
	  //   bias[col][v]->SetParameter(1, - bias[col][v]->GetParameter(0) * (x2+x3) );
	  //   bias[col][v]->SetParameter(2, r1 - bias[col][v]->GetParameter(0) * x1*x1 - bias[col][v]->GetParameter(1) * x1 );
	  //   bias[col][v]->SetLineColor(kOrange+10);
	  //   toysDipts[col][v]->SetMarkerColor(kOrange+10);
	  // }
	  if(m==0){ //power (linear in log-log)
	    
	    if(v==m){ 
	      cout<<"method "<<m<<endl;
	      cout<<"MC fit"<<endl;
	      //Exactly fit the MC
	      mcfit[col][m][0] = new TF1("MC_fit_meth"+(TString)(to_string(m)), "[0]*(x**[1])",_BcPtmin[0],_BcPtmax[0]);//(TF1*) bias[col][v]->Clone("MC_fit_meth"+(TString)(to_string(m)));
	      mcfit[col][m][0]->SetParameter(1, -4);//TMath::Log(r2/r1)/TMath::Log(x2/x1) );
	      mcfit[col][m][0]->SetParameter(0, 1e7);//r1*pow(x1, -bias[col][v]->GetParameter(1)));
	      ((col==0)?h_pTMC_cont_pp:h_pTMC_cont_PbPb)->Fit("MC_fit_meth"+(TString)(to_string(m)));

	      //"fit" the MC from the values of the integral on each bin range
	      cout<<"MC integral fit"<<endl;
	      mcfit[col][m][1] = (TF1*)mcfit[col][m][0]->Clone("MC_fitWithBinIntegral_meth"+(TString)(to_string(m)));
	      double parsMC[] = {mcfit[col][m][1]->GetParameter(0), mcfit[col][m][1]->GetParameter(1)};
	      double parLowLimsMC[] = {5e3, -8};
	      double parHiLimsMC[] = {1e9, -1.5};
	      int status = FindParamsPower(parsMC, yMC, xbins, parLowLimsMC, parHiLimsMC, 1);
	      if (status!=0) cout<<"ROOT FINDING FAILED!! status "<<status<<" (power)"<<endl;
	      mcfit[col][m][1]->SetParameters(parsMC[0],parsMC[1]);
	    }
	    else if(v==_biasNmeth+m*_biasNtoys) cout<<"Toys for method "<<m<<endl;

	    //basic function for biasing in a given variation
	    bias[col][v] = (TF1*)mcfit[col][m][0]->Clone("bias_var"+(TString)(to_string(v)));
	    color[v] = (v==m)?(kGreen+3):kGreen;

	    //"fit" the central value, from the integral on each bin range
	    double pars[] = {bias[col][v]->GetParameter(0), bias[col][v]->GetParameter(1)};
	    double parLowLims[] = {1e3, -8};
	    double parHiLims[] = {1e10, -1};
	    double y[] = {r1*ptwidth[0], r2*ptwidth[1]};
	    int status = FindParamsPower(pars, y, xbins, parLowLims, parHiLims, 0);
	    if (status!=0) {
	      cout<<"ROOT FINDING FAILED !! (power). RETRYING now (SUCCESS IF NO NEWS). y1, y2 = "<<r1<<" "<<r2<<endl;
	      pars[0] = 2*bias[col][v]->GetParameter(0); pars[1] = 0.5*bias[col][v]->GetParameter(1);
	      parLowLims[0] = 1e2; parLowLims[1] = -10;
	      parHiLims[0] = 1e12; parHiLims[1] = -0.2;
	      int status2 = FindParamsPower(pars, y, xbins, parLowLims, parHiLims, 0);
	      if(status2!=0){
		cout<<"FAILED AGAIN !! Giving up, resetting the variation"<<endl;
		if(v<_biasNmeth) {
		  bias[col][v]->Delete();
		  goto failedVar;}
	      }
	    }
	    bias[col][v]->SetParameters(pars[0],pars[1]);

	  }
	  if(m==1 || m==2){ //x^{n+m*log(x)} (quadratic in log-log scale)

	    double fix = 3-m;
	    double fixedPar;
	    
	    if(v==m){ 
	      cout<<"method "<<m<<endl;
	      cout<<"MC fit"<<endl;
	      //Exactly fit the MC //only needed once for the powerToLog function
	      if(m==1){
		mcfit[col][m][0] = new TF1("MC_fit_meth"+(TString)to_string(m), "[0]*(x**( [1]+[2]*TMath::Log(x) ))",_BcPtmin[0],_BcPtmax[0]);//(TF1*) bias[col][v]->Clone("MC_fit_meth"+(TString)(to_string(m)));
		mcfit[col][m][0]->SetParameter(1, 2);//TMath::Log(r2/r1)/TMath::Log(x2/x1) );
		mcfit[col][m][0]->SetParameter(0, 1e4);//r1*pow(x1, -bias[col][v]->GetParameter(1)));
		mcfit[col][m][0]->SetParameter(2, -1);
		((col==0)?h_pTMC_cont_pp:h_pTMC_cont_PbPb)->Fit("MC_fit_meth"+(TString)(to_string(m)));
	      } 
	      else mcfit[col][m][0] = (TF1*)mcfit[col][1][0]->Clone("MC_fit_meth"+(TString)(to_string(m)));
	      
	      //"fit" the MC from the values of the integral on each bin range
	      cout<<"MC integral fit"<<endl;
	      mcfit[col][m][1] = (TF1*)mcfit[col][m][0]->Clone("MC_fitWithBinIntegral_meth"+(TString)(to_string(m)));
	      fixedPar = mcfit[col][m][1]->GetParameter(fix);
	      double parsMC[] = {mcfit[col][m][1]->GetParameter(0), mcfit[col][m][1]->GetParameter(m)}; //parameter number m is varied
	      double parLowLimsMC[] = {1e1, (m==1)?(-0.6):(-5)};
	      double parHiLimsMC[] = {1e8, (m==1)?4.:(-0.01)};
	      int status = FindParamsPowerToLog(parsMC, yMC, xbins, fix, fixedPar, parLowLimsMC, parHiLimsMC, 1);
	      if(status!=0) cout<<"ROOT FINDING FAILED !! status "<<status<<" (power to log). y1, y2 = "<<r1<<" "<<r2<<endl;
	      mcfit[col][m][1]->SetParameters(parsMC[0],(m==1)?parsMC[1]:fixedPar,(m==2)?parsMC[1]:fixedPar);
	    }
	    else if(v==_biasNmeth+m*_biasNtoys) {
	      cout<<"Toys for method "<<m<<endl;
	      fixedPar = mcfit[col][m][1]->GetParameter(fix);
	    }

	    //basic function for biasing in a given variation
	    bias[col][v] = (TF1*) mcfit[col][m][0]->Clone("bias_var"+(TString)(to_string(v)));
	    color[v] = (v==m)?(kViolet+1):(kViolet-4);
	    if(m==2) color[v] = (v==m)?(kRed+2):(kOrange+6);
	    
	    //"fit" the central value, from the integral on each bin range
	    double pars[] = {bias[col][v]->GetParameter(0), bias[col][v]->GetParameter(m)}; //parameter number m is varied
	    double y[] = {r1*ptwidth[0], r2*ptwidth[1]};
	    double parLowLims[] = {1e1, (m==1)?(-0.9):(-5)};
	    double parHiLims[] = {1e8, (m==1)?4.5:(-0.01)};
	    int status = FindParamsPowerToLog(pars, y, xbins, fix, fixedPar, parLowLims, parHiLims, 0);
	    if (status!=0) {
	      cout<<"ROOT FINDING FAILED !! (power to log). RETRYING now (SUCCESS IF NO NEWS). y1, y2 = "<<r1<<" "<<r2<<endl;
	      pars[0] = 1.9*bias[col][v]->GetParameter(0); pars[1] = 0.5*bias[col][v]->GetParameter(1);
	      parLowLims[0] = 1e0; parLowLims[1] = (m==1)?(-3):(-8);
	      parHiLims[0] = 1e10; parHiLims[1] = (m==1)?7.:(0.5);
	      int status2 = FindParamsPowerToLog(pars, y, xbins, fix, fixedPar, parLowLims, parHiLims, 0);
	      if(status2!=0){
		cout<<"FAILED AGAIN !! Giving up, resetting the variation"<<endl;
		if(v<_biasNmeth) {
		  bias[col][v]->Delete();
		  goto failedVar;}
	      }
	    }
	    bias[col][v]->SetParameters(pars[0],(m==1)?pars[1]:fixedPar,(m==2)?pars[1]:fixedPar);

	  }
	  // else if(m==1){ //exp(-x)*x
	  //   double x3 = _BcPtmax[0];
	  //   bias[col][v] = new TF1("bias_var"+(TString)(to_string(v)),"[0]*exp(-[1]*x)/(x**[2])",_BcPtmin[0],_BcPtmax[0]);
	  //   float lnr21 = TMath::Log(r2/r1);
	  //   float lnx21 = TMath::Log(x2/x1);
	  //   float lnx23 = TMath::Log(x2/x3);
	  //   bias[col][v]->SetParameter(2, lnr21 / (lnx23*(x2-x1)/(x2-x3) - lnx21) );
	  //   bias[col][v]->SetParameter(1, - ( lnr21 + lnx21*bias[col][v]->GetParameter(2) ) / (x2-x1) );
	  //   bias[col][v]->SetParameter(0, pow(x1,bias[col][v]->GetParameter(2))*r1*TMath::Exp(bias[col][v]->GetParameter(1)*x1) );
	  //   bias[col][v]->SetLineColor(kOrange+1);
	  //   toysDipts[col][v]->SetMarkerColor(kOrange+1);
	  // }
	  // else if(m==3){ //exponential
	  //   bias[col][v] = new TF1("bias_var"+(TString)(to_string(v)),"[0]*exp([1]*x)",_BcPtmin[0],_BcPtmax[0]);
	  //   bias[col][v]->SetParameter(1, TMath::Log(r2/r1) / (x2-x1) );
	  //   bias[col][v]->SetParameter(0, r1*TMath::Exp(-x1 * bias[col][v]->GetParameter(1)) );
	  //   //cout<<"x1,x2,r1,r2,f(x1),f(x2)"<<x1<<" "<<x2<<" "<<r1<<" "<<r2<<" "<<bias[col][v]->Eval(x1)<<" "<<bias[col][v]->Eval(x2)<<endl;
	  //   bias[col][v]->SetLineColor(kRed+4);
	  //   toysDipts[col][v]->SetMarkerColor(kRed+4);
	  // }

	  if(v==m){
	    cout<<"values in 2 bins MC = "<<yMC[0]/ptwidth[0]<<" "<<yMC[1]/ptwidth[1]<<endl;
	    cout<<"integral/ptwidth of new function (MC \"fit\") on the first and 2nd bin = "<<mcfit[col][m][1]->Integral(xbins[0],xbins[1])/ptwidth[0]<<" "<<mcfit[col][m][1]->Integral(xbins[1],xbins[2])/ptwidth[1]<<endl;

	    mcfit[col][m][0]->SetLineColor(color[v]);
	    mcfit[col][m][1]->SetLineColor(color[v]);
	    mcfit[col][m][0]->SetLineStyle(3);
	    mcfit[col][m][1]->SetLineStyle(7);
	  }

	}
      }

      bias[col][v]->SetLineStyle((v<_biasNmeth)?1:4);
      bias[col][v]->SetLineColor(color[v]);
      toysDipts[col][v]->SetMarkerColor(color[v]);

      bias[col][v]->SetLineWidth((v<_biasNmeth)?2:1);
      //bias[col][v]->SetMinimum(0.05*min(ycorr[col][0],ycorr[col][1]));
    }

    //graph for nominal result
    double x[] = {dipoints[col][0][0],dipoints[col][0][1]};
    double y[] = {dipoints[col][0][2],dipoints[col][0][3]};
    double xelo[] = {x[0]-_BcPtmin[1] , x[1]-_BcPtmin[2]};
    double xehi[] = {_BcPtmax[1]-x[0] , _BcPtmax[2]-x[1]};
    double yelo[] = {yErrLo[col][0], yErrLo[col][1]};
    double yehi[] = {yErrHi[col][0], yErrHi[col][1]};
    TGraphAsymmErrors* nomi = new TGraphAsymmErrors(_NanaBins, x,y, xelo, xehi, yelo, yehi);
    nomi->SetTitle(";p_{T} [GeV];#frac{dN_{corr}}{dp_{T} dy}");
    nomi->SetMarkerStyle(20);
    nomi->SetMarkerSize(3);
    nomi->SetMarkerColor(kCyan+3);
    nomi->SetFillStyle(1001);
    nomi->SetFillColor(kCyan);
    nomi->SetFillColorAlpha(kCyan,0.3);//opacity 1=opaque
    nomi->GetHistogram()->GetYaxis()->SetRangeUser(0.1*min(y[0],y[1]) , 10*max(y[0],y[1]));
    nomi->GetHistogram()->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);

    //graph: ratio of nominal to MC, 2 points
    TH1F* h_pTMC_sc = (col==0)?h_pTMC_pp:h_pTMC_PbPb;
    TH1F* h_pTMC_cont_sc = (col==0)?h_pTMC_cont_pp:h_pTMC_cont_PbPb;
    TGraphAsymmErrors* nomi_r = CloneAndDivide(nomi, h_pTMC_sc, "nominal_ratio_"+(TString)((col==0)?"pp":"PbPb") ); 
    //    TGraph* 

    //Draw all spectra
    TCanvas *c1 = new TCanvas("c1","c1",2000,2000);
    c1->Divide(1,2);
    float rpad = 0.5;
    c1->cd(2)->SetPad(0., 0., 1., rpad);
    c1->cd(1)->SetPad(0., rpad, 1., 1);
    nomi_r->GetHistogram()->GetXaxis()->SetLabelSize(0.06);
    nomi_r->GetHistogram()->GetYaxis()->SetLabelSize(0.06);
    nomi->GetHistogram()->GetYaxis()->SetLabelSize(0.06);
    nomi->GetHistogram()->GetYaxis()->SetTitleSize(0.055);
    nomi->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);

    c1->cd(2)->SetLogx();
    c1->cd(2)->SetLogy();
    c1->cd(2)->SetLeftMargin(0.13);
    c1->cd(2)->SetRightMargin(0.04);
    c1->cd(2)->SetTopMargin(0.);
    c1->cd(2)->SetBottomMargin(0.15);
    gPad->SetTickx(1);

    c1->cd(1)->SetLogx();
    c1->cd(1)->SetLogy();
    c1->cd(1)->SetLeftMargin(0.13);
    c1->cd(1)->SetRightMargin(0.04);
    c1->cd(1)->SetTopMargin(0.04);
    c1->cd(1)->SetBottomMargin(0.);
    nomi->Draw("A2");
    nomi->Draw("Psame");

    for(int v=0;v<_biasNmeth+_biasNmeth*_biasNtoys;v++){
      bias[col][v]->Draw("Lsame");
      //multigraphs for varied (r1,r2) points
      toysDipts[col][v]->Draw("Psame");
    }

    TLatex colTx;
    colTx.SetNDC();
    colTx.SetTextFont(42);
    colTx.SetTextSize(0.06);
    colTx.DrawLatex(0.75,0.905,(TString)((col==0)?"pp":"PbPb"));

    //Draw MC
    if(col==0){
      h_pTMC_cont_pp->Draw("hist][same");
      h_pTMC_pp->Draw("hist][same");
    } else{
      h_pTMC_cont_PbPb->Draw("hist][same");
      h_pTMC_PbPb->Draw("hist][same");
    }
    for(int m=0;m<_biasNmeth;m++){
      if(m!=1) mcfit[col][m][0]->Draw("same");
      mcfit[col][m][1]->Draw("same");
    }

    //Re-draw data and nominal correction
    for(int v=0;v<_biasNmeth;v++){
      bias[col][v]->Draw("Lsame");
      //Get the LW abciss (x-position of points)
      x_LW_2ndstep[col][v] = GetLWabciss(bias[col][v], xbins);
      cout<<"x_LW_2ndstep method "<<v<<" bin1, bin2 = "<<x_LW_2ndstep[col][v][0]<<" "<<x_LW_2ndstep[col][v][1]<<endl;
    }
    nomi->Draw("P2same");

    TLatex fctname;
    fctname.SetNDC();
    fctname.SetTextFont(42);
    fctname.SetTextSize(0.05);

    TLegend *leg = new TLegend(0.64,0.74,0.95,0.88);
    leg->SetTextSize(0.047);
    leg->AddEntry(nomi, "1st step measurement");
    leg->AddEntry(((col==0)?h_pTMC_cont_pp:h_pTMC_cont_PbPb), "original MC");
    leg->AddEntry(((col==0)?h_pTMC_pp:h_pTMC_PbPb), "MC, integrated");
    leg->SetBorderSize(0);
    leg->Draw("same");

    TLegend *leg2 = new TLegend(0.79,0.61,1.,0.74);
    leg2->SetTextSize(0.045);
    leg2->AddEntry(mcfit[col][0][1], "MC fit");
    leg2->AddEntry(bias[col][0], "nominal");
    leg2->AddEntry(bias[col][_biasNmeth], "variations");
    leg2->SetBorderSize(0);
    leg2->Draw("same");
    fctname.DrawLatex(0.7,0.67, _biasMethName[0]);

    TLegend *leg3 = new TLegend(0.79,0.48,1.,0.59);
    leg3->SetTextSize(0.045);
    //leg3->AddEntry(mcfit[col][1][1], "MC fit");
    leg3->AddEntry(bias[col][1], "nominal");
    leg3->AddEntry(bias[col][_biasNmeth+_biasNtoys], "variations");
    leg3->SetBorderSize(0);
    leg3->Draw("same");
    fctname.DrawLatex(0.63,0.53, _biasMethName[1]);

    TLegend *leg4 = new TLegend(0.79,0.33,1.,0.46);
    leg4->SetTextSize(0.045);
    leg4->AddEntry(mcfit[col][2][1], "MC fit (free par.)");
    leg4->AddEntry(bias[col][2], "nominal");
    leg4->AddEntry(bias[col][_biasNmeth+2*_biasNtoys], "variations");
    leg4->SetBorderSize(0);
    leg4->Draw("same");
    fctname.DrawLatex(0.63,0.39, _biasMethName[2]);

    //Draw ratio pad
    c1->cd(2);
    nomi_r->GetYaxis()->SetRangeUser((col==0)?0.31:0.18,(col==0)?3.9:5.5);
    nomi_r->Draw("A2");
    nomi_r->Draw("Psame");

    for(int v=0;v<_biasNmeth+_biasNmeth*_biasNtoys;v++){
      bias_r[col][v] = CloneAndDivide(bias[col][v], mcfit[col][(v<_biasNmeth)?v:((int)(v-_biasNmeth)/_biasNtoys)][1], "biasFunction_"+(TString)((col==0)?"pp_":"PbPb_")+(TString)to_string(v));
      bias_r[col][v]->SetLineColor(color[v]);
      bias_r[col][v]->SetLineStyle((v<_biasNmeth)?1:4);
      bias_r[col][v]->SetLineWidth(1);
      bias_r[col][v]->Draw("L][same");
      //multigraphs for varied (r1,r2) points
      toysDipts_r[col][v] = CloneAndDivide(toysDipts[col][v], h_pTMC_sc, "toysDipts_ratioToMC_"+(TString)((col==0)?"pp_":"PbPb_")+(TString)to_string(v));
      toysDipts_r[col][v]->Draw("Psame");
    }

    //re-draw
    for(int v=0;v<_biasNmeth;v++){
      bias_r[col][v]->SetLineWidth(3);
      bias_r[col][v]->Draw("L0][same");
    }
    nomi_r->Draw("P2same");
    TLine one = TLine();
    one.SetLineStyle(7);
    one.DrawLine(_BcPtmin[0],1,_BcPtmax[0],1);

    c1->SaveAs("figs/ToyBiases_"+(TString)to_string(_biasNtoys)+"toysPerMethod_"+(TString)((col==0)?"pp":"PbPb")+".pdf");
  }

  //Store bias histograms
  if(!plotOnly){
    TFile *outfile = new TFile("pTBiases.root","recreate");
    outfile->WriteObject(&x_LW_MC,"x_LW_signalMC_pTspectrum");
    outfile->WriteObject(&x_LW_2ndstep,"x_LW_correctedpTspectrum");
    for(int col=0;col<2;col++){//pp or PbPb
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){  
	bias_r[col][v]->Write("pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v));
      }
    }
    outfile->Close();      
  }
  infile->Close();      

}
