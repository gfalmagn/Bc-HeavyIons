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
#include "../helpers/hub.cpp"

TF1* f1, f2;

//to find the x_LW_i value satisfying f(x_LW_i) = integralOnBin_i(f)
//h2 is the histogram containing the bin integrals
//Assumes that h is decreasing and monotonous
vector<float> GetLWabciss(TH1F* h, TH1F* h2){
  int nB = h2->GetNbinsX(); // =2
  int iB = 1;
  vector<float> res;
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
vector<float> GetLWabciss(TF1* f, vector<float> xbins){
  int nB = xbins.size()-1; // =2
  vector<float> integ;
  vector<float> res;
  int iB = 1;
  int nEval = 500;
  
  for(int i=1;i<=nB;i++){
    integ.push_back(f->Integral(xbins[i-1], xbins[i]));
    integ[i-1] /= xbins[i]-xbins[i-1];
  }

  for(int bin=0;bin<nEval;bin++){
    if(iB>nB) break;
    float xeval = xbins[0] + ((float)bin/nEval) * (xbins[nB]-xbins[0]);
    float ev = f->Eval( xeval );
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
  
  TH1F *res = new TH1F(name,name+";p_{T}^{#mu#mu#mu} [GeV];MC ratio biased/default",npts,_BcPtmin[0],_BcPtmax[0]);
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
  //Min function is to avoid the singularity n=1 (y=-1)
  TF2 * f1 = new TF2("f1",TString::Format("( (abs(y+1)<1e-4)*TMath::Log([0]/[1]) + (abs(y+1)>=1e-4)*( [0]**(y+1) - [1]**(y+1) )/(y+1) - [2]/x ) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) )",
					  parHiLim[0], parHiLim[0], parLoLim[0], parLoLim[0], parHiLim[1], parHiLim[1], parLoLim[1], parLoLim[1]));
  TF2 * f2 = new TF2("f2",TString::Format("( (abs(y+1)<1e-4)*TMath::Log([0]/[1]) + (abs(y+1)>=1e-4)*( [0]**(y+1) - [1]**(y+1) )/(y+1) - [2]/x ) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (x- %.3e )*1e3/abs(%.3e) ) * TMath::Max(1., (-x+ %.3e )*1e3/abs(%.3e) )",
					  parHiLim[0], parHiLim[0], parLoLim[0], parLoLim[0], parHiLim[1], parHiLim[1], parLoLim[1], parLoLim[1]));
  f1->SetParameters(xbins[1], xbins[0], y[0]);
  f2->SetParameters(xbins[2], xbins[1], y[1]);

  // wrap the functions
  ROOT::Math::WrappedMultiTF1 g1(*f1,2);
  ROOT::Math::WrappedMultiTF1 g2(*f2,2);
  r.AddFunction(g1);
  r.AddFunction(g2);

  r.Solve(pars,1000,1e-4,1e-5);
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

double KaplanIntegral(double *x, double *p) {
  TF1* integrand = new TF1("integrand","[0] / (1+ (x/[1])**2 )**[2]",_BcPtmin[0],_BcPtmax[0]);
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

int FindParamsKaplan(double *pars, double *y, vector<double> xbins, double fix, double fixedPar, double *parLoLim, double *parHiLim, int verbose=1){
  ROOT::Math::MultiRootFinder r(0);
  r.SetPrintLevel(verbose);

  //Parameters [0]-[4]: 3 bin limits , then 2 y values
  TF1 * f1 = new TF1("f1",KaplanIntegral,0,100,9); //0,100 are dummy borns
  TF1 * f2 = new TF1("f2",KaplanIntegral,0,100,9); //0,100 are dummy borns
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

void MakeToyBiases(bool secondStep=false, bool plotOnly=false, int MCclos=-1){

  gStyle->SetOptStat(0);
  auto h_test = new TH1F();
  h_test->SetDefaultSumw2(true);
  TVirtualFitter::SetMaxIterations( 10000 );//default is 5000
  gRandom = new TRandom3(234); //some seed give one or two doubly-failed fits

  bool useLWabciss = true; //Use Lafferty and Wyatt prescription for position of the bin centers
  bool preAEmeanPt = !useLWabciss; // use the mean pT in acceptance uncut sample rather than of observed events
  bool nominalOnly = !secondStep || (MCclos>-1);
  if(plotOnly) _biasNtoys = 40;
  if(nominalOnly || MCclos>-1) _biasNtoys = 0;

  //Create Hub gathering various data
  Hub H = Hub(secondStep,false);
  H.SetFit(false);
  H.SetMetafit();
  H.SetTnP();
  if(useLWabciss) H.SetxLW();
  H.SetpTbias();
  H.SetMCclosure();

  //pT distribution of gen signal MC (acceptance sample), before any cuts = original pT spectrum
  cout<<"Getting pT distribution from original MC..."<<endl;
  TFile *fileMC = TFile::Open(H.fAccName);
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
  vector<TH1F*> h_pTMC_cont;
  h_pTMC_cont.push_back(new TH1F("h_pTMC_cont_pp","h_pTMC_cont;p_{T}^{#mu#mu#mu} [GeV];#frac{dN_{MC}}{dp_{T}^{#mu#mu#mu} dy^{#mu#mu#mu}}",nbins,_BcPtmin[0],_BcPtmax[0]));
  h_pTMC_cont.push_back(new TH1F("h_pTMC_cont_PbPb","h_pTMC_cont;p_{T}^{#mu#mu#mu} [GeV];#frac{dN_{MC}}{dp_{T}^{#mu#mu#mu} dy^{#mu#mu#mu}}",nbins,_BcPtmin[0],_BcPtmax[0]));
  float bwpt_cont = (_BcPtmax[0]-_BcPtmin[0])/nbins;
  float anabins[] = {_BcPtmin[1], _BcPtmin[2],_BcPtmax[2]};
  vector<TH1F*> h_pTMC;
  h_pTMC.push_back(new TH1F("h_pTMC_pp","h_pTMC",_NanaBins,anabins));
  h_pTMC.push_back(new TH1F("h_pTMC_PbPb","h_pTMC",_NanaBins,anabins));

  vector<float> ntot=vector<float>(secondStep?2:1 , 0);
  vector<vector<float> > nbin = vector<vector<float> >(secondStep?2:1,vector<float>(_NanaBins, 0));
  vector<vector<float> > ptSum = vector<vector<float> >(secondStep?2:1,vector<float>(_NanaBins, 0));

  for(int j=0; j<nentries; j++){
    Gen_3mu_4mom->Clear();
    T->GetEntry(j);
    for(int iBc=0;iBc<Gen_Bc_size;iBc++){
      TLorentzVector *gen3mu = (TLorentzVector*) Gen_3mu_4mom->At(iBc);
      float pt = gen3mu->Pt();
      float rap = gen3mu->Rapidity();
      float mass = gen3mu->M();
      if(mass<_mBcMin || mass>_mBcMax) continue; //exclude the events excluded in the template fit

      vector<float> w;
      w.push_back((secondStep?( getBias(H.pTbias_step1[0][_nomMethVar],pt) ):1));
      if(secondStep) w.push_back( getBias((MCclos==-1)?(H.pTbias_step1[1][_nomMethVar]):(H.pTbias_MCclos_step1[MCclos]),pt) );

      //correct by pt and Y bin width
      float bwpt_2, bwY;
      
      for(int col=0;col<(int)(secondStep?2:1);col++){
	ntot[col]+=w[col];
	bool infid=false;
	for(int b=0;b<_NanaBins;b++){
	  if(inFidCuts(b+1, pt, rap)){
	    bwpt_2 = H.pTBinWidth[b+1];
	    bwY = H.YBinWidth[b+1];
	    //ptSum[col][b] += w[col]*pt;
	    nbin[col][b] += w[col];
	    infid=true;
	  }
	}
	if(!infid) continue;

	h_pTMC_cont[col]->Fill(pt, w[col] / (bwpt_cont*bwY));
	h_pTMC[col]->Fill(pt, w[col] / (bwpt_2*bwY));
      }

    }
  }

  if(!secondStep){
    h_pTMC[1] = (TH1F*)h_pTMC[0]->Clone("h_pTMC_PbPb");
    h_pTMC_cont[1] = (TH1F*)h_pTMC_cont[0]->Clone("h_pTMC_PbPb");
  }

  vector<vector<float> > avPt_preAE;// = vector<vector<float> >{{ptSum[0][0]/nbin[0][0], ptSum[0][1]/nbin[0][1]} , {ptSum[secondStep?1:0][0]/nbin[secondStep?1:0][0], ptSum[secondStep?1:0][1]/nbin[secondStep?1:0][1]}};
  vector<vector<float> > x_LW = vector<vector<float> >{GetLWabciss(h_pTMC_cont[0],h_pTMC[0]) ,
							 GetLWabciss(h_pTMC_cont[1],h_pTMC[1])};

  for(int col=0;col<(int)(secondStep?2:1);col++){
    cout<<(TString)((col==0)?"pp":"PbPb")<<endl;
    cout<<"in acc sample, ntot, nbin1, nbin2, fraction passing bin1 or bin2 = "<<ntot[col]<<" "<<nbin[col][0]<<" "<<nbin[col][1]<<" "<<nbin[col][0]/ntot[col]<<" "<<nbin[col][1]/ntot[col]<<" "<<endl;
    //    cout<<"average pT in acceptance MC sample, bin1 and bin2 = "<<avPt_preAE[col][0]<<" "<<avPt_preAE[col][1]<<endl;
    cout<<"LW abcisses from MC spectrum, common to pp and PbPb, bin1 and 2 = "<<x_LW[col][0]<<" "<<x_LW[col][1]<<endl;
  }

  //Plot original MC
  TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
  c2->SetLogy();
  c2->SetLogx();
  TLegend *leg7 = new TLegend(0.55,0.7,0.9,0.9);
  leg7->SetTextSize(0.033);
  leg7->SetBorderSize(0);

  for(int col=0;col<(int)(secondStep?2:1);col++){
    h_pTMC_cont[col]->SetLineWidth(2);
    h_pTMC[col]->SetLineWidth(2);
    h_pTMC[col]->SetLineColor((col==0)?kBlack:(kRed+2));
    h_pTMC_cont[col]->SetLineColor((col==0)?kBlue:kRed);
    h_pTMC_cont[col]->GetXaxis()->SetMoreLogLabels();
    h_pTMC_cont[col]->Draw((col==0)?"hist":"histsame");
    h_pTMC[col]->Draw("histsame");
    leg7->AddEntry(h_pTMC_cont[col], (TString)(secondStep?"2nd step":"original")+" MC, "+(TString)((col==0)?"pp":"PbPb"));
    leg7->AddEntry(h_pTMC[col], "integrated MC, "+(TString)((col==0)?"pp":"PbPb"));
  }
  if(secondStep) leg7->Draw("same");
  
  if(MCclos==-1) c2->SaveAs("figs/pT_originalMC"+(TString)(secondStep?"_2ndStep":"")+".pdf");  
  h_pTMC[1]->SetLineColor(kBlack);
  h_pTMC_cont[1]->SetLineColor(kBlue);

  //scaling to expected corrected yield in data (prefit)
  h_pTMC_cont[0]->Scale(_scaleMCsig[true] * 3000000 / nentries);
  h_pTMC[0]->Scale(_scaleMCsig[true] * 3000000 / nentries);
  h_pTMC_cont[1]->Scale(_scaleMCsig[false] * Ncoll_MB * 4200000 / nentries);
  h_pTMC[1]->Scale(_scaleMCsig[false] * Ncoll_MB * 4200000 / nentries);

  // //postfit (and prefit) yields

  // if(nominalOnly && !secondStep){
  // vector<vector<float> > *acc_oneBinned;
  // vector<vector<float> > *eff_oneBinned_pp;
  // vector<vector<float> > *eff_oneBinned_PbPb;
  // TFile *infile3, *infile2;
  // infile3 = new TFile("../acceptance/acceptanceMap.root","READ");
  // infile3->GetObject("acceptance_oneBinned", acc_oneBinned);
  
  // infile2 = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");
  // infile2->GetObject("efficiency_oneBinned_pp", eff_oneBinned_pp);
  // infile2->GetObject("efficiency_oneBinned_PbPb", eff_oneBinned_PbPb);

  //   vector<vector<vector<float> > > *Yields_postfit_pp;
  //   H.f_yields->GetObject("Yields_postfit_pp", Yields_postfit_pp);
  //   vector<vector<vector<float> > > *Yields_postfit_PbPb;
  //   H.f_yields->GetObject("Yields_postfit_PbPb", Yields_postfit_PbPb);

  //   y_nom_pp = new vector<float>;
  //   y_nom_PbPb = new vector<float>;
  //   for(int b=0;b<_NanaBins;b++){
  //     (H.Ycorr.pp_pt.Val).push_back((*Yields_postfit_pp)[0][b+1][0] / ((*acc_oneBinned)[b+1] * (*eff_oneBinned_pp)[b+1][0]));
  //     H.Ycorr.PbPb_pt.Val.push_back((*Yields_postfit_PbPb)[0][b+1][0] / ((*acc_oneBinned)[b+1] * (*eff_oneBinned_PbPb)[b+1][0]));
  //   }
  // }

  //pp and PbPb corrected yields
  //  vector<vector<float> > yield = vector<vector<float> >(2, vector<float>(_NanaBins,0));
  vector<vector<float> > ycorr = vector<vector<float> >(2, vector<float>(_NanaBins+1,0));
  for(int b=1;b<=_NanaBins;b++){
    //yield[0][b] = (*Yields_postfit_pp)[0][b+1][0];
    //yield[1][b] = (*Yields_postfit_PbPb)[0][b+1][0];

    if(MCclos==-1){
      ycorr[0][b] = H.Ycorr.pp_pt.Val[b] / ( H.pTBinWidth[b] * H.YBinWidth[b]); //yield[0][b] / ((*acc_oneBinned)[b] * (*eff_oneBinned_pp)[b][0]
      ycorr[1][b] = H.Ycorr.PbPb_pt.Val[b] / ( H.pTBinWidth[b] * H.YBinWidth[b]);//yield[1][b] / ((*acc_oneBinned)[b] * (*eff_oneBinned_PbPb)[b][0]
      ycorr[0][0] += H.Ycorr.pp_pt.Val[b];
      ycorr[1][0] += H.Ycorr.PbPb_pt.Val[b];
    } 
    else{
      ycorr[1][b] = H.Ycorr_MCclos[MCclos].Val[b] / ( H.pTBinWidth[b] * H.YBinWidth[b]);   
      ycorr[1][0] += H.Ycorr_MCclos[MCclos].Val[b];
    }
    cout<<"y_nom PbPb bin "<<b<<" = "<<ycorr[1][b]<<endl;
    cout<<"pt y width = "<<H.pTBinWidth[b] <<" "<< H.YBinWidth[b]<<endl;
  }

  cout<<"rrelerr_pp low bin1 bin2 = "<< H.fit.pp_pt.RelErrLo[1]<<" "<<H.fit.pp_pt.RelErrLo[2]<<endl;
  cout<<"rrelerr_PbPb low bin1 bin2 = "<< H.fit.PbPb_pt.RelErrLo[1]<<" "<<H.fit.PbPb_pt.RelErrLo[2]<<endl;
  cout<<"Metafit relerr_pp low bin1 bin2 = "<< H.metafit.pp_pt.RelErrLo[0]<<" "<<H.metafit.pp_pt.RelErrLo[1]<<endl;
  cout<<"Metafit relerr_PbPb low bin1 bin2 = "<< H.metafit.PbPb_pt.RelErrLo[0]<<" "<<H.metafit.PbPb_pt.RelErrLo[1]<<endl;

  vector<vector<float> > corr(2, vector<float>((int)(_NanaBins*(_NanaBins-1)/2), 0));
  //errors on corrected yield
  vector<vector<float> > yErrHi(2, vector<float>(_NanaBins,0));//[pp or PbPb][pt bins from 0]
  vector<vector<float> > yErrLo(2, vector<float>(_NanaBins,0));

  vector<Result> metaAndTnp(2);
  metaAndTnp[0] = SumResult({H.metafit.pp_pt,H.TnP.pp_pt} , ycorr[0],true);
  metaAndTnp[1] = SumResult({H.metafit.PbPb_pt,H.TnP.PbPb_pt} , ycorr[1],true);
  vector<Result> ResAll(2);
  ResAll[0] = SumResult({H.fit.pp_pt,H.metafit.pp_pt,H.TnP.pp_pt} , ycorr[0],true);
  ResAll[1] = SumResult({H.fit.PbPb_pt,H.metafit.PbPb_pt,H.TnP.PbPb_pt} , ycorr[1],true);
  
  for(int b=1;b<=_NanaBins;b++){
    cout<<"bin = "<<b<<endl;

    cout<<"fit POI pp_pt.Val = "<<H.fit.pp_pt.Val[b]<<endl;
    cout<<"pp metaAndTnp ErrLo = "<<metaAndTnp[0].ErrLo[b]<<endl;
    cout<<"relerrlo = "<<metaAndTnp[0].RelErrLo[b]<<endl;
    cout<<"pp ResultAllSources ErrLo = "<<ResAll[0].ErrLo[b]<<endl;
    cout<<"relerrlo = "<<ResAll[0].RelErrLo[b]<<endl;
    cout<<"corrected yield pp   bin "<<b<<" = "<<ycorr[0][b]<<endl;
    //cout<<"equivalent true number of accepted events (pp) bin "<<b<<" = "<<H.Ycorr.pp_pt.Val[b]* (*acc_oneBinned)[0][b] <<endl;
    cout<<"corrected yield pp (MC)   bin "<<b<<" = "<<h_pTMC[0]->GetBinContent(b)<<endl;
    cout<<"corrected yield pp (MC, true number)   bin "<<b<<" = "<<h_pTMC[0]->GetBinContent(b) * ( H.pTBinWidth[b] * H.YBinWidth[b]) <<endl;
    //cout<<"accepted yield pp (MC, true number)   bin "<<b<<" = "<<h_pTMC[0]->GetBinContent(b) * ( H.pTBinWidth[b] * H.YBinWidth[b]) * (*acc_oneBinned)[0][b] <<endl;

    cout<<"fit POI PbPb_pt.Val = "<<H.fit.PbPb_pt.Val[b]<<endl;
    cout<<"PbPb metaAndTnp ErrLo = "<<metaAndTnp[1].ErrLo[b]<<endl;
    cout<<"relerrlo = "<<metaAndTnp[1].RelErrLo[b]<<endl;
    cout<<"PbPb ResultAllSources ErrLo = "<<ResAll[1].ErrLo[b]<<endl;
    cout<<"relerrlo = "<<ResAll[1].RelErrLo[b]<<endl;
    cout<<"corrected yield PbPb bin "<<b<<" = "<<ycorr[1][b]<<endl;
    //cout<<"equivalent true number of accepted events (PbPb) bin "<<b<<" = "<<H.Ycorr.PbPb_pt.Val[b]* (*acc_oneBinned)[1][b]  <<endl;
    cout<<"corrected yield PbPb (MC) bin "<<b<<" = "<<h_pTMC[1]->GetBinContent(b)<<endl;
    cout<<"corrected yield PbPb (MC, true number)   bin "<<b<<" = "<<h_pTMC[1]->GetBinContent(b) * ( H.pTBinWidth[b] * H.YBinWidth[b]) <<endl;
    //cout<<"accepted yield PbPb (MC, true number)   bin "<<b<<" = "<<h_pTMC[1]->GetBinContent(b) * ( H.pTBinWidth[b] * H.YBinWidth[b]) * (*acc_oneBinned)[1][b] <<endl;
  }
  cout<<"correlation factor of bins 1-2 for metafit+TnP error, pp and PbPb: "<<metaAndTnp[0].Corr[0]<<" "<<metaAndTnp[1].Corr[0]<<endl;
  cout<<"correlation factor of bins 1-2 for fit+metafit+TnP error, pp and PbPb: "<<ResAll[0].Corr[0]<<" "<<ResAll[1].Corr[0]<<endl;

  //Get trees containing POIs distributed with the fit likelihood profile
  vector<TTree*> T_poi;
  TFile* f_poi_pp = new  TFile("~/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/higgsCombine_pp_POIfromNLL"+(TString)(secondStep?"_2ndStep":"")+".MultiDimFit.mH120.root","READ");
  TFile* f_poi_PbPb = new  TFile("~/Bc/templateFit/CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test/higgsCombine_PbPb_POIfromNLL"+(TString)(secondStep?"_2ndStep":"")+".MultiDimFit.mH120.root","READ");
  T_poi.push_back((TTree*)f_poi_pp->Get("limit"));
  T_poi.push_back((TTree*)f_poi_PbPb->Get("limit"));

  float r1[2],r1best[2];
  float r2[2],r2best[2];
  float probNLL[2];
  for(int col=0; col<2; col++){
    T_poi[col]->SetBranchAddress("r1", &r1[col]);
    T_poi[col]->SetBranchAddress("r2", &r2[col]);
    T_poi[col]->SetBranchAddress("quantileExpected", &probNLL[col]);

    T_poi[col]->GetEntry(0);
    if(probNLL[col]==-1) { //we expect best fit values to be in the first entry
      r1best[col] = r1[col];
      r2best[col] = r2[col];
    }
    else{
      cout<<"!!!!!! Did not find best-fit POI values in tree! Return now."<<endl;
      return;
    }

  }

  vector<vector<vector<float> > > dipoints = vector<vector<vector<float> > >(2, vector<vector<float> >(_biasNmeth*(_biasNtoys+1), vector<float>(4,0))); //[{x1,x2,r1,r2}] 
  vector<vector<vector<TF1*> > > mcfit = vector<vector<vector<TF1*> > >(2,vector<vector<TF1*> >(_biasNmeth, vector<TF1*>(2)));
  vector<vector<TF1*> > bias = vector<vector<TF1*> >(2,vector<TF1*>(_biasNmeth*(_biasNtoys+1)));
  vector<vector<TH1F*> > bias_r = vector<vector<TH1F*> >(2,vector<TH1F*>(_biasNmeth*(_biasNtoys+1)));
  vector<vector<TGraph*> > toysDipts = vector<vector<TGraph*> >(2 , vector<TGraph*>(_biasNmeth*(_biasNtoys+1)));
  vector<vector<TGraph*> > toysDipts_r = vector<vector<TGraph*> >(2 , vector<TGraph*>(_biasNmeth*(_biasNtoys+1)));
  vector<vector<vector<float> > > y_ratioToBest = vector<vector<vector<float> > >(2,  vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 1)));
  vector<Color_t> color = vector<Color_t>(_biasNmeth*(_biasNtoys+1));
  vector<vector<float> > x_LW_MC = vector<vector<float> >(2);
  vector<vector<vector<float> > > x_LW_2ndstep = vector<vector<vector<float> > >(2, vector<vector<float> >(_biasNmeth));
  vector<double> pTlims = {(double)H.pTLims[0],(double)H.pTLims[1],(double)H.pTLims[2]};

  //loop over variations - ycorr values
  cout<<"Running "<<_biasNtoys<<" toy biases for each of the "<<_biasNmeth<<" methods"<<endl;
  for(int col=0;col<2;col++){//pp or PbPb
    if(MCclos>-1 && col==0) continue; //PbPb only for closure test

    double yMC[] = {h_pTMC[col]->GetBinContent(1) * H.pTBinWidth[1],
		    h_pTMC[col]->GetBinContent(2) * H.pTBinWidth[2]};

    int jentry=0;

    for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
    failedVar: //if the fit failed twice for this varied measurement
      cout<<"variation v="<<v<<endl;

      if(MCclos>-1 && v!=_nomMethVar) continue;
      
      float r1v, r2v;
      for(int j=jentry; j<T_poi[col]->GetEntries(); j++){
	jentry++;
	T_poi[col]->GetEntry(j);

	float probtest = gRandom->Uniform();//between 0 and 1
	if(probtest<probNLL[col]){ //if... then keep these r1,r2 values
	  r1v = r1[col];
	  r2v = r2[col];
	  //cout<<"keep this entry, r1,r2 = "<<r1v<<" "<<r2v<<endl;
	  break;
	}

	//if we need more variations, start again looping on tree entries
	if(jentry==T_poi[col]->GetEntries()){
	  jentry=0; 
	  j=-1;
	}
      }

      //make toys for dipoints
      float z1 = gRandom->Gaus(),z2 = gRandom->Gaus();
      if(v<3) {
	z1=0; z2=0;
	r1v = r1best[col];
	r2v = r2best[col];
      }
      float z22 = (v<3)?0:( metaAndTnp[col].Corr[0]*z1+sqrt(1-pow(metaAndTnp[col].Corr[0],2))*z2 );
      float y1 = ycorr[col][1] * r1v/r1best[col] + z1 * ((z1>0) ? metaAndTnp[col].ErrHi[1] : metaAndTnp[col].ErrLo[1]);
      float y2 = ycorr[col][2] * r2v/r2best[col] + z22 * ((z22>0) ? metaAndTnp[col].ErrHi[2] : metaAndTnp[col].ErrLo[2]);
      if(y1<0.1*ycorr[col][1]) {y1=0.1*ycorr[col][1]; cout<<"using minimal y1"<<endl;}
      if(y2<0.1*ycorr[col][2]) {y2=0.1*ycorr[col][2]; cout<<"using minimal y2"<<endl;}
      if(v<_biasNmeth){
	x_LW_MC[col].push_back( useLWabciss?x_LW[col][0]:(preAEmeanPt?avPt_preAE[col][0]:ptAverage[col][0]) );
	x_LW_MC[col].push_back( useLWabciss?x_LW[col][1]:(preAEmeanPt?avPt_preAE[col][1]:ptAverage[col][1]) );
      }

      float xt[] = {x_LW_MC[col][0], x_LW_MC[col][1]};
      dipoints[col][v] = vector<float>{xt[0],xt[1],y1,y2};
      y_ratioToBest[col][1][v] = y1/ycorr[col][1];
      y_ratioToBest[col][2][v] = y2/ycorr[col][2];
      y_ratioToBest[col][0][v] = (y1*H.pTBinWidth[1]*H.YBinWidth[1] + y2*H.pTBinWidth[2]*H.YBinWidth[2] )/ycorr[col][0];

      float yt[] = {y1, y2};
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
	  //   bias[col][v]->SetParameter(0, (y2-y1) / (x2-x1) );
	  //   bias[col][v]->SetParameter(1, y2 - bias[col][v]->GetParameter(0) * x2 );
	  //   bias[col][v]->SetParameter(2, max(y1*(1+2.5*yErrHi[col][0]/ycorr[col][1]) , y2*(1+2.5*yErrHi[col][1]/ycorr[col][2]))); //2.5 times the relative error is allowed
	  //   bias[col][v]->SetParameter(3, min(y1*(1-2.5*yErrLo[col][0]/ycorr[col][1]) , y2*(1-2.5*yErrLo[col][1]/ycorr[col][2])) );
	  //   bias[col][v]->SetLineColor(kRed);
	  //   toysDipts[col][v]->SetMarkerColor(kRed);
	  // }
	  // else if(m==1){ //quadratic, passing through average pt x-points and highest pt edge
	  //   float x3 = _BcPtmax[0];
	  //   bias[col][v] = new TF1("bias_var"+(TString)(to_string(v)),"[0]*x*x+[1]*x+[2]",_BcPtmin[0],_BcPtmax[0]);
	  //   bias[col][v]->SetParameter(0, (y1-y2) / ((x1-x2)*(x1-x3)) );
	  //   bias[col][v]->SetParameter(1, - bias[col][v]->GetParameter(0) * (x2+x3) );
	  //   bias[col][v]->SetParameter(2, y1 - bias[col][v]->GetParameter(0) * x1*x1 - bias[col][v]->GetParameter(1) * x1 );
	  //   bias[col][v]->SetLineColor(kOrange+10);
	  //   toysDipts[col][v]->SetMarkerColor(kOrange+10);
	  // }
	  if(m==0){ //power (linear in log-log)
	    
	    if(v==m){ 
	      cout<<"method "<<m<<endl;
	      cout<<"MC fit"<<endl;
	      //Exactly fit the MC
	      mcfit[col][m][0] = new TF1("MC_fit_meth"+(TString)(to_string(m)), "[0]*(x**[1])",_BcPtmin[0],_BcPtmax[0]);//(TF1*) bias[col][v]->Clone("MC_fit_meth"+(TString)(to_string(m)));
	      mcfit[col][m][0]->SetParameter(1, -4);//TMath::Log(y2/y1)/TMath::Log(x2/x1) );
	      mcfit[col][m][0]->SetParameter(0, 1e7);//y1*pow(x1, -bias[col][v]->GetParameter(1)));
	      h_pTMC_cont[col]->Fit("MC_fit_meth"+(TString)(to_string(m)));

	      //"fit" the MC from the values of the integral on each bin range
	      cout<<"MC integral fit"<<endl;
	      mcfit[col][m][1] = (TF1*)mcfit[col][m][0]->Clone("MC_fitWithBinIntegral_meth"+(TString)(to_string(m)));
	      double parsMC[] = {(double)(mcfit[col][m][1]->GetParameter(0)), (double)(mcfit[col][m][1]->GetParameter(1))};
	      double parLowLimsMC[] = {1e4, -8};
	      double parHiLimsMC[] = {1e9, -1.5};
	      int status = FindParamsPower(parsMC, yMC, pTlims, parLowLimsMC, parHiLimsMC, 1);
	      if (status!=0) cout<<"ROOT FINDING FAILED!! status "<<status<<" (power)"<<endl;
	      mcfit[col][m][1]->SetParameters(parsMC[0],parsMC[1]);
	    }
	    else if(v==_biasNmeth+m*_biasNtoys) cout<<"Toys for method "<<m<<endl;

	    //basic function for biasing in a given variation
	    bias[col][v] = (TF1*)((v==m)?mcfit[col][m][0]:bias[col][m])->Clone("bias_var"+(TString)(to_string(v)));
	    color[v] = (v==m)?(kGreen+3):kGreen;

	    //"fit" the central value, from the integral on each bin range
	    double pars[] = {(double)(bias[col][v]->GetParameter(0)), (double)(bias[col][v]->GetParameter(1))};
            double parLowLims[] = {pars[0]*1e-2, pars[1]-2};
            double parHiLims[] = {pars[0]*1e2, pars[1]+2};
	    double y[] = {y1*H.pTBinWidth[1], y2*H.pTBinWidth[2]};
	    int status = FindParamsPower(pars, y, pTlims, parLowLims, parHiLims, 0);
	    if (status!=0) {
	      cout<<"ROOT FINDING FAILED !! (power). RETRYING now (SUCCESS IF NO NEWS). y1, y2 = "<<y1<<" "<<y2<<endl;
	      pars[0] = (double)(1.9*bias[col][m]->GetParameter(0)); pars[1] = (double)(0.3*bias[col][m]->GetParameter(1));
	      parLowLims[0] = (bias[col][m]->GetParameter(0))*1e-3; parLowLims[1] = bias[col][m]->GetParameter(1)-5;
	      parHiLims[0] = (bias[col][m]->GetParameter(0))*1e4; parHiLims[1] = bias[col][m]->GetParameter(1)+5;
	      int status2 = FindParamsPower(pars, y, pTlims, parLowLims, parHiLims, 3);
	      if(status2!=0){
		cout<<"FAILED AGAIN !! Giving up, resetting the variation"<<endl;
		if(v<_biasNmeth) {
		  bias[col][v]->Delete();
		  goto failedVar;}
	      }
	    }
	    bias[col][v]->SetParameters(pars[0],pars[1]);

	  }

	  if(m==1 || m==3){ //x^{n+m*log(x)} (quadratic in log-log scale)

	    double fix = (m==1)?2:1;
	    double fixedPar;
	    
	    if(v==m){ 
	      cout<<"method "<<m<<endl;
	      cout<<"MC fit"<<endl;
	      //Exactly fit the MC //only needed once for the powerToLog function
	      if(m==1){
		mcfit[col][m][0] = new TF1("MC_fit_meth"+(TString)to_string(m), "[0]*(x**( [1]+[2]*TMath::Log(x) ))",_BcPtmin[0],_BcPtmax[0]);//(TF1*) bias[col][v]->Clone("MC_fit_meth"+(TString)(to_string(m)));
		mcfit[col][m][0]->SetParameter(1, 2);//TMath::Log(y2/y1)/TMath::Log(x2/x1) );
		mcfit[col][m][0]->SetParameter(0, 1e4);//y1*pow(x1, -bias[col][v]->GetParameter(1)));
		mcfit[col][m][0]->SetParameter(2, -1);
		h_pTMC_cont[col]->Fit("MC_fit_meth"+(TString)(to_string(m)));
	      } 
	      else mcfit[col][m][0] = (TF1*)mcfit[col][1][0]->Clone("MC_fit_meth"+(TString)(to_string(m)));
	      
	      //"fit" the MC from the values of the integral on each bin range
	      cout<<"MC integral fit"<<endl;
	      mcfit[col][m][1] = (TF1*)mcfit[col][m][0]->Clone("MC_fitWithBinIntegral_meth"+(TString)(to_string(m)));
	      fixedPar = (double)(mcfit[col][m][1]->GetParameter(fix));
	      double parsMC[] = {(double)(mcfit[col][m][1]->GetParameter(0)), (double)(mcfit[col][m][1]->GetParameter(m))}; //parameter number m is varied
	      double parLowLimsMC[] = {1e1, (m==1)?(-0.6):(-5)};
	      double parHiLimsMC[] = {1e8, (m==1)?4.:(-0.01)};
	      int status = FindParamsPowerToLog(parsMC, yMC, pTlims, fix, fixedPar, parLowLimsMC, parHiLimsMC, 1);
	      if(status!=0) cout<<"ROOT FINDING FAILED in nominal MC !! status "<<status<<" (power to log). y1, y2 = "<<y1<<" "<<y2<<endl;
	      mcfit[col][m][1]->SetParameters(parsMC[0],(m==1)?parsMC[1]:fixedPar,(m==3)?parsMC[1]:fixedPar);
	    }
	    else if(v==_biasNmeth+m*_biasNtoys) {
	      cout<<"Toys for method "<<m<<endl;
	      fixedPar = (double)(mcfit[col][m][1]->GetParameter(fix));
	    }

	    //basic function for biasing in a given variation
	    if(v==m) cout<<"Fit on 2 data points"<<endl;
	    bias[col][v] = (TF1*) ((v==m)?mcfit[col][m][0]:bias[col][m])->Clone("bias_var"+(TString)(to_string(v)));
	    color[v] = (v==m)?(kMagenta+1):(kMagenta-4);
	    if(m==1) color[v] = (v==m)?(kRed+2):(kOrange+6);
	    
	    //"fit" the central value, from the integral on each bin range
	    double pars[] = {(double)(bias[col][v]->GetParameter(0)), (double)(bias[col][v]->GetParameter(m))}; //parameter number m is varied
	    double y[] = {y1*H.pTBinWidth[1], y2*H.pTBinWidth[2]};
	    double parLowLims[] = {pars[0]*1e-2, pars[1]-2};
	    double parHiLims[] = {pars[0]*1e2, pars[1]+2};
	    int status = FindParamsPowerToLog(pars, y, pTlims, fix, fixedPar, parLowLims, parHiLims, (MCclos==-1 && v!=m)?0:1);
	    if (status!=0) {
	      cout<<"ROOT FINDING FAILED !! (power to log). RETRYING now (SUCCESS IF NO NEWS). y1, y2 = "<<y1<<" "<<y2<<endl;
	      pars[0] = (double)(((m==1)?0.2:0.5)*bias[col][m]->GetParameter(0)); pars[1] = (double)(((m==1)?2.5:1.4)*bias[col][m]->GetParameter(1));
	      parLowLims[0] = (bias[col][m]->GetParameter(0))*1e-4-0.1; parLowLims[1] = bias[col][m]->GetParameter(1) - ((m==1)?8.:8);
	      parHiLims[0] = (bias[col][m]->GetParameter(0))*1e4; parHiLims[1] = bias[col][m]->GetParameter(1) + ((m==1)?8.:6);
	      int status2 = FindParamsPowerToLog(pars, y, pTlims, fix, fixedPar, parLowLims, parHiLims, 3);
	      if(status2!=0){
		cout<<"FAILED AGAIN !! Giving up, resetting the variation"<<endl;
		if(v<_biasNmeth) {
		  bias[col][v]->Delete();
		  goto failedVar;}
	      }
	    }
	    bias[col][v]->SetParameters(pars[0],(m==1)?pars[1]:fixedPar,(m==3)?pars[1]:fixedPar);

	  }

	  if(m==2){//Kaplan spectrum pow(1+pow(x/[0],2), -[1]) with [0] ~ 3-4 GeV and [1] ~ 5-6

	    double fix = 1;
	    double fixedPar;
	    
	    if(v==m){ 
	      cout<<"method "<<m<<endl;
	      cout<<"MC fit"<<endl;
	      //Exactly fit the MC //only needed once for the powerToLog function
	      mcfit[col][m][0] = new TF1("MC_fit_meth"+(TString)to_string(m), "[0] / (1+ (x/[1])**2 )**[2]",_BcPtmin[0],_BcPtmax[0]);
	      mcfit[col][m][0]->SetParameter(1, 3);
	      mcfit[col][m][0]->SetParameter(0, 1e4);//y1*pow(x1, -bias[col][v]->GetParameter(1)));
	      mcfit[col][m][0]->SetParameter(2, 5);
	      h_pTMC_cont[col]->Fit("MC_fit_meth"+(TString)(to_string(m)));
	      
	      //"fit" the MC from the values of the integral on each bin range
	      cout<<"MC integral fit"<<endl;
	      mcfit[col][m][1] = (TF1*)mcfit[col][m][0]->Clone("MC_fitWithBinIntegral_meth"+(TString)(to_string(m)));
	      fixedPar = (double)(mcfit[col][m][1]->GetParameter(fix));
	      double parsMC[] = {(double)(mcfit[col][m][1]->GetParameter(0)), (double)(mcfit[col][m][1]->GetParameter(3-fix))}; //parameter #m is varied
	      double parLowLimsMC[] = {1e1, (fix==1)?0.1:0.3};
	      double parHiLimsMC[] = {1e10, (fix==1)?10.:13.};
	      int status = FindParamsKaplan(parsMC, yMC, pTlims, fix, fixedPar, parLowLimsMC, parHiLimsMC, 1);
	      if(status!=0) cout<<"ROOT FINDING FAILED in nominal MC !! status "<<status<<" (power to log). y1, y2 = "<<y1<<" "<<y2<<endl;
	      mcfit[col][m][1]->SetParameters(parsMC[0],(fix==2)?parsMC[1]:fixedPar,(fix==1)?parsMC[1]:fixedPar);
	    }
	    else if(v==_biasNmeth+m*_biasNtoys) {
	      cout<<"Toys for method "<<m<<endl;
	      fixedPar = (double)(mcfit[col][m][1]->GetParameter(fix));
	    }

	    //basic function for biasing in a given variation
	    if(v==m) cout<<"Fit on 2 data points"<<endl;
	    bias[col][v] = (TF1*) ((v==m)?mcfit[col][m][0]:bias[col][m])->Clone("bias_var"+(TString)(to_string(v)));
	    color[v] = (v==m)?(kViolet+1):(kViolet-4);
	    if(m==3) color[v] = (v==m)?(kRed+2):(kOrange+6);
	    
	    //"fit" the central value, from the integral on each bin range
	    double pars[] = {(double)(bias[col][v]->GetParameter(0)), (double)(bias[col][v]->GetParameter(3-fix))}; //parameter number m is varied
	    double y[] = {y1*H.pTBinWidth[1], y2*H.pTBinWidth[2]};
	    double parLowLims[] = {pars[0]*1e-3, max(0.2,pars[1]-4)};
	    double parHiLims[] = {pars[0]*1e3, pars[1]+4};
	    int status = FindParamsKaplan(pars, y, pTlims, fix, fixedPar, parLowLims, parHiLims, (MCclos==-1 && v!=m)?0:1);
	    if (status!=0) {
	      cout<<"ROOT FINDING FAILED !! (power to log). RETRYING now (SUCCESS IF NO NEWS). y1, y2 = "<<y1<<" "<<y2<<endl;
	      pars[0] = (double)(((m==1)?0.2:0.5)*bias[col][m]->GetParameter(0)); pars[1] = (double)(((m==1)?2.5:1.4)*bias[col][m]->GetParameter(1));
	      parLowLims[0] = (bias[col][m]->GetParameter(0))*1e-4; parLowLims[1] = max(0.2, bias[col][m]->GetParameter(1) - 2);
	      parHiLims[0] = (bias[col][m]->GetParameter(0))*1e4; parHiLims[1] = bias[col][m]->GetParameter(1) + 3;
	      int status2 = FindParamsKaplan(pars, y, pTlims, fix, fixedPar, parLowLims, parHiLims, 3);
	      if(status2!=0){
		cout<<"FAILED AGAIN !! Giving up, resetting the variation"<<endl;
		if(v<_biasNmeth) {
		  bias[col][v]->Delete();
		  goto failedVar;}
	      }
	    }
	    bias[col][v]->SetParameters(pars[0],(fix==2)?pars[1]:fixedPar,(fix==1)?pars[1]:fixedPar);

	  }
	  // else if(m==1){ //exp(-x)*x
	  //   float x3 = _BcPtmax[0];
	  //   bias[col][v] = new TF1("bias_var"+(TString)(to_string(v)),"[0]*exp(-[1]*x)/(x**[2])",_BcPtmin[0],_BcPtmax[0]);
	  //   float lny21 = TMath::Log(y2/y1);
	  //   float lnx21 = TMath::Log(x2/x1);
	  //   float lnx23 = TMath::Log(x2/x3);
	  //   bias[col][v]->SetParameter(2, lny21 / (lnx23*(x2-x1)/(x2-x3) - lnx21) );
	  //   bias[col][v]->SetParameter(1, - ( lny21 + lnx21*bias[col][v]->GetParameter(2) ) / (x2-x1) );
	  //   bias[col][v]->SetParameter(0, pow(x1,bias[col][v]->GetParameter(2))*y1*TMath::Exp(bias[col][v]->GetParameter(1)*x1) );
	  //   bias[col][v]->SetLineColor(kOrange+1);
	  //   toysDipts[col][v]->SetMarkerColor(kOrange+1);
	  // }
	  // else if(m==3){ //exponential
	  //   bias[col][v] = new TF1("bias_var"+(TString)(to_string(v)),"[0]*exp([1]*x)",_BcPtmin[0],_BcPtmax[0]);
	  //   bias[col][v]->SetParameter(1, TMath::Log(y2/y1) / (x2-x1) );
	  //   bias[col][v]->SetParameter(0, y1*TMath::Exp(-x1 * bias[col][v]->GetParameter(1)) );
	  //   //cout<<"x1,x2,y1,y2,f(x1),f(x2)"<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<" "<<bias[col][v]->Eval(x1)<<" "<<bias[col][v]->Eval(x2)<<endl;
	  //   bias[col][v]->SetLineColor(kRed+4);
	  //   toysDipts[col][v]->SetMarkerColor(kRed+4);
	  // }

	  if(v==m){
	    //cout<<"values in 2 bins MC = "<<yMC[0]/H.pTBinWidth[1]<<" "<<yMC[1]/H.pTBinWidth[2]<<endl;
	    //cout<<"integral/pTBinWidth of new function (MC \"fit\") on the first and 2nd bin = "<<mcfit[col][m][1]->Integral(pTlims[0],pTlims[1])/H.pTBinWidth[1]<<" "<<mcfit[col][m][1]->Integral(pTlims[1],pTlims[2])/H.pTBinWidth[2]<<endl;

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
      //bias[col][v]->SetMinimum(0.05*min(ycorr[col][1],ycorr[col][2]));
    }

    //graph for nominal result
    float x[] = {dipoints[col][_nomMethVar][0],dipoints[col][_nomMethVar][1]};
    float y[] = {dipoints[col][_nomMethVar][2],dipoints[col][_nomMethVar][3]};
    float xelo[] = {x[0]-_BcPtmin[1] , x[1]-_BcPtmin[2]};
    float xehi[] = {_BcPtmax[1]-x[0] , _BcPtmax[2]-x[1]};
    float yelo[] = {ResAll[col].ErrLo[1], ResAll[col].ErrLo[2]};
    float yehi[] = {ResAll[col].ErrHi[1], ResAll[col].ErrHi[2]};
    if(MCclos>-1){
      yelo[0] = 0; yelo[1] = 0;
      yehi[0] = 0; yehi[1] = 0;
    }
    TGraphAsymmErrors* nomi = new TGraphAsymmErrors(_NanaBins, x,y, xelo, xehi, yelo, yehi);
    nomi->SetTitle(";p_{T}^{#mu#mu#mu} [GeV];#frac{dN_{corr}}{dp_{T}^{#mu#mu#mu} dy^{#mu#mu#mu}}");
    nomi->SetMarkerStyle(20);
    nomi->SetMarkerSize(3);
    nomi->SetMarkerColor(kCyan+3);
    nomi->SetFillStyle(1001);
    nomi->SetFillColor(kCyan);
    nomi->SetFillColorAlpha(kCyan,0.3);//opacity 1=opaque
    nomi->GetHistogram()->GetYaxis()->SetRangeUser(((MCclos==-1)?0.1:0.07)*min(y[0],y[1]) , ((MCclos==-1)?10:13)*max(y[0],y[1]));
    nomi->GetHistogram()->GetXaxis()->SetRangeUser(_BcPtmin[0],_BcPtmax[0]);

    //graph: ratio of nominal to MC, 2 points
    TH1F* h_pTMC_sc = h_pTMC[col];
    TH1F* h_pTMC_cont_sc = h_pTMC_cont[col];
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
      if(MCclos>-1 && v!=_nomMethVar) continue;
      bias[col][v]->Draw("Lsame");
      //multigraphs for varied (y1,y2) points
      toysDipts[col][v]->Draw("Psame");
    }

    TLatex colTx;
    colTx.SetNDC();
    colTx.SetTextFont(42);
    colTx.SetTextSize(0.06);
    colTx.DrawLatex(0.75,0.905,(TString)((col==0)?"pp":"PbPb"));

    //Draw MC
    h_pTMC_cont[col]->Draw("hist][same");
    h_pTMC[col]->Draw("hist][same");
    for(int m=0;m<_biasNmeth;m++){
      if(MCclos>-1 && m!=_nomMethVar) continue;
      //if(m==2) mcfit[col][m][0]->Draw("same"); //Exact fit of the MC, with 2 bins
      mcfit[col][m][1]->Draw("same"); //Exact fit of the MC, with continuous spectrum
    }

    //Re-draw data and nominal correction
    for(int v=0;v<_biasNmeth;v++){
      if(MCclos>-1 && v!=_nomMethVar) continue;
      bias[col][v]->Draw("Lsame");
      //Get the LW abciss (x-position of points)
      x_LW_2ndstep[col][v] = GetLWabciss(bias[col][v], H.pTLims);
      cout<<"x_LW_2ndstep method "<<v<<" bin1, bin2 = "<<x_LW_2ndstep[col][v][0]<<" "<<x_LW_2ndstep[col][v][1]<<endl;
    }
    nomi->Draw("P2same");

    TLatex fctname;
    fctname.SetNDC();
    fctname.SetTextFont(42);
    fctname.SetTextSize(0.05);

    TLegend *leg = new TLegend(0.64,0.74,0.95,0.88);
    leg->SetTextSize(0.047);
    leg->AddEntry(nomi, (TString)(secondStep?"2nd":"1st")+" step measurement");
    leg->AddEntry(h_pTMC_cont[col], (TString)(secondStep?"2nd step":"original")+" MC");
    leg->AddEntry(h_pTMC[col], "MC, integrated");
    leg->SetBorderSize(0);
    leg->Draw("same");

    TLegend *leg2 = new TLegend(0.79,0.61,1.,0.74);
    leg2->SetTextSize(0.045);
    leg2->AddEntry(mcfit[col][0][1], "MC fit");
    leg2->AddEntry(bias[col][0], "nominal");
    if(!nominalOnly) leg2->AddEntry(bias[col][_biasNmeth], "variations");
    leg2->SetBorderSize(0);
    if(MCclos==-1){
      leg2->Draw("same");
      fctname.DrawLatex(0.7,0.67, _biasMethName[0]);
    }

    TLegend *leg3 = new TLegend(0.79,0.48,1.,0.59);
    leg3->SetTextSize(0.045);
    leg3->AddEntry(mcfit[col][1][1], "MC fit (free par.)");
    leg3->AddEntry(bias[col][1], "nominal");
    if(!nominalOnly) leg3->AddEntry(bias[col][_biasNmeth+_biasNtoys], "variations");
    leg3->SetBorderSize(0);
    leg3->Draw("same");
    fctname.DrawLatex(0.63,0.53, _biasMethName[1]);

    TLegend *leg4 = new TLegend(0.79,0.33,1.,0.46);
    leg4->SetTextSize(0.045);
    leg4->AddEntry(mcfit[col][2][1], "MC fit (free par.)");
    leg4->AddEntry(bias[col][2], "nominal");
    if(!nominalOnly) leg4->AddEntry(bias[col][_biasNmeth+2*_biasNtoys], "variations");
    leg4->SetBorderSize(0);
    if(MCclos==-1) {
      leg4->Draw("same");
      fctname.SetTextSize(0.042);
      fctname.DrawLatex(0.63,0.39, _biasMethName[2]);
    }

    //Draw ratio pad
    c1->cd(2);
    nomi_r->GetYaxis()->SetRangeUser((col==0)?0.41:0.18,(col==0)?2.5:5.5);
    nomi_r->Draw("A2");
    nomi_r->Draw("Psame");

    for(int v=0;v<_biasNmeth+_biasNmeth*_biasNtoys;v++){
      if(MCclos>-1 && v!=_nomMethVar) continue;
      bias_r[col][v] = CloneAndDivide(bias[col][v], mcfit[col][(v<_biasNmeth)?v:((int)(v-_biasNmeth)/_biasNtoys)][1], "biasFunction_"+(TString)((MCclos==-1)?"":("MCclosureToy"+(TString)to_string(MCclos)+"_"))+(TString)((col==0)?"pp_":"PbPb_")+(TString)(secondStep?"_2ndStep":"")+(TString)to_string(v));
      bias_r[col][v]->SetLineColor(color[v]);
      bias_r[col][v]->SetLineStyle((v<_biasNmeth)?1:4);
      bias_r[col][v]->SetLineWidth(1);
      bias_r[col][v]->Draw("L][same");
      //multigraphs for varied (y1,y2) points
      toysDipts_r[col][v] = CloneAndDivide(toysDipts[col][v], h_pTMC_sc, "toysDipts_ratioToMC_"+(TString)((col==0)?"pp_":"PbPb_")+(TString)to_string(v));
      toysDipts_r[col][v]->Draw("Psame");
    }

    //re-draw
    for(int v=0;v<_biasNmeth;v++){
      if(MCclos>-1 && v!=_nomMethVar) continue;
      bias_r[col][v]->SetLineWidth(3);
      bias_r[col][v]->Draw("L0][same");
    }
    nomi_r->Draw("P2same");
    TLine one = TLine();
    one.SetLineStyle(7);
    one.DrawLine(_BcPtmin[0],1,_BcPtmax[0],1);

    c1->SaveAs("figs/ToyBiases_"+(TString)((MCclos==-1)?((TString)to_string(_biasNtoys)+"toysPerMethod_"):("MCclosureToy"+(TString)to_string(MCclos)))+(TString)((col==0)?"pp":"PbPb")+(TString)(secondStep?"_2ndStep":"")+".pdf");
  }

  //Store bias histograms
  if(!plotOnly){
    TFile *outfile = new TFile("pTBiases.root","update");
    if(MCclos==-1){
      outfile->WriteObject(&x_LW_MC,"x_LW_signalMC_pTspectrum"+(TString)(secondStep?"_2ndStep":""));
      outfile->WriteObject(&x_LW_2ndstep,"x_LW_correctedpTspectrum"+(TString)(secondStep?"_2ndStep":""));
      outfile->WriteObject(&y_ratioToBest,"VariedCorrYields_ratioToNominal"+(TString)(secondStep?"_2ndStep":""));
      outfile->WriteObject(&dipoints,"VariedCorrYields"+(TString)(secondStep?"_2ndStep":""));
      for(int col=0;col<2;col++){//pp or PbPb
	for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){  
	  bias_r[col][v]->Write("pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)+(TString)(secondStep?"_2ndStep":""));
	}
      }
    }
    else{
      bias_r[1][_nomMethVar]->Write("pTbias_PbPb_MCclosure_toy"+(TString)to_string(MCclos)+(TString)(secondStep?"_2ndStep":""));
    }
    outfile->Close();
  }

}
