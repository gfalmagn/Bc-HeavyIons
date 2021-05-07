#include <iostream>
#include "Definitions.h"
#include "Cuts.h"
#include "Tools.h"
using namespace std;

class Result {
private:
  int ispp, centDep; //0:PbPb, 1:pp, 2:RAA
public:
  vector<float> Val;
  vector<float> Err, ErrLo, ErrHi;
  vector<float> RelErr, RelErrLo, RelErrHi;
  vector<float> Corr; //correlation matrix is linearised
  void fillErr(bool relOnly=false);
  void print(TString name);
  void scale(vector<float> norm);

  vector<vector<float> > TripleErr(){
    return vector<vector<float> >{Err,ErrLo,ErrHi};
  }
  
  vector<vector<float> > TripleRelErr(){
    return vector<vector<float> >{RelErr,RelErrLo,RelErrHi};
  }
  
  Result(vector<float> value = vector<float>(_NanaBins+1,-1)) { 
    Val = value;
    Err = vector<float>(_NanaBins+1,-1);
    ErrLo = vector<float>(_NanaBins+1,-1);
    ErrHi = vector<float>(_NanaBins+1,-1);
    RelErr = vector<float>(_NanaBins+1,-1);
    RelErrLo = vector<float>(_NanaBins+1,-1);
    RelErrHi = vector<float>(_NanaBins+1,-1);
    Corr = vector<float>((int)(_NanaBins*(_NanaBins-1)/2),-2);
  }

};

void Result::scale(vector<float> norm){
  for(int b=0;b<=_NanaBins;b++){
    Val[b] /= norm[b];
    Err[b] /= norm[b];
    ErrLo[b] /= norm[b];
    ErrHi[b] /= norm[b];
  }
}

void Result::fillErr(bool relOnly=false){ //only non-covered case is when only the RelErr(symmetric) is already filled
  for(int b=0;b<=_NanaBins;b++){

    //set low and high errors to the symmetric one if needed
    if(RelErr[b]!=-1 && ErrLo[b]==-1 && RelErrLo[b]==-1)
      RelErrLo[b] = RelErr[b];
    if(RelErr[b]!=-1 && ErrHi[b]==-1 && RelErrHi[b]==-1)
      RelErrHi[b] = RelErr[b];

    //deduce absolute error from relative error
    if(ErrLo[b]==-1){
      if(RelErrLo[b]!=-1 && Val[b]!=-1) ErrLo[b] = RelErrLo[b]*Val[b];
      else if(Err[b]!=-1) ErrLo[b] = Err[b];
      //else cout<<"Result::fillErr: ErrLo["<<b<<"] cannot be defined"<<endl;
    } 
    if(ErrHi[b]==-1){
      if(RelErrHi[b]!=-1 && Val[b]!=-1) ErrHi[b] = RelErrHi[b]*Val[b];
      else if(Err[b]!=-1) ErrHi[b] = Err[b];
      //else cout<<"Result::fillErr: ErrHi["<<b<<"] cannot be defined"<<endl;
    } 
    //deduce relative error from absolute error
    if(RelErrLo[b]==-1){
      if(ErrLo[b]!=-1 && Val[b]!=-1 && Val[b]!=0) RelErrLo[b] = ErrLo[b]/Val[b];
      //else cout<<"Result::fillErr: RelErrLo["<<b<<"] cannot be defined"<<endl;
    } 
    if(RelErrHi[b]==-1){
      if(ErrHi[b]!=-1 && Val[b]!=-1 && Val[b]!=0) RelErrHi[b] = ErrHi[b]/Val[b];
      //else cout<<"Result::fillErr: RelErrHi["<<b<<"] cannot be defined"<<endl;
    } 

    //deduce symmetric error from asymmetric error
    if(Err[b]==-1) Err[b] = (ErrLo[b]+ErrHi[b])/2;
    if(relOnly && RelErr[b]==-1) RelErr[b] = (RelErrLo[b]+RelErrHi[b])/2;

    if(RelErr[b]==-1){
      if(Err[b]!=-1 && Val[b]!=-1 && Val[b]!=0) RelErr[b] = Err[b]/Val[b];
      //else cout<<"Result::fillErr: RelErr["<<b<<"] cannot be defined"<<endl;
    } 
  }
}

void Result::print(TString name){
  cout<<"\nContent of Result "<<name<<":"<<endl;
  for(int b=0;b<=_NanaBins;b++){
    cout<<(TString)((b==0)?"integrated":("     bin "+(TString)to_string(b)))<<": value = "<<TString::Format("%.4f", Val[b])<<endl;
    cout<<(TString)((b==0)?"integrated":("     bin "+(TString)to_string(b)))<<":          error, errLow, errHi = "<<TString::Format("%.4f %.4f %.4f", Err[b], ErrLo[b], ErrHi[b])<<endl;
    cout<<(TString)((b==0)?"integrated":("     bin "+(TString)to_string(b)))<<": relative error, errLow, errHi = "<<TString::Format("%.4f %.4f %.4f", RelErr[b], RelErrLo[b], RelErrHi[b])<<endl;
  }
  cout<<"Correlation between errors of bin 1 and 2 = "<<Corr[0]<<endl;
}

//Sum of two Result's, including errors and correlations
Result SumResult(vector<Result> vr, vector<float> value = {}, bool favourNewValue=true){
  Result res = Result();
  int nr = vr.size();

  //If no value is given, then it is the sum of input vector values
  if(value.size()!=0) res.Val = value;
  else{
    for(int b=0;b<=_NanaBins;b++){
      res.Val[b] = 0;
      for(int i=0;i<nr;i++)
	res.Val[b] += vr[i].Val[b];
    }
  }

  for(int b=0;b<=_NanaBins;b++){

    //sum the squares, when contents are !=-1
    for(int i=0;i<nr;i++){
      //absolute error
      if(!favourNewValue){ //keep absolute errors from the absolute errors of the components of the sum?
	if(vr[i].Err[b]!=-1)
	  res.Err[b] += pow(vr[i].Err[b] ,2);
	if(vr[i].ErrLo[b]!=-1)
	  res.ErrLo[b] += pow(vr[i].ErrLo[b] ,2);
	if(vr[i].ErrHi[b]!=-1)
	  res.ErrHi[b] += pow(vr[i].ErrHi[b] ,2);
      }

      //relative error
      if(value.size()!=0){ //only do relative error if the value of the sum is forcefully set
	if(vr[i].RelErr[b]!=-1)
	  res.RelErr[b] += pow(vr[i].RelErr[b] ,2);
	if(vr[i].RelErrLo[b]!=-1)
	  res.RelErrLo[b] += pow(vr[i].RelErrLo[b] ,2);
	if(vr[i].RelErrHi[b]!=-1)
	  res.RelErrHi[b] += pow(vr[i].RelErrHi[b] ,2);
      }
    }

    //sqrt for the final error
    if(res.Err[b]!=-1)
      res.Err[b] = sqrt(1+res.Err[b]); //the '1' is to compensate the default -1 value
    if(res.ErrLo[b]!=-1)
      res.ErrLo[b] = sqrt(1+res.ErrLo[b]);
    if(res.ErrHi[b]!=-1)
      res.ErrHi[b] = sqrt(1+res.ErrHi[b]);
    if(res.RelErr[b]!=-1)
      res.RelErr[b] = sqrt(1+res.RelErr[b]);
    if(res.RelErrLo[b]!=-1)
      res.RelErrLo[b] = sqrt(1+res.RelErrLo[b]);
    if(res.RelErrHi[b]!=-1)
      res.RelErrHi[b] = sqrt(1+res.RelErrHi[b]);
  }

  //if no value was given, recover relative error, assuming absolute error was filled. Or recover the absolute error when favourNewValue=true.
  res.fillErr();

  //correlation factor of the sum
  int idx = 0;
  for(int b=1;b<=_NanaBins;b++){
    for(int b2=b+1;b2<=_NanaBins;b2++){
      if(res.Err[b]!=-1 && res.Err[b2]!=-1 && res.Err[b]!=0 && res.Err[b2]!=0){ // if absolute errors are defined and not 0
	//loop over summed components
	float sumErr2b=0, sumErr2b2=0;
	for(int i=0;i<nr;i++){
	  if(vr[i].Corr[idx]!=-2 && (
				     (!favourNewValue && vr[i].Err[b]!=-1 && vr[i].Err[b2]!=-1)
				     || (favourNewValue && vr[i].RelErr[b]!=-1 && vr[i].RelErr[b2]!=-1)
				     ) ){ // if this component's correlation and errors are defined
	    if(res.Corr[idx]==-2) res.Corr[idx] = 0;
	    float errb = favourNewValue?( res.Val[b]*vr[i].RelErr[b] ):( vr[i].Err[b] );
	    float errb2 = favourNewValue?( res.Val[b2]*vr[i].RelErr[b2] ):( vr[i].Err[b2] );
	    res.Corr[idx] += vr[i].Corr[idx] * errb * errb2;
	    sumErr2b += pow(errb ,2);
	    sumErr2b2 += pow(errb2 ,2);
	  }
	}
	res.Corr[idx] /= favourNewValue?( sqrt(sumErr2b*sumErr2b2) ):( res.Err[b]*res.Err[b2] );
      }
      else if(res.Corr[idx]==-2) cout<<"cannot define correlation factor of the sum"<<endl;
      idx += 1;
    } 
  } 

  return res;
}


  
class Source {
private:
  int ispp;
public:
  Result pp_pt, PbPb_pt, RAA_pt; 
  Result PbPb_cent, RAA_cent; 
  void NormsAndRAA(vector<float> normpp, vector<float> normPbPb, vector<float> normPbPbcent);

  Source() {
    pp_pt = Result();
    PbPb_pt = Result();
    RAA_pt = Result();
    PbPb_cent = Result();
    RAA_cent = Result();
  }

  void uniformIntegRelErr(){//starts from well-defined relative error for pt dependence
    PbPb_cent.RelErrLo[0] = PbPb_pt.RelErrLo[0];
    PbPb_cent.RelErrHi[0] = PbPb_pt.RelErrHi[0];
    RAA_cent.RelErrLo[0] = RAA_pt.RelErrLo[0];
    RAA_cent.RelErrHi[0] = RAA_pt.RelErrHi[0];
  }

  void fillErrS(bool relOnly=false){
    pp_pt.fillErr();
    PbPb_pt.fillErr();
    RAA_pt.fillErr(true);
    PbPb_cent.fillErr();
    RAA_cent.fillErr(true);
  }

};

void Source::NormsAndRAA(vector<float> normpp, vector<float> normPbPb, vector<float> normPbPbcent){

  //normalise by luminosity or NMB*TAA, and bin widths
  pp_pt.scale(normpp);
  PbPb_pt.scale(normPbPb);
  PbPb_cent.scale(normPbPbcent);
    
  for(int b=0;b<=_NanaBins;b++){
    RAA_pt.Val[b] = PbPb_pt.Val[b] / pp_pt.Val[b];
    RAA_cent.Val[b] = PbPb_cent.Val[b] / pp_pt.Val[0];

    if(RAA_pt.Err[b]==-1 && RAA_pt.RelErr[b]==-1){
      if(PbPb_pt.RelErr[b]!=-1 && pp_pt.RelErr[b]!=-1)
	RAA_pt.RelErr[b] = quadSum(PbPb_pt.RelErr[b], pp_pt.RelErr[b] );
      if(PbPb_cent.RelErr[b]!=-1 && pp_pt.RelErr[0]!=-1)
	RAA_cent.RelErr[b] = quadSum(PbPb_cent.RelErr[b], pp_pt.RelErr[0] );
    }
    if(RAA_pt.ErrLo[b]==-1 && RAA_pt.RelErrLo[b]==-1){
      if(PbPb_pt.RelErrLo[b]!=-1 && pp_pt.RelErrLo[b]!=-1)
	RAA_pt.RelErrLo[b] = quadSum(PbPb_pt.RelErrLo[b], pp_pt.RelErrLo[b] );
      if(PbPb_cent.RelErrLo[b]!=-1 && pp_pt.RelErrLo[0]!=-1)
	RAA_cent.RelErrLo[b] = quadSum(PbPb_cent.RelErrLo[b], pp_pt.RelErrLo[0] );
    }
    if(RAA_pt.ErrHi[b]==-1 && RAA_pt.RelErrHi[b]==-1){
      if(PbPb_pt.RelErrHi[b]!=-1 && pp_pt.RelErrHi[b]!=-1)
	RAA_pt.RelErrHi[b] = quadSum(PbPb_pt.RelErrHi[b], pp_pt.RelErrHi[b] );
      if(PbPb_cent.RelErrHi[b]!=-1 && pp_pt.RelErrHi[0]!=-1)
	RAA_cent.RelErrHi[b] = quadSum(PbPb_cent.RelErrHi[b], pp_pt.RelErrHi[0] );
    }
  }

  //to fill absolute errors
  RAA_pt.fillErr();  
  RAA_cent.fillErr();  
}




//Hub of all Bc results  
class Hub {
private:
  bool secondStep, thirdStep;
  TString metaSyst, extSyst;
  TFile *f_yields, *f_yields3, *f_acceffToys, *f_fit, *f_pTbias, *f_accept, *f_effic;

public:
  Source Ycorr, fit, fitY, metafit, BcTau, AccEff, TnP, Full, Lumi;
  vector<Result> Ycorr_MCclos,Ycorr_MCclos_step1,Ycorr_MCclos_step2,Ycorr_MCclos_step3;
  vector<float> x_LW_pp,x_LW_PbPb, pTBinWidth, YBinWidth, centBinWidth, pTLims, centLims;
  vector<float> norm_pp, norm_PbPb, norm_PbPb_cent;
  vector<vector<TH1F*> > pTbias_step1,pTbias_step2;
  vector<TH1F*> pTbias_MCclos_step1,pTbias_MCclos_step2;
  TString fAccName;
  vector<vector<vector<float> > > acc_biased; //(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) )) 
  vector<vector<vector<float> > > eff_biased;
  vector<vector<vector<float> > > eff_biased_cent;
  vector<vector<float> > acc_oneBinned; //(2, vector<float>(_NanaBins+1, 0) )
  vector<vector<float> > eff_oneBinned;
  vector<vector<float> > eff_oneBinned_cent;
  vector<vector<float> > acc_oneBinned_2ndStep; //(2, vector<float>(_NanaBins+1, 0) )
  vector<vector<float> > eff_oneBinned_2ndStep;
  vector<vector<float> > eff_oneBinned_cent_2ndStep;
  vector<vector<vector<float> > > corrYr_varied; //corrected yields varied by all uncertainties except AccEff and BcTau //(2, vector<vector<float> >(_NanaBins+1, vector<float>(_biasNmeth*(_biasNtoys+1), 0) )) 

  void SetFit(bool corrRAA=true);
  void SetFullErr(bool corrRAA=true);
  void SetAccEffFinal(bool corrRAA=true);
  void SetMetafit();
  void SetBcTau();
  void SetTnP();
  void SetxLW();
  void SetLumi();
  void SetpTbias(bool allVarIn1stStep=false);
  void SetMCclosure();
  void SetAccEff();
  void ScaleByLumi();

  Hub(bool step2=false, bool step3=false, TString metafitSys="", TString extSys="") {
    secondStep = step2;
    thirdStep = step3;
    metaSyst = metafitSys;
    extSyst = extSys;

    //*** Small calculations
    for(int b=0;b<=_NanaBins;b++){
      pTLims.push_back((b==0)?_BcPtmin[1]:_BcPtmax[b]);
      centLims.push_back((b==0)?_Centmin[1]:_Centmax[b]);
      pTBinWidth.push_back(_BcPtmax[b]-_BcPtmin[b]);
      YBinWidth.push_back(_BcYmax[b]-_BcYmin[b]);
      centBinWidth.push_back(_Centmax[b]-_Centmin[b]);

      if(b>0){
	norm_pp.push_back( L_pp * pTBinWidth[b] * YBinWidth[b] );
	norm_PbPb.push_back( NMB_PbPb * TAA_090 * pTBinWidth[b] * YBinWidth[b]  * centBinWidth[0]/100 ); 
	norm_PbPb_cent.push_back( NMB_PbPb * ((b==1)?TAA_020:TAA_2090)  * centBinWidth[b]/100 ); 
      }
      else{
	norm_pp.push_back( L_pp );
	norm_PbPb.push_back(      NMB_PbPb * TAA_090 * centBinWidth[0]/100 ); 
	norm_PbPb_cent.push_back( NMB_PbPb * TAA_090 * centBinWidth[0]/100 ); 
      }
    }

    //*** FILES
    f_yields = new TFile("../AccEffCorr/corrected_yields"+(TString)(secondStep?"_2ndStep":"")+".root","READ");
    f_yields3 = new TFile("../AccEffCorr/corrected_yields"+(TString)(thirdStep?"_3rdStep":(secondStep?"_2ndStep":""))+".root","READ");
    f_acceffToys = new TFile("../twoSteps/AccEffFrom2ndStepToys.root","READ");
    f_pTbias = new TFile((step2 || step3)?"../twoSteps/pTBiases.root":"","READ");
    f_fit = new TFile();
    fAccName = "/data_CMS/cms/falmagne/tuples/pp17/Bc/TripleMu/acceptance/BcToJpsiMuNu_BCVEGPY_PYTHIA8_pp5TeV_RunIIpp5Spring18DR-00093_acceptance_14092020_ONIATREE_3640k.root";
    f_accept = new TFile("../acceptance/acceptanceMap.root","READ");
    f_effic = new TFile("../efficiency/AcceptanceEfficiencyMap.root","READ");

    //*** nominal corrected YIELD
    vector<float> *y_nom_pp;
    vector<float> *y_nom_PbPb;
    vector<float> *y_nom_PbPbcent;
    f_yields3->GetObject("FinalCorrectedYield_pp", y_nom_pp);
    f_yields3->GetObject("FinalCorrectedYield_PbPb", y_nom_PbPb);
    f_yields3->GetObject("FinalCorrectedYield_centralityDep_PbPb", y_nom_PbPbcent);

    Ycorr.pp_pt = Result(*y_nom_pp);
    Ycorr.PbPb_pt = Result(*y_nom_PbPb);
    Ycorr.PbPb_cent = Result(*y_nom_PbPbcent);

  }

  ~Hub(){
    if(f_yields!=NULL) f_yields->Close();    
    if(f_yields3!=NULL) f_yields3->Close();    
    if(f_acceffToys!=NULL) f_acceffToys->Close();    
    if(f_pTbias!=NULL) f_pTbias->Close();    
    if(f_fit!=NULL) f_fit->Close();    
    if(f_accept!=NULL) f_accept->Close();    
    if(f_effic!=NULL) f_effic->Close();    
  }
};

void Hub::ScaleByLumi(){
  Ycorr.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
  fitY.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
  metafit.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
  AccEff.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
  Full.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
  TnP.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
  BcTau.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
  Lumi.NormsAndRAA(norm_pp,norm_PbPb,norm_PbPb_cent);
}

void Hub::SetFit(bool corrRAA=true){
  //*** FIT POIs and errors
  vector<float> *rsig_pp;//[pt bin+1]
  vector<vector<float> > *rrelerr_pp; //[pt bin+1][sym err, lo err, hi err]
  vector<float> *r1r2Corr_pp;
  vector<float> * rsig_PbPb;
  vector<vector<float> > *rrelerr_PbPb; //[pt bin+1][sym err, lo err, hi err]
  vector<float> *r1r2Corr_PbPb;
  vector<float> * rsig_PbPbcent;
  vector<vector<float> > *rrelerr_PbPbcent; //[pt bin+1][sym err, lo err, hi err]
  vector<float> *r1r2Corr_PbPbcent;
  vector<float> *r1r2Corr_RAA;
  vector<float> *r1r2Corr_RAAcent;

  f_yields->GetObject("rsig_pp"+metaSyst+extSyst,rsig_pp);
  f_yields->GetObject("rsig_relerr_pp"+metaSyst+extSyst,rrelerr_pp);
  f_yields->GetObject("r1r2Correlation_pp"+metaSyst+extSyst, r1r2Corr_pp);
  f_yields->GetObject("rsig_PbPb"+metaSyst+extSyst,rsig_PbPb);
  f_yields->GetObject("rsig_relerr_PbPb"+metaSyst+extSyst,rrelerr_PbPb);
  f_yields->GetObject("r1r2Correlation_PbPb"+metaSyst+extSyst, r1r2Corr_PbPb);
  f_yields->GetObject("rsig_centralityDep_PbPb"+metaSyst+extSyst,rsig_PbPbcent);
  f_yields->GetObject("rsig_relerr_centralityDep_PbPb"+metaSyst+extSyst,rrelerr_PbPbcent);
  f_yields->GetObject("r1r2Correlation_centralityDep_PbPb"+metaSyst+extSyst, r1r2Corr_PbPbcent);

  fit.pp_pt = Result(*rsig_pp);
  fit.PbPb_pt = Result(*rsig_PbPb);
  fit.PbPb_cent = Result(*rsig_PbPbcent);

  for(int b=0;b<=_NanaBins;b++){
    fit.pp_pt.RelErrLo[b] = (*rrelerr_pp)[b][1];
    fit.pp_pt.RelErrHi[b] = (*rrelerr_pp)[b][2];      
    fit.PbPb_pt.RelErrLo[b] = (*rrelerr_PbPb)[b][1];
    fit.PbPb_pt.RelErrHi[b] = (*rrelerr_PbPb)[b][2];      
    fit.PbPb_cent.RelErrLo[b] = (*rrelerr_PbPbcent)[b][1];
    fit.PbPb_cent.RelErrHi[b] = (*rrelerr_PbPbcent)[b][2];      
  }

  //correlations
  fit.pp_pt.Corr = (*r1r2Corr_pp);
  fit.PbPb_pt.Corr = (*r1r2Corr_PbPb);
  fit.PbPb_cent.Corr = (*r1r2Corr_PbPbcent);
    
  if(corrRAA){
    f_yields3->GetObject("r1r2Correlation_RAA", r1r2Corr_RAA);
    f_yields3->GetObject("r1r2Correlation_RAA_centralityBins", r1r2Corr_RAAcent);
    fit.RAA_pt.Corr = (*r1r2Corr_RAA);
    fit.RAA_cent.Corr = (*r1r2Corr_RAAcent);
  }

  //with corr yields as values, instead of POI
  fitY = fit;
  fitY.pp_pt.Val = Ycorr.pp_pt.Val;
  fitY.PbPb_pt.Val = Ycorr.PbPb_pt.Val;
  fitY.PbPb_cent.Val = Ycorr.PbPb_cent.Val;

  fit.fillErrS();
  fitY.fillErrS();
}

void Hub::SetMetafit(){
  //*** METAFIT error
  metafit.pp_pt = Result(Ycorr.pp_pt.Val);
  metafit.PbPb_pt = Result(Ycorr.PbPb_pt.Val);
  metafit.PbPb_cent = Result(Ycorr.PbPb_cent.Val);

  vector<float> *y_metafitRelErrLo_pp;
  vector<float> *y_metafitRelErrHi_pp;
  vector<float> *metafitErrCorr_pp;
  vector<float> *y_metafitRelErrLo_PbPb;
  vector<float> *y_metafitRelErrHi_PbPb;
  vector<float> *metafitErrCorr_PbPb;
  vector<float> *y_metafitRelErrLo_integ;
  vector<float> *y_metafitRelErrHi_integ;
  vector<float> *y_metafitRelErrLo_PbPbcent;
  vector<float> *y_metafitRelErrHi_PbPbcent;
  vector<float> *metafitErrCorr_PbPbcent;
  vector<float> *y_metafitRelErrLo_RAA;
  vector<float> *y_metafitRelErrHi_RAA;
  vector<float> *metafitErrCorr_RAA;
  vector<float> *y_metafitRelErrLo_RAAcent;
  vector<float> *y_metafitRelErrHi_RAAcent;
  vector<float> *metafitErrCorr_RAAcent;
  vector<float> *y_metafitRelErrLo_RAAinteg;
  vector<float> *y_metafitRelErrHi_RAAinteg;

  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorLo_pp", y_metafitRelErrLo_pp);
  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorHi_pp", y_metafitRelErrHi_pp);
  f_yields->GetObject("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixpp", metafitErrCorr_pp);
  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorLo_PbPb", y_metafitRelErrLo_PbPb);
  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorHi_PbPb", y_metafitRelErrHi_PbPb);
  f_yields->GetObject("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixPbPb", metafitErrCorr_PbPb);
  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorLo__integratedPbPbOrpp", y_metafitRelErrLo_integ);
  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorHi__integratedPbPbOrpp", y_metafitRelErrHi_integ);
  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorLo__CentralityDep", y_metafitRelErrLo_PbPbcent);
  f_yields->GetObject("CorrectedYields_MetafitRelSystErrorHi__CentralityDep", y_metafitRelErrHi_PbPbcent);
  f_yields->GetObject("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrix_CentralityDep", metafitErrCorr_PbPbcent);

  f_yields->GetObject("RAA_MetafitRelSystErrorLo_", y_metafitRelErrLo_RAA);
  f_yields->GetObject("RAA_MetafitRelSystErrorHi_", y_metafitRelErrHi_RAA);
  f_yields->GetObject("RAA_MetafitSyst_LinearizedCorrelationMatrix", metafitErrCorr_RAA);
  f_yields->GetObject("RAA_MetafitRelSystErrorLo__CentralityDep", y_metafitRelErrLo_RAAcent);
  f_yields->GetObject("RAA_MetafitRelSystErrorHi__CentralityDep", y_metafitRelErrHi_RAAcent);
  f_yields->GetObject("RAA_MetafitSyst_LinearizedCorrelationMatrix_CentralityDep", metafitErrCorr_RAAcent);
  f_yields->GetObject("RAA_MetafitRelSystErrorLo__integrated", y_metafitRelErrLo_RAAinteg);
  f_yields->GetObject("RAA_MetafitRelSystErrorHi__integrated", y_metafitRelErrHi_RAAinteg);

  for(int b=1;b<=_NanaBins;b++){
    metafit.pp_pt.RelErrLo[b] = (*y_metafitRelErrLo_pp)[b-1];
    metafit.pp_pt.RelErrHi[b] = (*y_metafitRelErrHi_pp)[b-1];
    metafit.PbPb_pt.RelErrLo[b] = (*y_metafitRelErrLo_PbPb)[b-1];
    metafit.PbPb_pt.RelErrHi[b] = (*y_metafitRelErrHi_PbPb)[b-1];
    metafit.RAA_pt.RelErrLo[b] = (*y_metafitRelErrLo_RAA)[b-1];
    metafit.RAA_pt.RelErrHi[b] = (*y_metafitRelErrHi_RAA)[b-1];
    metafit.PbPb_cent.RelErrLo[b] = (*y_metafitRelErrLo_PbPbcent)[b-1];
    metafit.PbPb_cent.RelErrHi[b] = (*y_metafitRelErrHi_PbPbcent)[b-1];
    metafit.RAA_cent.RelErrLo[b] = (*y_metafitRelErrLo_RAAcent)[b-1];
    metafit.RAA_cent.RelErrHi[b] = (*y_metafitRelErrHi_RAAcent)[b-1];
  }

  //integrated
  metafit.pp_pt.RelErrLo[0] = (*y_metafitRelErrLo_integ)[1];
  metafit.pp_pt.RelErrHi[0] = (*y_metafitRelErrHi_integ)[1];
  metafit.PbPb_pt.RelErrLo[0] = (*y_metafitRelErrLo_integ)[0];
  metafit.PbPb_pt.RelErrHi[0] = (*y_metafitRelErrHi_integ)[0];
  metafit.RAA_pt.RelErrLo[0] = (*y_metafitRelErrLo_RAAinteg)[0];
  metafit.RAA_pt.RelErrHi[0] = (*y_metafitRelErrHi_RAAinteg)[0];
  metafit.uniformIntegRelErr();

  //correlations
  metafit.pp_pt.Corr = (*metafitErrCorr_pp);
  metafit.PbPb_pt.Corr = (*metafitErrCorr_PbPb);
  metafit.PbPb_cent.Corr = (*metafitErrCorr_PbPbcent);
  metafit.RAA_pt.Corr = (*metafitErrCorr_RAA);
  metafit.RAA_cent.Corr = (*metafitErrCorr_RAAcent);

  metafit.fillErrS(true); //relOnly if the central value is not defined
}

void Hub::SetAccEffFinal(bool corrRAA=true){
  //*** ACCEPTANCE and EFFICIENCY error
  AccEff.pp_pt = Result(Ycorr.pp_pt.Val);
  AccEff.PbPb_pt = Result(Ycorr.PbPb_pt.Val);
  AccEff.PbPb_cent = Result(Ycorr.PbPb_cent.Val);

  vector<vector<double> > *InvAccEff_pp;
  vector<vector<double> > *InvAccEff_PbPb;
  vector<vector<double> > *InvAccEff_PbPbcent;
  vector<float> *AEcorr12_pp;
  vector<float> *AEcorr12_PbPb;
  vector<float> *AEcorr12_PbPbcent;
  vector<float> *AcceffSystCorr_RAA;
  vector<float> *AcceffSystCorr_RAAcent;

  f_acceffToys->GetObject("InvAccEffFromCorrMC_LinearisedCorrelationFactor_pp"+(TString)(secondStep?"_2ndStep":""), AEcorr12_pp);
  f_acceffToys->GetObject("InvAccEffFromCorrMC_LinearisedCorrelationFactor_PbPb"+(TString)(secondStep?"_2ndStep":""), AEcorr12_PbPb);
  f_acceffToys->GetObject("InvAccEffFromCorrMC_LinearisedCorrelationFactor_PbPb_inCentBins"+(TString)(secondStep?"_2ndStep":""), AEcorr12_PbPbcent);
  f_acceffToys->GetObject("InvAccEffFromCorrMC_withSystErr_pp"+(TString)(secondStep?"_2ndStep":""), InvAccEff_pp);
  f_acceffToys->GetObject("InvAccEffFromCorrMC_withSystErr_PbPb"+(TString)(secondStep?"_2ndStep":""), InvAccEff_PbPb);
  f_acceffToys->GetObject("InvAccEffFromCorrMC_withSystErr_PbPb_inCentBins"+(TString)(secondStep?"_2ndStep":""), InvAccEff_PbPbcent);

  for(int b=0;b<=_NanaBins;b++){
    AccEff.pp_pt.RelErr[b] = (*InvAccEff_pp)[b][1] / (*InvAccEff_pp)[b][0];
    AccEff.pp_pt.RelErrLo[b] = (*InvAccEff_pp)[b][2] / (*InvAccEff_pp)[b][0];
    AccEff.pp_pt.RelErrHi[b] = (*InvAccEff_pp)[b][3] / (*InvAccEff_pp)[b][0];
    AccEff.PbPb_pt.RelErr[b] = (*InvAccEff_PbPb)[b][1] / (*InvAccEff_PbPb)[b][0];
    AccEff.PbPb_pt.RelErrLo[b] = (*InvAccEff_PbPb)[b][2] / (*InvAccEff_PbPb)[b][0];
    AccEff.PbPb_pt.RelErrHi[b] = (*InvAccEff_PbPb)[b][3] / (*InvAccEff_PbPb)[b][0];
    AccEff.PbPb_cent.RelErr[b] = (*InvAccEff_PbPbcent)[b][1] / (*InvAccEff_PbPbcent)[b][0];
    AccEff.PbPb_cent.RelErrLo[b] = (*InvAccEff_PbPbcent)[b][2] / (*InvAccEff_PbPbcent)[b][0];
    AccEff.PbPb_cent.RelErrHi[b] = (*InvAccEff_PbPbcent)[b][3] / (*InvAccEff_PbPbcent)[b][0];
  }

  AccEff.pp_pt.Corr = (*AEcorr12_pp);
  AccEff.PbPb_pt.Corr = (*AEcorr12_PbPb);
  AccEff.PbPb_cent.Corr = (*AEcorr12_PbPbcent);

  if(corrRAA){
    f_yields3->GetObject("AcceffSyst_Correlation_RAA", AcceffSystCorr_RAA);
    f_yields3->GetObject("AcceffSyst_Correlation_RAA_centralityBins", AcceffSystCorr_RAAcent);
    AccEff.RAA_pt.Corr = (*AcceffSystCorr_RAA);
    AccEff.RAA_cent.Corr = (*AcceffSystCorr_RAAcent);
  }

  AccEff.fillErrS();
}

void Hub::SetFullErr(bool corrRAA=true){
  //*** FULL error, excluding BcTau and Lumi
  Full.pp_pt = Result(Ycorr.pp_pt.Val);
  Full.PbPb_pt = Result(Ycorr.PbPb_pt.Val);
  Full.PbPb_cent = Result(Ycorr.PbPb_cent.Val);

  vector<vector<double> > *yFull_pp;
  vector<vector<double> > *yFull_PbPb;
  vector<vector<double> > *yFull_PbPbcent;
  vector<float> *fullcorr12_pp;
  vector<float> *fullcorr12_PbPb;
  vector<float> *fullcorr12_PbPbcent;
  vector<float> *FullCorr_RAA;
  vector<float> *FullCorr_RAAcent;

  f_acceffToys->GetObject("CorrYieldsFromCorrMC_LinearisedCorrelationFactor_pp"+(TString)(secondStep?"_2ndStep":""), fullcorr12_pp);
  f_acceffToys->GetObject("CorrYieldsFromCorrMC_LinearisedCorrelationFactor_PbPb"+(TString)(secondStep?"_2ndStep":""), fullcorr12_PbPb);
  f_acceffToys->GetObject("CorrYieldsFromCorrMC_LinearisedCorrelationFactor_PbPb_inCentBins"+(TString)(secondStep?"_2ndStep":""), fullcorr12_PbPbcent);
  f_acceffToys->GetObject("CorrYieldsFromCorrMC_withSystErr_pp"+(TString)(secondStep?"_2ndStep":""), yFull_pp);
  f_acceffToys->GetObject("CorrYieldsFromCorrMC_withSystErr_PbPb"+(TString)(secondStep?"_2ndStep":""), yFull_PbPb);
  f_acceffToys->GetObject("CorrYieldsFromCorrMC_withSystErr_PbPb_inCentBins"+(TString)(secondStep?"_2ndStep":""), yFull_PbPbcent);

  for(int b=0;b<=_NanaBins;b++){
    Full.pp_pt.RelErr[b] = (*yFull_pp)[b][1] / (*yFull_pp)[b][0];
    Full.pp_pt.RelErrLo[b] = (*yFull_pp)[b][2] / (*yFull_pp)[b][0];
    Full.pp_pt.RelErrHi[b] = (*yFull_pp)[b][3] / (*yFull_pp)[b][0];
    Full.PbPb_pt.RelErr[b] = (*yFull_PbPb)[b][1] / (*yFull_PbPb)[b][0];
    Full.PbPb_pt.RelErrLo[b] = (*yFull_PbPb)[b][2] / (*yFull_PbPb)[b][0];
    Full.PbPb_pt.RelErrHi[b] = (*yFull_PbPb)[b][3] / (*yFull_PbPb)[b][0];
    Full.PbPb_cent.RelErr[b] = (*yFull_PbPbcent)[b][1] / (*yFull_PbPbcent)[b][0];
    Full.PbPb_cent.RelErrLo[b] = (*yFull_PbPbcent)[b][2] / (*yFull_PbPbcent)[b][0];
    Full.PbPb_cent.RelErrHi[b] = (*yFull_PbPbcent)[b][3] / (*yFull_PbPbcent)[b][0];
  }

  Full.pp_pt.Corr = (*fullcorr12_pp);
  Full.PbPb_pt.Corr = (*fullcorr12_PbPb);
  Full.PbPb_cent.Corr = (*fullcorr12_PbPbcent);

  if(corrRAA){
    f_yields3->GetObject("FullErrorCorrelation_RAA", FullCorr_RAA);
    //f_yields3->GetObject("FullErrorCorrelation_RAA_centralityBins", FullCorr_RAAcent);
    Full.RAA_pt.Corr = (*FullCorr_RAA);
    //Full.RAA_cent.Corr = (*FullCorr_RAAcent);
  }

  Full.fillErrS();
}

void Hub::SetAccEff(){
  //*** ACCEPTANCE and EFFICIENCY variations and nominal values

  //Grab acc and eff for all toys (biased)
  vector<vector<vector<float> > > *acc_biased_tmp;
  f_accept->GetObject("acceptance_oneBinned_biased"+(TString)(secondStep?"_2ndStep":""), acc_biased_tmp);
  acc_biased = (*acc_biased_tmp);

  vector<vector<float> > *eff_biased_pp_tmp,*eff_biased_PbPb_tmp,*eff_biased_PbPb_cent_tmp;
  f_effic->GetObject("efficiency_oneBinned_biased_pp"+(TString)(secondStep?"_2ndStep":""), eff_biased_pp_tmp);
  f_effic->GetObject("efficiency_oneBinned_biased_PbPb"+(TString)(secondStep?"_2ndStep":""), eff_biased_PbPb_tmp);
  f_effic->GetObject("efficiency_oneBinned_centDep_biased_PbPb"+(TString)(secondStep?"_2ndStep":""), eff_biased_PbPb_cent_tmp);
  eff_biased.push_back(*eff_biased_pp_tmp);
  eff_biased_cent.push_back(vector<vector<float> >{});
  eff_biased.push_back(*eff_biased_PbPb_tmp);
  eff_biased_cent.push_back(*eff_biased_PbPb_cent_tmp);

  //Grab nominal acc and eff, with or without second-step pT biasing of MC 
  vector<vector<float> > *acc_oneBinned_tmp, *acc_oneBinned_2ndStep_tmp;
  f_accept->GetObject("acceptance_oneBinned", acc_oneBinned_tmp);
  acc_oneBinned = (*acc_oneBinned_tmp);
  if(secondStep){
    f_accept->GetObject("acceptance_oneBinned_2ndStep", acc_oneBinned_2ndStep_tmp);
    acc_oneBinned_2ndStep = (*acc_oneBinned_2ndStep_tmp);
  }

  vector<vector<float> > *eff_oneBinned_pp_tmp;
  vector<vector<float> > *eff_oneBinned_PbPb_tmp;
  vector<vector<float> > *eff_oneBinned_PbPb_cent_tmp;
  f_effic->GetObject("efficiency_oneBinned_pp", eff_oneBinned_pp_tmp);
  f_effic->GetObject("efficiency_oneBinned_PbPb", eff_oneBinned_PbPb_tmp);
  f_effic->GetObject("efficiency_oneBinned_centDep_PbPb", eff_oneBinned_PbPb_cent_tmp);
  eff_oneBinned = vector<vector<float> >(2);
  eff_oneBinned_cent = vector<vector<float> >(2);
  for(int b=0;b<=_NanaBins;b++){
    eff_oneBinned[0].push_back((*eff_oneBinned_pp_tmp)[b][0]);
    eff_oneBinned[1].push_back((*eff_oneBinned_PbPb_tmp)[b][0]);
    eff_oneBinned_cent[1].push_back((*eff_oneBinned_PbPb_cent_tmp)[b][0]);
  }

  if(secondStep){
    vector<vector<float> > *eff_oneBinned_pp_2ndStep_tmp;
    vector<vector<float> > *eff_oneBinned_PbPb_2ndStep_tmp;
    vector<vector<float> > *eff_oneBinned_PbPb_cent_2ndStep_tmp;
    f_effic->GetObject("efficiency_oneBinned_pp_2ndStep", eff_oneBinned_pp_2ndStep_tmp);
    f_effic->GetObject("efficiency_oneBinned_PbPb_2ndStep", eff_oneBinned_PbPb_2ndStep_tmp);
    f_effic->GetObject("efficiency_oneBinned_centDep_PbPb_2ndStep", eff_oneBinned_PbPb_cent_2ndStep_tmp);
    eff_oneBinned_2ndStep = vector<vector<float> >(2);
    eff_oneBinned_cent_2ndStep = vector<vector<float> >(2);
    for(int b=0;b<=_NanaBins;b++){
      eff_oneBinned_2ndStep[0].push_back((*eff_oneBinned_pp_2ndStep_tmp)[b][0]);
      eff_oneBinned_2ndStep[1].push_back((*eff_oneBinned_PbPb_2ndStep_tmp)[b][0]);
      eff_oneBinned_cent_2ndStep[1].push_back((*eff_oneBinned_PbPb_cent_2ndStep_tmp)[b][0]);
    }
  }

  //Grab variations of corrected yields associated to each pT bias
  vector<vector<vector<float> > > *corrYr_varied_tmp;  
  f_pTbias->GetObject("VariedCorrYields_ratioToNominal"+(TString)(secondStep?"_2ndStep":""),corrYr_varied_tmp);
  corrYr_varied = (*corrYr_varied_tmp);
}

void Hub::SetTnP(){
  //*** TNP error
  vector<float> *TnPrelErr_pp;
  vector<float> *TnPrelErr_PbPb;
  vector<float> *TnPrelErr_PbPbcent;
  f_yields3->GetObject("TagAndProbe_relError_pp", TnPrelErr_pp);
  f_yields3->GetObject("TagAndProbe_relError_PbPb", TnPrelErr_PbPb);
  f_yields3->GetObject("TagAndProbe_relError_centralityDep_PbPb", TnPrelErr_PbPbcent);

  TnP.pp_pt = Result(Ycorr.pp_pt.Val);
  TnP.PbPb_pt = Result(Ycorr.PbPb_pt.Val);
  TnP.PbPb_cent = Result(Ycorr.PbPb_cent.Val);
  TnP.pp_pt.RelErr = (*TnPrelErr_pp);
  TnP.PbPb_pt.RelErr = (*TnPrelErr_PbPb);
  TnP.PbPb_cent.RelErr = (*TnPrelErr_PbPbcent);
  TnP.pp_pt.Corr = vector<float>{_corrTnPerr};
  TnP.PbPb_pt.Corr = vector<float>{_corrTnPerr};
  TnP.PbPb_cent.Corr = vector<float>{_corrTnPerr};

  TnP.fillErrS(true);
}

void Hub::SetBcTau(){
  //*** BCTAU error
  BcTau.pp_pt = Result(Ycorr.pp_pt.Val);
  BcTau.PbPb_pt = Result(Ycorr.PbPb_pt.Val);
  BcTau.PbPb_cent = Result(Ycorr.PbPb_cent.Val);

  BcTau.pp_pt.RelErrLo = vector<float>(_NanaBins+1,_XS_BcTauRelSystLo); 
  BcTau.pp_pt.RelErrHi = vector<float>(_NanaBins+1,0); 
  BcTau.PbPb_pt.RelErrLo = vector<float>(_NanaBins+1,_XS_BcTauRelSystLo); 
  BcTau.PbPb_pt.RelErrHi = vector<float>(_NanaBins+1,0); 
  BcTau.PbPb_cent.RelErrLo = vector<float>(_NanaBins+1,_XS_BcTauRelSystLo); 
  BcTau.PbPb_cent.RelErrHi = vector<float>(_NanaBins+1,0); 
  BcTau.RAA_pt.RelErr = vector<float>(_NanaBins+1,_RAA_BcTauRelSyst); 
  BcTau.RAA_cent.RelErr = vector<float>(_NanaBins+1,_RAA_BcTauRelSyst); 

  //correlations    
  BcTau.pp_pt.Corr = vector<float>{_corr_BcTauSyst};
  BcTau.PbPb_pt.Corr = vector<float>{_corr_BcTauSyst};
  BcTau.PbPb_cent.Corr = vector<float>{_corr_BcTauSyst};
  BcTau.RAA_pt.Corr = vector<float>{_corr_BcTauSyst};
  BcTau.RAA_cent.Corr = vector<float>{_corr_BcTauSyst};

  BcTau.fillErrS(true);
}

void Hub::SetLumi(){
  //*** BCTAU error
  Lumi.pp_pt = Result(Ycorr.pp_pt.Val);
  Lumi.PbPb_pt = Result(Ycorr.PbPb_pt.Val);
  Lumi.PbPb_cent = Result(Ycorr.PbPb_cent.Val);

  Lumi.pp_pt.RelErr = vector<float>(_NanaBins+1,_lumiErr_pp); 
  Lumi.PbPb_pt.RelErr = vector<float>(_NanaBins+1,_lumiErr_PbPb); 
  Lumi.PbPb_cent.RelErr = vector<float>{_lumiErr_PbPb,_lumiErr_PbPb020,_lumiErr_PbPb2090};

  //correlations    
  Lumi.pp_pt.Corr = vector<float>{1.};
  Lumi.PbPb_pt.Corr = vector<float>{1.};
  Lumi.PbPb_cent.Corr = vector<float>{1.};
  Lumi.RAA_pt.Corr = vector<float>{1.};
  Lumi.RAA_cent.Corr = vector<float>{1.};

  Lumi.fillErrS(true);
}

//*** LW x-position for bins
void Hub::SetxLW(){
  if(!secondStep && !thirdStep) return;
  vector<vector<vector<float> > > *x_lw; //[col][method][bin]   
  f_pTbias->GetObject("x_LW_correctedpTspectrum"+(TString)(thirdStep?"_2ndStep":""), x_lw);
  x_LW_pp = vector<float>(_NanaBins+1,0);
  x_LW_PbPb = vector<float>(_NanaBins+1,0);

  for(int b=1;b<=_NanaBins;b++){
    x_LW_pp[b] = (*x_lw)[0][_nomMethVar][b-1];
    x_LW_PbPb[b] = (*x_lw)[1][_nomMethVar][b-1];
    //cout<<"bin "<<b<<" x_LW_pp, x_LW_PbPb = "<<x_LW_pp[b]<<" "<<x_LW_PbPb[b]<<endl;
  }
}

void Hub::SetpTbias(bool allVarIn1stStep=false){
  if(!secondStep && !thirdStep) return;
  pTbias_step1 = vector<vector<TH1F*> >(2);
  if(thirdStep) pTbias_step2 = vector<vector<TH1F*> >(2);

  for(int col=0;col<2;col++){
    for(int v=0;v<_biasNmeth*((allVarIn1stStep?_biasNtoys:0)+1);v++){
      //cout<<"pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)<<endl;
      pTbias_step1[col].push_back((TH1F*)f_pTbias->Get("pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)));
    }
    if(thirdStep){
      for(int v=0;v<_biasNmeth*(_biasNtoys+1);v++){
	//cout<<"pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)+"_2ndStep"<<endl;
	pTbias_step2[col].push_back((TH1F*)f_pTbias->Get("pTbias_"+(TString)((col==0)?"pp":"PbPb")+"_var"+(TString)to_string(v)+"_2ndStep"));
      }
    }
  }

}

void Hub::SetMCclosure(){
  pTbias_MCclos_step1 = vector<TH1F*>();
  if(thirdStep) pTbias_MCclos_step2 = vector<TH1F*>();
  for(int t=0;t<_nMCclos;t++){
    pTbias_MCclos_step1.push_back((TH1F*)f_pTbias->Get("pTbias_PbPb_MCclosure_toy"+(TString)to_string(t)));
    if(thirdStep) 
      pTbias_MCclos_step2.push_back((TH1F*)f_pTbias->Get("pTbias_PbPb_MCclosure_toy"+(TString)to_string(t)+"_2ndStep"));	
  }

  vector<vector<float> > *acc_oneBinned;                                                                                                                                                                                                   
  vector<vector<float> > *acc_MCclos_true;
  vector<vector<float> > *acc_MCclos_step2;
  vector<vector<float> > *acc_MCclos_step3;
  f_accept->GetObject("acceptance_oneBinned", acc_oneBinned);                                                                                                
  f_accept->GetObject("acceptance_MCclosure", acc_MCclos_true);   

  vector<vector<float> > *Nobs_MCclos;                                                                                                                                                                                                   
  vector<vector<float> > *eff_oneBinned_PbPb;                                                                                                                                                                                              
  vector<vector<float> > *eff_MCclos_PbPb_true;                                                                                                                                                                  
  vector<vector<float> > *eff_MCclos_PbPb_step2;
  vector<vector<float> > *eff_MCclos_PbPb_step3;
  f_effic->GetObject("Nobs_MCclosure_PbPb", Nobs_MCclos); //only this Nobs (step0) contains the toy MC biasing
  f_effic->GetObject("efficiency_oneBinned_PbPb", eff_oneBinned_PbPb);
  f_effic->GetObject("efficiency_MCclosure_PbPb", eff_MCclos_PbPb_true);                                                                                                    

  if(secondStep){
    f_accept->GetObject("acceptance_MCclosure_2ndStep", acc_MCclos_step2);   
    f_effic->GetObject("efficiency_MCclosure_PbPb_2ndStep", eff_MCclos_PbPb_step2);                                                                                                    
  }
  if(thirdStep){
    f_accept->GetObject("acceptance_MCclosure_3rdStep", acc_MCclos_step3);   
    f_effic->GetObject("efficiency_MCclosure_PbPb_3rdStep", eff_MCclos_PbPb_step3);                                                                                                    
  }

  Ycorr_MCclos = vector<Result>(_nMCclos);
  Ycorr_MCclos_step1 = vector<Result>(_nMCclos);
  Ycorr_MCclos_step2 = vector<Result>(_nMCclos);
  Ycorr_MCclos_step3 = vector<Result>(_nMCclos);

  for(int t=0;t<_nMCclos;t++){
    for(int b=0;b<=_NanaBins;b++){                            
      float acceff_true = (*acc_MCclos_true)[t][b] * (*eff_MCclos_PbPb_true)[b][t];
      float acceff_step1 = (*acc_oneBinned)[1][b] * (*eff_oneBinned_PbPb)[b][0];
      Ycorr_MCclos[t].Val[b] = (*Nobs_MCclos)[b][t] / acceff_true;
      Ycorr_MCclos_step1[t].Val[b] = (*Nobs_MCclos)[b][t] / acceff_step1;

      if(secondStep) {
	float acceff_step2 = (*acc_MCclos_step2)[t][b] * (*eff_MCclos_PbPb_step2)[b][t];
	Ycorr_MCclos_step2[t].Val[b] = (*Nobs_MCclos)[b][t] / acceff_step2;	
      }
      if(thirdStep){
	float acceff_step3 = (*acc_MCclos_step3)[t][b] * (*eff_MCclos_PbPb_step3)[b][t];
	Ycorr_MCclos_step3[t].Val[b] = (*Nobs_MCclos)[b][t] / acceff_step3;
      }
    }   
  }                                                                                                                                                                                                                                   

}



void hub(){
  Hub hu = Hub(true,true);
  hu.SetFit(true);
  hu.SetMetafit();
  hu.SetAccEffFinal(true);
  hu.SetAccEff();
  hu.SetTnP();
  hu.SetBcTau();
  hu.SetLumi();
  hu.SetxLW();
  hu.SetpTbias();
  hu.SetFullErr();
  hu.ScaleByLumi();
  hu.SetMCclosure();

  hu.fit.pp_pt.print("fit pt dep pp");
  hu.fit.PbPb_pt.print("fit pt dep PbPb");
  hu.fit.PbPb_cent.print("fit cent dep PbPb");
  hu.fit.RAA_pt.print("fit pt dep RAA");
  hu.fit.RAA_cent.print("fit cent dep RAA");
  hu.fitY.pp_pt.print("fitY pt dep pp");
  hu.fitY.PbPb_pt.print("fitY pt dep PbPb");
  hu.fitY.PbPb_cent.print("fitY cent dep PbPb");
  hu.fitY.RAA_pt.print("fitY pt dep RAA");
  hu.fitY.RAA_cent.print("fitY cent dep RAA");
  hu.Ycorr.pp_pt.print("Ycorr pt dep pp");
  hu.Ycorr.PbPb_pt.print("Ycorr pt dep PbPb");
  hu.Ycorr.PbPb_cent.print("Ycorr cent dep PbPb");
  hu.Ycorr.RAA_pt.print("Ycorr pt dep RAA");
  hu.Ycorr.RAA_cent.print("Ycorr cent dep RAA");
  hu.metafit.pp_pt.print("metafit pt dep pp");
  hu.metafit.PbPb_pt.print("metafit pt dep PbPb");
  hu.metafit.PbPb_cent.print("metafit cent dep PbPb");
  hu.metafit.RAA_pt.print("metafit pt dep RAA");
  hu.metafit.RAA_cent.print("metafit cent dep RAA");
  hu.AccEff.pp_pt.print("AccEff pt dep pp");
  hu.AccEff.PbPb_pt.print("AccEff pt dep PbPb");
  hu.AccEff.PbPb_cent.print("AccEff cent dep PbPb");
  hu.AccEff.RAA_pt.print("AccEff pt dep RAA");
  hu.AccEff.RAA_cent.print("AccEff cent dep RAA");
  hu.TnP.pp_pt.print("TnP pt dep pp");
  hu.TnP.PbPb_pt.print("TnP pt dep PbPb");
  hu.TnP.PbPb_cent.print("TnP cent dep PbPb");
  hu.BcTau.pp_pt.print("BcTau pt dep pp");
  hu.BcTau.PbPb_pt.print("BcTau pt dep PbPb");
  hu.BcTau.PbPb_cent.print("BcTau cent dep PbPb");
  hu.BcTau.RAA_pt.print("BcTau pt dep RAA");
  hu.BcTau.RAA_cent.print("BcTau cent dep RAA");
  hu.Lumi.pp_pt.print("Lumi pt dep pp");
  hu.Lumi.PbPb_pt.print("Lumi pt dep PbPb");
  hu.Lumi.PbPb_cent.print("Lumi cent dep PbPb");
  hu.Lumi.RAA_pt.print("Lumi pt dep RAA");
  hu.Lumi.RAA_cent.print("Lumi cent dep RAA");
  hu.Full.pp_pt.print("Full except BcTau, pt dep pp");
  hu.Full.PbPb_pt.print("Full except BcTau, pt dep PbPb");
  hu.Full.PbPb_cent.print("Full except BcTau, cent dep PbPb");
  hu.Full.RAA_pt.print("Full except BcTau, pt dep RAA");
  hu.Full.RAA_cent.print("Full except BcTau, cent dep RAA");

  // Result t = Result(vector<float>{7,8,7});
  // t.Corr = vector<float>{0.4};
  // t.Err = vector<float>{2,3,2};
  // t.RelErrLo = vector<float>{0.02,0.03,0.02};
  // Result t2 = Result(vector<float>{7,8,7});
  // t2.Corr = vector<float>{0.5};
  // t2.Err = vector<float>{2,3,2};
  // t2.RelErrLo = vector<float>{0.04,0.05,0.04};
  // Result t3 = SumResult(vector<Result>{t,t2}, t.Val);
  // t.print("1");
  // t2.print("2");
  // t3.print("3");

}
