//!!!! TGraphMultiErrors needs ROOTv20+, for example with cmsenv in ~/miniAOD/CMSSW_11_2_0_pre10

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphMultiErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "../helpers/hub.cpp"
#include "../helpers/nicePalette.h"

void XSandRAA(Hub H, bool centDep=false){
  gInterpreter->GenerateDictionary("vector<vector<vector<float> > >", "vector");

  Result NoAERes_pp,NoAERes_PbPb,NoAERes_RAA;
  Result FullRes_pp,FullRes_PbPb,FullRes_RAA;
  Result FullAltRes_pp,FullAltRes_PbPb,FullAltRes_RAA;
  Result FinRes_pp,FinRes_PbPb,FinRes_RAA;
  //All sources except AccEff
  NoAERes_pp = SumResult({H.fitY.pp_pt,H.metafit.pp_pt,H.TnP.pp_pt,H.BcTau.pp_pt},H.Ycorr.pp_pt.Val);

  //(quasi-)Full = All sources except BcTau and Lumi
  FullRes_pp = H.Full.pp_pt;
  FullAltRes_pp = SumResult({H.fitY.pp_pt,H.metafit.pp_pt,H.TnP.pp_pt,H.AccEff.pp_pt},H.Ycorr.pp_pt.Val);
  //Add quasi-full error and BcTau and Lumi error
  FinRes_pp = SumResult({FullRes_pp,H.BcTau.pp_pt,H.Lumi.pp_pt},H.Ycorr.pp_pt.Val);

  if(centDep){
    NoAERes_PbPb = SumResult({H.fitY.PbPb_cent,H.metafit.PbPb_cent,H.TnP.PbPb_cent,H.BcTau.PbPb_cent},H.Ycorr.PbPb_cent.Val);
    NoAERes_RAA = SumResult({H.fitY.RAA_cent,H.metafit.RAA_cent,H.TnP.RAA_cent,H.BcTau.RAA_cent},H.Ycorr.RAA_cent.Val);
    FullRes_PbPb = H.Full.PbPb_cent;
    FullRes_RAA = H.Full.RAA_cent;
    FullAltRes_PbPb = SumResult({H.fitY.PbPb_cent,H.metafit.PbPb_cent,H.TnP.PbPb_cent,H.AccEff.PbPb_cent},H.Ycorr.PbPb_cent.Val);
    FullAltRes_RAA = SumResult({H.fitY.RAA_cent,H.metafit.RAA_cent,H.TnP.RAA_cent,H.AccEff.RAA_cent},H.Ycorr.RAA_cent.Val);
    FinRes_PbPb = SumResult({FullAltRes_PbPb,H.BcTau.PbPb_cent,H.Lumi.PbPb_cent},H.Ycorr.PbPb_cent.Val); //take alternative Full for centrality
    FinRes_RAA = SumResult({FullAltRes_RAA,H.BcTau.RAA_cent,H.Lumi.RAA_cent},H.Ycorr.RAA_cent.Val);
  } else{
    NoAERes_PbPb = SumResult({H.fitY.PbPb_pt,H.metafit.PbPb_pt,H.TnP.PbPb_pt,H.BcTau.PbPb_pt},H.Ycorr.PbPb_pt.Val);
    NoAERes_RAA = SumResult({H.fitY.RAA_pt,H.metafit.RAA_pt,H.TnP.RAA_pt,H.BcTau.RAA_pt},H.Ycorr.RAA_pt.Val);
    FullRes_PbPb = H.Full.PbPb_pt;
    FullRes_RAA = H.Full.RAA_pt;
    FullAltRes_PbPb = SumResult({H.fitY.PbPb_pt,H.metafit.PbPb_pt,H.TnP.PbPb_pt,H.AccEff.PbPb_pt},H.Ycorr.PbPb_pt.Val);
    FullAltRes_RAA = SumResult({H.fitY.RAA_pt,H.metafit.RAA_pt,H.TnP.RAA_pt,H.AccEff.RAA_pt},H.Ycorr.RAA_pt.Val);
    FinRes_PbPb = SumResult({FullRes_PbPb,H.BcTau.PbPb_pt,H.Lumi.PbPb_pt},H.Ycorr.PbPb_pt.Val);
    FinRes_RAA = SumResult({FullRes_RAA,H.BcTau.RAA_pt,H.Lumi.RAA_pt},H.Ycorr.RAA_pt.Val);
  }

  const int nbins = _NanaBins;
  const int nYregions = 2;//hard-coded
  int nPts[nYregions+1] = {2,1,1};
  int nGr = 4; //pp,PbPb,RAA,RAA_integrated
  const int nMaxPts = 2;//hard-coded
  int nbins_test = 0;
  for(int b=1;b<=nYregions;b++) nbins_test += nPts[b];
  if(nbins_test!=nbins) cout<<"!!!!!!!!!!!! PROBLEM with number of analysis bins and/or number of Y regions ! Expected segfault"<<endl;

  double xErr[nGr][3][nYregions+1][nMaxPts], x[nGr][nYregions+1][nMaxPts];
  double y_metafitErr[nGr][3][nbins];
  double y_BcTauErr[nGr][3][nbins];
  double y_LumiErr[nGr][3][nbins];
  double yErrPart[nGr][3][nbins];
  double yfitErr[nGr][3][nYregions+1][nMaxPts];
  double yacceffErr[nGr][3][nbins];
  double yTnPErr[nGr][3][nbins];
  double yUncorrErr[nGr][3][nYregions+1][nMaxPts];
  double y[nGr][nYregions+1][nMaxPts], yErr[nGr][3][nYregions+1][nMaxPts];
  vector<vector<float> > corrtot(3, vector<float>((int)(nbins*(nbins-1)/2) , 0));

  for(int b=0;b<=nYregions;b++){ //b==0 is the graph with all points, without Y discrimination
    for(int m=nPts[b]-1;m>=0;m--){ //m decreasing is important (y[0][0][1] is needed for centrality dependence from the first bin on)

      int trueb = (b==0)?m:(b-1);
      //x values //take LW prescription for non-centered x
      if(centDep){
	x[0][b][m] = 0.5+trueb; //x=0.5: PbPb integrated, x=1.5: pp integrated
	x[1][b][m] = (_Centmax[trueb+1]+_Centmin[trueb+1])/2;
	x[2][b][m] = (_Centmax[trueb+1]+_Centmin[trueb+1])/2;
	x[3][b][m] = (_Centmax[0]+_Centmin[0])/2;
	xErr[0][1][b][m] = 0.4;
	xErr[0][2][b][m] = 0.4;
	for(int igr=1;igr<4;igr++){
	  xErr[igr][1][b][m] =  x[igr][b][m] - _Centmin[trueb+1];
	  xErr[igr][2][b][m] = -x[igr][b][m] + _Centmax[trueb+1];      
	}
      } else{
	x[0][b][m] = H.x_LW_pp[trueb+1];
	x[1][b][m] = H.x_LW_PbPb[trueb+1];
	x[2][b][m] = (x[0][b][m]+x[1][b][m])/2; //average of pp and PbPb bin centers
	x[3][b][m] = (H.x_LW_PbPb[1]+H.x_LW_PbPb[2])/2;
	for(int igr=0;igr<3;igr++){
	  xErr[igr][1][b][m] =  x[igr][b][m] - _BcPtmin[trueb+1];
	  xErr[igr][2][b][m] = -x[igr][b][m] + _BcPtmax[trueb+1];      
	}
      }
          
      //y values
      if(centDep){
	y[0][b][m] = ((trueb==0)? H.Ycorr.PbPb_pt : H.Ycorr.pp_pt).Val[0];
	y[1][b][m] = H.Ycorr.PbPb_cent.Val[trueb+1];
	y[2][b][m] = H.Ycorr.RAA_cent.Val[trueb+1];
      }
      else {
	y[0][b][m] = H.Ycorr.pp_pt.Val[trueb+1];//(*Yields_postfit_pp)[0][trueb+1][0] * (*InvAccEff_pp)[trueb][0] / normpp; //
	y[1][b][m] = H.Ycorr.PbPb_pt.Val[trueb+1];//(*Yields_postfit_PbPb)[0][trueb+1][0] * (*InvAccEff_PbPb)[trueb][0] / normPbPb;
       	y[2][b][m] = H.Ycorr.RAA_pt.Val[trueb+1];
      } 
      y[3][b][m] = (centDep?H.Ycorr.RAA_cent:H.Ycorr.RAA_pt).Val[0];  //integrated RAA
	
      //metafit error
      for(int pm=0;pm<3;pm++){
	if(centDep){
	  y_metafitErr[0][pm][trueb] = ((trueb==0)? H.metafit.PbPb_pt : H.metafit.pp_pt ).TripleErr()[pm][0];
	  y_metafitErr[1][pm][trueb] = H.metafit.PbPb_cent.TripleErr()[pm][trueb+1];
	  y_metafitErr[2][pm][trueb] = H.metafit.RAA_cent.TripleErr()[pm][trueb+1];
	}		    
	else{		    
	  y_metafitErr[0][pm][trueb] = H.metafit.pp_pt.TripleErr()[pm][trueb+1];
	  y_metafitErr[1][pm][trueb] = H.metafit.PbPb_pt.TripleErr()[pm][trueb+1];
	  y_metafitErr[2][pm][trueb] = H.metafit.RAA_pt.TripleErr()[pm][trueb+1];
	}
	y_metafitErr[3][pm][trueb] = (centDep?H.metafit.RAA_pt:H.metafit.RAA_pt).TripleErr()[pm][0];
      }

      //fit error
      for(int pm=0;pm<3;pm++){
	if(centDep){
	  yfitErr[0][pm][b][m] = ((trueb==0)? H.fitY.PbPb_pt : H.fitY.pp_pt ).TripleErr()[pm][0];
	  yfitErr[1][pm][b][m] = H.fitY.PbPb_cent.TripleErr()[pm][trueb+1];
	  yfitErr[2][pm][b][m] = H.fitY.RAA_cent.TripleErr()[pm][trueb+1];
	}	        
	else{	        
	  yfitErr[0][pm][b][m] = H.fitY.pp_pt.TripleErr()[pm][trueb+1];
	  yfitErr[1][pm][b][m] = H.fitY.PbPb_pt.TripleErr()[pm][trueb+1];
	  yfitErr[2][pm][b][m] = H.fitY.RAA_pt.TripleErr()[pm][trueb+1];
	}
	yfitErr[3][pm][b][m] = (centDep?H.fitY.RAA_cent:H.fitY.RAA_pt).TripleErr()[pm][0];
      }

      //acceff error //Special case, defined as final-total error (quadratic-)minus 'Full' error
      for(int pm=0;pm<3;pm++){
	yacceffErr[1][pm][trueb] = quadSubtract( FinRes_PbPb.TripleErr()[pm][trueb+1] , NoAERes_PbPb.TripleErr()[pm][trueb+1] );
	yacceffErr[2][pm][trueb] = quadSubtract( FinRes_RAA.TripleErr()[pm][trueb+1] , NoAERes_RAA.TripleErr()[pm][trueb+1] );
	if(centDep){
	  yacceffErr[0][pm][trueb] = ((trueb==0)? quadSubtract( FinRes_PbPb.TripleErr()[pm][0] , NoAERes_PbPb.TripleErr()[pm][0] ) :
				      quadSubtract( FinRes_pp.TripleErr()[pm][0] , NoAERes_pp.TripleErr()[pm][0] ) );
	}	        
	else{	        
	  //cout<<"b, pm, FinRes_pp.TripleErr()[pm][trueb+1]/y , NoAERes_pp.TripleErr()[pm][trueb+1]/y = "<<b<<" "<<pm<<" "<<FinRes_pp.TripleErr()[pm][trueb+1]/FinRes_pp.Val[trueb+1]<<" "<< NoAERes_pp.TripleErr()[pm][trueb+1]/FinRes_pp.Val[trueb+1]<<endl;
	  //cout<<"b, pm, FullAltRes_pp.TripleErr()[pm][trueb+1]/y , FullRes_pp.TripleErr()[pm][trueb+1]/y = "<<b<<" "<<pm<<" "<<FullAltRes_pp.TripleErr()[pm][trueb+1]/FullAltRes_pp.Val[trueb+1]<<" "<< FullRes_pp.TripleErr()[pm][trueb+1]/FullRes_pp.Val[trueb+1]<<endl;
	  yacceffErr[0][pm][trueb] = quadSubtract( FinRes_pp.TripleErr()[pm][trueb+1] , NoAERes_pp.TripleErr()[pm][trueb+1] );
	}
	yacceffErr[3][pm][trueb] = quadSubtract( FinRes_RAA.TripleErr()[pm][0] , NoAERes_RAA.TripleErr()[pm][0] );
      }

      //TnP error
      for(int pm=0;pm<3;pm++){
	if(centDep){
	  yTnPErr[0][pm][trueb] = ((trueb==0)? H.TnP.PbPb_pt : H.TnP.pp_pt ).TripleErr()[pm][0];
	  yTnPErr[1][pm][trueb] = H.TnP.PbPb_cent.TripleErr()[pm][trueb+1];
	  yTnPErr[2][pm][trueb] = H.TnP.RAA_cent.TripleErr()[pm][trueb+1];
	}	        
	else{	        
	  yTnPErr[0][pm][trueb] = H.TnP.pp_pt.TripleErr()[pm][trueb+1];
	  yTnPErr[1][pm][trueb] = H.TnP.PbPb_pt.TripleErr()[pm][trueb+1];
	  yTnPErr[2][pm][trueb] = H.TnP.RAA_pt.TripleErr()[pm][trueb+1];
	}
	yTnPErr[3][pm][trueb] = (centDep?H.TnP.RAA_cent:H.TnP.RAA_pt).TripleErr()[pm][0];
      }

      //Bc to tau channel systematic
      for(int pm=0;pm<3;pm++){
	if(centDep){
	  y_BcTauErr[0][pm][trueb] = ((trueb==0)? H.BcTau.PbPb_pt : H.BcTau.pp_pt ).TripleErr()[pm][0];
	  y_BcTauErr[1][pm][trueb] = H.BcTau.PbPb_cent.TripleErr()[pm][trueb+1];
	  y_BcTauErr[2][pm][trueb] = H.BcTau.RAA_cent.TripleErr()[pm][trueb+1];
	}	        
	else{	        
	  y_BcTauErr[0][pm][trueb] = H.BcTau.pp_pt.TripleErr()[pm][trueb+1];
	  y_BcTauErr[1][pm][trueb] = H.BcTau.PbPb_pt.TripleErr()[pm][trueb+1];
	  y_BcTauErr[2][pm][trueb] = H.BcTau.RAA_pt.TripleErr()[pm][trueb+1];
	}
	y_BcTauErr[3][pm][trueb] = (centDep?H.BcTau.RAA_cent:H.BcTau.RAA_pt).TripleErr()[pm][0];
      }

      //Lumi error
      for(int pm=0;pm<3;pm++){
	if(centDep){
	  y_LumiErr[0][pm][trueb] = ((trueb==0)? H.Lumi.PbPb_pt : H.Lumi.pp_pt ).TripleErr()[pm][0];
	  y_LumiErr[1][pm][trueb] = H.Lumi.PbPb_cent.TripleErr()[pm][trueb+1];
	  y_LumiErr[2][pm][trueb] = H.Lumi.RAA_cent.TripleErr()[pm][trueb+1];
	}	        
	else{	        
	  y_LumiErr[0][pm][trueb] = H.Lumi.pp_pt.TripleErr()[pm][trueb+1];
	  y_LumiErr[1][pm][trueb] = H.Lumi.PbPb_pt.TripleErr()[pm][trueb+1];
	  y_LumiErr[2][pm][trueb] = H.Lumi.RAA_pt.TripleErr()[pm][trueb+1];
	}
	y_LumiErr[3][pm][trueb] = (centDep?H.Lumi.RAA_cent:H.Lumi.RAA_pt).TripleErr()[pm][0];
      }

      //total error
      for(int pm=0;pm<3;pm++){
	yErr[1][pm][b][m] = FinRes_PbPb.TripleErr()[pm][trueb+1];
	yErr[2][pm][b][m] = FinRes_RAA.TripleErr()[pm][trueb+1];
	if(centDep){
	  yErr[0][pm][b][m] = (trueb==0)? FinRes_PbPb.TripleErr()[pm][0] : FinRes_pp.TripleErr()[pm][0] ;
	}	        
	else{	        
	  yErr[0][pm][b][m] = FinRes_pp.TripleErr()[pm][trueb+1];
	}
	yErr[3][pm][b][m] = FinRes_RAA.TripleErr()[pm][0];
      }

      if(b==0){
	cout<<"point "<<m<<": x_LW pp, PbPb, RAA = "<<x[0][b][m]<<" "<<x[1][b][m]<<" "<<x[2][b][m]<<endl;
	cout<<"point "<<m<<": y_pp, y_PbPb, y_RAA = "<<y[0][b][m]<<" "<<y[1][b][m]<<" "<<y[2][b][m]<<endl;
	cout<<"point "<<m<<": low  errors on y_pp, y_PbPb, y_RAA = "<<yErr[0][1][b][m]<<" "<<yErr[1][1][b][m]<<" "<<yErr[2][1][b][m]<<endl;
	cout<<"point "<<m<<": high errors on y_pp, y_PbPb, y_RAA = "<<yErr[0][2][b][m]<<" "<<yErr[1][2][b][m]<<" "<<yErr[2][2][b][m]<<endl;
	cout<<"fit     rel errors on pp, PbPb, RAA = "<<yfitErr[0][0][b][m]/y[0][b][m]<<" "<<yfitErr[1][0][b][m]/y[1][b][m]<<" "<<yfitErr[2][0][b][m]/y[2][b][m]<<endl;
	cout<<"metafit rel errors on pp, PbPb, RAA = "<<y_metafitErr[0][0][trueb]/y[0][b][m]<<" "<<y_metafitErr[1][0][trueb]/y[1][b][m]<<" "<<y_metafitErr[2][0][trueb]/y[2][b][m]<<endl;
	cout<<"acceff  rel errors on pp, PbPb, RAA = "<<yacceffErr[0][0][trueb]/y[0][b][m]<<" "<<yacceffErr[1][0][trueb]/y[1][b][m]<<" "<<yacceffErr[2][0][trueb]/y[2][b][m]<<endl;
	cout<<"TnP  rel errors on pp, PbPb, RAA = "<<yTnPErr[0][0][trueb]/y[0][b][m]<<" "<<yTnPErr[1][0][trueb]/y[1][b][m]<<" "<<yTnPErr[2][0][trueb]/y[2][b][m]<<endl;
      }
      
    }
  }
  
  if(!centDep){
    corrtot[0] = FinRes_pp.Corr;// ( (*r1r2Corr_pp)[0] * yfitErr[0][0][0][0] * yfitErr[0][0][0][1]  
		    //   +(*metafitErrCorr_pp)[0] * y_metafitErr[0][0][0] * y_metafitErr[0][0][1] 
		    //   +(*AEcorr12_pp)[0] * yacceffErr[0][0][0] * yacceffErr[0][0][1]
		    //   +_corr_BcTauSyst * y_BcTauErr[0][0][0] * y_BcTauErr[0][0][1]
		    //   +_corrTnPerr * yTnPErr[0][0][0] * yTnPErr[0][0][1]
		    //   )/( yErr[0][0][0][0]*yErr[0][0][0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
    cout<<"correlation factor between 2 bins for pp, from fit, metafit, acceff = "<<H.fitY.pp_pt.Corr[0]<<" "<<H.metafit.pp_pt.Corr[0]<<" "<< H.AccEff.pp_pt.Corr[0]<<endl;
  }
  corrtot[1] = FinRes_PbPb.Corr;// ( (*r1r2Corr_PbPb)[0] * yfitErr[1][0][0][0] * yfitErr[1][0][0][1]  
		  //   +(*metafitErrCorr_PbPb)[0] * y_metafitErr[1][0][0] * y_metafitErr[1][0][1]
		  //   +(*AEcorr12_PbPb)[0] * yacceffErr[1][0][0] * yacceffErr[1][0][1]
		  //   +_corr_BcTauSyst * y_BcTauErr[1][0][0] * y_BcTauErr[1][0][1]
		  //   +_corrTnPerr * yTnPErr[1][0][0] * yTnPErr[1][0][1]
		  //   )/( yErr[1][0][0][0]*yErr[1][0][0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  corrtot[2] = FinRes_RAA.Corr;//  ( (*r1r2Corr_RAA)[0] * yfitErr[2][0][0][0] * yfitErr[2][0][0][1]  
		    // +(*metafitErrCorr_RAA)[0] * y_metafitErr[2][0][0] * y_metafitErr[2][0][1] 
		    // +(*AcceffSystCorr_RAA)[0] * yacceffErr[2][0][0] * yacceffErr[2][0][1]
		    // +_corr_BcTauSyst * y_BcTauErr[2][0][0] * y_BcTauErr[2][0][1]
		    // +_corrTnPerr * yTnPErr[2][0][0] * yTnPErr[2][0][1]
		    // )/( yErr[2][0][0][0]*yErr[2][0][0][1] ); //(Cov(1,2)+Cov'(1,2)) / (sigma_tot(1)*sigma_tot(2))
  if(centDep)  cout<<"correlation factor between 2 bins for PbPb, from fit, metafit, acceff = "<<H.fitY.PbPb_cent.Corr[0]<<" "<<H.metafit.PbPb_cent.Corr[0]<<" "<< H.AccEff.PbPb_cent.Corr[0]<<endl;
  else cout<<"correlation factor between 2 bins for PbPb, from fit, metafit, acceff = "<<H.fitY.PbPb_pt.Corr[0]<<" "<<H.metafit.PbPb_pt.Corr[0]<<" "<< H.AccEff.PbPb_pt.Corr[0]<<endl;
  if(centDep)  cout<<"correlation factor between 2 bins for RAA, from fit, metafit, acceff = "<<H.fitY.RAA_cent.Corr[0]<<" "<<H.metafit.RAA_cent.Corr[0]<<" "<< H.AccEff.RAA_cent.Corr[0]<<endl;
  else cout<<"correlation factor between 2 bins for RAA, from fit, metafit, acceff = "<<H.fitY.RAA_pt.Corr[0]<<" "<<H.metafit.RAA_pt.Corr[0]<<" "<< H.AccEff.RAA_pt.Corr[0]<<endl;
  cout<<"correlation factor between 2 bins, for pp, PbPb, RAA = "<<corrtot[0][0]<<" "<<corrtot[1][0]<<" "<<corrtot[2][0]<<endl;

  //fully uncorrelated part of the error
  for(int b=0;b<=nYregions;b++){
    for(int m=nPts[b]-1;m>=0;m--){
      int trueb = (b==0)?m:(b-1);

      for(int pm=0;pm<3;pm++){
	if(!centDep) yUncorrErr[0][pm][b][m] = UncorrelatedError(corrtot[0][0] , yErr[0][pm][0][trueb] , yErr[0][pm][0][1-trueb]);
	yUncorrErr[1][pm][b][m] = UncorrelatedError(corrtot[1][0] , yErr[1][pm][0][trueb] , yErr[1][pm][0][1-trueb]);
	yUncorrErr[2][pm][b][m] = UncorrelatedError(corrtot[2][0] , yErr[2][pm][0][trueb] , yErr[2][pm][0][1-trueb]);
      }

    }
  }
  
  if(!centDep){
    float RAApTsubtr = H.Ycorr.RAA_pt.Val[1] - H.Ycorr.RAA_pt.Val[2];
    float RAApTsubtrErr = sqrt( pow(FinRes_RAA.ErrLo[1] ,2) + pow(FinRes_RAA.ErrHi[2] ,2) - 2 * FinRes_RAA.Corr[0] * FinRes_RAA.ErrLo[1] * FinRes_RAA.ErrHi[2] );
    cout<<"difference between two pT bins of RAA, value, error, significance of 'difference!=0' : "<<RAApTsubtr<<" "<<RAApTsubtrErr<<" "<<RAApTsubtr/RAApTsubtrErr<<endl;
  }

  //Record final results
  vector<vector<vector<float> > > XSRAA(nGr, vector<vector<float> >(4, vector<float>(nbins, 0))); //dimension1: pp,PbPb,RAA. dimension2: value, errlo, errhi, errsym. dimension3: pt bins
  vector<vector<vector<float> > > metafitRelerr(4, vector<vector<float> >(3, vector<float>(nbins, 0)));//same, without value for dimension2
  for(int b=0;b<nbins;b++){
    for(int p=0;p<nGr;p++){
      XSRAA[p][0][b] = y[p][0][b];
      XSRAA[p][1][b] = yErr[p][1][0][b];
      XSRAA[p][2][b] = yErr[p][2][0][b];
      XSRAA[p][3][b] = yErr[p][0][0][b];
      cout<<"XS or RAA for ngraph,b = "<<p<<" "<<b<<" is = "<<XSRAA[p][0][b]<<" - "<<XSRAA[p][1][b]<<" + "<<XSRAA[p][2][b]<<endl;
    }
    for(int pm=0;pm<3;pm++){
      metafitRelerr[0][pm][b] = centDep?((b==0)?(H.metafit.PbPb_pt.TripleRelErr()[pm][0]):(H.metafit.pp_pt.TripleRelErr()[pm][0])):(H.metafit.pp_pt.TripleRelErr()[pm][b+1]);
      metafitRelerr[1][pm][b] = (centDep?(H.metafit.PbPb_cent):(H.metafit.PbPb_pt)).TripleRelErr()[pm][b+1];
      metafitRelerr[2][pm][b] = (centDep?(H.metafit.RAA_cent):(H.metafit.RAA_pt)).TripleRelErr()[pm][b+1];
      metafitRelerr[3][pm][b] = (centDep?(H.metafit.RAA_cent):(H.metafit.RAA_pt)).TripleRelErr()[pm][0];
    }
  }
  TFile *infile = new TFile("../AccEffCorr/corrected_yields_3rdStep.root","UPDATE");
  infile->WriteObject(&XSRAA,"XSRAA_final"+(TString)(centDep?"_centralityBins":""));
  infile->WriteObject(&metafitRelerr,"MetaFit_RelErr"+(TString)(centDep?"_centralityBins":""));
  infile->WriteObject(&corrtot,"LinearizedCorrelationMatrix_total"+(TString)(centDep?"_centralityBins":""));
  infile->Close();
    
  //graph styles
  gStyle->SetEndErrorSize(15);//size in pixels
  vector<vector<Style_t> > Mstyle = {{20,20,(Style_t)(centDep?20:89)} , {21,21,(Style_t)(centDep?21:90)} , {20,20,(Style_t)(centDep?20:89)}}; //89,90: empty markers, line width 4 (go to 107,108 for width 5)
  // vector<Color_t> Mcol = {(Color_t)(kblue)  , (Color_t)(kspring) , (Color_t)(kazure+2)};
  // vector<Color_t> Lcol = {(Color_t)(kazure) , (Color_t)(kspring) , (Color_t)(kazure)};
  // vector<Color_t> Fcol = {(Color_t)(kcyanLight)  , (Color_t)(kspringLight) , (Color_t)kcyanLight};
  vector<Color_t> Mcol = {(Color_t)(centDep?korange:(kBlue+2))  , kSpring-7 , (Color_t)(kazure)};
  vector<Color_t> Lcol = {(Color_t)(centDep?(kOrange-7):(kAzure-2)) , kSpring-6 , (Color_t)(kazure+2)};
  vector<Color_t> Fcol = {(Color_t)(centDep?korangeLight:(kCyan-6))  , kSpring+9 , (Color_t)kcyanLight};
  // vector<Color_t> Mcol = {kBlue+2  , kSpring-7 , kAzure+4};
  // vector<Color_t> Lcol = {kAzure-2 , kSpring-6 , kAzure+5};
  // vector<Color_t> Fcol = {kCyan-6  , kSpring+9 , kCyan};
  vector<Width_t> Lwidth = {3,3,4};

  vector<vector<TGraphMultiErrors*> > gr(nGr);
  for(int igr=0;igr<3;igr++){
    bool useFitErr = (igr!=2);
    for(int b=0;b<=nYregions;b++){ 
      gr[igr].push_back( new TGraphMultiErrors(nPts[b], x[igr][b],y[igr][b], xErr[igr][1][b], xErr[igr][2][b],yErr[igr][1][b], yErr[igr][2][b]) );
      gr[igr][b]->AddYError(nPts[b], (useFitErr?yfitErr:yUncorrErr)[igr][1][b], (useFitErr?yfitErr:yUncorrErr)[igr][2][b]);
      gr[igr][b]->AddYError(nPts[b], yErr[igr][1][b], yErr[igr][2][b]);
      gr[igr][b]->AddYError(nPts[b], (useFitErr?yfitErr:yUncorrErr)[igr][1][b], (useFitErr?yfitErr:yUncorrErr)[igr][2][b]);

      gr[igr][b]->SetMarkerColor(Mcol[igr]);
      gr[igr][b]->SetMarkerSize(5);
      gr[igr][b]->SetMarkerStyle(Mstyle[igr][b]);
      gr[igr][b]->SetLineColor(Lcol[igr]);
      gr[igr][b]->SetLineWidth(Lwidth[igr]);
      gr[igr][b]->SetFillColor(Fcol[igr]); 
      gr[igr][b]->SetFillStyle(1001);
      for(int p=0;p<4;p++){
	gr[igr][b]->GetAttLine(p)->SetLineColor(Lcol[igr]);
	gr[igr][b]->GetAttLine(p)->SetLineWidth(Lwidth[igr]);
      }

      gr[igr][b]->GetAttFill(0)->SetFillStyle(3001); //full error, draw boxes
      gr[igr][b]->GetAttFill(0)->SetFillColor(Fcol[igr]); 
      gr[igr][b]->GetAttFill(1)->SetFillStyle(1001); //fit/uncorrelated error only, draw boxes
      gr[igr][b]->GetAttFill(1)->SetFillColor(Fcol[igr]); 
      gr[igr][b]->GetAttLine(2)->SetLineStyle(2);
      gr[igr][b]->GetAttLine(3)->SetLineStyle(1);
      
    }
  }
  
  //Draw XS and RAA
  TCanvas *c1 = new TCanvas("c1","XS",2000,2000);
  c1->SetBottomMargin(0.13);
  c1->SetLeftMargin(0.16);
  c1->SetTopMargin(0.08);
  c1->SetRightMargin(0.05);
  gr[1][0]->SetTitle("");//Cross-sections
  gr[1][0]->GetYaxis()->SetRangeUser(0.5*y[1][0][1],5.7*y[1][0][0]);
  if(!centDep) gr[1][0]->GetXaxis()->SetLimits(0,_BcPtmax[0]+0.95);
  if(centDep){
    c1->SetBottomMargin(0.09);
    gr[1][0]->GetYaxis()->SetTitleOffset(1.9);
    gr[1][0]->GetYaxis()->SetRangeUser(0.5*y[0][0][0],2.6*y[0][0][1]);
    gr[1][0]->GetXaxis()->SetRangeUser(0,1.9);
    gr[1][0]->GetXaxis()->SetBinLabel(1,"PbPb");
    gr[1][0]->GetXaxis()->SetBinLabel(2,"pp");
    gr[1][0]->GetXaxis()->SetLabelSize(0.045);
  }
  gr[1][0]->GetYaxis()->SetTitleOffset(1.7);
  gr[1][0]->GetXaxis()->SetTitleOffset(1.2);
  gr[1][0]->GetXaxis()->SetTitle(centDep?"":"p_{T}^{#mu#mu#mu} [GeV]");
  gr[1][0]->GetYaxis()->SetTitle(centDep?"BF  #times  #left(#sigma_{pp}    or    #frac{1}{N_{MB} T_{PbPb}} N_{PbPb}^{corr}#right)  [pb]":"BF  #times  #left(#frac{d#sigma_{pp}}{dp_{T}^{#mu#mu#mu} dy^{#mu#mu#mu}}    or    #frac{1}{N_{MB} T_{PbPb}} #frac{dN_{PbPb}^{corr}}{dp_{T}^{#mu#mu#mu} dy^{#mu#mu#mu}}#right)  [pb/GeV]");
  gr[1][0]->SetMarkerColor(kWhite);
  gr[1][0]->SetLineColor(kWhite);
  gr[1][0]->Draw("APX");
  if(!centDep){
    gr[1][1]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
    gr[1][2]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
  }
  gr[0][1]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
  gr[0][2]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
    
  vector<TLegend*> leg;
  for(int b=0;b<nYregions;b++){
    leg.push_back(new TLegend((b==0)?0.53:0.76,0.6,(b==0)?0.87:0.99,0.75));
    leg[b]->AddEntry(gr[0][b+1],(b==0)?"pp":"pp");
    leg[b]->AddEntry(gr[1][b+1],(b==0)?"PbPb":"PbPb");
    leg[b]->SetTextSize(0.039);
    leg[b]->SetBorderSize(0);
    leg[b]->SetFillStyle(0);
    if(!centDep) leg[b]->Draw("same");
  }
  TLatex legEntry;
  legEntry.SetNDC();
  legEntry.SetTextFont(42);
  legEntry.SetTextSize(0.035);
  if(centDep){
    legEntry.DrawLatex(0.49,0.8,"6<p_{T}^{#mu#mu#mu}<11 GeV & 1.3<|y^{#mu#mu#mu}|<2.3");    
    legEntry.DrawLatex(0.45,0.72,"#font[52]{or} 11<p_{T}^{#mu#mu#mu}<35 GeV & 0<|y^{#mu#mu#mu}|<2.3");
  } else {
    legEntry.DrawLatex(0.53,0.76,"1.3<|y^{#mu#mu#mu}|<2.3");
    legEntry.DrawLatex(0.76,0.76,"0<|y^{#mu#mu#mu}|<2.3");
    legEntry.DrawLatex(0.7,0.46,Form("#rho_{1-2}^{pp} = %.2f",corrtot[0][0]));
    legEntry.DrawLatex(0.7,0.39,Form("#rho_{1-2}^{PbPb} = %.2f",corrtot[1][0]));
    legEntry.DrawLatex(0.58,0.81,"Centrality 0-90%");
  }
  
  if(!centDep) gPad->SetLogy();

  //DRAW CMS Preliminary                                                                                                                                                                                                 
  TLatex CMStag;
  CMStag.SetNDC();
  CMStag.SetTextFont(42);
  CMStag.SetTextSize(0.045);
  CMStag.DrawLatex(0.19,0.875,"#font[61]{CMS}");
  CMStag.DrawLatex(0.19,0.825,"#font[52]{Preliminary}");
  TLatex lumitag;
  lumitag.SetNDC();
  lumitag.SetTextFont(42);
  lumitag.SetTextSize(0.035);
  lumitag.DrawLatex(0.35,0.93, Form("pp (%.0f pb^{-1}) + PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));

  //DRAW CHANNEL                                                                                                                                                                                                 
  TLatex channel;
  channel.SetNDC();
  channel.SetTextFont(42);
  channel.SetTextSize(0.043);
  channel.DrawLatex(0.47,0.875,"#font[61]{B_{c}^{#pm} #rightarrow (J/#psi #rightarrow #mu^{-} #mu^{+}) #mu^{#pm} #nu_{#mu} }");

  c1->SaveAs("CrossSections"+(TString)(centDep?"_integrated":"")+".pdf");

  TCanvas *c2 = new TCanvas("c2","RAA",2000,2000);
  c2->SetLeftMargin(0.11);
  c2->SetBottomMargin(0.13);
  c2->SetTopMargin(0.08);
  c2->SetRightMargin(0.05);
  gr[2][0]->SetTitle("");//R_{PbPb}
  gr[2][0]->SetMarkerColor(kWhite);
  gr[2][0]->SetLineColor(kWhite);
  gr[2][0]->GetYaxis()->SetTitleOffset(1.3);
  gr[2][0]->GetXaxis()->SetTitleOffset(1.2);
  gr[2][0]->GetYaxis()->SetRangeUser(0,3.2);
  gr[2][0]->GetXaxis()->SetLimits(0,_BcPtmax[0]+0.95);
  if(centDep) gr[2][0]->GetXaxis()->SetLimits(_Centmin[0]-3.,_Centmax[0]+3.);
  gr[2][0]->GetXaxis()->SetTitle(centDep?"Centrality [%]":"p_{T}^{#mu#mu#mu} [GeV]");
  gr[2][0]->GetYaxis()->SetTitle("R_{PbPb}(B_{c})");
  gr[2][0]->Draw("AP");
  gr[2][1]->Draw("PX s same; 2; 2; P s=0; P s=0"); //"general; yfullerr; yfiterr;"
  gr[2][2]->Draw("PX s same; 2; 2; P s=0; P s=0");

  TLegend *leg2 = new TLegend(0.55,0.62,0.94,0.78);
  leg2->AddEntry(gr[2][1], "1.3<|y^{#mu#mu#mu}|<2.3");
  leg2->AddEntry(gr[2][2], "0<|y^{#mu#mu#mu}|<2.3");
  leg2->SetTextSize(0.039);
  leg2->SetBorderSize(0);

  if(centDep){
    legEntry.DrawLatex(0.49,0.79,"6<p_{T}^{#mu#mu#mu}<11 GeV & 1.3<|y^{#mu#mu#mu}|<2.3");    
    legEntry.DrawLatex(0.45,0.71,"#font[52]{or} 11<p_{T}^{#mu#mu#mu}<35 GeV & 0<|y^{#mu#mu#mu}|<2.3");
  } else {
    leg2->Draw("same");
    legEntry.DrawLatex(0.55,0.81,"Centrality 0-90%");
  }
  legEntry.DrawLatex(0.7,centDep?0.23:0.48,Form("#rho_{1-2} = %.2f",corrtot[2][0]));

  CMStag.DrawLatex(0.14,0.87,"#font[61]{CMS}");
  CMStag.DrawLatex(0.14,0.82,"#font[52]{Preliminary}");
  lumitag.DrawLatex(0.35,0.93, Form("pp (%.0f pb^{-1}) + PbPb (%.2f nb^{-1}), 5.02 TeV",L_pp,L_PbPb*1e3));
  channel.DrawLatex(0.47,0.87,"#font[61]{B_{c}^{#pm} #rightarrow (J/#psi #rightarrow #mu^{-} #mu^{+}) #mu^{#pm} #nu_{#mu} }");

  TLine *one = new TLine((centDep?(_Centmin[0]-3):0),1,(centDep?(_Centmax[0]+3):(_BcPtmax[0]+0.95)),1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->Draw("same");

  c2->SaveAs("RAA"+(TString)(centDep?"_centralityBins":"")+".pdf");
  
  ///////////////////////////////////////
  // Draw all systematics
  ///////////////////////////////////////
  int nerr = 7;
  vector<TH1F*> S(nerr);
  vector<TString> name = {"fit","metafit","acceff","TnP","BcTau","Lumi","total"};
  vector<TString> prettyname = {"fit","fit method","Acc #times Eff","Tag-and-probe","B_{c}#rightarrow J/#psi #tau","Luminosity+Glauber","total"};
  vector<int> Sorder = {5,4,3,2,1,0,6};
  vector<Color_t> col = {kRed,kBlue,kGreen+1,kOrange-6,kMagenta-4, kYellow+1, kBlack};
  
  for(int i=0;i<nerr;i++)
    S[i] = new TH1F(name[i],"",centDep?7:6,0,centDep?7:6);

  S[0]->SetBinContent(centDep?2:1,yfitErr[0][0][0][0]/y[0][0][0]);//take average of hi and lo errors
  S[0]->SetBinContent(centDep?1:2,yfitErr[0][0][0][1]/y[0][0][1]);
  S[0]->SetBinContent(3,yfitErr[1][0][0][0]/y[1][0][0]);
  S[0]->SetBinContent(4,yfitErr[1][0][0][1]/y[1][0][1]);
  S[0]->SetBinContent(5,yfitErr[2][0][0][0]/y[2][0][0]);
  S[0]->SetBinContent(6,yfitErr[2][0][0][1]/y[2][0][1]);

  S[1]->SetBinContent(centDep?2:1,y_metafitErr[0][0][0]/y[0][0][0]);//take average of hi and lo errors
  S[1]->SetBinContent(centDep?1:2,y_metafitErr[0][0][1]/y[0][0][1]);
  S[1]->SetBinContent(3,y_metafitErr[1][0][0]/y[1][0][0]);
  S[1]->SetBinContent(4,y_metafitErr[1][0][1]/y[1][0][1]);
  S[1]->SetBinContent(5,y_metafitErr[2][0][0]/y[2][0][0]);
  S[1]->SetBinContent(6,y_metafitErr[2][0][1]/y[2][0][1]);

  S[2]->SetBinContent(centDep?2:1,yacceffErr[0][0][0]/y[0][0][0]);
  S[2]->SetBinContent(centDep?1:2,yacceffErr[0][0][1]/y[0][0][1]);
  S[2]->SetBinContent(3,yacceffErr[1][0][0]/y[1][0][0]);
  S[2]->SetBinContent(4,yacceffErr[1][0][1]/y[1][0][1]);
  S[2]->SetBinContent(5,yacceffErr[2][0][0]/y[2][0][0]);
  S[2]->SetBinContent(6,yacceffErr[2][0][1]/y[2][0][1]);

  S[3]->SetBinContent(centDep?2:1,yTnPErr[0][0][0]/y[0][0][0]);
  S[3]->SetBinContent(centDep?1:2,yTnPErr[0][0][1]/y[0][0][1]);
  S[3]->SetBinContent(3,yTnPErr[1][0][0]/y[1][0][0]);
  S[3]->SetBinContent(4,yTnPErr[1][0][1]/y[1][0][1]);
  S[3]->SetBinContent(5,yTnPErr[2][0][0]/y[2][0][0]);
  S[3]->SetBinContent(6,yTnPErr[2][0][1]/y[2][0][1]);

  S[4]->SetBinContent(centDep?2:1,y_BcTauErr[0][1][0]/y[0][0][0]); //take the lower error here
  S[4]->SetBinContent(centDep?1:2,y_BcTauErr[0][1][1]/y[0][0][1]);
  S[4]->SetBinContent(3,y_BcTauErr[1][1][0]/y[1][0][0]);
  S[4]->SetBinContent(4,y_BcTauErr[1][1][1]/y[1][0][1]);
  S[4]->SetBinContent(5,y_BcTauErr[2][1][0]/y[2][0][0]);
  S[4]->SetBinContent(6,y_BcTauErr[2][1][1]/y[2][0][1]);

  S[5]->SetBinContent(centDep?2:1,y_LumiErr[0][0][0]/y[0][0][0]);
  S[5]->SetBinContent(centDep?1:2,y_LumiErr[0][0][1]/y[0][0][1]);
  S[5]->SetBinContent(3,y_LumiErr[1][0][0]/y[1][0][0]);
  S[5]->SetBinContent(4,y_LumiErr[1][0][1]/y[1][0][1]);
  S[5]->SetBinContent(5,y_LumiErr[2][0][0]/y[2][0][0]);
  S[5]->SetBinContent(6,y_LumiErr[2][0][1]/y[2][0][1]);

  S[6]->SetBinContent(centDep?2:1,yErr[0][0][0][0]/y[0][0][0]);//take average of hi and lo errors
  S[6]->SetBinContent(centDep?1:2,yErr[0][0][0][1]/y[0][0][1]);
  S[6]->SetBinContent(3,yErr[1][0][0][0]/y[1][0][0]);
  S[6]->SetBinContent(4,yErr[1][0][0][1]/y[1][0][1]);
  S[6]->SetBinContent(5,yErr[2][0][0][0]/y[2][0][0]);
  S[6]->SetBinContent(6,yErr[2][0][0][1]/y[2][0][1]);

  if(centDep){
    S[0]->SetBinContent(7,yfitErr[3][0][0][0]/y[3][0][0]);
    S[1]->SetBinContent(7,y_metafitErr[3][0][0]/y[3][0][0]);
    S[2]->SetBinContent(7,yacceffErr[3][0][0]/y[3][0][0]);
    S[3]->SetBinContent(7,yTnPErr[3][0][0]/y[3][0][0]);
    S[4]->SetBinContent(7,y_BcTauErr[3][1][0]/y[3][0][0]);
    S[5]->SetBinContent(7,y_LumiErr[3][0][0]/y[3][0][0]);
    S[6]->SetBinContent(7,yErr[3][0][0][0]/y[3][0][0]);
  }

  for(int i=0;i<nerr;i++){
    int idx = Sorder[i];
    S[idx]->SetLineColor(col[idx]);
    S[idx]->SetFillColor(col[idx]);
    S[idx]->SetLineWidth(2);
    S[idx]->GetXaxis()->SetBinLabel(centDep?2:1,centDep?"PbPb integrated":"pp bin1");
    S[idx]->GetXaxis()->SetBinLabel(centDep?1:2,centDep?"pp integrated":"pp bin2");
    S[idx]->GetXaxis()->SetBinLabel(3,"PbPb bin1");
    S[idx]->GetXaxis()->SetBinLabel(4,"PbPb bin2");
    S[idx]->GetXaxis()->SetBinLabel(5,"R_{PbPb} bin1");
    S[idx]->GetXaxis()->SetBinLabel(6,"R_{PbPb} bin2");
    if(centDep)
      S[idx]->GetXaxis()->SetBinLabel(7,"R_{PbPb} integrated");      
    S[idx]->GetXaxis()->SetLabelSize(0.05);
    S[idx]->GetYaxis()->SetTitle("relative error");
    S[idx]->GetYaxis()->SetTitleSize(0.055);
    S[idx]->GetYaxis()->SetLabelSize(0.04);
    S[idx]->GetYaxis()->SetTitleOffset(1.);
    S[idx]->GetYaxis()->SetRangeUser(0,1.15*S[nerr-1]->GetMaximum());
    S[idx]->GetXaxis()->SetRangeUser(-1,centDep?7:6);
    S[idx]->SetBarWidth(1/(nerr+2.));
    S[idx]->SetBarOffset((i+1)/(nerr+2.));
  }
  
  gStyle->SetOptStat(0);

  TCanvas *c3 = new TCanvas("c3","Syst",2000,1500);
  c3->SetBottomMargin(centDep?0.13:0.1);
  c3->SetTopMargin(0.04);
  c3->SetRightMargin(centDep?0.11:0.07);
  c3->SetLeftMargin(0.12);
  S[0]->SetTitle("");
  for(int i=0;i<nerr;i++)
    S[Sorder[i]]->Draw((TString)((i==0)?"bar":"barsame"));

  TLegend *leg3 = new TLegend(0.15,0.5,0.25,0.95);
  for(int i=0;i<nerr;i++)
    leg3->AddEntry(S[Sorder[i]], prettyname[Sorder[i]]);
  leg3->SetHeader("Uncertainty sources");
  leg3->SetTextSize(0.043);
  leg3->SetBorderSize(0);
  leg3->Draw("same");

  c3->SaveAs("SystematicsSummary"+(TString)(centDep?"_centralityBins":"")+".pdf");

  //Output table for all uncertainties
  cout<<Form("%25s %13s %13s %13s %13s %13s %13s","Uncertainty source",centDep?"PbPb integrated":"pp bin 1",centDep?"pp integrated":"pp bin 2","PbPb bin 1","PbPb bin 2","R_PbPb bin 1","R_PbPb bin 2")<<endl;
  cout<<Form("%25s %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f","Fit",S[0]->GetBinContent(1),S[0]->GetBinContent(2),S[0]->GetBinContent(3),S[0]->GetBinContent(4),S[0]->GetBinContent(5),S[0]->GetBinContent(6))<<endl;
  cout<<Form("%25s %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f","Fit method",S[1]->GetBinContent(1),S[1]->GetBinContent(2),S[1]->GetBinContent(3),S[1]->GetBinContent(4),S[1]->GetBinContent(5),S[1]->GetBinContent(6))<<endl;
  cout<<Form("%25s %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f","Acceptance and efficiency",S[2]->GetBinContent(1),S[2]->GetBinContent(2),S[2]->GetBinContent(3),S[2]->GetBinContent(4),S[2]->GetBinContent(5),S[2]->GetBinContent(6))<<endl;
  cout<<Form("%25s %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f","Tag and probe",S[3]->GetBinContent(1),S[3]->GetBinContent(2),S[3]->GetBinContent(3),S[3]->GetBinContent(4),S[3]->GetBinContent(5),S[3]->GetBinContent(6))<<endl;
  cout<<Form("%25s %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f","Bc -> Jpsi tau nu",S[4]->GetBinContent(1),S[4]->GetBinContent(2),S[4]->GetBinContent(3),S[4]->GetBinContent(4),S[4]->GetBinContent(5),S[4]->GetBinContent(6))<<endl;
  cout<<Form("%25s %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f","Lumi",S[5]->GetBinContent(1),S[5]->GetBinContent(2),S[5]->GetBinContent(3),S[5]->GetBinContent(4),S[5]->GetBinContent(5),S[5]->GetBinContent(6))<<endl;
  cout<<Form("%25s %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f","Total",S[6]->GetBinContent(1),S[6]->GetBinContent(2),S[6]->GetBinContent(3),S[6]->GetBinContent(4),S[6]->GetBinContent(5),S[6]->GetBinContent(6))<<endl;

}

void Draw_XSandRAA(){

  //Create Hub gathering various data 
  Hub H = Hub(true,true);
  H.SetFit(true);
  H.SetMetafit();
  H.SetAccEffFinal(true);
  H.SetAccEff();
  H.SetTnP();
  H.SetBcTau();
  H.SetLumi();
  H.SetxLW();
  H.SetFullErr();
  H.ScaleByLumi();

  cout<<"\npT dependence\n\n"<<endl;
  XSandRAA(H,false);
  cout<<"\ncentrality dependence\n\n"<<endl;
  XSandRAA(H,true);
}
