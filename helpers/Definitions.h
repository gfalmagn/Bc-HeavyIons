bool _withTM = false; //Whether to allow one muon to be non-global in the pre-selection
bool _preFitCorrAE = false;
bool _keepFirstMassBin = false;

float m_Jpsi = 3.096916;
float m_Bc = 6.2749;
float m_mu = 0.10566;
float m_K = 0.493677;
float m_Pi = 0.13957;

float L_pp = 302.260; //+-1.5% //pb-1 //brilcalc lumi -b "STABLE BEAMS" -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/5TeV/ReReco/Cert_306546-306826_5TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json --hltpath "HLT_HIL1DoubleMu0_v1" --begin 306546 --end 306826
float L_PbPb = 1.6054e-3; //+-1.9% //pb-1 //brilcalc lumi --normtag  /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt --hltpath "HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1"
float Leq_PbPb = 208*208*L_PbPb; 
float NMB_PbPb = 1.1194e+10; //+-1.26% https://indico.cern.ch/event/935265/contributions/3930641/attachments/2068596/3472117/PAG_200703_NMB.pdf
float TAA_090 = 6.27e-9; //pb-1 //6.27+-0.14 mb-1 (2.2%), HIN-19-007-pas-v5

float temp_corr_AcceffSyst_pp = 0.8;
float temp_corr_AcceffSyst_PbPb = 0.8;

float _corr_BcTauSyst = 1.;
float _XS_BcTauRelSystLo = 0.03;
float _RAA_BcTauRelSyst = 0.01;

//hard-coded, temporary, should be measured on signal MC and signal-enriched data (average of the two)
vector<vector<float> > ptAverage = {{8.75,17.8} , {9.2,16.9} , {9.0,17.3}}; 

vector<int> _nbinMSR(bool ispp){
  if(_withTM) return (ispp?vector<int>{25,23,20}:vector<int>{16,15,13});
  else return (ispp?(_preFitCorrAE?vector<int>{10,8,7}:vector<int>{16,15,13}):(
		     _preFitCorrAE?vector<int>{7,6,6}:vector<int>{10,8,7}));
}
vector<int> _nbinMCR(bool ispp){
  if(_withTM) return (ispp?vector<int>{5,5,4}:vector<int>{4,3,2});
  else return (ispp?(_preFitCorrAE?vector<int>{10,8,7}:vector<int>{4,3,2}):(
		     _preFitCorrAE?vector<int>{2,1,1}:vector<int>{2,1,1}));
}
vector<int> _nbinM(bool ispp){
  return vector<int>{_nbinMSR(ispp)[0] + (int)_keepFirstMassBin + _nbinMCR(ispp)[0],
                     _nbinMSR(ispp)[1] + (int)_keepFirstMassBin + _nbinMCR(ispp)[1],
                     _nbinMSR(ispp)[2] + (int)_keepFirstMassBin + _nbinMCR(ispp)[2]};}
float _mBcMax = 6.2;
float _mBcMin = 3.5;
float _mMax = 7.3;

float Mbins[100];//max number of mass bins = 100
float* _Mbinning(bool ispp, int BDTbin){
  int nbinsCR = _nbinMCR(ispp)[BDTbin];
  int nbinsSR = _nbinMSR(ispp)[BDTbin];
  float Mstep = (_mBcMax-_mBcMin)/(float)nbinsSR;
  if(_keepFirstMassBin) Mbins[0] = _mBcMin - Mstep;
  for(int b=0; b<=nbinsSR; b++)
    Mbins[b+(int)_keepFirstMassBin] = _mBcMin + b*Mstep;
  for(int b=(int)_keepFirstMassBin; b<=nbinsCR; b++)
    Mbins[b+(int)_keepFirstMassBin+nbinsSR] = _mBcMax + b*(_mMax-_mBcMax)/(float)nbinsCR;

  return Mbins;
}

bool isSoft(bool includeTMOST, bool isGlb, bool isTrk, bool TMOneStaTight, bool highP, float dxy, float dz, int nPix, int nTrk){

  if(includeTMOST){//soft
    return TMOneStaTight && highP && fabs(dxy)<0.3 && fabs(dz)<20 && nPix>0 && nTrk>5;
  } else{ //hybrid-soft, with or without global
    if(isGlb) {return isTrk && dxy<0.3 && dz<20 && nPix>0 && nTrk>5;}
    else {return TMOneStaTight && highP && fabs(dxy)<0.3 && fabs(dz)<20 && nPix>0 &&nTrk>5;} //if not global, return to standard softID
  }

}
