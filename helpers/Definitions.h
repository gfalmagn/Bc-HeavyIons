bool _withTM = false; //Whether to allow one muon to be non-global in the pre-selection

float m_Jpsi = 3.096916;
float m_Bc = 6.2749;
float m_mu = 0.10566;
float m_K = 0.493677;
float m_Pi = 0.13957;

int _nbinMSR(bool ispp){return _withTM?(ispp?25:16):(ispp?16:10);}
int _nbinMCR(bool ispp){return _withTM?(ispp?5:4):(ispp?4:2);}
int _nbinM(bool ispp){return _nbinMSR(ispp) + 1 + _nbinMCR(ispp);}

float Mbins[100];//max number of mass bins = 100
float* _Mbinning(bool ispp){
  int nbinsCR = _nbinMCR(ispp);
  int nbinsSR = _nbinMSR(ispp);
  float Mstep = (m_Bc-3.3)/(float)nbinsSR;
  Mbins[0] = 3.3 - Mstep;
  for(int b=0;b<=nbinsSR;b++)
    Mbins[b+1] = 3.3 + b*Mstep;
  for(int b=1;b<=nbinsCR;b++)
    Mbins[b+1+nbinsSR] = m_Bc + b*1./(float)nbinsCR; //CR is 1GeV wide

  return Mbins;
}

std::vector<float> _BDTcuts(bool ispp){
  vector<float> res;
  res.push_back(_withTM?-0.5:-1);
  if(_withTM) res.push_back(-0.2);
  res.push_back(0.);
  res.push_back(_withTM?0.15:0.4);//ispp?0.15:0.2);
  res.push_back(_withTM?0.4:1);
  return res;
}

std::vector<float> _corrBDTcuts(bool ispp){
  vector<float> res;
  res.push_back(_withTM?(ispp?-0.35:-0.26):(ispp?-1:-0.8));
  if(_withTM) res.push_back(ispp?-0.05:0.05);
  res.push_back(_withTM?(ispp?0.1:0.23):(ispp?0.:0.32));
  res.push_back(_withTM?(ispp?0.23:0.42):(ispp?0.3:0.73));
  res.push_back(_withTM?(ispp?0.45:0.64):(ispp?1.:1.4));
  return res;
}

//String manipulation for output                                                                                                                                                                                                            
int _nChan(bool ispp){return _BDTcuts(ispp).size()-1; }

std::vector<std::string> _BDTcut_s(bool ispp, bool UncorrBDTfromM=false){
  std::vector<std::string> res2;
  char res[_nChan(ispp)+1][10];
  for(int k=0;k<_nChan(ispp)+1;k++){
    sprintf(res[k], "%.2f", (UncorrBDTfromM?_corrBDTcuts:_BDTcuts)(ispp)[k]);
    res2.push_back((std::string)res[k]);
  }
  return res2;
}

std::vector<std::string> _BDTcut_s2(bool ispp){
  std::vector<std::string> res;
  for(int k=0;k<_nChan(ispp)+1;k++){
    res.push_back((std::string)_BDTcut_s(ispp)[k]);
    if(_BDTcuts(ispp)[k]<-0.001){
      res[k].replace(2,1,"p");
      res[k].replace(0,1,"Minus");}
    else{
      res[k].replace(1,1,"p");}
  }
  return res;
}

bool isSoft(bool includeTM, bool isGlb, bool isTrk, bool TMOneStaTight, bool highP, float dxy, float dz, int nPix, int nTrk){

  if(includeTM){//pp
    return TMOneStaTight && highP && dxy<0.3 && dz<20 && nPix>0 && nTrk>5;
  } else{ //PbPb
    if(isGlb) {return isTrk && dxy<0.3 && dz<20 && nPix>0 && nTrk>5;}
    else {return TMOneStaTight && highP && dxy<0.3 && dz<20 && nPix>0 &&nTrk>5;} //if not global, return to standard softID
  }

}
