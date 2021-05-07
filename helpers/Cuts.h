//FIDUCIAL cuts
int _NanaBins = 2;
int _NcentBins = 2;
vector<float> _BcYmin{0., 1.3,  0.}; //don't modify spacing and dots, for subsequent automatic tables
vector<float> _BcYmax{2.3, 2.3, 2.3};
vector<float> _BcPtmin{6., 6., 11.};
vector<float> _BcPtmax{35., 11., 35.};
vector<float> _Centmin{0., 0., 20.};
vector<float> _Centmax{90., 20., 90.};

bool inFidCuts(int kinBin, float pt, float rap){
  if(kinBin>0)
    return (pt>_BcPtmin[kinBin] && pt<_BcPtmax[kinBin] && fabs(rap)>_BcYmin[kinBin] && fabs(rap)<_BcYmax[kinBin]);
  if(kinBin==0){//integrated bin
    bool inFidCuts = false;
    for(int B=1;B<=_NanaBins;B++){
      inFidCuts = inFidCuts || (pt>_BcPtmin[B] && pt<_BcPtmax[B] && fabs(rap)>_BcYmin[B] && fabs(rap)<_BcYmax[B]);}
    return inFidCuts;
  }
  return false; 
}

//BEGIN Pre-selection //Don't touch a character of this except the values, for subsequent automatic tables
float _ctauSignif_cut = 1.2;//1.5;//1.0;
float _ctauSignif3D_cut = 1.2;//1.5;//1.0;
float _alpha_cut(bool ispp=true){return ispp?0.6:0.6;}//0.6:0.5
float _alpha3D_cut(bool ispp=true){return ispp?0.3:0.3;}//0.35:0.25
float _vtxProb_cut = 0.008;//0.01;//0.007;
float _vtxProb_cutLoose = 0.001;
float _QQvtxProb_cut = 0.005;
float _QQvtxProb_cutTight = 0.02;
float _QQdca_cut = 10.0; //basically no cuts
float _dRsum_cut = 4.5;//7.5;//5.3;
float _BcCorrM_cut(bool ispp=true){return ispp?20:20;}//20:18 //no explicit floating point here, for subsequent automatic tables
//END Pre-selection

//Jpsi mass peak and sidebands regions
float JpeakLo = 0.15, JpeakHi = 0.11;
float JpeakLoT = 0.10, JpeakHiT = 0.08;
float JpeakHiBuf = 0.04, JpeakLoBuf = 0.05;
float JpeakWid = JpeakLo+JpeakHi;
float JpeakWidT = JpeakLoT+JpeakHiT;

bool inJpsiMassRange(float M, bool tight=false){
  return (m_Jpsi-(tight?JpeakLoT:JpeakLo) < M 
	  && M < m_Jpsi+(tight?JpeakHiT:JpeakHi));
}

bool inLooseMassRange(float M){
  return (m_Jpsi-0.4 < M && M < m_Jpsi+0.4);
}

bool inJpsiMassSB(float M, bool tight=false){
  return ( (m_Jpsi-(tight?(JpeakLoT+JpeakLoBuf+JpeakWidT/2):(JpeakLo+JpeakLoBuf+JpeakWid/2)) < M 
	    && M < m_Jpsi-(tight?(JpeakLoT+JpeakLoBuf):(JpeakLo+JpeakLoBuf)))
	   || (m_Jpsi+(tight?(JpeakHiT+JpeakHiBuf):(JpeakHi+JpeakHiBuf)) < M 
	       && M < m_Jpsi+(tight?(JpeakHiT+JpeakHiBuf+JpeakWidT/2):(JpeakHi+JpeakHiBuf+JpeakWid/2))) );
}
