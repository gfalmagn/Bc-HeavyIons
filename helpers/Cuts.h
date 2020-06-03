//FIDUCIAL cuts
int _NanaBins = 2;
vector<float> _BcYmin{1.3,  0.};
vector<float> _BcYmax{2.3, 2.3};
vector<float> _BcPtmin{6,  11};
vector<float> _BcPtmax{11, 50};

//BEGIN Pre-selection
float _ctauSignif_cut = 1.5;
float _ctauSignif3D_cut = 1.5;
float _alpha_cut = 0.8;
float _alpha3D_cut = 0.8;
float _vtxProb_cut = 0.01;
float _vtxProb_cutLoose = 0.001;
float _QQvtxProb_cut = 0.005;
float _QQvtxProb_cutTight = 0.02;
float _QQdca_cut = 0.3;
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
