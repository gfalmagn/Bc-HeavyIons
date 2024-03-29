
float quadSum(float f1, float f2, float f3=0, float f4=0){
  return sqrt(f1*f1+f2*f2+f3*f3+f4*f4);
}

float quadSubtract(float f1, float f2, float f3=0, float f4=0){
  float tmp = f1*f1-f2*f2-f3*f3-f4*f4;
  return ((tmp>=0)?1:(-1)) * sqrt(fabs(f1*f1-f2*f2-f3*f3-f4*f4));
}

std::vector<float> dRcorrectedForQQM(float dRjpsi,float dRmuWmi, float dRmuWpl, float QQM_correction, float mumiPt, float muplPt){
  float dRjpsi_corr = (dRjpsi<=0)?1:( TMath::ACos(1- QQM_correction*QQM_correction* (1-TMath::Cos(dRjpsi)) ) / dRjpsi );
  dRjpsi = dRjpsi_corr * dRjpsi;
  if(dRjpsi<=0) dRjpsi = 1e-8; if(dRjpsi>4) dRjpsi = 4; //check that deltaR is in [0,4]

  float mumiCorr = 0.5*(dRjpsi_corr-1) * muplPt/(muplPt+mumiPt);
  float muplCorr = mumiCorr*mumiPt/muplPt;
  float dRmuWmi_new = dRmuWmi*( 1+ mumiCorr*(1+ (pow(dRjpsi/dRjpsi_corr ,2) - pow( dRmuWpl,2))/pow(dRmuWmi ,2) ) );
  float dRmuWpl_new = dRmuWpl*( 1+ muplCorr*(1+ (pow(dRjpsi/dRjpsi_corr ,2) - pow( dRmuWmi,2))/pow(dRmuWpl ,2) ) );
  if(dRmuWmi_new<=0) dRmuWmi_new = 1e-8; if(dRmuWmi_new>4) dRmuWmi_new = 4; if(std::isnan(dRmuWmi_new)) dRmuWmi_new = dRmuWmi; //check that deltaR is in [0,4]
  if(dRmuWpl_new<=0) dRmuWpl_new = 1e-8; if(dRmuWpl_new>4) dRmuWpl_new = 4; if(std::isnan(dRmuWpl_new)) dRmuWpl_new = dRmuWpl; //check that deltaR is in [0,4]

  return std::vector<float>{dRjpsi,dRmuWmi_new,dRmuWpl_new};
}

void AddTH1(TH1F* h, TH1F* h2){ //add h2 to h, even with different binnings /complexity ~nb*nb2
  int nb = h->GetNbinsX();
  int nb2 = h2->GetNbinsX();
  if(h->GetBinLowEdge(nb+1)!=h2->GetBinLowEdge(nb2+1)) cout<<"AddTH1 WARNING! The two histograms to be added have different ranges! Upper range difference: "<<h->GetBinLowEdge(nb+1)-h2->GetBinLowEdge(nb+1)<<endl;
  for(int b2=0; b2<=nb2+1;b2++){
    if(b2==0 || b2==nb2+1){
      int b = (b2==0)?0:(nb+1); 
      h->SetBinContent(b, h->GetBinContent(b)+h2->GetBinContent(b2));
      h->SetBinError(b, quadSum(h->GetBinError(b),h2->GetBinError(b2)));
    }
    else{
      float b2l = h2->GetBinLowEdge(b2), b2h = h2->GetBinLowEdge(b2+1), b2w = h2->GetBinWidth(b2);
      
      for(int b=1; b<=nb;b++){
	float bl = h->GetBinLowEdge(b);
	float bh = h->GetBinLowEdge(b+1);
	float widthPart = -1;
	if((bh>=b2l && bh<b2h) || (bl>=b2l && bl<b2h)){//lower or higher bin edge of h is within the bin of h2
	  if(bh>=b2l && bh<b2h && bl<b2l) //bin b is across the lower limit b2l
	    widthPart = bh-b2l;
	  else if(bl>=b2l && bl<b2h && bh>=b2h) //bin b is across the upper limit b2h
	    widthPart = b2h-bl;
	  else if(bl>=b2l && bh<b2h) //bin b is fully within b2
	    widthPart = bh-bl;
	}
	else if(bl<b2l && bh>=b2h) //lower AND higher bin edge of h2 is within the bin of h
	  widthPart = b2w;
	
	if(widthPart>=0){
	  h->SetBinContent(b, h->GetBinContent(b) + h2->GetBinContent(b2) * widthPart/b2w );
	  h->SetBinError(b, quadSum(h->GetBinError(b), h2->GetBinError(b2) * sqrt(widthPart/b2w)) );
	}
      
      }
    }
  }
}

double getBias(TH1F* h, double pt){
  int n = h->GetNbinsX();
  if(pt<_BcPtmin[0]){
    float y1 = h->GetBinContent(1);
    float y4 = h->GetBinContent(4);
    float x = (int)(n * (pt-_BcPtmin[0])/(_BcPtmax[0]-_BcPtmin[0]) );
    return y1 + (y4-y1)*x/3;
  }
  else
    return h->GetBinContent( 1 + (int)(n * (pt-_BcPtmin[0])/(_BcPtmax[0]-_BcPtmin[0])) );
}

double MaxVec(vector<float> v, float lastitem, int maxsize  = 1e9, vector<float> v2 = vector<float>(0)){
  float res = -1e20;
  if (lastitem>res) res = lastitem;
  for(int i=0;i<(int)v.size()-1;i++){
    if(i>=maxsize) break;
    double denom = 1.;
    if (v2.size()!=0) denom = v2[i];
    if(v[i]/denom > res) res = v[i]/denom;
  }
  return res;
}

double MinVec(vector<float> v, float lastitem, int maxsize = 1e9, vector<float> v2 = vector<float>(0)){
  float res = 1e20;
  if (lastitem<res) res = lastitem;
  for(int i=0;i<(int)v.size()-1;i++){
    if(i>=maxsize) break;
    double denom = 1.;
    if (v2.size()!=0) denom = v2[i];
    if(v[i]/denom < res) res = v[i]/denom;
  }
  return res;
}

double SimpleMean(vector<float> v, int firstit, int lastit){
  int n=lastit-firstit+1;

  double res=0;
  for(int i=firstit; i<=lastit;i++){
    res += v[i];
  }
  return res/n;
}

vector<double> DoubleSidedRMS(vector<float> v, int firstit, int lastit, double mean){
  int nLo=0, nHi=0, nTot=0;
  double resLo=0, resHi=0, resTot=0;
  for(int i=firstit; i<=lastit;i++){
    resTot += pow(v[i] - mean,2);
    nTot += 1;
    if(v[i]>mean){
      resHi += pow(v[i] - mean,2);
      nHi += 1;
    }
    else{
      resLo += pow(v[i] - mean,2);
      nLo += 1;
    }
  }

  vector<double> res;
  res.push_back(mean);
  res.push_back(sqrt(resTot/(nTot-1))); //we do not give it exactly the mean of the elements, so it could be n and not n-1, but whatever
  res.push_back(sqrt(resLo/(nLo-1)));
  res.push_back(sqrt(resHi/(nHi-1)));

  return res;
}

double UncorrelatedError(double corrfact, double errMain, double errOther){
  double correrr2 = errMain*errOther* 2*fabs(corrfact) / (1+pow(errOther/errMain,2));
  //cout<<"rho, sigmaMain, sigmaOther, fraction of correlated error (squared, simple) = "<<corrfact<<" "<<errMain<<" "<<errOther<<" "<<correrr2/pow(errMain,2)<<" "<<sqrt(correrr2)/errMain<<endl;
  return sqrt(pow(errMain,2) - correrr2);
}

int _nMCclos=2;
float MCclosurePTw(float pt, int itoy){
  if (itoy==0) return pow(pt/11. , 1.7);
  else return pow(pt/11. , -1.7);
}
