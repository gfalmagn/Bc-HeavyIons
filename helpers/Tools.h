
float quadSum(float f1, float f2, float f3=0, float f4=0){
  return sqrt(f1*f1+f2*f2+f3*f3+f4*f4);
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
  return h->GetBinContent( 1 + (int)(n * (pt-_BcPtmin[0])/(_BcPtmax[0]-_BcPtmin[0])) );
}
