//All functions are from true pp TnP results, in: ~/TnP/CMSSW_9_4_12/src/MuonAnalysis/TagAndProbe/test/results
//pp track->Global-muon eff
double Efficiency1(double pt, double eta, bool isMC){
  if(isMC){
    if(fabs(eta)<1.2){
      return 0.90292*(TMath::Erf((pt-3.01346)/0.97673)*TMath::Exp(-0.00000*pt))+0.07733;
    }else if (fabs(eta)>=1.2 && fabs(eta)<1.8){
      return 0.86142*(TMath::Erf((pt-1.61130)/1.07755)*TMath::Exp(-0.00000*pt))+0.10642;
    }else if (fabs(eta)>=1.8 && fabs(eta)<2.4){
      return 0.69379*(TMath::Erf((pt-0.96622)/1.19922)*TMath::Exp(-0.00000*pt))+0.28553;
    } else return 0;
  }
  else{
    if(fabs(eta)<1.2){
      return 3.09240*(TMath::Erf((pt-1.97227)/1.41718)*TMath::Exp(-0.00042*pt))-2.10974;
    }else if (fabs(eta)>=1.2 && fabs(eta)<1.8){
      return 3.04183*(TMath::Erf((pt-0.74461)/1.33635)*TMath::Exp(-0.00000*pt))-2.08233;
    }else if (fabs(eta)>=1.8 && fabs(eta)<2.4){
      return 0.85986*(TMath::Erf((pt-0.85652)/1.20367)*TMath::Exp(-0.00182*pt))+0.10307;
    }else return 0;
  }
}

//pp Global muon->HybridSoftID eff
double Efficiency2(double pt, double eta, bool isMC){
  if(isMC){
    if(fabs(eta)<1.2){
      return 0.90205*(TMath::Erf((pt-1.36986)/1.68108)*TMath::Exp(-0.00000*pt))+0.08958;
    }else if (fabs(eta)>=1.2 && fabs(eta)<1.8){
      return 0.01674*(TMath::Erf((pt-8.49396)/39.86501)*TMath::Exp(-0.00218*pt))+0.99197;
    }else if (fabs(eta)>=1.8 && fabs(eta)<2.1){
      return 0.76933*(TMath::Erf((pt-0.77297)/0.47811)*TMath::Exp(-0.00000*pt))+0.22857;
    }else if (fabs(eta)>=2.1 && fabs(eta)<2.4){
      return 0.78151*(TMath::Erf((pt-0.65113)/0.56473)*TMath::Exp(-0.00000*pt))+0.21622;
    }else return 0;
  }
  else{
    if(fabs(eta)<1.2){
      return 1.43717*(TMath::Erf((pt-0.82453)/1.93104)*TMath::Exp(-0.00000*pt))-0.44607;
    }else if (fabs(eta)>=1.2 && fabs(eta)<1.8){
      return 0.93441*(TMath::Erf((pt+5.19818)/4.65100)*TMath::Exp(-0.00000*pt))+0.05988;
    }else if (fabs(eta)>=1.8 && fabs(eta)<2.1){
      return 0.90888*(TMath::Erf((pt+4.42062)/3.42231)*TMath::Exp(-0.00000*pt))+0.08874;
    }else if (fabs(eta)>=2.1 && fabs(eta)<2.4){
      return 1.47600*(TMath::Erf((pt-0.67408)/0.51438)*TMath::Exp(-0.00000*pt))-0.47932;
    }else return 0;
  }
}

//pp global-> hybridsoft+trigger eff
double Efficiency2times3(double pt, double eta, bool isMC){
  if(isMC){
    if(fabs(eta)<1.2){
      return 0.88451*(TMath::Erf((pt-2.17252)/1.74171)*TMath::Exp(-0.00000*pt))+0.06476;
    }else if (fabs(eta)>=1.2 && fabs(eta)<1.8){
      return 3.50849*(TMath::Erf((pt+5.15056)/6.02503)*TMath::Exp(-0.00000*pt))-2.57487;
    }else if (fabs(eta)>=1.8 && fabs(eta)<2.1){
      return 0.80644*(TMath::Erf((pt-1.11158)/1.31148)*TMath::Exp(-0.01102*pt))+0.17642;
    }else if (fabs(eta)>=2.1 && fabs(eta)<2.4){
      return 0.77663*(TMath::Erf((pt-0.69756)/1.25528)*TMath::Exp(-0.00515*pt))+0.18916;
    }else return 0;
  }
  else{
    if(fabs(eta)<1.2){
      return 1.57388*(TMath::Erf((pt-1.75360)/1.68781)*TMath::Exp(-0.00105*pt))-0.61811;
    }else if (fabs(eta)>=1.2 && fabs(eta)<1.8){
      return 3.10710*(TMath::Erf((pt+8.02355)/8.30784)*TMath::Exp(-0.00001*pt))-2.17027;
    }else if (fabs(eta)>=1.8 && fabs(eta)<2.1){
      return 4.54910*(TMath::Erf((pt+0.75206)/1.89159)*TMath::Exp(-0.00180*pt))-3.58225;
    }else if (fabs(eta)>=2.1 && fabs(eta)<2.4){
      return 0.92282*(TMath::Erf((pt-0.55830)/1.07141)*TMath::Exp(-0.00463*pt))+0.04455;
    }else return 0;
  }
}

/* //pp Trigger-only */
/* double Efficiency3(double pt, double eta, bool isMC){ */
/*   return Efficiency2times3(pt,eta,isMC)/Efficiency2(pt,eta,isMC); */
/* } */
