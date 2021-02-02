#include "../helpers/SystematicsList.h"

void runBuildEffMap(bool runAEtoys=true){
  gROOT->ProcessLine(".L BuildEffMap.C++");

  if(!runAEtoys){
    cout<<"\n        *********     BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",true,true)"<<endl;
    gROOT->ProcessLine("BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",true,true)");
    cout<<"\n        *********     BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",true,false)"<<endl;
    gROOT->ProcessLine("BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",true,false)");
    cout<<"\n        *********     BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",false,true)"<<endl;
    gROOT->ProcessLine("BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",false,true)");
  }
  cout<<"\n        *********     BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",false,false)"<<endl;
  gROOT->ProcessLine("BuildEffMap(true,"+(TString)(runAEtoys?"true":"false")+",false,false)");

  if(!runAEtoys){
    cout<<"\n        *********     BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",true,true)"<<endl;
    gROOT->ProcessLine("BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",true,true)");
    cout<<"\n        *********     BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",true,false)"<<endl;
    gROOT->ProcessLine("BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",true,false)");
    cout<<"\n        *********     BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",false,true)"<<endl;
    gROOT->ProcessLine("BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",false,true)");
  }
  cout<<"\n        *********     BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",false,false)"<<endl;
  gROOT->ProcessLine("BuildEffMap(false,"+(TString)(runAEtoys?"true":"false")+",false,false)");

}
