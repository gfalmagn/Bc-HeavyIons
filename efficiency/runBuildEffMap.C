#include "../helpers/SystematicsList.h"

void runBuildEffMap(){
  gROOT->ProcessLine(".L BuildEffMap.C++");

  cout<<"\n        *********     BuildEffMap(true,true,true)"<<endl;
  gROOT->ProcessLine("BuildEffMap(true,true,true)");
  cout<<"\n        *********     BuildEffMap(true,true,false)"<<endl;
  gROOT->ProcessLine("BuildEffMap(true,true,false)");
  cout<<"\n        *********     BuildEffMap(true,false,true)"<<endl;
  gROOT->ProcessLine("BuildEffMap(true,false,true)");
  cout<<"\n        *********     BuildEffMap(true,false,false)"<<endl;
  gROOT->ProcessLine("BuildEffMap(true,false,false)");

  cout<<"\n        *********     BuildEffMap(false,true,true)"<<endl;
  gROOT->ProcessLine("BuildEffMap(false,true,true)");
  cout<<"\n        *********     BuildEffMap(false,true,false)"<<endl;
  gROOT->ProcessLine("BuildEffMap(false,true,false)");
  cout<<"\n        *********     BuildEffMap(false,false,true)"<<endl;
  gROOT->ProcessLine("BuildEffMap(false,false,true)");
  cout<<"\n        *********     BuildEffMap(false,false,false)"<<endl;
  gROOT->ProcessLine("BuildEffMap(false,false,false)");

}
