#include "../helpers/SystematicsList.h"

void runMetafitSysts(){
  gROOT->ProcessLine(".L CorrectPostfitYields.C++");

  for(int isys=0;isys<_nMetafitSyst;isys++){
    //pp
    cout<<"\n        *********     CorrectPostfitYields(true,\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)"<<endl;
    gROOT->ProcessLine("CorrectPostfitYields(true,\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)");
    //PbPb
    cout<<"\n        *********     CorrectPostfitYields(false,\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)"<<endl;
    gROOT->ProcessLine("CorrectPostfitYields(false,\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)");    
  }

}
