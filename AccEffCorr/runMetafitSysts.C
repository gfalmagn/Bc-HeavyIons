#include "../helpers/SystematicsList.h"

void runMetafitSysts(bool secondStep=false){
  gROOT->ProcessLine(".L CorrectPostfitYields.C++");
  gSystem->Exec("rm corrected_yields"+(TString)(secondStep?"_2ndStep":"")+".root");

  for(int isys=0;isys<_nMetafitSyst;isys++){
    //pp
    cout<<"\n        *********     CorrectPostfitYields(true,"+(TString)(secondStep?"true":"false")+",\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)"<<endl;
    gROOT->ProcessLine("CorrectPostfitYields(true,"+(TString)(secondStep?"true":"false")+",\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)");
    //PbPb
    cout<<"\n        *********     CorrectPostfitYields(false,"+(TString)(secondStep?"true":"false")+",\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)"<<endl;
    gROOT->ProcessLine("CorrectPostfitYields(false,"+(TString)(secondStep?"true":"false")+",\""+systName[isys]+"\",\""+systExtName[isys]+"\",false)");    
  }

}
