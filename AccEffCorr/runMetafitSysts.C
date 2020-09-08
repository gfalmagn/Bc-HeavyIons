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

  // //PbPb
  // cout<<"\n        *********     CorrectPostfitYields(false,\"\",\"\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"\",\"\",false)");
  // cout<<"\n        *********     CorrectPostfitYields(false,\"\",\"_noBDT1\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"\",\"_noBDT1\",false)");
  // cout<<"\n        *********     CorrectPostfitYields(false,\"_BDTuncorrFromM\",\"\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"_BDTuncorrFromM\",\"\",false)");
  // cout<<"\n        *********     CorrectPostfitYields(false,\"_BDTuncorrFromM\",\"_noBDT1\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"_BDTuncorrFromM\",\"_noBDT1\",false)");
  // cout<<"\n        *********     CorrectPostfitYields(false,\"_regulLowStatShapes\",\"_autoMCstatsNoBDT23\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"_regulLowStatShapes\",\"_autoMCstatsNoBDT23\",false)");
  // cout<<"\n        *********     CorrectPostfitYields(false,\"_regulLowStatShapes\",\"_autoMCstatsNoBDT3\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"_regulLowStatShapes\",\"_autoMCstatsNoBDT3\",false)");
  // cout<<"\n        *********     CorrectPostfitYields(false,\"_scaleSystBDTintegrated_regulLowStatShapes\",\"_autoMCstatsNoBDT23\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"_scaleSystBDTintegrated_regulLowStatShapes\",\"_autoMCstatsNoBDT23\",false)");
  // cout<<"\n        *********     CorrectPostfitYields(false,\"_scaleSystBDTintegrated_regulLowStatShapes\",\"_autoMCstatsNoBDT3\",false)"<<endl;
  // gROOT->ProcessLine("CorrectPostfitYields(false,\"_scaleSystBDTintegrated_regulLowStatShapes\",\"_autoMCstatsNoBDT3\",false)");

}
