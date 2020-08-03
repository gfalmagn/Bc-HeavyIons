#include "TString.h"
#include "TColor.h"
#include "TStyle.h"

const int _nMetafitSyst = 8;//!!!hard-coded
TString systName[] = {"","","_BDTuncorrFromM","_BDTuncorrFromM","_regulLowStatShapes","_regulLowStatShapes","_scaleSystBDTintegrated_regulLowStatShapes","_scaleSystBDTintegrated_regulLowStatShapes"};
TString systExtName[] = {"","_noBDT1","","_noBDT1","_autoMCstatsNoBDT23","_autoMCstatsNoBDT3","_autoMCstatsNoBDT23","_autoMCstatsNoBDT3"};
TString systPrettyName[] = {"nominal","w/o BDT bin1","BDT uncorrelated from mass","BDT uncorrelated from mass, w/o BDT bin1","regularized low-stats shapes, no MC stats BDT2-3","regularized low-stats shapes, no MC stats BDT3","BDT-integrated shape scaling, no MC stats BDT2-3","BDT-integrated shape scaling, no MC stats BDT3"};
Color_t systCol[] = {kRed,kMagenta, kGreen,kGreen+3, kBlue,kBlue+2, kCyan,kCyan-8};
Style_t systStyle[] = {20,20, 21,25, 22,26, 23,32};
