#include "TString.h"
#include "TColor.h"
#include "TStyle.h"

const int _nMetafitSyst = 12;//!!!hard-coded //without the average nominal
TString systName[] = {"","","_BDTuncorrFromM","_BDTuncorrFromM","_regulLowStatShapes","_regulLowStatShapes","_scaleSystBDTintegrated_regulLowStatShapes","_scaleSystBDTintegrated_regulLowStatShapes","_MbinsVar1","_MbinsVar2","_BDTbinsUp","_BDTbinsDown"};
TString systExtName[] = {"","_noBDT1","","_noBDT1","_autoMCstatsNoBDT23","_autoMCstatsNoBDT3","_autoMCstatsNoBDT23","_autoMCstatsNoBDT3","","","",""};
TString systPrettyName[] = {"nominal (old)","w/o BDT bin1","BDT uncorrelated from mass","BDT uncorrelated from mass, w/o BDT bin1","Shape regularisation, no MCstats BDT2-3","Shape regularisation, no MCstats BDT3","BDT-integrated shape norm, no MCstats BDT2-3","BDT-integrated shape norm, no MCstats BDT3","coarser mass binning","finer mass binning","thinner signal-enriched BDT bin","wider signal-enriched BDT bin","nominal (average over good methods)"}; //last one for averaged nominal
Color_t systCol[] = {kRed,kMagenta, kGreen,kGreen+3, kBlue,kBlue+2, kCyan,kCyan-8,kOrange+2,kOrange-8,kViolet+1,kViolet-8  ,kBlack}; //last one for averaged nominal
Style_t systStyle[] = {20,33, 21,25, 22,26, 23,32,34,47,39,41  ,20}; //last one for averaged nominal

bool usedInRMS[] = {0,1,1,1,1,1,1,1,0,0,0,0}; //false for nominal
int _nRMSblocks = 3;
int RMSblock[] = {0,1,2,2,3,3,3,3,0,0,0,0}; //needs to be in increasing order, then zero's
float usedInNominalAverage[] = {1,0,0,0,0,0,0,0,0,0,0,0};//{1,0,1,0,0.5,0.5,0.5,0.5,0,0,0,0};
