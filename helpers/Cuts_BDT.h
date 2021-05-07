#include "Cuts_BDT_preliminary.h"

//String manipulation for output  
int _nChan(bool ispp){return _BDTcuts(ispp).size()-1; }

std::vector<std::string> _BDTcut_s(bool ispp, int kinBin=0, int centBin=0, bool secondStep=false, bool UncorrBDTfromM=false, int varyBDTbin=0){
  std::vector<std::string> res2;
  char res[_nChan(ispp)+1][10];
  for(int k=0;k<_nChan(ispp)+1;k++){
    sprintf(res[k], "%.2f", _BDTcuts(ispp,kinBin,centBin,secondStep,UncorrBDTfromM,varyBDTbin)[k]);
    res2.push_back((std::string)res[k]);
  }
  return res2;
}
