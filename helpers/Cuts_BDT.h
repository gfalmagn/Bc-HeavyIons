#include "Cuts_BDT_preliminary.h"

//String manipulation for output  
int _nChan(bool ispp){return _BDTcuts(ispp).size()-1; }

std::vector<std::string> _BDTcut_s(bool ispp, int kinBin=0, bool UncorrBDTfromM=false){
  std::vector<std::string> res2;
  char res[_nChan(ispp)+1][10];
  for(int k=0;k<_nChan(ispp)+1;k++){
    sprintf(res[k], "%.2f", _BDTcuts(ispp,kinBin,UncorrBDTfromM)[k]);
    res2.push_back((std::string)res[k]);
  }
  return res2;
}

std::vector<std::string> _BDTcut_s2(bool ispp, int kinBin=0, bool UncorrBDTfromM=false){
  std::vector<std::string> res;
  for(int k=0;k<_nChan(ispp)+1;k++){
    res.push_back((std::string)_BDTcut_s(ispp,kinBin,UncorrBDTfromM)[k]);
    if(_BDTcuts(ispp)[k]<-0.001){
      res[k].replace(2,1,"p");
      res[k].replace(0,1,"Minus");}
    else{
      res[k].replace(1,1,"p");}
  }
  return res;
}
