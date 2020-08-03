#include "Definitions.h"

std::vector<float> _BDTcuts(bool ispp, int kinBin=0, bool BDTuncorrFromM=false){
  if(ispp){
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
          return std::vector<float>{-0.54,0.00,0.19,0.51};}
        if(kinBin==1){
          return std::vector<float>{-0.41,0.05,0.21,0.53};}
        if(kinBin==2){
          return std::vector<float>{-0.49,0.07,0.23,0.54};}
      }
      else{
        if(kinBin==0){
          return std::vector<float>{-0.54,0.00,0.19,0.51};}
        if(kinBin==1){
          return std::vector<float>{-0.50,-0.04,0.13,0.41};}
        if(kinBin==2){
          return std::vector<float>{-0.54,0.04,0.22,0.51};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          return std::vector<float>{-0.54,0.00,0.19,0.51};}
        if(kinBin==1){
          return std::vector<float>{-0.41,0.05,0.21,0.53};}
        if(kinBin==2){
          return std::vector<float>{-0.49,0.07,0.23,0.54};}
      }
      else{
        if(kinBin==0){
          return std::vector<float>{-0.54,0.00,0.19,0.51};}
        if(kinBin==1){
          return std::vector<float>{-0.50,-0.04,0.13,0.41};}
        if(kinBin==2){
          return std::vector<float>{-0.54,0.04,0.22,0.51};}
      }
    }
  }
  else{
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
          return std::vector<float>{-0.58,0.02,0.23,0.59};}
        if(kinBin==1){
          return std::vector<float>{-0.37,0.16,0.35,0.70};}
        if(kinBin==2){
          return std::vector<float>{-0.30,0.21,0.40,0.79};}
      }
      else{
        if(kinBin==0){
          return std::vector<float>{-0.58,0.02,0.23,0.59};}
        if(kinBin==1){
          return std::vector<float>{-0.58,-0.09,0.10,0.42};}
        if(kinBin==2){
          return std::vector<float>{-0.46,0.09,0.26,0.59};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          return std::vector<float>{-0.58,0.02,0.23,0.59};}
        if(kinBin==1){
          return std::vector<float>{-0.37,0.16,0.35,0.70};}
        if(kinBin==2){
          return std::vector<float>{-0.30,0.21,0.40,0.79};}
      }
      else{
        if(kinBin==0){
          return std::vector<float>{-0.58,0.02,0.23,0.59};}
        if(kinBin==1){
          return std::vector<float>{-0.58,-0.09,0.10,0.42};}
        if(kinBin==2){
          return std::vector<float>{-0.46,0.09,0.26,0.59};}
      }
    }
  }
  return std::vector<float>{};
}
