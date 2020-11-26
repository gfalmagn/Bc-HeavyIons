#include "Definitions.h"

std::vector<float> _BDTcuts(bool ispp, int kinBin=0, bool BDTuncorrFromM=false, int varyBDTbin=0){
  if(ispp){
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
            return std::vector<float>{-0.55,0.00,0.17,0.49};}
        if(kinBin==1){
            return std::vector<float>{-0.35,0.06,0.23,0.64};}
        if(kinBin==2){
            return std::vector<float>{-0.47,0.05,0.20,0.51};}
      }
      else{
        if(kinBin==0){
            return std::vector<float>{-0.55,0.00,0.17,0.49};}
        if(kinBin==1){
            return std::vector<float>{-0.50,-0.03,0.14,0.49};}
        if(kinBin==2){
            return std::vector<float>{-0.55,0.03,0.19,0.48};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.55,0.00,0.14,0.49};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.55,0.00,0.17,0.49};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.55,0.00,0.21,0.49};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.35,0.06,0.19,0.64};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.35,0.06,0.23,0.64};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.35,0.06,0.27,0.64};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.47,0.05,0.17,0.51};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.47,0.05,0.20,0.51};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.47,0.05,0.23,0.51};}
        }
      }
      else{
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.55,0.00,0.14,0.49};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.55,0.00,0.17,0.49};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.55,0.00,0.21,0.49};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.50,-0.03,0.10,0.49};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.50,-0.03,0.14,0.49};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.50,-0.03,0.19,0.49};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.55,0.03,0.16,0.48};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.55,0.03,0.19,0.48};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.55,0.03,0.22,0.48};}
        }
      }
    }
  }
  else{
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
            return std::vector<float>{-0.64,0.01,0.22,0.54};}
        if(kinBin==1){
            return std::vector<float>{-0.45,0.11,0.31,0.78};}
        if(kinBin==2){
            return std::vector<float>{-0.49,0.22,0.41,0.80};}
      }
      else{
        if(kinBin==0){
            return std::vector<float>{-0.64,0.01,0.22,0.54};}
        if(kinBin==1){
            return std::vector<float>{-0.64,-0.08,0.14,0.56};}
        if(kinBin==2){
            return std::vector<float>{-0.62,0.05,0.24,0.53};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,0.01,0.17,0.54};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,0.01,0.22,0.54};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,0.01,0.26,0.54};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.45,0.11,0.26,0.78};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.45,0.11,0.31,0.78};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.45,0.11,0.36,0.78};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.49,0.22,0.37,0.80};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.49,0.22,0.41,0.80};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.49,0.22,0.46,0.80};}
        }
      }
      else{
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,0.01,0.17,0.54};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,0.01,0.22,0.54};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,0.01,0.26,0.54};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,-0.08,0.08,0.56};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,-0.08,0.14,0.56};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,-0.08,0.19,0.56};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.62,0.05,0.20,0.53};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.62,0.05,0.24,0.53};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.62,0.05,0.28,0.53};}
        }
      }
    }
  }
  return std::vector<float>{};
}
