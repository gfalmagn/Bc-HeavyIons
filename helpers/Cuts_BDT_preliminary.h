#include "Definitions.h"

std::vector<float> _BDTcuts(bool ispp, int kinBin=0, bool secondStep=false, bool BDTuncorrFromM=false, int varyBDTbin=0){
  if(secondStep){
  if(ispp){
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
            return std::vector<float>{-0.59,0.01,0.17,0.49};}
        if(kinBin==1){
            return std::vector<float>{-0.37,0.04,0.21,0.57};}
        if(kinBin==2){
            return std::vector<float>{-0.47,0.08,0.24,0.54};}
      }
      else{
        if(kinBin==0){
            return std::vector<float>{-0.59,0.01,0.17,0.49};}
        if(kinBin==1){
            return std::vector<float>{-0.55,-0.05,0.13,0.46};}
        if(kinBin==2){
            return std::vector<float>{-0.59,0.03,0.19,0.49};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.59,0.01,0.14,0.49};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.59,0.01,0.17,0.49};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.59,0.01,0.21,0.49};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.37,0.04,0.17,0.57};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.37,0.04,0.21,0.57};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.37,0.04,0.24,0.57};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.47,0.08,0.20,0.54};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.47,0.08,0.24,0.54};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.47,0.08,0.27,0.54};}
        }
      }
      else{
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.59,0.01,0.14,0.49};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.59,0.01,0.17,0.49};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.59,0.01,0.21,0.49};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.55,-0.05,0.09,0.46};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.55,-0.05,0.13,0.46};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.55,-0.05,0.17,0.46};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.59,0.03,0.15,0.49};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.59,0.03,0.19,0.49};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.59,0.03,0.23,0.49};}
        }
      }
    }
  }
  else{
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
            return std::vector<float>{-0.64,-0.01,0.19,0.55};}
        if(kinBin==1){
            return std::vector<float>{-0.41,0.13,0.33,0.84};}
        if(kinBin==2){
            return std::vector<float>{-0.44,0.25,0.44,0.81};}
      }
      else{
        if(kinBin==0){
            return std::vector<float>{-0.64,-0.01,0.19,0.55};}
        if(kinBin==1){
            return std::vector<float>{-0.60,-0.09,0.12,0.57};}
        if(kinBin==2){
            return std::vector<float>{-0.64,0.03,0.22,0.53};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,-0.01,0.14,0.55};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,-0.01,0.19,0.55};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,-0.01,0.24,0.55};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.41,0.13,0.28,0.84};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.41,0.13,0.33,0.84};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.41,0.13,0.39,0.84};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.44,0.25,0.40,0.81};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.44,0.25,0.44,0.81};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.44,0.25,0.49,0.81};}
        }
      }
      else{
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,-0.01,0.14,0.55};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,-0.01,0.19,0.55};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,-0.01,0.24,0.55};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.60,-0.09,0.07,0.57};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.60,-0.09,0.12,0.57};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.60,-0.09,0.17,0.57};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,0.03,0.18,0.53};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,0.03,0.22,0.53};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,0.03,0.26,0.53};}
        }
      }
    }
  }
  }
  else{
  if(ispp){
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
            return std::vector<float>{-0.56,0.00,0.17,0.47};}
        if(kinBin==1){
            return std::vector<float>{-0.37,0.04,0.21,0.56};}
        if(kinBin==2){
            return std::vector<float>{-0.47,0.08,0.23,0.52};}
      }
      else{
        if(kinBin==0){
            return std::vector<float>{-0.56,0.00,0.17,0.47};}
        if(kinBin==1){
            return std::vector<float>{-0.56,-0.04,0.13,0.46};}
        if(kinBin==2){
            return std::vector<float>{-0.56,0.03,0.19,0.47};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.56,0.00,0.13,0.47};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.56,0.00,0.17,0.47};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.56,0.00,0.20,0.47};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.37,0.04,0.17,0.56};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.37,0.04,0.21,0.56};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.37,0.04,0.24,0.56};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.47,0.08,0.20,0.52};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.47,0.08,0.23,0.52};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.47,0.08,0.26,0.52};}
        }
      }
      else{
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.56,0.00,0.13,0.47};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.56,0.00,0.17,0.47};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.56,0.00,0.20,0.47};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.56,-0.04,0.09,0.46};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.56,-0.04,0.13,0.46};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.56,-0.04,0.17,0.46};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.56,0.03,0.15,0.47};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.56,0.03,0.19,0.47};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.56,0.03,0.22,0.47};}
        }
      }
    }
  }
  else{
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
            return std::vector<float>{-0.64,0.00,0.21,0.55};}
        if(kinBin==1){
            return std::vector<float>{-0.41,0.13,0.34,0.84};}
        if(kinBin==2){
            return std::vector<float>{-0.46,0.25,0.45,0.87};}
      }
      else{
        if(kinBin==0){
            return std::vector<float>{-0.64,0.00,0.21,0.55};}
        if(kinBin==1){
            return std::vector<float>{-0.61,-0.09,0.12,0.57};}
        if(kinBin==2){
            return std::vector<float>{-0.64,0.03,0.23,0.55};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,0.00,0.16,0.55};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,0.00,0.21,0.55};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,0.00,0.25,0.55};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.41,0.13,0.29,0.84};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.41,0.13,0.34,0.84};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.41,0.13,0.39,0.84};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.46,0.25,0.41,0.87};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.46,0.25,0.45,0.87};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.46,0.25,0.50,0.87};}
        }
      }
      else{
        if(kinBin==0){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,0.00,0.16,0.55};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,0.00,0.21,0.55};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,0.00,0.25,0.55};}
        }
        if(kinBin==1){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.61,-0.09,0.07,0.57};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.61,-0.09,0.12,0.57};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.61,-0.09,0.17,0.57};}
        }
        if(kinBin==2){
          if(varyBDTbin==-1){
            return std::vector<float>{-0.64,0.03,0.19,0.55};}
          if(varyBDTbin==0){
            return std::vector<float>{-0.64,0.03,0.23,0.55};}
          if(varyBDTbin==1){
            return std::vector<float>{-0.64,0.03,0.27,0.55};}
        }
      }
    }
  }
  }
  return std::vector<float>{};
}
