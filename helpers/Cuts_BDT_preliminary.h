#ifndef DEFINITIONS_H 
#define DEFINITIONS_H 
#include "Definitions.h"

std::vector<float> _BDTcuts(bool ispp, int kinBin=0, int centBin=0, bool secondStep=false, bool BDTuncorrFromM=false, int varyBDTbin=0){
  if(secondStep){
  if(ispp){
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
              return std::vector<float>{-1.73,0.70,1.75,3.97};}
        if(kinBin==1){
              return std::vector<float>{-1.88,0.52,1.72,4.27};}
        if(kinBin==2){
              return std::vector<float>{-1.73,0.70,1.63,4.06};}
      }
      else{
        if(kinBin==0){
              return std::vector<float>{-0.38,0.03,0.20,0.50};}
        if(kinBin==1){
              return std::vector<float>{-0.40,-0.03,0.15,0.47};}
        if(kinBin==2){
              return std::vector<float>{-0.38,0.05,0.21,0.50};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.73,0.55,1.51,3.97};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.73,0.70,1.75,3.97};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.73,0.85,1.96,3.97};}
          }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.88,0.34,1.42,4.27};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.88,0.52,1.72,4.27};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.88,0.67,2.02,4.27};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.73,0.55,1.42,4.06};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.73,0.70,1.63,4.06};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.73,0.85,1.84,4.06};}
          }
      }
      else{
        if(kinBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.38,0.00,0.16,0.50};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.38,0.03,0.20,0.50};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.38,0.05,0.24,0.50};}
          }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.40,-0.06,0.11,0.47};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.40,-0.03,0.15,0.47};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.40,-0.01,0.20,0.47};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.38,0.03,0.18,0.50};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.38,0.05,0.21,0.50};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.38,0.07,0.25,0.50};}
          }
      }
    }
  }
  else{
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(centBin==0){
              return std::vector<float>{-1.19,1.60,2.83,4.90};}
          if(centBin==1){
              return std::vector<float>{-1.19,1.66,2.95,4.87};}
          if(centBin==2){
              return std::vector<float>{-1.73,1.42,2.65,4.69};}
        }
        if(kinBin==1){
              return std::vector<float>{-2.18,1.21,2.68,5.23};}
        if(kinBin==2){
              return std::vector<float>{-0.83,1.48,2.38,4.54};}
      }
      else{
        if(kinBin==0){
          if(centBin==0){
              return std::vector<float>{-0.46,0.02,0.22,0.54};}
          if(centBin==1){
              return std::vector<float>{-0.46,0.03,0.23,0.50};}
          if(centBin==2){
              return std::vector<float>{-0.50,0.01,0.22,0.54};}
        }
        if(kinBin==1){
              return std::vector<float>{-0.60,-0.05,0.16,0.53};}
        if(kinBin==2){
              return std::vector<float>{-0.39,0.07,0.25,0.54};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(centBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.19,1.39,2.56,4.90};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.19,1.60,2.83,4.90};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.19,1.78,3.13,4.90};}
          }
          if(centBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.19,1.45,2.68,4.87};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.19,1.66,2.95,4.87};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.19,1.84,3.22,4.87};}
          }
          if(centBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.73,1.24,2.38,4.69};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.73,1.42,2.65,4.69};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.73,1.66,2.92,4.69};}
          }
        }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-2.18,1.00,2.38,5.23};}
            if(varyBDTbin==0){
              return std::vector<float>{-2.18,1.21,2.68,5.23};}
            if(varyBDTbin==1){
              return std::vector<float>{-2.18,1.45,3.01,5.23};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.83,1.30,2.17,4.54};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.83,1.48,2.38,4.54};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.83,1.63,2.59,4.54};}
          }
      }
      else{
        if(kinBin==0){
          if(centBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.46,-0.02,0.18,0.54};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.46,0.02,0.22,0.54};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.46,0.05,0.27,0.54};}
          }
          if(centBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.46,-0.02,0.18,0.50};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.46,0.03,0.23,0.50};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.46,0.06,0.27,0.50};}
          }
          if(centBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.50,-0.03,0.18,0.54};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.50,0.01,0.22,0.54};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.50,0.05,0.27,0.54};}
          }
        }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.60,-0.10,0.12,0.53};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.60,-0.05,0.16,0.53};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.60,-0.03,0.21,0.53};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.39,0.03,0.21,0.54};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.39,0.07,0.25,0.54};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.39,0.09,0.29,0.54};}
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
              return std::vector<float>{-1.85,0.64,1.66,3.88};}
        if(kinBin==1){
              return std::vector<float>{-1.88,0.49,1.72,4.09};}
        if(kinBin==2){
              return std::vector<float>{-1.73,0.64,1.57,3.76};}
      }
      else{
        if(kinBin==0){
              return std::vector<float>{-0.39,0.02,0.19,0.47};}
        if(kinBin==1){
              return std::vector<float>{-0.40,-0.03,0.15,0.47};}
        if(kinBin==2){
              return std::vector<float>{-0.37,0.05,0.20,0.47};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.85,0.46,1.42,3.88};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.85,0.64,1.66,3.88};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.85,0.79,1.90,3.88};}
          }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.88,0.34,1.42,4.09};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.88,0.49,1.72,4.09};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.88,0.67,2.02,4.09};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.73,0.49,1.36,3.76};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.73,0.64,1.57,3.76};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.73,0.79,1.78,3.76};}
          }
      }
      else{
        if(kinBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.39,-0.01,0.15,0.47};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.39,0.02,0.19,0.47};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.39,0.04,0.23,0.47};}
          }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.40,-0.06,0.11,0.47};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.40,-0.03,0.15,0.47};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.40,-0.01,0.20,0.47};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.37,0.02,0.17,0.47};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.37,0.05,0.20,0.47};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.37,0.07,0.24,0.47};}
          }
      }
    }
  }
  else{
    if(_withTM){
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(centBin==0){
              return std::vector<float>{-1.10,1.72,2.92,5.11};}
          if(centBin==1){
              return std::vector<float>{-0.98,1.81,3.01,5.35};}
          if(centBin==2){
              return std::vector<float>{-1.58,1.60,2.74,4.87};}
        }
        if(kinBin==1){
              return std::vector<float>{-2.06,1.27,2.74,5.11};}
        if(kinBin==2){
              return std::vector<float>{-0.92,1.51,2.44,4.63};}
      }
      else{
        if(kinBin==0){
          if(centBin==0){
              return std::vector<float>{-0.44,0.04,0.23,0.56};}
          if(centBin==1){
              return std::vector<float>{-0.42,0.04,0.23,0.54};}
          if(centBin==2){
              return std::vector<float>{-0.48,0.03,0.23,0.56};}
        }
        if(kinBin==1){
              return std::vector<float>{-0.59,-0.05,0.16,0.52};}
        if(kinBin==2){
              return std::vector<float>{-0.39,0.07,0.25,0.56};}
      }
    }
    else{
      if(BDTuncorrFromM){
        if(kinBin==0){
          if(centBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.10,1.51,2.65,5.11};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.10,1.72,2.92,5.11};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.10,1.90,3.19,5.11};}
          }
          if(centBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.98,1.60,2.77,5.35};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.98,1.81,3.01,5.35};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.98,1.99,3.31,5.35};}
          }
          if(centBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-1.58,1.36,2.50,4.87};}
            if(varyBDTbin==0){
              return std::vector<float>{-1.58,1.60,2.74,4.87};}
            if(varyBDTbin==1){
              return std::vector<float>{-1.58,1.78,3.01,4.87};}
          }
        }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-2.06,1.00,2.41,5.11};}
            if(varyBDTbin==0){
              return std::vector<float>{-2.06,1.27,2.74,5.11};}
            if(varyBDTbin==1){
              return std::vector<float>{-2.06,1.51,3.07,5.11};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.92,1.36,2.23,4.63};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.92,1.51,2.44,4.63};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.92,1.66,2.65,4.63};}
          }
      }
      else{
        if(kinBin==0){
          if(centBin==0){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.44,0.00,0.19,0.56};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.44,0.04,0.23,0.56};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.44,0.07,0.28,0.56};}
          }
          if(centBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.42,0.00,0.19,0.54};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.42,0.04,0.23,0.54};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.42,0.07,0.28,0.54};}
          }
          if(centBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.48,-0.01,0.19,0.56};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.48,0.03,0.23,0.56};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.48,0.06,0.28,0.56};}
          }
        }
        if(kinBin==1){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.59,-0.09,0.12,0.52};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.59,-0.05,0.16,0.52};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.59,-0.02,0.22,0.52};}
          }
        if(kinBin==2){
            if(varyBDTbin==-1){
              return std::vector<float>{-0.39,0.04,0.21,0.56};}
            if(varyBDTbin==0){
              return std::vector<float>{-0.39,0.07,0.25,0.56};}
            if(varyBDTbin==1){
              return std::vector<float>{-0.39,0.10,0.30,0.56};}
          }
      }
    }
  }
  }
  return std::vector<float>{};
}
#endif
