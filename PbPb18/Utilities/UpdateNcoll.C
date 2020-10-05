#include <iostream>
#include <fstream>
#include "EVENTUTILS.h"
using namespace std;

void UpdateNcoll(){

  const int nbins = 200;

  //Ncoll from DOI: 10.1103/PhysRevC.97.054910 (improved MC Glauber Loizides Kamin d'Enterria)
  double NcollMCG[20] = {1762,1380,1088,855.3,667.6,515.7,392.9,294.5,216.4,155.5,109.2,74.73,49.88,32.38,20.54,12.85,8.006,5.084,3.27,2.035};
  double NcollOld[20];
  double NcollMB = 382.2;

  int nav = 10;
  for(int wbin=0;wbin<20;wbin++){
    NcollOld[wbin] = 0;
    for(int bin=0;bin<10;bin++){
      NcollOld[wbin] += findNcoll(wbin*nav + bin);
    }
    NcollOld[wbin] /= (float)nav;
    cout<<"NcollOld for wbin "<<wbin<<" = "<<NcollOld[wbin]<<endl;
    cout<<"NcollMCG for wbin "<<wbin<<" = "<<NcollMCG[wbin]<<endl;
  }

  ofstream f;
  f.open ("EVENTUTILS_new.h");
  f << "\nDouble_t findNcoll(int hiBin) {"<<
    "\n  const int nbins = "<<nbins<<";"<<
    "\n  const Double_t Ncoll[nbins] = {";

  for(int bin=0;bin<nbins;bin++){
    std::string s(Form((bin<nbins-1)?"%.2f, ":"%.2f",  findNcoll(bin) * NcollMCG[bin/nav] / NcollOld[bin/nav]));
    f << s;
  }  

  f << "};"<<
    "\n  return Ncoll[hiBin];"<<
    "\n};";
  
  f.close();
}
