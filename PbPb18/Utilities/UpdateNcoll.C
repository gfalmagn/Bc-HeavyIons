#include <iostream>
#include <fstream>
#include "EVENTUTILS.h"
using namespace std;

void UpdateNcoll(){

  const int nbins = 200;
  ofstream f;
  f.open ("EVENTUTILS_new.h");

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

  //Npart from DOI: 10.1103/PhysRevC.97.054910 (improved MC Glauber Loizides Kamin d'Enterria)
  double NpartMCG[20] = {384.3,331.2,283,240.9,204,171.6,143.2,118.3,96.49,77.48,61.19,47.31,35.74,26.26,18.75,13.09,9.038,6.304,4.452,3.103};
  double NpartOld[20];
  double NpartMB = 382.2;

  for(int wbin=0;wbin<20;wbin++){
    NpartOld[wbin] = 0;
    for(int bin=0;bin<10;bin++){
      NpartOld[wbin] += findNpart(wbin*nav + bin);
    }
    NpartOld[wbin] /= (float)nav;
    cout<<"NpartOld for wbin "<<wbin<<" = "<<NpartOld[wbin]<<endl;
    cout<<"NpartMCG for wbin "<<wbin<<" = "<<NpartMCG[wbin]<<endl;
  }

  f << "\nDouble_t findNpart(int hiBin) {"<<
    "\n  const int nbins = "<<nbins<<";"<<
    "\n  const Double_t Npart[nbins] = {";

  for(int bin=0;bin<nbins;bin++){
    std::string s(Form((bin<nbins-1)?"%.2f, ":"%.2f",  findNpart(bin) * NpartMCG[bin/nav] / NpartOld[bin/nav]));
    f << s;
  }  

  f << "};"<<
    "\n  return Npart[hiBin];"<<
    "\n};";
  
  f.close();
}
