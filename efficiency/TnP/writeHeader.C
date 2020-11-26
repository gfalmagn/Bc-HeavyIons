#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

void writeHeader (bool looseAcc=false) {
  char line[1024];
  ifstream* in = NULL;
  ofstream file_header("tnp_weight_lowptPbPb"+(TString)(looseAcc?"_looseAccPart":"")+".h");
  //set the header
  if(!looseAcc){
  file_header<<"#ifndef tnp_weight_lowptPbPb_h"<<endl;
  file_header<<"#define tnp_weight_lowptPbPb_h"<<endl;
  file_header<<"#include \"TMath.h\""<<endl;
  file_header<<endl;
  file_header<<"// IN THIS FILE YOU WILL FIND:"<<endl;
  file_header<<"// +++++++++++++++++++++++++++++++++++++++"<<endl;
  file_header<<"// - Trigger: (tnp_weight_trg_pbpb)"<<endl;
  file_header<<"//   * filterId = 0: Jpsi L2 filter"<<endl;
  file_header<<"//   * filterId = 1: Jpsi L3 filter"<<endl;
  file_header<<"//   * filterId = 2: Upsi L2 filter"<<endl;
  file_header<<"//   * filterId = 3: Upsi L3 filter"<<endl;
  file_header<<"//   * idx = 0:  nominal"<<endl;
  file_header<<"//   * idx = -1: syst variation, +1 sigma"<<endl;
  file_header<<"//   * idx = -2: syst variation, -1 sigma"<<endl;
  file_header<<"//   * idx = +1: stat variation, +1 sigma"<<endl;
  file_header<<"//   * idx = +2: stat variation, -1 sigma"<<endl;
  file_header<<"//   * idx = 99: syst variation, tag selection"<<endl;
  file_header<<"// - MuID: (tnp_weight_muid_pbpb)"<<endl;
  file_header<<"//   * idx = 0:  nominal"<<endl;
  file_header<<"//   * idx = -1: syst variation, +1 sigma"<<endl;
  file_header<<"//   * idx = -2: syst variation, -1 sigma"<<endl;
  file_header<<"//   * idx = +1: stat variation, +1 sigma"<<endl;
  file_header<<"//   * idx = +2: stat variation, -1 sigma"<<endl;
  file_header<<"//   * idx = 99: syst variation, tag selection"<<endl;
  file_header<<"// - Inner tracking: (tnp_weight_trk_pbpb)"<<endl;
  file_header<<"//   * idx = 0:  nominal"<<endl;
  file_header<<"//   * idx = -1: syst variation, +1 sigma"<<endl;
  file_header<<"//   * idx = -2: syst variation, -1 sigma"<<endl;
  file_header<<"//   * idx = +1: stat variation, +1 sigma"<<endl;
  file_header<<"//   * idx = +2: stat variation, -1 sigma"<<endl;
  file_header<<"//   * idx = 99: syst variation, tag selection"<<endl;
  file_header<<"// +++++++++++++++++++++++++++++++++++++++"<<endl;
  file_header<<endl;
  file_header<<"double tnp_weight_muid_pbpb(double pt, double eta, int idx=0);"<<endl;
  file_header<<"double tnp_weight_trk_pbpb(double eta, int idx=0);"<<endl;
  file_header<<"double tnp_weight_trg_pbpb(double pt, double eta, int filterId=0,int idx=0);"<<endl;
  file_header<<endl;
  file_header<<"///////////////////////////////////////////////////"<<endl;
  file_header<<"//             M u I D    P b P b                //"<<endl;
  file_header<<"///////////////////////////////////////////////////"<<endl;
  file_header<<endl;
  file_header<<"double tnp_weight_muid_pbpb(double pt, double eta, int idx) {"<<endl;
  file_header<<"double x = pt;"<<endl;
  file_header<<"double num=1, den=1, syst=0, statUp=0, statDown=0;"<<endl;

  file_header<<"if (idx != 99) {"<<endl; // nominal tag selection

  in = new ifstream("Effplots/MuonID/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;
  file_header<<"if (fabs(eta) >= 0 && fabs(eta) < 1.2) {"<<endl;
  in = new ifstream("Effplots/MuonID/syst_data_pt_1.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {"<<endl;
  in = new ifstream("Effplots/MuonID/syst_data_pt_2.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {"<<endl;
  in = new ifstream("Effplots/MuonID/syst_data_pt_3.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {"<<endl;
  in = new ifstream("Effplots/MuonID/syst_data_pt_4.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;
  file_header<<"}"<<endl; // end of nominal tag selection

  file_header<<"else {"<<endl; // syst tag selection
  in = new ifstream("Effplots/MuonIDTag3/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl; // end of syst tag selection
  file_header<<endl;

  file_header<<"double syst_factor = 0; double stat_factor = 0;"<<endl;
  file_header<<"if (idx == 3) return den;"<<endl; 
  file_header<<"else if (idx == -1) syst_factor = syst;"<<endl;
  file_header<<"else if (idx == -2) syst_factor = -1*syst;"<<endl;
  file_header<<"else if (idx == +1) stat_factor = statUp;"<<endl;
  file_header<<"else if (idx == +2) stat_factor = -1*statDown;"<<endl;
  file_header<<"return ((num+syst_factor+stat_factor)/den);"<<endl;
  file_header<<"}"<<endl; //end of tnp_weight_muid_pbpb function
  file_header<<endl;

  file_header<<"///////////////////////////////////////////////////"<<endl;
  file_header<<"//              T R G     P b P b                //"<<endl;
  file_header<<"///////////////////////////////////////////////////"<<endl;
  file_header<<endl;
  file_header<<"double tnp_weight_trg_pbpb(double pt, double eta, int filterId,int idx) {"<<endl;
  file_header<<"double x = pt;"<<endl;
  file_header<<"double num=1, den=1, syst=0, statUp=0, statDown=0;"<<endl;

  file_header<<"if (filterId==0) { //L2 Jpsi"<<endl;
  file_header<<"if (idx != 99) {"<<endl; // nominal tag selection
  in = new ifstream("Effplots/Trg/L2Jpsi/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;
  file_header<<"if (fabs(eta) >= 0 && fabs(eta) < 1.2) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_1.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_2.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_3.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_4.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;
  file_header<<"}"<<endl; //end of nominal tag selection

  file_header<<"else {"<<endl; // syst tag selection
  in = new ifstream("Effplots/TrgTag3/L2Jpsi/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;
  file_header<<"}"<<endl; //end of syst tag selection 
  file_header<<"}"<<endl; //end of L2 Jpsi if

  file_header<<"else if (filterId==1) { //L3 Jpsi"<<endl;
  file_header<<"if (idx != 99) {"<<endl; // nominal tag selection
  in = new ifstream("Effplots/Trg/L3Jpsi/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;
  file_header<<"if (fabs(eta) >= 0 && fabs(eta) < 1.2) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_1.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_2.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_3.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_4.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;
  file_header<<"}"<<endl; //end of nominal tag selection

  file_header<<"else {"<<endl; // syst tag selection
  in = new ifstream("Effplots/TrgTag3/L3Jpsi/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;
  file_header<<"}"<<endl; //end of syst tag selection 
  file_header<<"}"<<endl; //end of L3 Jpsi if

  //////////////////////////////// for now use Jpsi for Upsi
  file_header<<"else if (filterId==2) { //L2 Jpsi"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;
  file_header<<"if (fabs(eta) >= 0 && fabs(eta) < 1.2) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_1.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_2.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_3.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {"<<endl;
  in = new ifstream("Effplots/Trg/L2Jpsi/syst_data_pt_4.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;
  file_header<<"}"<<endl; //end of L2 Upsi if //////////////////////////////// for now use Jpsi for Upsi

  file_header<<"else if (filterId==3) { //L3 Jpsi"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/correction_binned.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;
  file_header<<"if (fabs(eta) >= 0 && fabs(eta) < 1.2) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_1.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_2.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_3.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;

  file_header<<"else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {"<<endl;
  in = new ifstream("Effplots/Trg/L3Jpsi/syst_data_pt_4.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl;
  file_header<<"}"<<endl; //end of L3 Upsi if


  file_header<<"double syst_factor = 0; double stat_factor = 0;"<<endl;
  file_header<<"if (idx == 3) return den;"<<endl; 
  file_header<<"else if (idx == -1) syst_factor = syst;"<<endl;
  file_header<<"else if (idx == -2) syst_factor = -1*syst;"<<endl;
  file_header<<"else if (idx == +1) stat_factor = statUp;"<<endl;
  file_header<<"else if (idx == +2) stat_factor = -1*statDown;"<<endl;
  file_header<<"return ((num+syst_factor+stat_factor)/den);"<<endl;

  file_header<<"}"<<endl; //end of tnp_weight_trg_pbpb function 
  file_header<<endl;

  file_header<<"///////////////////////////////////////////////////"<<endl;
  file_header<<"//              T R K     P b P b                //"<<endl;
  file_header<<"///////////////////////////////////////////////////"<<endl;
  file_header<<endl;
  file_header<<"double tnp_weight_trk_pbpb(double eta, int idx) {"<<endl;
  file_header<<"double x = eta;"<<endl;
  file_header<<"double num=1, den=1, syst=0, statUp=0, statDown=0;"<<endl;

  file_header<<"if (idx != 99) {"<<endl; // nominal tag selection
  in = new ifstream("Effplots/Trk/EtaValues_Trk.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<endl;

  in = new ifstream("Effplots/Trk/syst_data_eta.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl; // end of nominal tag selection
  file_header<<"else {"<<endl; // syst tag selection
  in = new ifstream("Effplots/TrkTag3/EtaValues_Trk.txt");
  while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
  file_header<<"}"<<endl; // end of syst tag selection

  file_header<<endl;
  file_header<<"double syst_factor = 0; double stat_factor = 0;"<<endl;
  file_header<<"if (idx == 3) return den;"<<endl;
  file_header<<"else if (idx == -1) syst_factor = syst;"<<endl;
  file_header<<"else if (idx == -2) syst_factor = -1*syst;"<<endl;
  file_header<<"else if (idx == +1) stat_factor = statUp;"<<endl;
  file_header<<"else if (idx == +2) stat_factor = -1*statDown;"<<endl;
  file_header<<"return ((num+syst_factor+stat_factor)/den);"<<endl;

  file_header<<"}"<<endl; //end of tnp_weight_trk_pbpb function 
  file_header<<endl;

  /////////
  file_header<<"#endif"<<endl;
  }
  
  //Loose Acceptance part only here
  else{
    file_header<<"///////////////////////////////////////////////////"<<endl;
    file_header<<"//      M u I D  LOOSE ACCEPTANCE  P b P b       //"<<endl;
    file_header<<"///////////////////////////////////////////////////"<<endl;
    file_header<<endl;
    file_header<<"double tnp_weight_muid_looseacceptance_pbpb(double pt, double eta, int idx) {"<<endl;
    file_header<<"double x = pt;"<<endl;
    file_header<<"double num=1, den=1, syst=0, statUp=0, statDown=0;"<<endl;

    file_header<<"if (idx != 99) {"<<endl; // nominal tag selection

    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/MuonIDLooseAcc/correction_binned.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<endl;
    file_header<<"if (fabs(eta) >= 0 && fabs(eta) < 1.1) {"<<endl;
    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/MuonIDLooseAcc/syst_data_pt_1.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<"}"<<endl;

    file_header<<"else if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {"<<endl;
    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/MuonIDLooseAcc/syst_data_pt_2.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<"}"<<endl;

    file_header<<"else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {"<<endl;
    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/MuonIDLooseAcc/syst_data_pt_3.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<"}"<<endl;

    file_header<<"else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {"<<endl;
    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/MuonIDLooseAcc/syst_data_pt_4.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<"}"<<endl;
    file_header<<"}"<<endl; // end of nominal tag selection

    file_header<<"else {"<<endl; // syst tag selection
    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/MuonIDLooseAccTag3/correction_binned.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<"}"<<endl; // end of syst tag selection
    file_header<<endl;

    file_header<<"double syst_factor = 0; double stat_factor = 0;"<<endl;
    file_header<<"if (idx == 3) return den;"<<endl;
    file_header<<"else if (idx == -1) syst_factor = syst;"<<endl;
    file_header<<"else if (idx == -2) syst_factor = -1*syst;"<<endl;
    file_header<<"else if (idx == +1) stat_factor = statUp;"<<endl;
    file_header<<"else if (idx == +2) stat_factor = -1*statDown;"<<endl;
    file_header<<"return ((num+syst_factor+stat_factor)/den);"<<endl;
    file_header<<"}"<<endl; //end of tnp_weight_muid_pbpb function
    file_header<<endl;


    file_header<<"///////////////////////////////////////////////////"<<endl;
    file_header<<"//   T R K  LOOSE ACCEPTANCE     P b P b         //"<<endl;
    file_header<<"///////////////////////////////////////////////////"<<endl;
    file_header<<endl;
    file_header<<"double tnp_weight_trk_looseacceptance_pbpb(double eta, int idx) {"<<endl;
    file_header<<"double x = eta;"<<endl;
    file_header<<"double num=1, den=1, syst=0, statUp=0, statDown=0;"<<endl;

    file_header<<"if (idx != 99) {"<<endl; // nominal tag selection
    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/TrkLooseAcc/EtaValues_Trk.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<endl;

    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/TrkLooseAcc/syst_data_eta.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<"}"<<endl; // end of nominal tag selection
    file_header<<"else {"<<endl; // syst tag selection
    in = new ifstream("/home/llr/cms/falmagne/TnP/CMSSW_10_3_2/src/MuonAnalysis/TagAndProbe/test/jpsiHI/Effplots/TrkLooseAccTag3/EtaValues_Trk.txt");
    while (in->getline(line,1024,'\n')){file_header<<line<<endl;}
    file_header<<"}"<<endl; // end of syst tag selection

    file_header<<endl;
    file_header<<"double syst_factor = 0; double stat_factor = 0;"<<endl;
    file_header<<"if (idx == 3) return den;"<<endl;
    file_header<<"else if (idx == -1) syst_factor = syst;"<<endl;
    file_header<<"else if (idx == -2) syst_factor = -1*syst;"<<endl;
    file_header<<"else if (idx == +1) stat_factor = statUp;"<<endl;
    file_header<<"else if (idx == +2) stat_factor = -1*statDown;"<<endl;
    file_header<<"return ((num+syst_factor+stat_factor)/den);"<<endl;

    file_header<<"}"<<endl; //end of tnp_weight_trk_pbpb function 
    file_header<<endl;

  }

  file_header.close();

}
