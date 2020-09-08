#ifndef tnp_weight_pp_efficiencies_h
#define tnp_weight_pp_efficiencies_h
#include "TMath.h"

// IN THIS FILE YOU WILL FIND:
// +++++++++++++++++++++++++++++++++++++++
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
//   * idx = 99: syst variation, tag selection
//   * idx = 3: MC efficiency 
// +++++++++++++++++++++++++++++++++++++++


double tnp_weight_glb_looseacceptance_pp(double pt, double eta, int idx=0);
double tnp_weight_muid_looseacceptance_pp(double pt, double eta, int idx=0);
double tnp_weight_glb_tightacceptance_pp(double pt, double eta, int idx=0);
double tnp_weight_muidtrg_tightacceptance_pp(double pt, double eta, int idx=0);

///////////////////////////////////////////////////
//             GLB LooseAcceptance               //
///////////////////////////////////////////////////

double tnp_weight_glb_looseacceptance_pp(double pt, double eta, int idx=0) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (idx != 99) {
    
    // SF for 0 < |eta| < 1.1
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      if (x >= 3.3 && x <3.7) {num = 0.545674; den = 0.479689; statUp = 0.00279995; statDown = 0.00278558;}
      else if (x >= 3.7 && x <3.9) {num = 0.765766; den = 0.750925; statUp = 0.00457248; statDown = 0.00455396;}
      else if (x >= 3.9 && x <4.2) {num = 0.855967; den = 0.850531; statUp = 0.00398287; statDown = 0;}
      else if (x >= 4.2 && x <4.6) {num = 0.913168; den = 0.920223; statUp = 0.00362134; statDown = 0.00360798;}
      else if (x >= 4.6 && x <5.2) {num = 0.95689; den = 0.966469; statUp = 0.00313401; statDown = 0.00313982;}
      else if (x >= 5.2 && x <7) {num = 0.980718; den = 0.98453; statUp = 0.00217621; statDown = 0.0021778;}
      else if (x >= 7 && x <10.5) {num = 0.977835; den = 0.98558; statUp = 0.00244144; statDown = 0.00243946;}
      else {num = 0.963895; den = 0.98272; statUp = 0.00370609; statDown = 0.00371316;}
    }
    // SF for 1.1 < |eta| < 1.8
    if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      if (x >= 1.75 && x <2.4) {num = 0.358092; den = 0.368398; statUp = 0.00387365; statDown = 0.00376561;}
      else if (x >= 2.4 && x <2.8) {num = 0.696959; den = 0.727652; statUp = 0.00548732; statDown = 0.00542193;}
      else if (x >= 2.8 && x <3.2) {num = 0.826358; den = 0.82493; statUp = 0.00555076; statDown = 0.00549516;}
      else if (x >= 3.2 && x <3.9) {num = 0.909702; den = 0.910116; statUp = 0.004166; statDown = 0.00415827;}
      else if (x >= 3.9 && x <5) {num = 0.96088; den = 0.968761; statUp = 0.00369013; statDown = 0.00362025;}
      else if (x >= 5 && x <8) {num = 0.974763; den = 0.976658; statUp = 0.00310129; statDown = 0.00320536;}
      else {num = 0.957222; den = 0.972933; statUp = 0.00458492; statDown = 0.00456555;}
    }
    // SF for 1.8 < |eta| < 2.4
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.4) {
      if (x >= 1.2 && x <1.8) {num = 0.470546; den = 0.493486; statUp = 0.0079263; statDown = 0.00614067;}
      else if (x >= 1.8 && x <2.3) {num = 0.779572; den = 0.80093; statUp = 0.00813203; statDown = 0.00796026;}
      else if (x >= 2.3 && x <2.9) {num = 0.917865; den = 0.939274; statUp = 0.00770755; statDown = 0.00759797;}
      else if (x >= 2.9 && x <3.7) {num = 0.948355; den = 0.969628; statUp = 0.00664304; statDown = 0.00657473;}
      else if (x >= 3.7 && x <6) {num = 0.962889; den = 0.982337; statUp = 0.0046804; statDown = 0.00469451;}
      else {num = 0.948709; den = 0.980496; statUp = 0.00543928; statDown = 0.0053934;}
    }
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      // syst uncertainties
      if (x >= 3.3 && x < 3.7) syst = 0.00219018;
      else if (x >= 3.7 && x < 3.9) syst = 0.00267876;
      else if (x >= 3.9 && x < 4.2) syst = 0.000886909;
      else if (x >= 4.2 && x < 4.6) syst = 0.0010186;
      else if (x >= 4.6 && x < 5.2) syst = 0.00058522;
      else if (x >= 5.2 && x < 7) syst = 0.00129457;
      else if (x >= 7 && x < 10.5) syst = 0.000990161;
      else syst = 0.00211349;
    }  
   
    else if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 1.75 && x < 2.4) syst = 0.0016153;
      else if (x >= 2.4 && x < 2.8) syst = 0.00177856;
      else if (x >= 2.8 && x < 3.2) syst = 0.00147102;
      else if (x >= 3.2 && x < 3.9) syst = 0.00201313;
      else if (x >= 3.9 && x < 5) syst = 0.00129165;
      else if (x >= 5 && x < 8) syst = 0.00199723;
      else syst = 0.000925577;
    }

    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.2 && x < 1.8) syst = 0.00649426;
      else if (x >= 1.8 && x < 2.3) syst = 0.00481814;
      else if (x >= 2.3 && x < 2.9) syst = 0.00191262;
      else if (x >= 2.9 && x < 3.7) syst = 0.00237759;
      else if (x >= 3.7 && x < 6) syst = 0.00150681;
      else syst = 0.0054402;
    }
  }
  else {
    // SF for 0 < |eta| < 1.1
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      if (x >= 3.3 && x <3.7) {num = 0.545674; den = 0.479689; statUp = 0.00279995; statDown = 0.00278558;}
      else if (x >= 3.7 && x <3.9) {num = 0.765766; den = 0.750925; statUp = 0.00457248; statDown = 0.00455396;}
      else if (x >= 3.9 && x <4.2) {num = 0.855967; den = 0.850531; statUp = 0.00398287; statDown = 0;}
      else if (x >= 4.2 && x <4.6) {num = 0.913168; den = 0.920223; statUp = 0.00362134; statDown = 0.00360798;}
      else if (x >= 4.6 && x <5.2) {num = 0.95689; den = 0.966469; statUp = 0.00313401; statDown = 0.00313982;}
      else if (x >= 5.2 && x <7) {num = 0.980718; den = 0.98453; statUp = 0.00217621; statDown = 0.0021778;}
      else if (x >= 7 && x <10.5) {num = 0.977835; den = 0.98558; statUp = 0.00244144; statDown = 0.00243946;}
      else {num = 0.963895; den = 0.98272; statUp = 0.00370609; statDown = 0.00371316;}
    }
    // SF for 1.1 < |eta| < 1.8
    if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      if (x >= 1.75 && x <2.4) {num = 0.358092; den = 0.368398; statUp = 0.00387365; statDown = 0.00376561;}
      else if (x >= 2.4 && x <2.8) {num = 0.696959; den = 0.727652; statUp = 0.00548732; statDown = 0.00542193;}
      else if (x >= 2.8 && x <3.2) {num = 0.826358; den = 0.82493; statUp = 0.00555076; statDown = 0.00549516;}
      else if (x >= 3.2 && x <3.9) {num = 0.909702; den = 0.910116; statUp = 0.004166; statDown = 0.00415827;}
      else if (x >= 3.9 && x <5) {num = 0.96088; den = 0.968761; statUp = 0.00369013; statDown = 0.00362025;}
      else if (x >= 5 && x <8) {num = 0.974763; den = 0.976658; statUp = 0.00310129; statDown = 0.00320536;}
      else {num = 0.957222; den = 0.972933; statUp = 0.00458492; statDown = 0.00456555;}
    }
    // SF for 1.8 < |eta| < 2.4
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.4) {
      if (x >= 1.2 && x <1.8) {num = 0.470546; den = 0.493486; statUp = 0.0079263; statDown = 0.00614067;}
      else if (x >= 1.8 && x <2.3) {num = 0.779572; den = 0.80093; statUp = 0.00813203; statDown = 0.00796026;}
      else if (x >= 2.3 && x <2.9) {num = 0.917865; den = 0.939274; statUp = 0.00770755; statDown = 0.00759797;}
      else if (x >= 2.9 && x <3.7) {num = 0.948355; den = 0.969628; statUp = 0.00664304; statDown = 0.00657473;}
      else if (x >= 3.7 && x <6) {num = 0.962889; den = 0.982337; statUp = 0.0046804; statDown = 0.00469451;}
      else {num = 0.948709; den = 0.980496; statUp = 0.00543928; statDown = 0.0053934;}
    }
  }

  double syst_factor = 0; double stat_factor = 0;
  if (idx ==3) {return den;}
  else{
    if (idx == -1) syst_factor = syst;
    else if (idx == -2) syst_factor = -1*syst;
    else if (idx == +1) stat_factor = statUp;
    else if (idx == +2) stat_factor = -1*statDown;
    return ((num+syst_factor+stat_factor)/den);
  }
}


///////////////////////////////////////////////////
//            M u I D         LooseAcceptance    //
///////////////////////////////////////////////////
double tnp_weight_muid_looseacceptance_pp(double pt, double eta, int idx=0) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (idx != 99) {
    // SF for 0 < |eta| < 1.1
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      if (x >= 3.3 && x <3.7) {num = 0.919213; den = 0.927524; statUp = 0.00130712; statDown = 0.0013169;}
      else if (x >= 3.7 && x <4) {num = 0.958336; den = 0.964014; statUp = 0.00098851; statDown = 0.00100085;}
      else if (x >= 4 && x <4.5) {num = 0.972815; den = 0.976775; statUp = 0.000658825; statDown = 0.000667667;}
      else if (x >= 4.5 && x <5) {num = 0.981845; den = 0.984312; statUp = 0.000607121; statDown = 0.000618402;}
      else if (x >= 5 && x <5.5) {num = 0.987175; den = 0.988898; statUp = 0.000585655; statDown = 0.000600501;}
      else if (x >= 5.5 && x <6) {num = 0.989389; den = 0.991812; statUp = 0.000599164; statDown = 0.000618297;}
      else if (x >= 6 && x <7) {num = 0.990073; den = 0.991969; statUp = 0.000494077; statDown = 0.000507297;}
      else if (x >= 7 && x <8) {num = 0.990831; den = 0.992648; statUp = 0.000585538; statDown = 0.00060653;}
      else if (x >= 8 && x <10.5) {num = 0.992668; den = 0.99256; statUp = 0.00046987; statDown = 0.000485967;}
      else if (x >= 10.5 && x <14) {num = 0.99413; den = 0.993753; statUp = 0.000594339; statDown = 0.000624891;}
      else if (x >= 14 && x <18) {num = 0.991855; den = 0.992597; statUp = 0.00107351; statDown = 0.00114532;}
      else {num = 0.992834; den = 0.992808; statUp = 0.00716634; statDown = 0.00138491;}
    }
    // SF for 1.1 < |eta| < 1.8
    if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      if (x >= 1.75 && x <2.5) {num = 0.975024; den = 0.976703; statUp = 0.000976144; statDown = 0.000992569;}
      else if (x >= 2.5 && x <3) {num = 0.97956; den = 0.982742; statUp = 0.0006431; statDown = 0.000652864;}
      else if (x >= 3 && x <3.5) {num = 0.980467; den = 0.983917; statUp = 0.000598349; statDown = 0.000607471;}
      else if (x >= 3.5 && x <4) {num = 0.984111; den = 0.987936; statUp = 0.00059104; statDown = 0.00060266;}
      else if (x >= 4 && x <4.5) {num = 0.98905; den = 0.991363; statUp = 0.000495048; statDown = 0.000495048;}
      else if (x >= 4.5 && x <5) {num = 0.989798; den = 0.992516; statUp = 0.000627955; statDown = 0.000648085;}
      else if (x >= 5 && x <6) {num = 0.992447; den = 0.994514; statUp = 0.000475251; statDown = 0.000490007;}
      else if (x >= 6 && x <7) {num = 0.992888; den = 0.996081; statUp = 0.0005844; statDown = 0.000609488;}
      else if (x >= 7 && x <9) {num = 0.993985; den = 0.996469; statUp = 0.000539269; statDown = 0.00056272;}
      else if (x >= 9 && x <14) {num = 0.994103; den = 0.99704; statUp = 0.000598822; statDown = 0.000629629;}
      else if (x >= 14 && x <18) {num = 0.995286; den = 0.996088; statUp = 0.00128415; statDown = 0.00143594;}
      else {num = 0.997023; den = 0.997445; statUp = 0.00144995; statDown = 0.0017078;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.33 && x <1.75) {num = 0.987365; den = 0.992475; statUp = 0.00181843; statDown = 0.0018746;}
      else if (x >= 1.75 && x <2.1) {num = 0.992106; den = 0.996675; statUp = 0.000896418; statDown = 0.000928494;}
      else if (x >= 2.1 && x <2.5) {num = 0.996755; den = 0.998413; statUp = 0.00053316; statDown = 0.000554546;}
      else if (x >= 2.5 && x <3) {num = 0.997319; den = 0.998958; statUp = 0.000447122; statDown = 0.000466366;}
      else if (x >= 3 && x <3.5) {num = 0.99755; den = 0.999041; statUp = 0.0024505; statDown = 0.000484142;}
      else if (x >= 3.5 && x <4) {num = 0.997469; den = 0.999548; statUp = 0.000507592; statDown = 0.000545444;}
      else if (x >= 4 && x <4.5) {num = 0.997685; den = 0.999644; statUp = 0.000527844; statDown = 0.000581495;}
      else if (x >= 4.5 && x <5.25) {num = 0.998079; den = 0.999339; statUp = 0.00048508; statDown = 0.000531115;}
      else if (x >= 5.25 && x <6.5) {num = 0.998826; den = 0.999537; statUp = 0.000412939; statDown = 0.000460615;}
      else if (x >= 6.5 && x <9) {num = 0.997839; den = 0.99955; statUp = 0.000493282; statDown = 0.000546225;}
      else if (x >= 9 && x <13) {num = 0.998093; den = 0.999791; statUp = 0.000690858; statDown = 0.000799614;}
      else {num = 0.997906; den = 0.999422; statUp = 0.00209421; statDown = 0.00138463;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.2 && x <1.9) {num = 0.979743; den = 0.986272; statUp = 0.00109581; statDown = 0.00111105;}
      else if (x >= 1.9 && x <2.5) {num = 0.994461; den = 0.997793; statUp = 0.00553899; statDown = 0.000611357;}
      else if (x >= 2.5 && x <3) {num = 0.996774; den = 0.999019; statUp = 0.000569912; statDown = 0.000597349;}
      else if (x >= 3 && x <3.5) {num = 0.997965; den = 0.998333; statUp = 0.000546018; statDown = 0.000584973;}
      else if (x >= 3.5 && x <4) {num = 0.997536; den = 0.998867; statUp = 0.000637923; statDown = 0.000692719;}
      else if (x >= 4 && x <4.75) {num = 0.999282; den = 0.999551; statUp = 0.000442867; statDown = 0.000497266;}
      else if (x >= 4.75 && x <5.5) {num = 0.998631; den = 0.999213; statUp = 0.000596208; statDown = 0.00067702;}
      else if (x >= 5.5 && x <8) {num = 0.997838; den = 0.999565; statUp = 0.000516111; statDown = 0.000568641;}
      else if (x >= 8 && x <12) {num = 0.99848; den = 0.998524; statUp = 0.000744772; statDown = 0.000875927;}
      else {num = 1; den = 0.999295; statUp = 2.88878e-11; statDown = 0.000509368;}
    }
    
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      // syst uncertainties
      if (x >= 3.3 && x < 3.7) syst = 7.47299e-05;
      else if (x >= 3.7 && x < 4) syst = 0.000122414;
      else if (x >= 4 && x < 4.5) syst = 0.000143316;
      else if (x >= 4.5 && x < 5) syst = 0.000136017;
      else if (x >= 5 && x < 5.5) syst = 4.67405e-05;
      else if (x >= 5.5 && x < 6) syst = 9.11727e-05;
      else if (x >= 6 && x < 7) syst = 9.55021e-05;
      else if (x >= 7 && x < 8) syst = 0.000233307;
      else if (x >= 8 && x < 10.5) syst = 0.000229503;
      else if (x >= 10.5 && x < 14) syst = 6.08627e-05;
      else if (x >= 14 && x < 18) syst = 0.000402988;
      else syst = 0.000519493;
    }
    else if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 1.75 && x < 2.5) syst = 0.00014918;
      else if (x >= 2.5 && x < 3) syst = 0.000348243;
      else if (x >= 3 && x < 3.5) syst = 0.000174388;
      else if (x >= 3.5 && x < 4) syst = 0.000366704;
      else if (x >= 4 && x < 4.5) syst = 3.12523e-05;
      else if (x >= 4.5 && x < 5) syst = 0.00026829;
      else if (x >= 5 && x < 6) syst = 0.000187476;
      else if (x >= 6 && x < 7) syst = 6.21438e-05;
      else if (x >= 7 && x < 9) syst = 0.000172745;
      else if (x >= 9 && x < 14) syst = 2.47709e-05;
      else if (x >= 14 && x < 18) syst = 0.000148418;
      else syst = 0.000222137;
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.33 && x < 1.75) syst = 0.000383119;
      else if (x >= 1.75 && x < 2.1) syst = 0.000856907;
      else if (x >= 2.1 && x < 2.5) syst = 0.000414077;
      else if (x >= 2.5 && x < 3) syst = 0.000286409;
      else if (x >= 3 && x < 3.5) syst = 6.91918e-05;
      else if (x >= 3.5 && x < 4) syst = 0.000155472;
      else if (x >= 4 && x < 4.5) syst = 0.000341682;
      else if (x >= 4.5 && x < 5.25) syst = 0.00015471;
      else if (x >= 5.25 && x < 6.5) syst = 3.73735e-05;
      else if (x >= 6.5 && x < 9) syst = 0.000171324;
      else if (x >= 9 && x < 13) syst = 0.000569324;
      else syst = 0.000158047;
    }
    else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.2 && x < 1.9) syst = 0.000156846;
      else if (x >= 1.9 && x < 2.5) syst = 0.000249898;
      else if (x >= 2.5 && x < 3) syst = 0.000332898;
      else if (x >= 3 && x < 3.5) syst = 0.000277506;
      else if (x >= 3.5 && x < 4) syst = 0.000530462;
      else if (x >= 4 && x < 4.75) syst = 9.82921e-05;
      else if (x >= 4.75 && x < 5.5) syst = 0.000436148;
      else if (x >= 5.5 && x < 8) syst = 0.000182531;
      else if (x >= 8 && x < 12) syst = 0.000374469;
      else syst = 0.000559234;
    }
  }
  else {
    // SF for 0 < |eta| < 1.1
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      if (x >= 3.3 && x <3.7) {num = 0.919213; den = 0.927524; statUp = 0.00130712; statDown = 0.0013169;}
      else if (x >= 3.7 && x <4) {num = 0.958336; den = 0.964014; statUp = 0.00098851; statDown = 0.00100085;}
      else if (x >= 4 && x <4.5) {num = 0.972815; den = 0.976775; statUp = 0.000658825; statDown = 0.000667667;}
      else if (x >= 4.5 && x <5) {num = 0.981845; den = 0.984312; statUp = 0.000607121; statDown = 0.000618402;}
      else if (x >= 5 && x <5.5) {num = 0.987175; den = 0.988898; statUp = 0.000585655; statDown = 0.000600501;}
      else if (x >= 5.5 && x <6) {num = 0.989389; den = 0.991812; statUp = 0.000599164; statDown = 0.000618297;}
      else if (x >= 6 && x <7) {num = 0.990073; den = 0.991969; statUp = 0.000494077; statDown = 0.000507297;}
      else if (x >= 7 && x <8) {num = 0.990831; den = 0.992648; statUp = 0.000585538; statDown = 0.00060653;}
      else if (x >= 8 && x <10.5) {num = 0.992668; den = 0.99256; statUp = 0.00046987; statDown = 0.000485967;}
      else if (x >= 10.5 && x <14) {num = 0.99413; den = 0.993753; statUp = 0.000594339; statDown = 0.000624891;}
      else if (x >= 14 && x <18) {num = 0.991855; den = 0.992597; statUp = 0.00107351; statDown = 0.00114532;}
      else {num = 0.992834; den = 0.992808; statUp = 0.00716634; statDown = 0.00138491;}
    }
      // SF for 1.1 < |eta| < 1.8
    if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      if (x >= 1.75 && x <2.5) {num = 0.975024; den = 0.976703; statUp = 0.000976144; statDown = 0.000992569;}
      else if (x >= 2.5 && x <3) {num = 0.97956; den = 0.982742; statUp = 0.0006431; statDown = 0.000652864;}
      else if (x >= 3 && x <3.5) {num = 0.980467; den = 0.983917; statUp = 0.000598349; statDown = 0.000607471;}
      else if (x >= 3.5 && x <4) {num = 0.984111; den = 0.987936; statUp = 0.00059104; statDown = 0.00060266;}
      else if (x >= 4 && x <4.5) {num = 0.98905; den = 0.991363; statUp = 0.000495048; statDown = 0.000495048;}
      else if (x >= 4.5 && x <5) {num = 0.989798; den = 0.992516; statUp = 0.000627955; statDown = 0.000648085;}
      else if (x >= 5 && x <6) {num = 0.992447; den = 0.994514; statUp = 0.000475251; statDown = 0.000490007;}
      else if (x >= 6 && x <7) {num = 0.992888; den = 0.996081; statUp = 0.0005844; statDown = 0.000609488;}
      else if (x >= 7 && x <9) {num = 0.993985; den = 0.996469; statUp = 0.000539269; statDown = 0.00056272;}
      else if (x >= 9 && x <14) {num = 0.994103; den = 0.99704; statUp = 0.000598822; statDown = 0.000629629;}
      else if (x >= 14 && x <18) {num = 0.995286; den = 0.996088; statUp = 0.00128415; statDown = 0.00143594;}
      else {num = 0.997023; den = 0.997445; statUp = 0.00144995; statDown = 0.0017078;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.33 && x <1.75) {num = 0.987365; den = 0.992475; statUp = 0.00181843; statDown = 0.0018746;}
      else if (x >= 1.75 && x <2.1) {num = 0.992106; den = 0.996675; statUp = 0.000896418; statDown = 0.000928494;}
      else if (x >= 2.1 && x <2.5) {num = 0.996755; den = 0.998413; statUp = 0.00053316; statDown = 0.000554546;}
      else if (x >= 2.5 && x <3) {num = 0.997319; den = 0.998958; statUp = 0.000447122; statDown = 0.000466366;}
      else if (x >= 3 && x <3.5) {num = 0.99755; den = 0.999041; statUp = 0.0024505; statDown = 0.000484142;}
      else if (x >= 3.5 && x <4) {num = 0.997469; den = 0.999548; statUp = 0.000507592; statDown = 0.000545444;}
      else if (x >= 4 && x <4.5) {num = 0.997685; den = 0.999644; statUp = 0.000527844; statDown = 0.000581495;}
      else if (x >= 4.5 && x <5.25) {num = 0.998079; den = 0.999339; statUp = 0.00048508; statDown = 0.000531115;}
      else if (x >= 5.25 && x <6.5) {num = 0.998826; den = 0.999537; statUp = 0.000412939; statDown = 0.000460615;}
      else if (x >= 6.5 && x <9) {num = 0.997839; den = 0.99955; statUp = 0.000493282; statDown = 0.000546225;}
      else if (x >= 9 && x <13) {num = 0.998093; den = 0.999791; statUp = 0.000690858; statDown = 0.000799614;}
      else {num = 0.997906; den = 0.999422; statUp = 0.00209421; statDown = 0.00138463;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.2 && x <1.9) {num = 0.979743; den = 0.986272; statUp = 0.00109581; statDown = 0.00111105;}
      else if (x >= 1.9 && x <2.5) {num = 0.994461; den = 0.997793; statUp = 0.00553899; statDown = 0.000611357;}
      else if (x >= 2.5 && x <3) {num = 0.996774; den = 0.999019; statUp = 0.000569912; statDown = 0.000597349;}
      else if (x >= 3 && x <3.5) {num = 0.997965; den = 0.998333; statUp = 0.000546018; statDown = 0.000584973;}
      else if (x >= 3.5 && x <4) {num = 0.997536; den = 0.998867; statUp = 0.000637923; statDown = 0.000692719;}
      else if (x >= 4 && x <4.75) {num = 0.999282; den = 0.999551; statUp = 0.000442867; statDown = 0.000497266;}
      else if (x >= 4.75 && x <5.5) {num = 0.998631; den = 0.999213; statUp = 0.000596208; statDown = 0.00067702;}
      else if (x >= 5.5 && x <8) {num = 0.997838; den = 0.999565; statUp = 0.000516111; statDown = 0.000568641;}
      else if (x >= 8 && x <12) {num = 0.99848; den = 0.998524; statUp = 0.000744772; statDown = 0.000875927;}
      else {num = 1; den = 0.999295; statUp = 2.88878e-11; statDown = 0.000509368;}
    }
  }
  
  double syst_factor = 0; double stat_factor = 0;
  if (idx == 3){return den;}
  else{
    if (idx == -1) syst_factor = syst;
    else if (idx == -2) syst_factor = -1*syst;
    else if (idx == +1) stat_factor = statUp;
    else if (idx == +2) stat_factor = -1*statDown;
    return ((num+syst_factor+stat_factor)/den);
  }
}
  


///////////////////////////////////////////////////                                                                                                                                                
//             GLB TightAcceptance               //                                                                                                                                                
///////////////////////////////////////////////////                                                                                                                                                

double tnp_weight_glb_tightacceptance_pp(double pt, double eta, int idx=0) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (idx != 99) {
    
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <3.7) {num = 0.647134; den = 0.61234; statUp = 0.00402448; statDown = 0.00400457;}
      else if (x >= 3.7 && x <3.9) {num = 0.776324; den = 0.763515; statUp = 0.00447844; statDown = 0.00444373;}
      else if (x >= 3.9 && x <4.2) {num = 0.861661; den = 0.858798; statUp = 0.00385527; statDown = 0.00386838;}
      else if (x >= 4.2 && x <4.6) {num = 0.917883; den = 0.925158; statUp = 0.00351346; statDown = 0.00350028;}
      else if (x >= 4.6 && x <5.2) {num = 0.959048; den = 0.96808; statUp = 0.00305373; statDown = 0.00306687;}
      else if (x >= 5.2 && x <7) {num = 0.98193; den = 0.984917; statUp = 0.00212598; statDown = 0.00212453;}
      else if (x >= 7 && x <10.5) {num = 0.976915; den = 0.985588; statUp = 0.00238157; statDown = 0.00238884;}
      else {num = 0.963402; den = 0.982844; statUp = 0.00364871; statDown = 0.00363991;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.06 && x <2.4) {num = 0.641029; den = 0.639099; statUp = 0.0158524; statDown = 0.0152448;}
      else if (x >= 2.4 && x <2.8) {num = 0.814069; den = 0.813375; statUp = 0.00876633; statDown = 0.00862426;}
      else if (x >= 2.8 && x <3.2) {num = 0.90978; den = 0.907716; statUp = 0.00672604; statDown = 0.00674953;}
      else if (x >= 3.2 && x <3.9) {num = 0.937477; den = 0.942342; statUp = 0.00461469; statDown = 0.00456798;}
      else if (x >= 3.9 && x <5) {num = 0.958227; den = 0.969385; statUp = 0.00400184; statDown = 0.0039995;}
      else if (x >= 5 && x <8) {num = 0.970878; den = 0.974695; statUp = 0.00345924; statDown = 0.00344718;}
      else {num = 0.957445; den = 0.971303; statUp = 0.00506534; statDown = 0.0050402;}
    }
    // SF for 1.8 < |eta| < 2.4
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.4) {
      if (x >= 1.5 && x <1.8) {num = 0.655698; den = 0.688643; statUp = 0.0132701; statDown = 0.0131456;}
      else if (x >= 1.8 && x <2.3) {num = 0.816705; den = 0.836733; statUp = 0.00889631; statDown = 0.00867543;}
      else if (x >= 2.3 && x <2.9) {num = 0.917865; den = 0.939274; statUp = 0.00770755; statDown = 0.00759797;}
      else if (x >= 2.9 && x <3.7) {num = 0.948355; den = 0.969628; statUp = 0.00664304; statDown = 0.00657473;}
      else if (x >= 3.7 && x <6) {num = 0.962889; den = 0.982337; statUp = 0.0046804; statDown = 0.00469451;}
      else {num = 0.948709; den = 0.980496; statUp = 0.00543928; statDown = 0.0053934;}
    }
    
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties                 
      if (x >= 3.5 && x < 3.7) syst = 0.00108251;
      else if (x >= 3.7 && x < 3.9) syst = 0.000211988;
      else if (x >= 3.9 && x < 4.2) syst = 0.000743086;
      else if (x >= 4.2 && x < 4.6) syst = 0.0014736;
      else if (x >= 4.6 && x < 5.2) syst = 9.585e-05;
      else if (x >= 5.2 && x < 7) syst = 0.00149636;
      else if (x >= 7 && x < 10.5) syst = 0.000896526;
      else syst = 0.00096564;
      }
    
    else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties   
      if (x >= 2.06 && x < 2.4) syst = 0.0104505;
      else if (x >= 2.4 && x < 2.8) syst = 0.003583;
      else if (x >= 2.8 && x < 3.2) syst = 0.00125896;
      else if (x >= 3.2 && x < 3.9) syst = 0.00233562;
      else if (x >= 3.9 && x < 5) syst = 0.00213652;
      else if (x >= 5 && x < 8) syst = 0.00213652;
      else syst = 0.000619272;
      
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.4) {
      // syst uncertainties                        
      if (x >= 1.5 && x < 1.8) syst = 0.0103348;
      else if (x >= 1.8 && x < 2.3) syst = 0.00791429;
      else if (x >= 2.3 && x < 2.9) syst = 0.00155826;
      else if (x >= 2.9 && x < 3.7) syst = 0.00197487;
      else if (x >= 3.7 && x < 6) syst = 0.000700885;
      else syst = 0.00538024;
    }
  }
  else {
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <3.7) {num = 0.647134; den = 0.61234; statUp = 0.00402448; statDown = 0.00400457;}
      else if (x >= 3.7 && x <3.9) {num = 0.776324; den = 0.763515; statUp = 0.00447844; statDown = 0.00444373;}
      else if (x >= 3.9 && x <4.2) {num = 0.861661; den = 0.858798; statUp = 0.00385527; statDown = 0.00386838;}
      else if (x >= 4.2 && x <4.6) {num = 0.917883; den = 0.925158; statUp = 0.00351346; statDown = 0.00350028;}
      else if (x >= 4.6 && x <5.2) {num = 0.959048; den = 0.96808; statUp = 0.00305373; statDown = 0.00306687;}
      else if (x >= 5.2 && x <7) {num = 0.98193; den = 0.984917; statUp = 0.00212598; statDown = 0.00212453;}
      else if (x >= 7 && x <10.5) {num = 0.976915; den = 0.985588; statUp = 0.00238157; statDown = 0.00238884;}
      else {num = 0.963402; den = 0.982844; statUp = 0.00364871; statDown = 0.00363991;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.06 && x <2.4) {num = 0.641029; den = 0.639099; statUp = 0.0158524; statDown = 0.0152448;}
      else if (x >= 2.4 && x <2.8) {num = 0.814069; den = 0.813375; statUp = 0.00876633; statDown = 0.00862426;}
      else if (x >= 2.8 && x <3.2) {num = 0.90978; den = 0.907716; statUp = 0.00672604; statDown = 0.00674953;}
      else if (x >= 3.2 && x <3.9) {num = 0.937477; den = 0.942342; statUp = 0.00461469; statDown = 0.00456798;}
      else if (x >= 3.9 && x <5) {num = 0.958227; den = 0.969385; statUp = 0.00400184; statDown = 0.0039995;}
      else if (x >= 5 && x <8) {num = 0.970878; den = 0.974695; statUp = 0.00345924; statDown = 0.00344718;}
      else {num = 0.957445; den = 0.971303; statUp = 0.00506534; statDown = 0.0050402;}
    }
    // SF for 1.8 < |eta| < 2.4
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.4) {
      if (x >= 1.5 && x <1.8) {num = 0.655698; den = 0.688643; statUp = 0.0132701; statDown = 0.0131456;}
      else if (x >= 1.8 && x <2.3) {num = 0.816705; den = 0.836733; statUp = 0.00889631; statDown = 0.00867543;}
      else if (x >= 2.3 && x <2.9) {num = 0.917865; den = 0.939274; statUp = 0.00770755; statDown = 0.00759797;}
      else if (x >= 2.9 && x <3.7) {num = 0.948355; den = 0.969628; statUp = 0.00664304; statDown = 0.00657473;}
      else if (x >= 3.7 && x <6) {num = 0.962889; den = 0.982337; statUp = 0.0046804; statDown = 0.00469451;}
      else {num = 0.948709; den = 0.980496; statUp = 0.00543928; statDown = 0.0053934;}
    }    
  }

  double syst_factor = 0; double stat_factor = 0;
  if (idx ==3) {return den;}
  else{
    if (idx == -1) syst_factor = syst;
    else if (idx == -2) syst_factor = -1*syst;
    else if (idx == +1) stat_factor = statUp;
    else if (idx == +2) stat_factor = -1*statDown;
    return ((num+syst_factor+stat_factor)/den);
  }
}
///////////////////////////////////////////////////                                                                                                                                                
//         MuIdTrg TightAcceptance               //                                                                                                                                              
///////////////////////////////////////////////////                                                                                                                                                

double tnp_weight_muidtrg_tightacceptance_pp(double pt, double eta, int idx=0) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (idx != 99) {
    
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <3.75) {num = 0.760025; den = 0.726008; statUp = 0.00223153; statDown = 0.00224102;}
      else if (x >= 3.75 && x <4) {num = 0.840932; den = 0.810266; statUp = 0.00187357; statDown = 0.00188621;}
      else if (x >= 4 && x <4.25) {num = 0.881662; den = 0.863344; statUp = 0.00169811; statDown = 0.00171232;}
      else if (x >= 4.25 && x <4.5) {num = 0.903897; den = 0.892805; statUp = 0.00163168; statDown = 0.00164881;}
      else if (x >= 4.5 && x <4.75) {num = 0.917579; den = 0.905505; statUp = 0.00160761; statDown = 0.00162638;}
      else if (x >= 4.75 && x <5) {num = 0.925012; den = 0.913046; statUp = 0.00163619; statDown = 0.00165509;}
      else if (x >= 5 && x <5.25) {num = 0.933373; den = 0.923486; statUp = 0.00164321; statDown = 0.00166401;}
      else if (x >= 5.25 && x <5.5) {num = 0.941896; den = 0.928054; statUp = 0.00164337; statDown = 0.00166749;}
      else if (x >= 5.5 && x <5.75) {num = 0.940442; den = 0.929635; statUp = 0.00175782; statDown = 0.00178401;}
      else if (x >= 5.75 && x <6) {num = 0.940325; den = 0.932643; statUp = 0.00186084; statDown = 0.001894;}
      else if (x >= 6 && x <6.5) {num = 0.942744; den = 0.943386; statUp = 0.00141057; statDown = 0.00142927;}
      else if (x >= 6.5 && x <7) {num = 0.942266; den = 0.95177; statUp = 0.00159219; statDown = 0.00161629;}
      else if (x >= 7 && x <7.5) {num = 0.943744; den = 0.953602; statUp = 0.00174535; statDown = 0.00177521;}
      else if (x >= 7.5 && x <8) {num = 0.946452; den = 0.954787; statUp = 0.00188743; statDown = 0.00192476;}
      else if (x >= 8 && x <9.25) {num = 0.949518; den = 0.955381; statUp = 0.00137582; statDown = 0.00139715;}
      else if (x >= 9.25 && x <10.5) {num = 0.948893; den = 0.959961; statUp = 0.00172295; statDown = 0.00175923;}
      else if (x >= 10.5 && x <12.25) {num = 0.947729; den = 0.958987; statUp = 0.00190692; statDown = 0.00194425;}
      else if (x >= 12.25 && x <14) {num = 0.943215; den = 0.95612; statUp = 0.00254603; statDown = 0.00260929;}
      else if (x >= 14 && x <16) {num = 0.927823; den = 0.958966; statUp = 0.00341274; statDown = 0.00350551;}
      else if (x >= 16 && x <18) {num = 0.924853; den = 0.949437; statUp = 0.00447626; statDown = 0.00460933;}
      else if (x >= 18 && x <24) {num = 0.920994; den = 0.947671; statUp = 0.0040602; statDown = 0.00417625;}
      else {num = 0.907235; den = 0.947741; statUp = 0.00726321; statDown = 0.00757945;}
    }
    
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.06 && x <2.5) {num = 0.663024; den = 0.629104; statUp = 0.00465948; statDown = 0.0046752;}
      else if (x >= 2.5 && x <2.75) {num = 0.716054; den = 0.67925; statUp = 0.00346314; statDown = 0.00347783;}
      else if (x >= 2.75 && x <3) {num = 0.77954; den = 0.755625; statUp = 0.00273296; statDown = 0.00274738;}
      else if (x >= 3 && x <3.25) {num = 0.778439; den = 0.777372; statUp = 0.00251347; statDown = 0.0025254;}
      else if (x >= 3.25 && x <3.5) {num = 0.778047; den = 0.78947; statUp = 0.00256436; statDown = 0.00257709;}
      else if (x >= 3.5 && x <3.75) {num = 0.790273; den = 0.807203; statUp = 0.00268526; statDown = 0.00270314;}
      else if (x >= 3.75 && x <4) {num = 0.79865; den = 0.81617; statUp = 0.00284142; statDown = 0.00285024;}
      else if (x >= 4 && x <4.5) {num = 0.811593; den = 0.826727; statUp = 0.00216103; statDown = 0.00217519;}
      else if (x >= 4.5 && x <5) {num = 0.833137; den = 0.842921; statUp = 0.00237555; statDown = 0.00239234;}
      else if (x >= 5 && x <5.5) {num = 0.856332; den = 0.866264; statUp = 0.00255814; statDown = 0.00257031;}
      else if (x >= 5.5 && x <6) {num = 0.879643; den = 0.883357; statUp = 0.00268117; statDown = 0.00271087;}
      else if (x >= 6 && x <7) {num = 0.894543; den = 0.906841; statUp = 0.00216179; statDown = 0.00217536;}
      else if (x >= 7 && x <8.5) {num = 0.911399; den = 0.91876; statUp = 0.00215268; statDown = 0.00217989;}
      else if (x >= 8.5 && x <10) {num = 0.922513; den = 0.928345; statUp = 0.00272674; statDown = 0.00278094;}
      else if (x >= 10 && x <12) {num = 0.936125; den = 0.945127; statUp = 0.00295207; statDown = 0.0030282;}
      else if (x >= 12 && x <14.5) {num = 0.941952; den = 0.942103; statUp = 0.00362431; statDown = 0.00374608;}
      else if (x >= 14.5 && x <19) {num = 0.936863; den = 0.937644; statUp = 0.0631369; statDown = 0.00469025;}
      else {num = 0.934587; den = 0.936006; statUp = 0.00658549; statDown = 0.00693284;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.5 && x <1.9) {num = 0.69191; den = 0.595152; statUp = 0.0046497; statDown = 0.00457911;}
      else if (x >= 1.9 && x <2.2) {num = 0.763801; den = 0.705953; statUp = 0.00306162; statDown = 0.0030754;}
      else if (x >= 2.2 && x <2.5) {num = 0.857503; den = 0.818846; statUp = 0.0024947; statDown = 0.00251812;}
      else if (x >= 2.5 && x <2.8) {num = 0.906732; den = 0.896244; statUp = 0.00222889; statDown = 0.00225497;}
      else if (x >= 2.8 && x <3.1) {num = 0.927348; den = 0.92319; statUp = 0.00221453; statDown = 0.00223841;}
      else if (x >= 3.1 && x <3.4) {num = 0.929258; den = 0.930977; statUp = 0.0023626; statDown = 0.00240264;}
      else if (x >= 3.4 && x <4) {num = 0.927743; den = 0.935697; statUp = 0.00192051; statDown = 0.00193394;}
      else if (x >= 4 && x <5) {num = 0.921443; den = 0.938186; statUp = 0.00189124; statDown = 0.00192131;}
      else if (x >= 5 && x <6) {num = 0.917552; den = 0.936965; statUp = 0.00250432; statDown = 0.00254367;}
      else if (x >= 6 && x <8) {num = 0.909679; den = 0.929456; statUp = 0.00262321; statDown = 0.00266266;}
      else if (x >= 8 && x <11) {num = 0.897524; den = 0.913692; statUp = 0.00375973; statDown = 0.00383064;}
      else if (x >= 11 && x <15) {num = 0.868831; den = 0.885758; statUp = 0.00645408; statDown = 0.0066088;}
      else {num = 0.8119; den = 0.823906; statUp = 0.0107466; statDown = 0.0110123;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.5 && x <1.9) {num = 0.836604; den = 0.754014; statUp = 0.00245146; statDown = 0.0024644;}
      else if (x >= 1.9 && x <2.3) {num = 0.926632; den = 0.884537; statUp = 0.00181535; statDown = 0.00183481;}
      else if (x >= 2.3 && x <2.7) {num = 0.940324; den = 0.916059; statUp = 0.00180535; statDown = 0.00183144;}
      else if (x >= 2.7 && x <3) {num = 0.940721; den = 0.930185; statUp = 0.00229731; statDown = 0.00234184;}
      else if (x >= 3 && x <3.5) {num = 0.94583; den = 0.932965; statUp = 0.00192795; statDown = 0.00196199;}
      else if (x >= 3.5 && x <4.25) {num = 0.945673; den = 0.942458; statUp = 0.00188936; statDown = 0.00192219;}
      else if (x >= 4.25 && x <5) {num = 0.945181; den = 0.947329; statUp = 0.00238284; statDown = 0.00243443;}
      else if (x >= 5 && x <6) {num = 0.94597; den = 0.953674; statUp = 0.00259303; statDown = 0.0026595;}
      else if (x >= 6 && x <7) {num = 0.946883; den = 0.949923; statUp = 0.00344602; statDown = 0.00344602;}
      else if (x >= 7 && x <8.5) {num = 0.955594; den = 0.955484; statUp = 0.00337735; statDown = 0.00353239;}
      else if (x >= 8.5 && x <13) {num = 0.943172; den = 0.936958; statUp = 0.00378355; statDown = 0.00391933;}
      else {num = 0.879442; den = 0.885817; statUp = 0.00978302; statDown = 0.00971173;}
    }
    

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 3.75) syst = 0.000900168;
      else if (x >= 3.75 && x < 4) syst =  0.000865001;
      else if (x >= 4 && x < 4.25) syst = 0.000129247;
      else if (x >= 4.25 && x < 4.5) syst = 0.000452042;
      else if (x >= 4.5 && x < 4.75) syst = 0.000501881;
      else if (x >= 4.75 && x < 5) syst = 0.00038986;
      else if (x >= 5 && x < 5.25) syst = 0.000558043;
      else if (x >= 5.25 && x < 5.5) syst = 0.000621342;
      else if (x >= 5.5 && x < 5.75) syst = 0.000708795;
      else if (x >= 5.75 && x < 6) syst = 0.00144451;
      else if (x >= 6 && x < 6.5) syst = 0.000387211;
      else if (x >= 6.5 && x < 7) syst = 0.000334128;
      else if (x >= 7 && x < 7.5) syst = 0.000377601;
      else if (x >= 7.5 && x < 8) syst = 0.000397838;
      else if (x >= 8 && x < 9.25) syst = 0.000219424;
      else if (x >= 9.25 && x < 10.5) syst = 0.00046839;
      else if (x >= 10.5 && x < 12.25) syst = 0.000342797;
      else if (x >= 12.25 && x < 14) syst = 0.000312191;
      else if (x >= 14 && x < 16) syst = 0.000816031;
      else if (x >= 16 && x < 18) syst = 0.00173171;
      else if (x >= 18 && x < 24) syst = 0.000750224;
      else syst = 0.00260996;
      
    }
    
    else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.06 && x < 2.5) syst = 0.0086659;
      else if (x >= 2.5 && x < 2.75) syst = 0.000493877;
      else if (x >= 2.75 && x < 3) syst = 0.00051552;
      else if (x >= 3 && x < 3.25) syst = 0.000650532;
      else if (x >= 3.25 && x < 3.5) syst = 0.00101988;
      else if (x >= 3.5 && x < 3.75) syst = 0.00080565;
      else if (x >= 3.75 && x < 4) syst = 0.00128796;
      else if (x >= 4 && x < 4.5) syst = 0.00138638;
      else if (x >= 4.5 && x < 5) syst = 0.000574534;
      else if (x >= 5 && x < 5.5) syst = 0.0015193;
      else if (x >= 5.5 && x < 6) syst = 0.000849741;
      else if (x >= 6 && x < 7) syst = 0.00132843;
      else if (x >= 7 && x < 8.5) syst = 0.00141202;
      else if (x >= 8.5 && x < 10) syst = 0.00125018;
      else if (x >= 10 && x < 12) syst =  0.00082619;
      else if (x >= 12 && x < 14.5) syst = 0.00226697;
      else if (x >= 14.5 && x < 19) syst = 0.00146068;
      else syst = 0.000958602;
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties                                                                                                                                                                              
      if (x >= 1.5 && x < 1.9) syst = 0.00518218;                                                                                                                                             
      else if (x >= 1.9 && x < 2.2) syst = 0.00165702;                                                                                                                                        
      else if (x >= 2.2 && x < 2.5) syst = 0.00198365;                                                                                                                                          
      else if (x >= 2.5 && x < 2.8) syst = 0.00103727;                                                                                                                                            
      else if (x >= 2.8 && x < 3.1) syst = 0.000301947;                                                                                                                                      
      else if (x >= 3.1 && x < 3.4) syst = 0.00107157;                                                                                                                                         
      else if (x >= 3.4 && x < 4) syst = 0.000162552;                                                                                                                                             
      else if (x >= 4 && x < 5) syst = 0.00109106;                                                                                                                                               
      else if (x >= 5 && x < 6) syst =  0.00139839;                                                                                                                                               
      else if (x >= 6 && x < 8) syst = 0.00105435;                                                                                                                                            
      else if (x >= 8 && x < 11) syst = 0.000812637;                                                                                                                                             
      else if (x >= 11 && x < 15) syst = 0.000877749;                                                                                                                                             
      else syst = 0.012101;  
    }
    else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.5 && x < 1.9) syst = 0.000888127;
      else if (x >= 1.9 && x < 2.3) syst = 0.00049553;
      else if (x >= 2.3 && x < 2.7) syst = 0.00156344;
      else if (x >= 2.7 && x < 3) syst = 0.000958998;
      else if (x >= 3 && x < 3.5) syst = 0.00104966;
      else if (x >= 3.5 && x < 4.25) syst = 0.000761361;
      else if (x >= 4.25 && x < 5) syst = 0.00299566;
      else if (x >= 5 && x < 6) syst = 0.000262191;
      else if (x >= 6 && x < 7) syst =  0.00239638;
      else if (x >= 7 && x < 8.5) syst = 0.00168078;
      else if (x >= 8.5 && x < 13) syst = 0.00220707;
      else syst = 0.0109342;
    }
  }
  else {
    
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <3.75) {num = 0.760025; den = 0.726008; statUp = 0.00223153; statDown = 0.00224102;}
      else if (x >= 3.75 && x <4) {num = 0.840932; den = 0.810266; statUp = 0.00187357; statDown = 0.00188621;}
      else if (x >= 4 && x <4.25) {num = 0.881662; den = 0.863344; statUp = 0.00169811; statDown = 0.00171232;}
      else if (x >= 4.25 && x <4.5) {num = 0.903897; den = 0.892805; statUp = 0.00163168; statDown = 0.00164881;}
      else if (x >= 4.5 && x <4.75) {num = 0.917579; den = 0.905505; statUp = 0.00160761; statDown = 0.00162638;}
      else if (x >= 4.75 && x <5) {num = 0.925012; den = 0.913046; statUp = 0.00163619; statDown = 0.00165509;}
      else if (x >= 5 && x <5.25) {num = 0.933373; den = 0.923486; statUp = 0.00164321; statDown = 0.00166401;}
      else if (x >= 5.25 && x <5.5) {num = 0.941896; den = 0.928054; statUp = 0.00164337; statDown = 0.00166749;}
      else if (x >= 5.5 && x <5.75) {num = 0.940442; den = 0.929635; statUp = 0.00175782; statDown = 0.00178401;}
      else if (x >= 5.75 && x <6) {num = 0.940325; den = 0.932643; statUp = 0.00186084; statDown = 0.001894;}
      else if (x >= 6 && x <6.5) {num = 0.942744; den = 0.943386; statUp = 0.00141057; statDown = 0.00142927;}
      else if (x >= 6.5 && x <7) {num = 0.942266; den = 0.95177; statUp = 0.00159219; statDown = 0.00161629;}
      else if (x >= 7 && x <7.5) {num = 0.943744; den = 0.953602; statUp = 0.00174535; statDown = 0.00177521;}
      else if (x >= 7.5 && x <8) {num = 0.946452; den = 0.954787; statUp = 0.00188743; statDown = 0.00192476;}
      else if (x >= 8 && x <9.25) {num = 0.949518; den = 0.955381; statUp = 0.00137582; statDown = 0.00139715;}
      else if (x >= 9.25 && x <10.5) {num = 0.948893; den = 0.959961; statUp = 0.00172295; statDown = 0.00175923;}
      else if (x >= 10.5 && x <12.25) {num = 0.947729; den = 0.958987; statUp = 0.00190692; statDown = 0.00194425;}
      else if (x >= 12.25 && x <14) {num = 0.943215; den = 0.95612; statUp = 0.00254603; statDown = 0.00260929;}
      else if (x >= 14 && x <16) {num = 0.927823; den = 0.958966; statUp = 0.00341274; statDown = 0.00350551;}
      else if (x >= 16 && x <18) {num = 0.924853; den = 0.949437; statUp = 0.00447626; statDown = 0.00460933;}
      else if (x >= 18 && x <24) {num = 0.920994; den = 0.947671; statUp = 0.0040602; statDown = 0.00417625;}
      else {num = 0.907235; den = 0.947741; statUp = 0.00726321; statDown = 0.00757945;}
    }

    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.06 && x <2.5) {num = 0.663024; den = 0.629104; statUp = 0.00465948; statDown = 0.0046752;}
      else if (x >= 2.5 && x <2.75) {num = 0.716054; den = 0.67925; statUp = 0.00346314; statDown = 0.00347783;}
      else if (x >= 2.75 && x <3) {num = 0.77954; den = 0.755625; statUp = 0.00273296; statDown = 0.00274738;}
      else if (x >= 3 && x <3.25) {num = 0.778439; den = 0.777372; statUp = 0.00251347; statDown = 0.0025254;}
      else if (x >= 3.25 && x <3.5) {num = 0.778047; den = 0.78947; statUp = 0.00256436; statDown = 0.00257709;}
      else if (x >= 3.5 && x <3.75) {num = 0.790273; den = 0.807203; statUp = 0.00268526; statDown = 0.00270314;}
      else if (x >= 3.75 && x <4) {num = 0.79865; den = 0.81617; statUp = 0.00284142; statDown = 0.00285024;}
      else if (x >= 4 && x <4.5) {num = 0.811593; den = 0.826727; statUp = 0.00216103; statDown = 0.00217519;}
      else if (x >= 4.5 && x <5) {num = 0.833137; den = 0.842921; statUp = 0.00237555; statDown = 0.00239234;}
      else if (x >= 5 && x <5.5) {num = 0.856332; den = 0.866264; statUp = 0.00255814; statDown = 0.00257031;}
      else if (x >= 5.5 && x <6) {num = 0.879643; den = 0.883357; statUp = 0.00268117; statDown = 0.00271087;}
      else if (x >= 6 && x <7) {num = 0.894543; den = 0.906841; statUp = 0.00216179; statDown = 0.00217536;}
      else if (x >= 7 && x <8.5) {num = 0.911399; den = 0.91876; statUp = 0.00215268; statDown = 0.00217989;}
      else if (x >= 8.5 && x <10) {num = 0.922513; den = 0.928345; statUp = 0.00272674; statDown = 0.00278094;}
      else if (x >= 10 && x <12) {num = 0.936125; den = 0.945127; statUp = 0.00295207; statDown = 0.0030282;}
      else if (x >= 12 && x <14.5) {num = 0.941952; den = 0.942103; statUp = 0.00362431; statDown = 0.00374608;}
      else if (x >= 14.5 && x <19) {num = 0.936863; den = 0.937644; statUp = 0.0631369; statDown = 0.00469025;}
      else {num = 0.934587; den = 0.936006; statUp = 0.00658549; statDown = 0.00693284;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.5 && x <1.9) {num = 0.69191; den = 0.595152; statUp = 0.0046497; statDown = 0.00457911;}
      else if (x >= 1.9 && x <2.2) {num = 0.763801; den = 0.705953; statUp = 0.00306162; statDown = 0.0030754;}
      else if (x >= 2.2 && x <2.5) {num = 0.857503; den = 0.818846; statUp = 0.0024947; statDown = 0.00251812;}
      else if (x >= 2.5 && x <2.8) {num = 0.906732; den = 0.896244; statUp = 0.00222889; statDown = 0.00225497;}
      else if (x >= 2.8 && x <3.1) {num = 0.927348; den = 0.92319; statUp = 0.00221453; statDown = 0.00223841;}
      else if (x >= 3.1 && x <3.4) {num = 0.929258; den = 0.930977; statUp = 0.0023626; statDown = 0.00240264;}
      else if (x >= 3.4 && x <4) {num = 0.927743; den = 0.935697; statUp = 0.00192051; statDown = 0.00193394;}
      else if (x >= 4 && x <5) {num = 0.921443; den = 0.938186; statUp = 0.00189124; statDown = 0.00192131;}
      else if (x >= 5 && x <6) {num = 0.917552; den = 0.936965; statUp = 0.00250432; statDown = 0.00254367;}
      else if (x >= 6 && x <8) {num = 0.909679; den = 0.929456; statUp = 0.00262321; statDown = 0.00266266;}
      else if (x >= 8 && x <11) {num = 0.897524; den = 0.913692; statUp = 0.00375973; statDown = 0.00383064;}
      else if (x >= 11 && x <15) {num = 0.868831; den = 0.885758; statUp = 0.00645408; statDown = 0.0066088;}
      else {num = 0.8119; den = 0.823906; statUp = 0.0107466; statDown = 0.0110123;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.5 && x <1.9) {num = 0.836604; den = 0.754014; statUp = 0.00245146; statDown = 0.0024644;}
      else if (x >= 1.9 && x <2.3) {num = 0.926632; den = 0.884537; statUp = 0.00181535; statDown = 0.00183481;}
      else if (x >= 2.3 && x <2.7) {num = 0.940324; den = 0.916059; statUp = 0.00180535; statDown = 0.00183144;}
      else if (x >= 2.7 && x <3) {num = 0.940721; den = 0.930185; statUp = 0.00229731; statDown = 0.00234184;}
      else if (x >= 3 && x <3.5) {num = 0.94583; den = 0.932965; statUp = 0.00192795; statDown = 0.00196199;}
      else if (x >= 3.5 && x <4.25) {num = 0.945673; den = 0.942458; statUp = 0.00188936; statDown = 0.00192219;}
      else if (x >= 4.25 && x <5) {num = 0.945181; den = 0.947329; statUp = 0.00238284; statDown = 0.00243443;}
      else if (x >= 5 && x <6) {num = 0.94597; den = 0.953674; statUp = 0.00259303; statDown = 0.0026595;}
      else if (x >= 6 && x <7) {num = 0.946883; den = 0.949923; statUp = 0.00344602; statDown = 0.00344602;}
      else if (x >= 7 && x <8.5) {num = 0.955594; den = 0.955484; statUp = 0.00337735; statDown = 0.00353239;}
      else if (x >= 8.5 && x <13) {num = 0.943172; den = 0.936958; statUp = 0.00378355; statDown = 0.00391933;}
      else {num = 0.879442; den = 0.885817; statUp = 0.00978302; statDown = 0.00971173;}
    }
  }
  double syst_factor = 0; double stat_factor = 0;
  if (idx ==3) {return den;}
  else{
    if (idx == -1) syst_factor = syst;
    else if (idx == -2) syst_factor = -1*syst;
    else if (idx == +1) stat_factor = statUp;
    else if (idx == +2) stat_factor = -1*statDown;
    return ((num+syst_factor+stat_factor)/den);
  }
}
#endif
