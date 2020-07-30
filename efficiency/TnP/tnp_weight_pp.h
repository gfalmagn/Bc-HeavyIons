#ifndef tnp_weight_h
#define tnp_weight_h

#include "TMath.h"
#include <tuple>

// IN THIS FILE YOU WILL FIND:
// ++++++++++++++

//BINNED SCALE FACTORS for pp 2017

//To extract the scale factor for a muon of PT=pt and ETA=eta, do:

//#include "tnp_weight.h"
//void testtnpweight(){
// auto SFwithError = tnp_weight_theEfficiencyYouWant_pp( pt, eta );
// double central_value = std::get<0>(SFwithError);
// double stat_error = std::get<1>(SFwithError);
// double syst_error = std::get<2>(SFwithError);
// double tot_error = std::get<3>(SFwithError);
//}

// ++++++++++++++

///////////////////////////////////////////////////
//              G l o b a l    p p               //
///////////////////////////////////////////////////

//To be used for global muon efficiency, with tight muon kinematic acceptance cuts                                                                                                                   
//Central scale-factor value, statistical error, systematic error, total error                                                                                                                        

std::tuple<double,double,double,double> tnp_weight_GlobalMuon_TightAcceptance_pp(double pt, double eta){

  if( (fabs(eta) > 0) && (fabs(eta) < 1.2) ){
    if (pt<3.7) return std::make_tuple(1.05682, 0.0065723, 0.00108251, 0.00666085);
    else if (pt<3.9) return std::make_tuple(1.01678, 0.00586556, 0.000211988, 0.00586939);
    else if (pt<4.2) return std::make_tuple(1.00333, 0.00448914, 0.000743086, 0.00455023);
    else if (pt<4.6) return std::make_tuple(0.992136, 0.00379769, 0.0014736, 0.00407357);
    else if (pt<5.2) return std::make_tuple(0.99067, 0.00315441, 9.585e-05, 0.00315587);
    else if (pt<7) return std::make_tuple(0.996967, 0.00215854, 0.00149636, 0.00262647);
    else if (pt<10.5) return std::make_tuple(0.991201, 0.0024164, 0.000896526, 0.00257735);
    else if (pt<30) return std::make_tuple(0.980219, 0.0037124, 0.00096564, 0.00383593);
  }

  if( (fabs(eta) > 1.2) && (fabs(eta) < 1.8) ){
    if (pt<2.4) return std::make_tuple(1.00302, 0.0248043, 0.0104505, 0.0269159);
    else if (pt<2.8) return std::make_tuple(1.00085, 0.0107777, 0.003583, 0.0113577);
    else if (pt<3.2) return std::make_tuple(1.00227, 0.00740985, 0.00125896, 0.00751604);
    else if (pt<3.9) return std::make_tuple(0.994838, 0.00489704, 0.00233562, 0.00542551);
    else if (pt<5) return std::make_tuple(0.98849, 0.00412823, 0.00213652, 0.00464833);
    else if (pt<8) return std::make_tuple(0.996084, 0.00354905, 0.00058147, 0.00359637);
    else if (pt<30) return std::make_tuple(0.985733, 0.00521499, 0.000619272, 0.00525163);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.4) ){
    if (pt<1.8) return std::make_tuple(0.95216, 0.0192699, 0.0103348, 0.0218664);
    else if (pt<2.3) return std::make_tuple(0.976063, 0.0106322, 0.00791429, 0.0132544);
    else if (pt<2.9) return std::make_tuple(0.977207, 0.00820585, 0.00155826, 0.0083525);
    else if (pt<3.7) return std::make_tuple(0.978061, 0.00685112, 0.00197487, 0.00713008);
    else if (pt<6) return std::make_tuple(0.980203, 0.00476456, 0.000700885, 0.00481583);
    else if (pt<30) return std::make_tuple(0.967581, 0.00554748, 0.00538024, 0.00772797);
  }

  return std::make_tuple(1,0,0,0);
}

//To be used for global muon efficiency, with loose muon kinematic acceptance cuts
//Central scale-factor value, statistical error, systematic error, total error

std::tuple<double,double,double,double> tnp_weight_GlobalMuon_LooseAcceptance_pp(double pt, double eta){

  if( (fabs(eta) > 0) && (fabs(eta) < 1.1) ){
    if (pt<3.7) return std::make_tuple(1.13756, 0.00582119, 0.00380783, 0.00695599);
    else if (pt<3.9) return std::make_tuple(1.01976, 0.00608097, 0.00336339, 0.00694914);
    else if (pt<4.2) return std::make_tuple(1.00639, 0.00467267, 0.00121721, 0.00482861);
    else if (pt<4.6) return std::make_tuple(0.992333, 0.00392665, 0.00114831, 0.00409112);
    else if (pt<5.2) return std::make_tuple(0.990089, 0.00324541, 0.000429282, 0.00327367);
    else if (pt<7) return std::make_tuple(0.996127, 0.00221153, 0.00147196, 0.0026566);
    else if (pt<10.5) return std::make_tuple(0.992142, 0.00247593, 0.000986955, 0.00266539);
    else if (pt<30) return std::make_tuple(0.980844, 0.00377456, 0.00214214, 0.00434006);
  }

  if( (fabs(eta) > 1.1) && (fabs(eta) < 1.8) ){
    if (pt<2.4) return std::make_tuple(0.972026, 0.0103711, 0.00567522, 0.0118224);
    else if (pt<2.8) return std::make_tuple(0.957819, 0.00749527, 0.00308897, 0.00810684);
    else if (pt<3.2) return std::make_tuple(1.00173, 0.00669556, 0.00168225, 0.00690366);
    else if (pt<3.9) return std::make_tuple(0.999546, 0.00457292, 0.00294545, 0.00543942);
    else if (pt<5) return std::make_tuple(0.991865, 0.00377255, 0.00139078, 0.00402075);
    else if (pt<8) return std::make_tuple(0.99806, 0.00322876, 0.00187509, 0.00373375);
    else if (pt<30) return std::make_tuple(0.983852, 0.00469324, 0.00102468, 0.0048038);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.4) ){
    if (pt<1.8) return std::make_tuple(0.953515, 0.013659, 0.014648, 0.0200282);
    else if (pt<2.3) return std::make_tuple(0.973334, 0.0100595, 0.00573586, 0.0115799);
    else if (pt<2.9) return std::make_tuple(0.977207, 0.00815023, 0.00220486, 0.0084432);
    else if (pt<3.7) return std::make_tuple(0.978061, 0.00682009, 0.00248013, 0.00725704);
    else if (pt<6) return std::make_tuple(0.980203, 0.00476855, 0.0014999, 0.00499888);
    else if (pt<30) return std::make_tuple(0.967581, 0.0055242, 0.00540491, 0.0077285);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.4) ){
    if (pt<1.8) return std::make_tuple(0.953515, 0.013659, 0.014648, 0.0200282);
    else if (pt<2.3) return std::make_tuple(0.973334, 0.0100595, 0.00573586, 0.0115799);
    else if (pt<2.9) return std::make_tuple(0.977207, 0.00815023, 0.00220486, 0.0084432);
    else if (pt<3.7) return std::make_tuple(0.978061, 0.00682009, 0.00248013, 0.00725704);
    else if (pt<6) return std::make_tuple(0.980203, 0.00476855, 0.0014999, 0.00499888);
    else if (pt<30) return std::make_tuple(0.967581, 0.0055242, 0.00540491, 0.0077285);
  }
  return std::make_tuple(1,0,0,0);
}

//To be used for HybridSoftID efficiency, with loose muon kinematic acceptance cuts
//Central scale-factor value, statistical error, systematic error, total error 
std::tuple<double,double,double,double> tnp_weight_HybridSoftID_LooseAcceptance_pp(double pt, double eta){ 

  if( (fabs(eta) > 0) && (fabs(eta) < 1.1) ){
    if (pt<3.7) return std::make_tuple(0.991039, 0.00141451, 0.000509766, 0.00150356);
    else if (pt<4) return std::make_tuple(0.99411, 0.0010317, 0.000204188, 0.00105171);
    else if (pt<4.5) return std::make_tuple(0.995946, 0.000679001, 6.97228e-05, 0.000682571);
    else if (pt<5) return std::make_tuple(0.997494, 0.000622336, 0.000383673, 0.000731099);
    else if (pt<5.5) return std::make_tuple(0.998257, 0.000599561, 4.79224e-05, 0.000601474);
    else if (pt<6) return std::make_tuple(0.997557, 0.000613637, 5.19586e-05, 0.000615833);
    else if (pt<7) return std::make_tuple(0.998089, 0.000504536, 9.82613e-05, 0.000514016);
    else if (pt<8) return std::make_tuple(0.99817, 0.000600291, 0.000129404, 0.00061408);
    else if (pt<10.5) return std::make_tuple(1.00011, 0.000481444, 0.000263544, 0.000548857);
    else if (pt<14) return std::make_tuple(1.00038, 0.00061336, 0.000527764, 0.000809164);
    else if (pt<18) return std::make_tuple(0.999253, 0.00111745, 0.000554964, 0.00124767);
    else if (pt<30) return std::make_tuple(1.00003, 0.00133897, 0.000341248, 0.00138177);
  }

  if( (fabs(eta) > 1.1) && (fabs(eta) < 1.8) ){
    if (pt<2.5) return std::make_tuple(0.998281, 0.00100785, 0.000295922, 0.00105039);
    else if (pt<3) return std::make_tuple(0.996762, 0.000659372, 0.000166904, 0.000680168);
    else if (pt<3.5) return std::make_tuple(0.996493, 0.000612724, 0.000136018, 0.000627639);
    else if (pt<4) return std::make_tuple(0.996128, 0.000604109, 0.000232836, 0.000647426);
    else if (pt<4.5) return std::make_tuple(0.997668, 0.000499361, 0.000208141, 0.000541003);
    else if (pt<5) return std::make_tuple(0.997261, 0.000642767, 4.43061e-05, 0.000644293);
    else if (pt<6) return std::make_tuple(0.997922, 0.000485236, 6.43679e-05, 0.000489487);
    else if (pt<7) return std::make_tuple(0.996794, 0.000599223, 0.000133639, 0.000613944);
    else if (pt<9) return std::make_tuple(0.997508, 0.000552848, 0.000234433, 0.0006005);
    else if (pt<14) return std::make_tuple(0.997054, 0.000615961, 9.69372e-05, 0.000623542);
    else if (pt<18) return std::make_tuple(0.999195, 0.00136452, 7.8333e-05, 0.00136677);
    else if (pt<30) return std::make_tuple(0.999577, 0.00157859, 0.000303196, 0.00160744);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.1) ){
    if (pt<1.75) return std::make_tuple(0.994851, 0.00186009, 0.000212624, 0.0018722);
    else if (pt<2.1) return std::make_tuple(0.995416, 0.000915212, 0.000641968, 0.00111792);
    else if (pt<2.5) return std::make_tuple(0.998339, 0.000544649, 0.000506432, 0.000743718);
    else if (pt<3) return std::make_tuple(0.998359, 0.000457187, 0.000227392, 0.000510615);
    else if (pt<3.5) return std::make_tuple(0.998507, 0.000469804, 6.26562e-05, 0.000473963);
    else if (pt<4) return std::make_tuple(0.997919, 0.00052654, 0.000175154, 0.000554908);
    else if (pt<4.5) return std::make_tuple(0.99804, 0.000554418, 0.000332236, 0.000646344);
    else if (pt<5.25) return std::make_tuple(0.998739, 0.000508105, 0.000140418, 0.00052715);
    else if (pt<6.5) return std::make_tuple(0.999288, 0.000436503, 8.59923e-05, 0.000444893);
    else if (pt<9) return std::make_tuple(0.998288, 0.000519603, 0.000285149, 0.000592703);
    else if (pt<13) return std::make_tuple(0.998302, 0.000743805, 0.000569406, 0.000936733);
    else if (pt<30) return std::make_tuple(0.998483, 0.00124804, 0.000142791, 0.00125618);
  }
 
  if( (fabs(eta) > 2.1) && (fabs(eta) < 2.4) ){
    if (pt<1.9) return std::make_tuple(0.99338, 0.00111883, 0.000402815, 0.00118914);
    else if (pt<2.5) return std::make_tuple(0.99666, 0.000604072, 0.000339534, 0.000692955);
    else if (pt<3) return std::make_tuple(0.997752, 0.000584138, 0.000383278, 0.000698655);
    else if (pt<3.5) return std::make_tuple(0.999632, 0.000566087, 0.000461662, 0.00073047);
    else if (pt<4) return std::make_tuple(0.998667, 0.000665826, 0.000433435, 0.000794474);
    else if (pt<4.75) return std::make_tuple(0.999731, 0.000469729, 0.000167938, 0.000498847);
    else if (pt<5.5) return std::make_tuple(0.999417, 0.000636667, 0.000601941, 0.000876172);
    else if (pt<8) return std::make_tuple(0.998273, 0.000542231, 0.000172143, 0.000568901);
    else if (pt<12) return std::make_tuple(0.999956, 0.000809011, 0.000373559, 0.000891093);
    else if (pt<30) return std::make_tuple(1.00071, 0.000542843, 0.000541856, 0.000766998);
  }
       
  return std::make_tuple(1,0,0,0);                                                                                                                                                                   }                

//To be used for HybridSoftID+(J/psi trigger) efficiency, with tight muon kinematic acceptance cuts (the loose cuts cannot be used with trigger efficiency)
//Central scale-factor value, statistical error, systematic error, total error

std::tuple<double,double,double,double> tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(double pt, double eta){
  if( (fabs(eta) > 0) && (fabs(eta) < 1.2) ){
    if (pt<3.75) return std::make_tuple(1.04685, 0.00307983, 0.000900168, 0.00320868);
    else if (pt<4) return std::make_tuple(1.03785, 0.00231993, 0.000865001, 0.00247594);
    else if (pt<4.25) return std::make_tuple(1.02122, 0.0019751, 0.000129247, 0.00197932);
    else if (pt<4.5) return std::make_tuple(1.01242, 0.00183717, 0.000452042, 0.00189196);
    else if (pt<4.75) return std::make_tuple(1.01333, 0.00178571, 0.000501881, 0.0018549);
    else if (pt<5) return std::make_tuple(1.01311, 0.00180214, 0.00038986, 0.00184382);
    else if (pt<5.25) return std::make_tuple(1.01071, 0.00179057, 0.000558043, 0.00187552);
    else if (pt<5.5) return std::make_tuple(1.01492, 0.00178375, 0.000621342, 0.00188887);
    else if (pt<5.75) return std::make_tuple(1.01162, 0.00190487, 0.000708795, 0.00203246);
    else if (pt<6) return std::make_tuple(1.00824, 0.00201297, 0.00144451, 0.00247763);
    else if (pt<6.5) return std::make_tuple(0.99932, 0.00150513, 0.000387211, 0.00155414);
    else if (pt<7) return std::make_tuple(0.990015, 0.00168554, 0.000334128, 0.00171833);
    else if (pt<7.5) return std::make_tuple(0.989662, 0.00184611, 0.000377601, 0.00188433);
    else if (pt<8) return std::make_tuple(0.991271, 0.00199633, 0.000397838, 0.00203559);
    else if (pt<9.25) return std::make_tuple(0.993863, 0.00145115, 0.000219424, 0.00146765);
    else if (pt<10.5) return std::make_tuple(0.98847, 0.00181365, 0.00046839, 0.00187316);
    else if (pt<12.25) return std::make_tuple(0.988261, 0.0020079, 0.000342797, 0.00203695);
    else if (pt<14) return std::make_tuple(0.986502, 0.00269587, 0.000312191, 0.00271388);
    else if (pt<16) return std::make_tuple(0.967524, 0.00360689, 0.000816031, 0.00369805);
    else if (pt<18) return std::make_tuple(0.974107, 0.0047847, 0.00173171, 0.00508843);
    else if (pt<24) return std::make_tuple(0.971851, 0.00434543, 0.000750224, 0.00440971);
    else if (pt<30) return std::make_tuple(0.95726, 0.00782987, 0.00260996, 0.0082534);
  }

  if( (fabs(eta) > 1.2) && (fabs(eta) < 1.8) ){
    if (pt<2.5) return std::make_tuple(1.05392, 0.00741907, 0.0086659, 0.0114079);
    else if (pt<2.75) return std::make_tuple(1.05418, 0.00510863, 0.000493877, 0.00513245);
    else if (pt<3) return std::make_tuple(1.03165, 0.00362626, 0.00051552, 0.00366272);
    else if (pt<3.25) return std::make_tuple(1.00137, 0.00324096, 0.000650532, 0.00330561);
    else if (pt<3.5) return std::make_tuple(0.985531, 0.00325619, 0.00101988, 0.00341217);
    else if (pt<3.75) return std::make_tuple(0.979026, 0.00333719, 0.00080565, 0.00343306);
    else if (pt<4) return std::make_tuple(0.978533, 0.00348675, 0.00128796, 0.00371702);
    else if (pt<4.5) return std::make_tuple(0.981693, 0.00262222, 0.00138638, 0.00296615);
    else if (pt<5) return std::make_tuple(0.988393, 0.00282795, 0.000574534, 0.00288572);
    else if (pt<5.5) return std::make_tuple(0.988535, 0.00296014, 0.0015193, 0.00332727);
    else if (pt<6) return std::make_tuple(0.995795, 0.003052, 0.000849741, 0.00316809);
    else if (pt<7) return std::make_tuple(0.986439, 0.00239136, 0.00132843, 0.00273557);
    else if (pt<8.5) return std::make_tuple(0.991989, 0.0023577, 0.00141202, 0.00274819);
    else if (pt<10) return std::make_tuple(0.993719, 0.00296632, 0.00125018, 0.003219);
    else if (pt<12) return std::make_tuple(0.990475, 0.0031633, 0.00082619, 0.00326941);
    else if (pt<14.5) return std::make_tuple(0.99984, 0.00391145, 0.00226697, 0.00452091);
    else if (pt<19) return std::make_tuple(0.999167, 0.00490761, 0.00146068, 0.00512038);
    else if (pt<30) return std::make_tuple(0.998484, 0.00721995, 0.000958602, 0.00728331);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.1) ){
    if (pt<1.9) return std::make_tuple(1.16258, 0.00663765, 0.00518218, 0.008421);
    else if (pt<2.2) return std::make_tuple(1.08194, 0.00434659, 0.00165702, 0.00465173);
    else if (pt<2.5) return std::make_tuple(1.04721, 0.00279163, 0.00198365, 0.00342462);
    else if (pt<2.8) return std::make_tuple(1.0117, 0.00250143, 0.00103727, 0.00270797);
    else if (pt<3.1) return std::make_tuple(1.0045, 0.00241179, 0.000301947, 0.00243062);
    else if (pt<3.4) return std::make_tuple(0.998153, 0.00255819, 0.00107157, 0.00277355);
    else if (pt<4) return std::make_tuple(0.991499, 0.00205979, 0.000162552, 0.00206619);
    else if (pt<5) return std::make_tuple(0.982154, 0.00203163, 0.00109106, 0.00230607);
    else if (pt<6) return std::make_tuple(0.979281, 0.00269371, 0.00139839, 0.00303506);
    else if (pt<8) return std::make_tuple(0.978722, 0.00284354, 0.00105435, 0.00303272);
    else if (pt<11) return std::make_tuple(0.982305, 0.00415391, 0.000812637, 0.00423265);
    else if (pt<15) return std::make_tuple(0.98089, 0.00737357, 0.000877749, 0.00742563);
    else if (pt<30) return std::make_tuple(0.985428, 0.0132041, 0.012101, 0.0179104);
  }

  if( (fabs(eta) > 2.1) && (fabs(eta) < 2.4) ){
    if (pt<1.9) return std::make_tuple(1.10953, 0.00267574, 0.000888127, 0.00281928);
    else if (pt<2.3) return std::make_tuple(1.04759, 0.00206322, 0.00049553, 0.0021219);
    else if (pt<2.7) return std::make_tuple(1.02649, 0.00198497, 0.00156344, 0.00252675);
    else if (pt<3) return std::make_tuple(1.01133, 0.00249575, 0.000958998, 0.00267366);
    else if (pt<3.5) return std::make_tuple(1.01379, 0.00208413, 0.00104966, 0.00233354);
    else if (pt<4.25) return std::make_tuple(1.00341, 0.00202201, 0.000761361, 0.0021606);
    else if (pt<5) return std::make_tuple(0.997733, 0.00254249, 0.00299566, 0.00392915);
    else if (pt<6) return std::make_tuple(0.991922, 0.00275362, 0.000262191, 0.00276608);
    else if (pt<7) return std::make_tuple(0.9968, 0.00362768, 0.00239638, 0.00434772);
    else if (pt<8.5) return std::make_tuple(1.00011, 0.00361481, 0.00168078, 0.00398646);
    else if (pt<13) return std::make_tuple(1.00663, 0.00411005, 0.00220707, 0.00466516);
    else if (pt<30) return std::make_tuple(0.992803, 0.010354, 0.0109342, 0.0150586);
  }
  return std::make_tuple(1,0,0,0);

}

//To be used for tracker muon reconstruction efficiency, with loose muon kinematic acceptance cuts                                          
//Central scale-factor value, statistical error, systematic error, total error
std::tuple<double,double,double,double> tnp_weight_TM_LooseAcceptance_pp(double pt, double eta){

  if( (fabs(eta) > 0) && (fabs(eta) < 1.1) ){
    if (pt<3.8) return std::make_tuple(1.02525, 0.00394571, 0.00224281, 0.0045386);
    else if (pt<4.4) return std::make_tuple(1.00546, 0.00309167, 0.000778025, 0.00318806);
    else if (pt<5.2) return std::make_tuple(0.997276, 0.00275521, 0.000185849, 0.00276147);
    else if (pt<7) return std::make_tuple(1.00185, 0.00219137, 0.00134363, 0.0025705);
    else if (pt<30) return std::make_tuple(0.996559, 0.00203593, 0.000884415, 0.00221973);
  }

  if( (fabs(eta) > 1.1) && (fabs(eta) < 1.8) ){
    if (pt<2.8) return std::make_tuple(0.980744, 0.00554272, 0.00199386, 0.00589043);
    else if (pt<4) return std::make_tuple(0.999097, 0.00352956, 0.00185843, 0.00398893);
    else if (pt<8) return std::make_tuple(1.00099, 0.00247998, 0.0052363, 0.00579389);
    else if (pt<30) return std::make_tuple(1.00119, 0.00451636, 0.000920447, 0.0046092);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.4) ){
    if (pt<2.5) return std::make_tuple(1.04151, 0.029129, 0.0442224, 0.052954);
    else if (pt<3.7) return std::make_tuple(0.996334, 0.00508262, 0.00563347, 0.00758742);
    else if (pt<6) return std::make_tuple(1.00072, 0.00209385, 0.00196219, 0.00286957);
    else if (pt<30) return std::make_tuple(0.998784, 0.00200216, 0.00681056, 0.00709876);
  }

  return std::make_tuple(1,0,0,0);
}

//To be used for SoftID efficiency, with loose muon kinematic acceptance cuts                                                                                                                  
//Central scale-factor value, statistical error, systematic error, total error                                                                                                                       
std::tuple<double,double,double,double> tnp_weight_SoftID_LooseAcceptance_pp(double pt, double eta){

  if( (fabs(eta) > 0) && (fabs(eta) < 1.1) ){
    if (pt<3.8) return std::make_tuple(0.999324, 0.000142399, 3.98581e-05, 0.000147872);
    else if (pt<4.4) return std::make_tuple(0.999053, 0.000133731, 3.58144e-05, 0.000138444);
    else if (pt<5.2) return std::make_tuple(0.999506, 0.000105829, 2.87498e-05, 0.000109665);
    else if (pt<7) return std::make_tuple(0.999614, 8.10057e-05, 4.25893e-06, 8.11176e-05);
    else if (pt<30) return std::make_tuple(0.999995, 3.57268e-05, 4.66432e-06, 3.603e-05);
  }
  
  if( (fabs(eta) > 1.1) && (fabs(eta) < 1.8) ){
    if (pt<2.8) return std::make_tuple(0.999011, 0.000320735, 0.000245804, 0.000404093);
    else if (pt<4) return std::make_tuple(0.999208, 0.000160307, 0.000129319, 0.000205965);
    else if (pt<8) return std::make_tuple(0.999498, 0.000103203, 2.62778e-05, 0.000106496);
    else if (pt<30) return std::make_tuple(1, 1.77276e-05, 2.36063e-07, 1.77292e-05);
  }
  
  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.4) ){
    if (pt<2.5) return std::make_tuple(0.997624, 0.000434151, 0.000260878, 0.000506502);
    else if (pt<3.7) return std::make_tuple(0.999134, 0.000274953, 8.50146e-05, 0.000287796);
    else if (pt<6) return std::make_tuple(0.999693, 0.000200237, 5.08814e-05, 0.000206601);
    else if (pt<30) return std::make_tuple(0.999584, 0.000228825, 9.01656e-05, 0.000245948);
  }
  return std::make_tuple(1,0,0,0);
}

#endif
