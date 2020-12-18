///////////////////////////////////////////////////
//      M u I D  LOOSE ACCEPTANCE  P b P b       //
///////////////////////////////////////////////////

double tnp_weight_muid_looseacceptance_pbpb(double pt, double eta, int idx) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (idx != 99) {
    // SF for 0 < |eta| < 1.1
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      if (x >= 3.3 && x <4.5) {num = 0.948733; den = 0.955024; statUp = 0.00380327; statDown = 0.00391521;}
      else if (x >= 4.5 && x <5) {num = 0.9866; den = 0.985274; statUp = 0.0041967; statDown = 0.0045537;}
      else if (x >= 5 && x <6.5) {num = 0.987123; den = 0.99045; statUp = 0.00263437; statDown = 0.00279898;}
      else if (x >= 6.5 && x <8) {num = 0.989174; den = 0.992717; statUp = 0.00297485; statDown = 0.00330039;}
      else if (x >= 8 && x <10.5) {num = 0.989445; den = 0.988638; statUp = 0.00331444; statDown = 0.00372826;}
      else {num = 0.982816; den = 0.97592; statUp = 0.00439878; statDown = 0.00484155;}
    }
    // SF for 1.1 < |eta| < 1.8
    if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      if (x >= 1.75 && x <4) {num = 0.982791; den = 0.971937; statUp = 0.00479601; statDown = 0.00490521;}
      else if (x >= 4 && x <5) {num = 0.987852; den = 0.989626; statUp = 0.00347727; statDown = 0.00372541;}
      else if (x >= 5 && x <6) {num = 0.985747; den = 0.993259; statUp = 0.00436325; statDown = 0.0048202;}
      else if (x >= 6 && x <7.5) {num = 0.988915; den = 0.995815; statUp = 0.00435158; statDown = 0.00490696;}
      else if (x >= 7.5 && x <15) {num = 0.982434; den = 0.984478; statUp = 0.00456413; statDown = 0.00502111;}
      else {num = 0.972038; den = 0.981322; statUp = 0.0140267; statDown = 0.0166183;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.33 && x <3) {num = 0.967731; den = 0.956475; statUp = 0.0133546; statDown = 0.0133268;}
      else if (x >= 3 && x <3.5) {num = 0.988271; den = 0.991548; statUp = 0.00658562; statDown = 0.00748651;}
      else if (x >= 3.5 && x <4.5) {num = 0.991137; den = 0.997078; statUp = 0.00406749; statDown = 0.00471666;}
      else if (x >= 4.5 && x <5.5) {num = 0.996091; den = 0.998014; statUp = 0.00334903; statDown = 0.00440572;}
      else if (x >= 5.5 && x <7) {num = 1; den = 0.998654; statUp = 2.35228e-07; statDown = 0.00316994;}
      else if (x >= 7 && x <9) {num = 0.991205; den = 0.999535; statUp = 0.00595637; statDown = 0.00734183;}
      else {num = 1; den = 0.997536; statUp = 7.62371e-09; statDown = 0.00194303;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.2 && x <2.7) {num = 0.90688; den = 0.912775; statUp = 0.0233261; statDown = 0.0226667;}
      else if (x >= 2.7 && x <3.7) {num = 0.977493; den = 0.979099; statUp = 0.00763304; statDown = 0.00812459;}
      else if (x >= 3.7 && x <8) {num = 0.998764; den = 0.996218; statUp = 0.00123579; statDown = 0.00376785;}
      else if (x >= 8 && x <11) {num = 1; den = 0.998724; statUp = 2.73793e-09; statDown = 0.00535402;}
      else {num = 1; den = 0.996093; statUp = 1.56969e-08; statDown = 0.0071709;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      // syst uncertainties
      if (x >= 3.3 && x < 4.5) syst = 0.0013846;
      else if (x >= 4.5 && x < 5) syst = 0.000744057;
      else if (x >= 5 && x < 6.5) syst = 0.000224688;
      else if (x >= 6.5 && x < 8) syst = 0.000860391;
      else if (x >= 8 && x < 10.5) syst = 0.00209347;
      else syst = 0.00269482;
    }
    else if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 1.75 && x < 4) syst = 0.000881792;
      else if (x >= 4 && x < 5) syst = 0.00158573;
      else if (x >= 5 && x < 6) syst = 0.000957972;
      else if (x >= 6 && x < 7.5) syst = 0.00208154;
      else if (x >= 7.5 && x < 15) syst = 0.000923833;
      else syst = 0.0044031;
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.33 && x < 3) syst = 0.00955292;
      else if (x >= 3 && x < 3.5) syst = 0.000328874;
      else if (x >= 3.5 && x < 4.5) syst = 0.00355421;
      else if (x >= 4.5 && x < 5.5) syst = 0.000508896;
      else if (x >= 5.5 && x < 7) syst = 2.35114e-07;
      else if (x >= 7 && x < 9) syst = 0.00487386;
      else syst = 7.58011e-09;
    }
    else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.2 && x < 2.7) syst = 0.0100094;
      else if (x >= 2.7 && x < 3.7) syst = 0.00625813;
      else if (x >= 3.7 && x < 8) syst = 0.00283863;
      else if (x >= 8 && x < 11) syst = 2.7243e-09;
      else syst = 1.55679e-08;
    }
  }
  else {
    // SF for 0 < |eta| < 1.1
    if (fabs(eta) >= 0 && fabs(eta) < 1.1) {
      if (x >= 3.3 && x <4.5) {num = 0.938695; den = 0.955514; statUp = 0.00664236; statDown = 0.00687315;}
      else if (x >= 4.5 && x <5) {num = 0.98844; den = 0.984727; statUp = 0.00616829; statDown = 0.00703819;}
      else if (x >= 5 && x <6.5) {num = 0.993642; den = 0.990815; statUp = 0.00423012; statDown = 0.00474381;}
      else if (x >= 6.5 && x <8) {num = 0.989747; den = 0.99268; statUp = 0.00532569; statDown = 0.00630142;}
      else if (x >= 8 && x <10.5) {num = 0.985607; den = 0.989021; statUp = 0.00606928; statDown = 0.00732847;}
      else {num = 0.97083; den = 0.975779; statUp = 0.00989657; statDown = 0.0113238;}
    }
    // SF for 1.1 < |eta| < 1.8
    if (fabs(eta) >= 1.1 && fabs(eta) < 1.8) {
      if (x >= 1.75 && x <4) {num = 0.980611; den = 0.972609; statUp = 0.0073829; statDown = 0.007479;}
      else if (x >= 4 && x <5) {num = 0.9777; den = 0.989743; statUp = 0.00642669; statDown = 0.00683214;}
      else if (x >= 5 && x <6) {num = 0.997809; den = 0.993825; statUp = 0.00656549; statDown = 0.00656549;}
      else if (x >= 6 && x <7.5) {num = 0.990085; den = 0.995743; statUp = 0.00684176; statDown = 0.00790556;}
      else if (x >= 7.5 && x <15) {num = 0.999376; den = 0.985708; statUp = 0.000624266; statDown = 0.00742711;}
      else {num = 0.973238; den = 0.983101; statUp = 0.0267622; statDown = 0.0407917;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.33 && x <3) {num = 0.967203; den = 0.957168; statUp = 0.0159714; statDown = 0.0157842;}
      else if (x >= 3 && x <3.5) {num = 0.99004; den = 0.991599; statUp = 0.00996034; statDown = 0.011883;}
      else if (x >= 3.5 && x <4.5) {num = 0.993569; den = 0.996469; statUp = 0.00594978; statDown = 0.00693082;}
      else if (x >= 4.5 && x <5.5) {num = 0.995885; den = 0.998046; statUp = 0.00411469; statDown = 0.00664728;}
      else if (x >= 5.5 && x <7) {num = 1; den = 0.998865; statUp = 1.09426e-07; statDown = 0.006463;}
      else if (x >= 7 && x <9) {num = 0.993975; den = 0.999619; statUp = 0.00602541; statDown = 0.0140919;}
      else {num = 0.991651; den = 0.997466; statUp = 0.00834881; statDown = 0.016729;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.2 && x <2.7) {num = 0.894759; den = 0.910137; statUp = 0.0267882; statDown = 0.0256294;}
      else if (x >= 2.7 && x <3.7) {num = 0.97625; den = 0.983037; statUp = 0.011189; statDown = 0.0116927;}
      else if (x >= 3.7 && x <8) {num = 0.991373; den = 0.99574; statUp = 0.00628986; statDown = 0.00699593;}
      else if (x >= 8 && x <11) {num = 1; den = 0.998833; statUp = 9.63229e-12; statDown = 0.00513902;}
      else {num = 1; den = 0.996783; statUp = 3.178e-08; statDown = 0.0187675;}
    }
  }

  double syst_factor = 0; double stat_factor = 0;
  if (idx == 3) return den;
  else if (idx == -1) syst_factor = syst;
  else if (idx == -2) syst_factor = -1*syst;
  else if (idx == +1) stat_factor = statUp;
  else if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

///////////////////////////////////////////////////
//   T R K  LOOSE ACCEPTANCE     P b P b         //
///////////////////////////////////////////////////

double tnp_weight_trk_looseacceptance_pbpb(double eta, int idx) {
  double x = eta;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (idx != 99) {
    //SF in eta bins
    if (x >= -2.4 && x < -1.6) {num = 0.994686; den = 0.998801; statUp = 0.00454507; statDown = 0.00462482;}
    if (x >= -1.6 && x < -1.1) {num = 0.975586; den = 0.968435; statUp = 0.0060169; statDown = 0.00621157;}
    if (x >= -1.1 && x < -0.9) {num = 0.97513; den = 0.97459; statUp = 0.0157419; statDown = 0.0159012;}
    if (x >= -0.9 && x < -0.6) {num = 0.965031; den = 0.973323; statUp = 0.0108792; statDown = 0.011148;}
    if (x >= -0.6 && x < -0.3) {num = 0.953727; den = 0.976661; statUp = 0.0102523; statDown = 0.010391;}
    if (x >= -0.3 && x < 0.3) {num = 0.965483; den = 0.967508; statUp = 0.00794441; statDown = 0.00810298;}
    if (x >= 0.3 && x < 0.6) {num = 0.960845; den = 0.967106; statUp = 0.0114093; statDown = 0.0116539;}
    if (x >= 0.6 && x < 0.9) {num = 0.972216; den = 0.965041; statUp = 0.0113936; statDown = 0.0117501;}
    if (x >= 0.9 && x < 1.1) {num = 0.953803; den = 0.963747; statUp = 0.0196302; statDown = 0.0197941;}
    if (x >= 1.1 && x < 1.6) {num = 0.951295; den = 0.9643; statUp = 0.00684901; statDown = 0.00704341;}
    if (x >= 1.6 && x < 2.4) {num = 0.998071; den = 0.998608; statUp = 0.00192877; statDown = 0.00508664;}

    // syst uncertainties
    if (x >= -2.4 && x < -1.6) syst = 0.0236203;
    else if (x >= -1.6 && x < -1.1) syst = 0.003712;
    else if (x >= -1.1 && x < -0.9) syst = 0.0145127;
    else if (x >= -0.9 && x < -0.6) syst = 0.00364391;
    else if (x >= -0.6 && x < -0.3) syst = 0.0046084;
    else if (x >= -0.3 && x < 0.3) syst = 0.00482188;
    else if (x >= 0.3 && x < 0.6) syst = 0.0111176;
    else if (x >= 0.6 && x < 0.9) syst = 0.00135623;
    else if (x >= 0.9 && x < 1.1) syst = 0.0143408;
    else if (x >= 1.1 && x < 1.6) syst = 0.00675803;
    else syst = 0.00192877;
  }
  else {
    //SF in eta bins
    if (x >= -2.4 && x < -1.6) {num = 1; den = 0.998369; statUp = 9.56358e-08; statDown = 0.00749361;}
    if (x >= -1.6 && x < -1.1) {num = 0.967693; den = 0.966985; statUp = 0.00994236; statDown = 0.0101406;}
    if (x >= -1.1 && x < -0.9) {num = 0.928328; den = 0.973725; statUp = 0.0329219; statDown = 0.0328953;}
    if (x >= -0.9 && x < -0.6) {num = 0.936238; den = 0.973701; statUp = 0.0383609; statDown = 0.036511;}
    if (x >= -0.6 && x < -0.3) {num = 0.971953; den = 0.978004; statUp = 0.0177767; statDown = 0.0183502;}
    if (x >= -0.3 && x < 0.3) {num = 0.957823; den = 0.9665; statUp = 0.0155799; statDown = 0.0159396;}
    if (x >= 0.3 && x < 0.6) {num = 0.952916; den = 0.967395; statUp = 0.0183144; statDown = 0.0191085;}
    if (x >= 0.6 && x < 0.9) {num = 0.942845; den = 0.966915; statUp = 0.0244783; statDown = 0.0254036;}
    if (x >= 0.9 && x < 1.1) {num = 1; den = 0.961867; statUp = 1.80899e-07; statDown = 0.0211132;}
    if (x >= 1.1 && x < 1.6) {num = 0.972972; den = 0.962446; statUp = 0.0104528; statDown = 0.0106489;}
    if (x >= 1.6 && x < 2.4) {num = 0.999999; den = 0.998751; statUp = 5.60152e-07; statDown = 0.00472985;}
  }

  double syst_factor = 0; double stat_factor = 0;
  if (idx == 3) return den;
  else if (idx == -1) syst_factor = syst;
  else if (idx == -2) syst_factor = -1*syst;
  else if (idx == +1) stat_factor = statUp;
  else if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

