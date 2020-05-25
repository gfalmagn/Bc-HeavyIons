
// bool trkAcc(float pt, float eta){
//   return (fabs(eta) < 2.4 &&
// 	  ((fabs(eta) < 0.8 && pt >= 3.3) ||
// 	   (0.8 <= fabs(eta) && fabs(eta) < 1.5 && pt >= 5.81-3.14*fabs(eta) ) ||
// 	   (1.5 <= fabs(eta) && pt >= 0.8 && pt >= 1.89-0.526*fabs(eta) )));
// }

bool looseAcc(float pt, float eta, bool TM = false){
  if(TM) return (fabs(eta) < 2.4 &&
		 ((fabs(eta) < 1.1 && pt >= 3.3) ||
		  (1.1 <= fabs(eta) && fabs(eta) < 1.3 && pt >= 13.2-9.*fabs(eta) ) ||
		  (1.3 <= fabs(eta) && pt >= 0.8 && pt >= 3.02-1.17*fabs(eta) )));
  else return (fabs(eta) < 2.4 &&
	       ((fabs(eta) < 0.3 && pt >= 3.4) ||
		(fabs(eta) > 0.3 && fabs(eta) < 1.1 && pt >= 3.3) ||
		(fabs(eta) > 1.1 && fabs(eta) < 1.55 && pt >= 2.1 && pt >= 7.7-4.0*fabs(eta) ) ||
		(fabs(eta) > 1.55 && pt >= 1.2 && pt >= 4.25-1.39*fabs(eta)) ));
}

bool tightAcc(float pt, float eta){
  return (fabs(eta) < 2.4 &&
	 ((fabs(eta) < 1.2 && pt >= 3.5) ||
	  (1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 5.47-1.89*fabs(eta)) ||
	  (2.1 <= fabs(eta) && pt >= 1.5)));
}

bool InAcc(TLorentzVector mu1, TLorentzVector mu2, TLorentzVector mu3, bool TM = false){
  return (   (looseAcc(mu1.Pt(),mu1.Eta(),TM) && tightAcc(mu2.Pt(),mu2.Eta()) && tightAcc(mu3.Pt(),mu3.Eta()) )
	  || (tightAcc(mu1.Pt(),mu1.Eta()) && looseAcc(mu2.Pt(),mu2.Eta(),TM) && tightAcc(mu3.Pt(),mu3.Eta()) )
	  || (tightAcc(mu1.Pt(),mu1.Eta()) && tightAcc(mu2.Pt(),mu2.Eta()) && looseAcc(mu3.Pt(),mu3.Eta(),TM) )
	  );
}
