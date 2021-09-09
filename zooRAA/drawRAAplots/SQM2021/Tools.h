#include "TGraphAsymmErrors.h"

void SetEx(TGraphAsymmErrors* gae, Double_t Ex)
{
  Int_t np = gae->GetN();
  for (Int_t i=0; i<np; i++) {
    gae->SetPointEXhigh(i,Ex);
    gae->SetPointEXlow(i,Ex);
  }
}
