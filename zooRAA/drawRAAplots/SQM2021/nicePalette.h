#include <iostream>

TColor *pal = new TColor();     

Int_t kmagenta = pal->GetColor(124,  0,124);
Int_t kviolet  = pal->GetColor( 72,  0,190);
Int_t kblue    = pal->GetColor(  9,  0,200);
Int_t kazure   = pal->GetColor(  0, 48, 97);
Int_t kcyan    = pal->GetColor(  0, 83, 98);
Int_t kteal    = pal->GetColor(  0, 92, 46);
Int_t kgreen   = pal->GetColor( 15, 85, 15);
Int_t kspring  = pal->GetColor( 75, 97, 53);
Int_t kyellow  = pal->GetColor(117,118,  0);
Int_t korange  = pal->GetColor(101, 42,  0);
Int_t kred     = pal->GetColor(190,  0,  3);
Int_t kpink    = pal->GetColor(180, 35,145);
// good for systematic band fill
Int_t kmagentaLight = pal->GetColor(215,165,215);
Int_t kvioletLight  = pal->GetColor(200,160,255);
Int_t kblueLight    = pal->GetColor(178,185,254);
Int_t kazureLight   = pal->GetColor(153,195,225);
Int_t kcyanLight    = pal->GetColor(140,209,224);
Int_t ktealLight    = pal->GetColor( 92,217,141);
Int_t kgreenLight   = pal->GetColor(135,222,135);
Int_t kspringLight  = pal->GetColor(151,207,116);
Int_t kyellowLight  = pal->GetColor(225,225,100);
Int_t korangeLight  = pal->GetColor(255,168,104);
Int_t kredLight     = pal->GetColor(253,169,179);
Int_t kpinkLight    = pal->GetColor(255,192,224);

void nicePalette(){
  std::cout <<kmagenta<<std::endl;
}
