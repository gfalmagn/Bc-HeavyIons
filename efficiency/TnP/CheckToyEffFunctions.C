#include "ToyEffFunctions.h"
#include "../../acceptance/SgMuonAcceptanceCuts.h" //acceptance cuts here

double efficiency(double *x, double *p){
  //x = {pt,eta}
  //p[0] = 0 (data) or 1 (MC)
  //p[1] = 0 (loose efficiency all->l), 1 (tight efficiency l->t), 2 (tighter efficiency t->t'), 3 = 2*3 (l->t'), 4 = 1*2 (all->t), 5 = 1*2*3 (all->t')
  if(!looseAcc(x[0],x[1])) return 0;

  if(p[0]<0.5){ //DATA
    if(p[1]<0.5) return Efficiency1(x[0],x[1],false);
    if(p[1]>0.5 && p[1]<1.5) return Efficiency2(x[0],x[1],false);
    if(p[1]>1.5 && p[1]<2.5) return (Efficiency2(x[0],x[1],false)==0)?0:( Efficiency2times3(x[0],x[1],false)/Efficiency2(x[0],x[1],false) );
    if(p[1]>2.5 && p[1]<3.5) return Efficiency2times3(x[0],x[1],false);
    if(p[1]>3.5 && p[1]<4.5) return (Efficiency2(x[0],x[1],false)==0)?0:( Efficiency1(x[0],x[1],false) * Efficiency2times3(x[0],x[1],false)/Efficiency2(x[0],x[1],false) );
    else return Efficiency1(x[0],x[1],false)*Efficiency2times3(x[0],x[1],false);
  }
  else { //MC
    if(p[1]<0.5) return Efficiency1(x[0],x[1],true);
    if(p[1]>0.5 && p[1]<1.5) return Efficiency2(x[0],x[1],true);
    if(p[1]>1.5 && p[1]<2.5) return (Efficiency2(x[0],x[1],true)==0)?0:( Efficiency2times3(x[0],x[1],true)/Efficiency2(x[0],x[1],true) );
    if(p[1]>2.5 && p[1]<3.5) return Efficiency2times3(x[0],x[1],true);
    if(p[1]>3.5 && p[1]<4.5) return (Efficiency2(x[0],x[1],true)==0)?0:( Efficiency1(x[0],x[1],true) * Efficiency2times3(x[0],x[1],true)/Efficiency2(x[0],x[1],true) );
    else return Efficiency1(x[0],x[1],true)*Efficiency2times3(x[0],x[1],true);
  }
}

void CheckToyEffFunctions(){

  TF2 *eff = new TF2("eff",efficiency, 0,15, -2.4,2.4 ,2); //2 parameters
  TF2 *eff2 = new TF2("eff2",efficiency, 0,15, -2.4,2.4, 2);
  TF2 *eff3 = new TF2("eff3",efficiency, 0,15, -2.4,2.4, 2);
  eff->SetNpx(80);
  eff2->SetNpx(80);
  eff3->SetNpx(80);
  
  TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
  eff->SetParameters(0, 0); //data l
  eff->Draw("surf");  
  eff2->SetParameters(0, 1); //data l->t
  eff2->SetLineColor(kBlue);
  eff2->Draw("surfsame");  
  eff3->SetParameters(0, 2); //data t->t'
  eff3->SetLineColor(kGreen+2);
  eff3->Draw("surfsame");  

  TF2 *eff4 = new TF2("eff4",efficiency, 0,15, -2.4,2.4 ,2); //2 parameters
  TF2 *eff5 = new TF2("eff5",efficiency, 0,15, -2.4,2.4, 2);
  TF2 *eff6 = new TF2("eff6",efficiency, 0,15, -2.4,2.4, 2);
  eff4->SetNpx(80);
  eff5->SetNpx(80);
  eff6->SetNpx(80);

  TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
  eff4->SetParameters(1, 0); //mc l
  eff4->Draw("surf");  
  eff5->SetParameters(1, 1); //mc l->t
  eff5->SetLineColor(kBlue);
  eff5->Draw("surfsame");  
  eff6->SetParameters(1, 2); //mc t->t'
  eff6->SetLineColor(kGreen+2);
  eff6->Draw("surfsame");  

  //"loose" efficiency should be:
  //eff->SetParameters(dataOrMC, 0);
  //"normal" efficiency (product of the first 2 efficiencies) should be:
  //eff->SetParameters(dataOrMC, 4);
  //"tight" efficiency (product of the 3 efficiencies) should be:
  //eff->SetParameters(dataOrMC, 5);
  
}
