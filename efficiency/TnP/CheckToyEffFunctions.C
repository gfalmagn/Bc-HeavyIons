#include "ToyEffFunctions.h"
#include "../../helpers/SgMuonAcceptanceCuts.h" //acceptance cuts here

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
    if(p[1]>3.5 && p[1]<4.5) return Efficiency1(x[0],x[1],false) * Efficiency2(x[0],x[1],false);
    else return Efficiency1(x[0],x[1],false)*Efficiency2times3(x[0],x[1],false);
  }
  else { //MC
    if(p[1]<0.5) return Efficiency1(x[0],x[1],true);
    if(p[1]>0.5 && p[1]<1.5) return Efficiency2(x[0],x[1],true);
    if(p[1]>1.5 && p[1]<2.5) return (Efficiency2(x[0],x[1],true)==0)?0:( Efficiency2times3(x[0],x[1],true)/Efficiency2(x[0],x[1],true) );
    if(p[1]>2.5 && p[1]<3.5) return Efficiency2times3(x[0],x[1],true);
    if(p[1]>3.5 && p[1]<4.5) return Efficiency1(x[0],x[1],true) * Efficiency2(x[0],x[1],true);
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
  eff->GetXaxis()->SetTitle("p_{T}");
  eff->GetYaxis()->SetTitle("#eta");
  eff->Draw("surf");  
  eff2->SetParameters(0, 4); //data t
  eff2->SetLineColor(kBlue);
  eff2->Draw("surfsame");  
  eff3->SetParameters(0, 5); //data t'
  eff3->SetLineColor(kGreen+2);
  eff3->Draw("surfsame");  

  c1->SaveAs("ToyEffFunctions_data_2D.pdf");

  TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
  vector<TF1*> effpT,eff2pT,eff3pT;
  float etalim[] = {0,1.2,1.8,2.1,2.4};
  Color_t cols[] = {kBlack,kBlue,kRed,kGreen};

  effpT.push_back( new TF1("effpT_eta0","[&](double *x, double *p){ return eff->Eval(x[0],0.1);}",0,15,0) );
  effpT.push_back( new TF1("effpT_eta1","[&](double *x, double *p){ return eff->Eval(x[0],1.3);}",0,15,0) );
  effpT.push_back( new TF1("effpT_eta2","[&](double *x, double *p){ return eff->Eval(x[0],1.9);}",0,15,0) );
  effpT.push_back( new TF1("effpT_eta3","[&](double *x, double *p){ return eff->Eval(x[0],2.2);}",0,15,0) );
  eff2pT.push_back( new TF1("eff2pT_eta0","[&](double *x, double *p){ return eff2->Eval(x[0],0.1);}",0,15,0) );
  eff2pT.push_back( new TF1("eff2pT_eta1","[&](double *x, double *p){ return eff2->Eval(x[0],1.3);}",0,15,0) );
  eff2pT.push_back( new TF1("eff2pT_eta2","[&](double *x, double *p){ return eff2->Eval(x[0],1.9);}",0,15,0) );
  eff2pT.push_back( new TF1("eff2pT_eta3","[&](double *x, double *p){ return eff2->Eval(x[0],2.2);}",0,15,0) );
  eff3pT.push_back( new TF1("eff3pT_eta0","[&](double *x, double *p){ return eff3->Eval(x[0],0.1);}",0,15,0) );
  eff3pT.push_back( new TF1("eff3pT_eta1","[&](double *x, double *p){ return eff3->Eval(x[0],1.3);}",0,15,0) );
  eff3pT.push_back( new TF1("eff3pT_eta2","[&](double *x, double *p){ return eff3->Eval(x[0],1.9);}",0,15,0) );
  eff3pT.push_back( new TF1("eff3pT_eta3","[&](double *x, double *p){ return eff3->Eval(x[0],2.2);}",0,15,0) );

  TLegend *leg = new TLegend(0.6,0.1,0.9,0.3);
  for(int i=0;i<4;i++){
    effpT[i]->SetLineColor(cols[i]);
    eff2pT[i]->SetLineColor(cols[i]);
    eff3pT[i]->SetLineColor(cols[i]);
    leg->AddEntry(effpT[i] , Form("%.1f<|#eta|<%.1f",etalim[i],etalim[i+1]), "l");
    effpT[i]->GetXaxis()->SetTitle("p_{T}");
    effpT[i]->GetYaxis()->SetTitle("single-muon efficiency");
    effpT[i]->SetTitle("Loose, medium, tight efficiencies (data)");
    effpT[i]->Draw((i==0)?"":"same");
    eff2pT[i]->Draw("same");
    eff3pT[i]->Draw("same");
  }
  leg->SetBorderSize(0);
  leg->Draw("same");

  c3->SaveAs("ToyEffFunctions_data_1D.pdf");

  TF2 *eff4 = new TF2("eff4",efficiency, 0,15, -2.4,2.4 ,2); //2 parameters
  TF2 *eff5 = new TF2("eff5",efficiency, 0,15, -2.4,2.4, 2);
  TF2 *eff6 = new TF2("eff6",efficiency, 0,15, -2.4,2.4, 2);
  eff4->SetNpx(80);
  eff5->SetNpx(80);
  eff6->SetNpx(80);

  TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
  eff4->SetParameters(1, 0); //mc l
  eff4->GetXaxis()->SetTitle("p_{T}");
  eff4->GetYaxis()->SetTitle("#eta");
  eff4->Draw("surf");  
  eff5->SetParameters(1, 4); //mc t
  eff5->SetLineColor(kBlue);
  eff5->Draw("surfsame");  
  eff6->SetParameters(1, 5); //mc t'
  eff6->SetLineColor(kGreen+2);
  eff6->Draw("surfsame");  

  c2->SaveAs("ToyEffFunctions_MC_2D.pdf");


  TCanvas *c4 = new TCanvas("c4","c4",1500,1500);
  vector<TF1*> eff4pT,eff5pT,eff6pT;

  eff4pT.push_back( new TF1("eff4pT_eta0","[&](double *x, double *p){ return eff4->Eval(x[0],0.1);}",0,15,0) );
  eff4pT.push_back( new TF1("eff4pT_eta1","[&](double *x, double *p){ return eff4->Eval(x[0],1.3);}",0,15,0) );
  eff4pT.push_back( new TF1("eff4pT_eta2","[&](double *x, double *p){ return eff4->Eval(x[0],1.9);}",0,15,0) );
  eff4pT.push_back( new TF1("eff4pT_eta3","[&](double *x, double *p){ return eff4->Eval(x[0],2.2);}",0,15,0) );
  eff5pT.push_back( new TF1("eff5pT_eta0","[&](double *x, double *p){ return eff5->Eval(x[0],0.1);}",0,15,0) );
  eff5pT.push_back( new TF1("eff5pT_eta1","[&](double *x, double *p){ return eff5->Eval(x[0],1.3);}",0,15,0) );
  eff5pT.push_back( new TF1("eff5pT_eta2","[&](double *x, double *p){ return eff5->Eval(x[0],1.9);}",0,15,0) );
  eff5pT.push_back( new TF1("eff5pT_eta3","[&](double *x, double *p){ return eff5->Eval(x[0],2.2);}",0,15,0) );
  eff6pT.push_back( new TF1("eff6pT_eta0","[&](double *x, double *p){ return eff6->Eval(x[0],0.1);}",0,15,0) );
  eff6pT.push_back( new TF1("eff6pT_eta1","[&](double *x, double *p){ return eff6->Eval(x[0],1.3);}",0,15,0) );
  eff6pT.push_back( new TF1("eff6pT_eta2","[&](double *x, double *p){ return eff6->Eval(x[0],1.9);}",0,15,0) );
  eff6pT.push_back( new TF1("eff6pT_eta3","[&](double *x, double *p){ return eff6->Eval(x[0],2.2);}",0,15,0) );

  TLegend *leg2 = new TLegend(0.6,0.1,0.9,0.3);
  for(int i=0;i<4;i++){
    eff4pT[i]->SetLineColor(cols[i]);
    eff5pT[i]->SetLineColor(cols[i]);
    eff6pT[i]->SetLineColor(cols[i]);
    leg2->AddEntry(eff4pT[i] , Form("%.1f<|#eta|<%.1f",etalim[i],etalim[i+1]), "l");
    eff4pT[i]->GetXaxis()->SetTitle("p_{T}");
    eff4pT[i]->GetYaxis()->SetTitle("single-muon efficiency");
    eff4pT[i]->SetTitle("Loose, medium, tight efficiencies (MC)");
    eff4pT[i]->Draw((i==0)?"":"same");
    eff5pT[i]->Draw("same");
    eff6pT[i]->Draw("same");
  }
  leg2->SetBorderSize(0);
  leg2->Draw("same");

  c4->SaveAs("ToyEffFunctions_MC_1D.pdf");








  int n=100;
  for(int i=0;i<n;i++){  
    for(int j=0;j<n;j++){  
      float x = i*20./n;
      float y = -2.4+i*4.8/n; 
      // cout<<"x,y,eff,eff2,eff3 = "<<x<<" "<<y<<" "<<eff->Eval(x,y)<<" "<<eff2->Eval(x,y)<<" "<<eff3->Eval(x,y)<<endl;
      // cout<<"x,y,eff4,eff5,eff6 = "<<x<<" "<<y<<" "<<eff4->Eval(x,y)<<" "<<eff5->Eval(x,y)<<" "<<eff6->Eval(x,y)<<endl;
      if(eff->Eval(x,y)<eff2->Eval(x,y)) cout<<"WARNING: data l<t for pt,eta = "<<x<<" "<<y<<endl;
      if(eff2->Eval(x,y)<eff3->Eval(x,y)) cout<<"WARNING: data t<t' for pt,eta = "<<x<<" "<<y<<endl;
      if(eff4->Eval(x,y)<eff5->Eval(x,y)) cout<<"WARNING: MC l<t for pt,eta = "<<x<<" "<<y<<endl;
      if(eff5->Eval(x,y)<eff6->Eval(x,y)) cout<<"WARNING: MC t<t' for pt,eta = "<<x<<" "<<y<<endl;
    }
  }
  //"loose" efficiency should be:
  //eff->SetParameters(dataOrMC, 0);
  //"normal" efficiency (product of the first 2 efficiencies) should be:
  //eff->SetParameters(dataOrMC, 4);
  //"tight" efficiency (product of the 3 efficiencies) should be:
  //eff->SetParameters(dataOrMC, 5);
  
}
