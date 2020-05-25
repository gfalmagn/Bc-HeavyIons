
//TH2Poly's for sophisticated shape
TH2Poly *_hp(){ 

  TH2Poly *res = new TH2Poly();
  //pt 6-10, y 0-1.
  res->AddBin(0,7,0.4,8); res->AddBin(0.4,7,0.6,8); res->AddBin(0.6,7,0.8,8); res->AddBin(0.8,7,1.,8); res->AddBin(0.8,6,1.,7);
  res->AddBin(0,8,0.5,10); res->AddBin(0.5,8,0.8,10); res->AddBin(0.8,8,1.,9); res->AddBin(0.8,9,1.,10);  
  //pt 6-10, y 1-1.2
  res->AddBin(1.,6,1.2,8); res->AddBin(1.,8,1.2,9); res->AddBin(1.,9,1.2,10);
  //pt 0-2
  res->AddBin(1.6,0,1.7,2); res->AddBin(1.7,0,1.85,2); res->AddBin(1.85,0,1.95,2); res->AddBin(1.95,0,2.05,2); res->AddBin(2.05,0,2.15,2); res->AddBin(2.15,0,2.3,2);
  //pt 2-4
  res->AddBin(1.2,2,1.4,4); res->AddBin(1.4,2,1.6,4); for(int i=0;i<=6;i++) res->AddBin( 1.6+i*0.1 ,2, 1.7+i*0.1 ,4);
  //pt 4-6
  res->AddBin(1.,4,1.25,6); res->AddBin(1.2,4,1.3,6); res->AddBin(1.3,4,1.45,5); res->AddBin(1.3,5,1.45,6);  res->AddBin(1.45,4,1.6,5); res->AddBin(1.45,5,1.6,6); 
  for(int i=0;i<=6;i++) {
    res->AddBin( 1.6+i*0.1 ,4, 1.7+i*0.1 ,5);
    res->AddBin( 1.6+i*0.1 ,5, 1.7+i*0.1 ,6);
  }

  //pt 6-8
  res->AddBin(1.2,6,1.3,7); res->AddBin(1.2,7,1.3,8); res->AddBin(1.3,6,1.45,7); res->AddBin(1.3,7,1.45,8); res->AddBin(1.45,6,1.6,7); res->AddBin(1.45,7,1.6,8);
  for(int i=0;i<=3;i++) {
    res->AddBin( 1.6+i*0.1 ,7, 1.7+i*0.1 ,8);
    res->AddBin( 1.6+i*0.1 ,6, 1.7+i*0.1 ,7);
  }
  res->AddBin(2.,6,2.15,7); res->AddBin(2.,7,2.15,8); res->AddBin(2.15,6,2.3,7); res->AddBin(2.15,7,2.3,8);
  //pt 8-10
  res->AddBin(1.2,8,1.3,9); res->AddBin(1.2,9,1.3,10); res->AddBin(1.3,8,1.45,9); res->AddBin(1.3,9,1.45,10); res->AddBin(1.45,8,1.6,9); res->AddBin(1.45,9,1.6,10);
  for(int i=0;i<=3;i++) {
    res->AddBin( 1.6+i*0.1 ,9, 1.7+i*0.1 ,10);
    res->AddBin( 1.6+i*0.1 ,8, 1.7+i*0.1 ,9);
  }
  res->AddBin(2.,8,2.15,9); res->AddBin(2.,9,2.15,10); res->AddBin(2.15,8,2.3,9); res->AddBin(2.15,9,2.3,10);
  //pt 10-11
  res->AddBin(0,10,0.5,11.); res->AddBin(0.5,10,0.8,11.); res->AddBin(0.8,10,1.,11.); res->AddBin(1.,10,1.2,11.); res->AddBin(1.2,10,1.3,11.); res->AddBin(1.3,10,1.45,11.); res->AddBin(1.45,10,1.6,11.);
  for(int i=0;i<=3;i++) 
    res->AddBin( 1.6+i*0.1 ,10, 1.7+i*0.1 ,11);
  res->AddBin(2.,10,2.15,11.); res->AddBin(2.15,10,2.3,11.);

  //pt 11-13
  res->AddBin(0,11,0.5,13.); res->AddBin(0.5,11,0.7,13.); res->AddBin(0.7,11,0.95,13.);
  for(int i=0;i<=8;i++) 
    res->AddBin( 0.95+i*0.15 ,11, 1.1+i*0.15 ,13);
  //pt 13-15
  res->AddBin(0,13,0.5,15.); res->AddBin(0.5,13,0.7,15.); res->AddBin(0.7,13,0.95,15.);
  for(int i=0;i<=8;i++)
    res->AddBin( 0.95+i*0.15 ,13, 1.1+i*0.15 ,15);
  
  //pt 15-18
  res->AddBin(0,15,0.5,18.); res->AddBin(0.5,15,0.9,18.); 
  for(int i=0;i<=3;i++){
    res->AddBin( 0.9+i*0.3 ,15, 1.2+i*0.3 ,16.5);
    res->AddBin( 0.9+i*0.3 ,16.5, 1.2+i*0.3 ,18);
  }
  res->AddBin(2.1,15,2.3,18.);
  //pt 18-22
  res->AddBin(0,18,0.5,20.);   res->AddBin(0,20,0.5,22.); 
  for(int i=0;i<=5;i++) 
    res->AddBin(0.5+i*0.3,18,0.8+i*0.3,22.); 
  //pt 22-30
  res->AddBin(0,22,0.5,30.); 
  for(int i=0;i<=5;i++) res->AddBin(0.5+i*0.3,22,0.8+i*0.3,30.); 
  //pt 30-50
  res->AddBin(0,30,0.3,50.); 
  for(int i=0;i<=4;i++) res->AddBin(0.3+i*0.4,30,0.7+i*0.4,50.); 

  //y 2.3-2.4
  res->AddBin(2.3,0,2.4,4); res->AddBin(2.3,4,2.4,7); res->AddBin(2.3,7,2.4,10); res->AddBin(2.3,10,2.4,14); res->AddBin(2.3,14,2.4,22); res->AddBin(2.3,22,2.4,50);

  return res;
}
