//This code takes as input two histograms, h1 and h2, containing sidebands of a signal, a title for the pdf you will save plots in, an option for geometric or arithmetic averaging in the parameters averages, arrays of initial parameters to begin fits with, and the expected signal region, xlow to xhigh.  It refurns a TF1 object containing the final fit function for the signal region and a pdf with several plots which compare the individual sidebands to the averaged parameter fit and the averaged parameter fit to the sum of the sidebands.

//Author: Natalie Blot, February 2020

#include <iostream>
#include <cmath>
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"


TF1* SideBandFitter(TH1 * h1,TH1 * h2, std::string title, char opt,float xlow,float xhigh, std::vector<double> guess1, std::vector<double> guess2, bool crystalball) {
  //Define the function we use to fit the histograms.  In this case, the product of a error function and an exponential decay.  The function has 4 parameters,
  //[0] = normalization
  // [1] = Dispalcement of the error function
  // [2] = steepness of the error function
  // [3] = exponential decay parameter
  // [4] = displacement of the exponential.  By default, this is fixed to 0 but making it free may improve your fit.  You simply have to do some commenting/uncommenting in this code

auto h_test = new TH1F();
h_test->SetDefaultSumw2(true);
//Function starts
float low = h1->GetXaxis()->GetXmin();
float high = h1->GetXaxis()->GetXmax();

 std::cout << "Before entering any crystalball if block" << std::endl;
 if (!crystalball) {
   cout << "Entering the not crystal ball if block" << endl;
   float xred = xlow;
   //float xred = xlow+0.1;
   TF1 * leftfit = new TF1("leftfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);

   TF1 * rightfit = new TF1("rightfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xred,xhigh);

   //Initialize fits to guess1
   leftfit->SetParameters(guess1[0],guess1[1],guess1[2],guess1[3],guess1[4]);

   //Fix parameters to low BDT cut values to try and mitigate drop for tighter BDT cuts
   // leftfit->FixParameter(2,0.61);
   // leftfit->FixParameter(3,1.16);

   leftfit->FixParameter(4,0);
   //Now fit h1
   h1->Fit(leftfit,"","",xlow,xhigh);

   //initialize to guess2
   rightfit->SetParameters(guess2[0],guess2[1],guess2[2],guess2[3],guess2[4]);
   //Fix exp shift to the same value as leftfit
   rightfit->FixParameter(4,leftfit->GetParameter(4));

   //Fix x intercept to same value as leftfit
   //rightfit->FixParameter(1,leftfit->GetParameter(1));

   //Fix parameters to low BDT cut values to try and mitigate drop for tighter BDT cuts
   // rightfit->FixParameter(2,1.69);
   // rightfit->FixParameter(3,1.16);
   //Fit h2
   h2->Fit(rightfit,"","",xred,xhigh);

   //Average the parameters
   //Weight average by number of data points
   double leftN = h1->Integral(h1->FindBin(low),h1->FindBin(high));
   double rightN = h2->Integral(h2->FindBin(low),h2->FindBin(high));
   double totalN = leftN + rightN;

   double avgp0,avgp1,avgp2,avgp3,avgp4,w1,w2;

   if (opt=='A'){
     avgp0 = (leftN*leftfit->GetParameter(0) + rightN*rightfit->GetParameter(0))/totalN;
     avgp1 = (leftN*leftfit->GetParameter(1) + rightN*rightfit->GetParameter(1))/totalN;
     avgp2 = (leftN*leftfit->GetParameter(2) + rightN*rightfit->GetParameter(2))/totalN;
     avgp3 = (leftN*leftfit->GetParameter(3) + rightN*rightfit->GetParameter(3))/totalN;
     avgp4 = (leftN*leftfit->GetParameter(4) + rightN*rightfit->GetParameter(4))/totalN;}

   else {
     w1 = leftN/totalN;
     w2 = rightN/totalN;
     avgp0 = exp(w1*log(leftfit->GetParameter(0))+w2*log(rightfit->GetParameter(0)));
     avgp1 = exp(w1*log(leftfit->GetParameter(1))+w2*log(rightfit->GetParameter(1)));
     avgp2 = exp(w1*log(leftfit->GetParameter(2))+w2*log(rightfit->GetParameter(2)));
     avgp3 = exp(w1*log(leftfit->GetParameter(3))+w2*log(rightfit->GetParameter(3)));
     avgp4 = exp(w1*log(leftfit->GetParameter(4))+w2*log(rightfit->GetParameter(4)));}

   TF1 * signalfit = new TF1("signalfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);
   signalfit->SetParameters(avgp0,avgp1,avgp2,avgp3,avgp4);

   //normalise the fits
   double leftint = leftfit->Integral(xlow,xhigh);
   double rightint = rightfit->Integral(xlow,xhigh);
   double sigint = signalfit->Integral(xlow,xhigh);


   TF1 * normleftfit = new TF1("normleftfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);
   normleftfit->SetParameters(leftfit->GetParameter(0)/leftint,leftfit->GetParameter(1), leftfit->GetParameter(2),leftfit->GetParameter(3),leftfit->GetParameter(4));
   normleftfit->SetTitle("Normalized fits-Crystal Ball");

   TF1 * normrightfit = new TF1("normrightfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);
   normrightfit->SetParameters(rightfit->GetParameter(0)/rightint,rightfit->GetParameter(1), rightfit->GetParameter(2),rightfit->GetParameter(3),rightfit->GetParameter(4));
   TF1 * normsigfit =  new TF1("normsigfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);
   normsigfit->SetParameters(signalfit->GetParameter(0)/sigint,signalfit->GetParameter(1), signalfit->GetParameter(2),signalfit->GetParameter(3),signalfit->GetParameter(4));

   //Sum of histograms and fits, also fit the sum independently to compare them all.  sumfit is the fit of hsum, fitsum is the sum rightfit and leftfit
   TH1 * hsum = new TH1F("hsum","hsum",h1->GetNbinsX(),low,high);
   hsum->Add(h1);
   hsum->Add(h2);
   //Make new function to fit hsum and then fit it
   //TF1 * sumfit = new TF1("sumfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);
   //Guess
   //sumfit->SetParameters(250,0,2,4,0);
   //Fix exponential shift
   // sumfit->FixParameter(4,leftfit->GetParameter(4));
   //Fix the x intercept to be the same as rightfit
   //sumfit->FixParameter(1,rightfit->GetParameter(1));
   //hsum->Fit(sumfit,"","",low,high);

   //Renormalization constants
   double Nleft = h1->Integral(h1->FindBin(xlow),h1->FindBin(xhigh));
   double Nright = h2->Integral(h2->FindBin(xlow),h2->FindBin(xhigh));

   //Sum leftfit and rightfit
   TF1 * fitsum = new TF1("fitsum","[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3]) + [5]*TMath::Erf((x-[6])/[7])*exp(-(x-[8])/[9])",xlow,xhigh);
   fitsum->SetParameter(0,normleftfit->GetParameter(0)*Nleft*h1->GetBinWidth(3));
   fitsum->SetParameter(1,leftfit->GetParameter(1));
   fitsum->SetParameter(2,leftfit->GetParameter(2));
   fitsum->SetParameter(3,leftfit->GetParameter(3));
   fitsum->SetParameter(4,leftfit->GetParameter(4));
   fitsum->SetParameter(5,normrightfit->GetParameter(0)*Nright*h2->GetBinWidth(3));
   fitsum->SetParameter(6,rightfit->GetParameter(1));
   fitsum->SetParameter(7,rightfit->GetParameter(2));
   fitsum->SetParameter(8,rightfit->GetParameter(4));
   fitsum->SetParameter(9,rightfit->GetParameter(3));

   //Renormalize the signal fit to the data
   cout << "Nleft+Nright=" << Nleft+Nright <<endl;
   TF1 * datanormsigfit = new TF1("datanormsigfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);
   datanormsigfit->SetParameters(normsigfit->GetParameter(0)*(Nleft+Nright)*h1->GetBinWidth(3),normsigfit->GetParameter(1),normsigfit->GetParameter(2),normsigfit->GetParameter(3),normsigfit->GetParameter(4));

   //Aesthetics- set all the colors
   signalfit->SetLineColor(1);
   normsigfit->SetLineColor(1);
   h1->SetMarkerColor(3);
   h1->SetMarkerSize(1.5);
   h2->SetMarkerSize(1.5);
   h2->SetMarkerColor(2);
   leftfit->SetLineColor(3);
   normleftfit->SetLineColor(3);
   rightfit->SetLineColor(2);
   normrightfit->SetLineColor(2);
   fitsum->SetLineColor(kGreen+3);
   //sumfit->SetLineColor(4);
   h1->SetMarkerStyle(23);
   h2->SetMarkerStyle(8);
   datanormsigfit->SetLineColor(1);
   hsum->SetMarkerColor(4);
   hsum->SetMarkerStyle(8);

   //Make legend for pads 1 and 2
   TLegend * legend = new TLegend(0.6,0.4,1,0.7);
   legend->AddEntry(h1, "Left Sideband", "lp");
   legend->AddEntry(h2,"Right Sideband", "lp");
   legend->AddEntry(signalfit, "Signal Region Fit","l");

   //Make a legend for pad3
   TLegend * legend3 =new  TLegend(0.6,0.4,1,0.7);
   legend3->AddEntry(fitsum,"Sum of left and right fits", "l");
   //legend3->AddEntry(sumfit,"Fit of the summed histogram", "lp");
   legend3->AddEntry(datanormsigfit,"Average parameter fit", "l");

   gStyle->SetOptStat(0);

   //Draw and save things
   TCanvas * c1 = new TCanvas("c1","c1",1500,650);
   //c1->Divide(2,2);
   c1->cd();
   h1->SetTitle("Unnormalized fits and data");
   h1->Draw("P");
   h2->Draw("SAMEP");
   signalfit->GetYaxis()->SetRangeUser(0,signalfit->GetMaximum()*1.1);
   signalfit->Draw("SAME");
   legend->Draw("SAME");
   leftfit->Draw("SAME");
   rightfit->Draw("SAME");

   TCanvas * c2 = new TCanvas("c2","c2",1500,650);
   c2->cd();
   normleftfit->GetYaxis()->SetRangeUser(0,normsigfit->GetMaximum()*1.3);
   normleftfit->Draw();
   normrightfit->Draw("SAME");
   normsigfit->Draw("SAME");
   legend->Draw("SAME");

   TCanvas * c3 = new TCanvas("c3","c3",1500,650);
   c3->cd();
   hsum->SetTitle("Comparison between new and old methods");
   fitsum->GetYaxis()->SetRangeUser(0,fitsum->GetMaximum()*1.1);
   hsum->Draw("P");
   //sumfit->Draw("SAME");
   fitsum->Draw("SAME");
   datanormsigfit->Draw("SAME");
   legend3->Draw("SAME");

   //Create the titles we want for each plot
   std::string title1 = title + "Unnormalized.pdf";
   const char * titlepointer1=title1.c_str();

   std::string title2 = title + "Normalized.pdf";
   const char * titlepointer2=title2.c_str();

   std::string title3 = title + "Comparison.pdf";
   const char * titlepointer3=title3.c_str();

   c1->Modified();
   c1->Update();
   //c1->SaveAs(titlepointer1);

   c2->Modified();
   c2->Update();
   //c2->SaveAs(titlepointer2);

   c3->Modified();
   c3->Update();
   //c3->SaveAs(titlepointer3);

   cout << "Left fit chisqr/NDF " << leftfit->GetChisquare()/leftfit->GetNDF()<< endl <<"Probability of left sideband fit " << leftfit->GetProb() << endl << "right fit chisqr/NDF " << rightfit->GetChisquare()/rightfit->GetNDF()<< endl << "Probability of right sideband fit " << rightfit->GetProb() << endl;

   return normsigfit;
 }//end if !crystalball

 //***************CRYSTALBALL VERSION****************************************************************************************
 else {
   cout << "#########In the crystal ball block#############" << endl;
   float xred = xlow;

   //parameters of the flipped crystal ball (gaussian to the left, power law tail to the right)
   //[0] = normalization, N
   //[1] = alpha, crystal ball cut off
   //[2] = mean of gaussian core
   //[3]= std dev of gaussian core
   //[4] = exponential of the power law
   
   //This definition finallly works TF1 * testcb = new TF1("testcb","[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4])",-2,2);

   TF1 * leftfit = new TF1("leftfit", "[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4])",xlow,xhigh);

   TF1 * rightfit = new TF1("rightfit", "[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4])",xlow,xhigh);

   //Initialize fits to guess1
   leftfit->SetParameters(guess1[0],guess1[1],guess1[2],guess1[3],guess1[4]);

   //Now fit h1
   cout << "Fitting left sideband" << endl;
   h1->Fit(leftfit,"","",xlow,xhigh);

   //initialize to guess2
   rightfit->SetParameters(guess2[0],guess2[1],guess2[2],guess2[3],guess2[4]);
   // rightfit->FixParameter(0,guess2[0]);
   // rightfit->FixParameter(1,guess2[1]);
   // rightfit->FixParameter(2,guess2[2]);
   // rightfit->FixParameter(3,guess2[3]);
   // rightfit->FixParameter(4,guess2[4]);
   //Fit h2
   cout << "Fitting right sideband" << endl;
   h2->Fit(rightfit,"","",xred,xhigh);

   //Average the parameters
   //Weight average by number of data points
   double leftN = h1->Integral(h1->FindBin(low),h1->FindBin(high));
   double rightN = h2->Integral(h2->FindBin(low),h2->FindBin(high));
   double totalN = leftN + rightN;

   double avgp0,avgp1,avgp2,avgp3,avgp4,w1,w2;

   if (opt=='A'){
     avgp0 = (leftN*leftfit->GetParameter(0) + rightN*rightfit->GetParameter(0))/totalN;
     avgp1 = (leftN*leftfit->GetParameter(1) + rightN*rightfit->GetParameter(1))/totalN;
     avgp2 = (leftN*leftfit->GetParameter(2) + rightN*rightfit->GetParameter(2))/totalN;
     avgp3 = (leftN*leftfit->GetParameter(3) + rightN*rightfit->GetParameter(3))/totalN;
     avgp4 = (leftN*leftfit->GetParameter(4) + rightN*rightfit->GetParameter(4))/totalN;}

   else {
     w1 = leftN/totalN;
     w2 = rightN/totalN;
     avgp0 = exp(w1*log(leftfit->GetParameter(0))+w2*log(rightfit->GetParameter(0)));
     avgp1 = exp(w1*log(leftfit->GetParameter(1))+w2*log(rightfit->GetParameter(1)));
     avgp2 = exp(w1*log(leftfit->GetParameter(2))+w2*log(rightfit->GetParameter(2)));
     avgp3 = exp(w1*log(leftfit->GetParameter(3))+w2*log(rightfit->GetParameter(3)));
     avgp4 = exp(w1*log(leftfit->GetParameter(4))+w2*log(rightfit->GetParameter(4)));}

   TF1 * signalfit = new TF1("signalfit", "[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4])",xlow,xhigh);
   signalfit->SetParameters(avgp0,avgp1,avgp2,avgp3,avgp4);

   //normalise the fits
   double leftint = leftfit->Integral(xlow,xhigh);
   double rightint = rightfit->Integral(xlow,xhigh);
   double sigint = signalfit->Integral(xlow,xhigh);

   TF1 * normleftfit = new TF1("normleftfit", "[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4])",xlow,xhigh);
   normleftfit->SetParameters(leftfit->GetParameter(0)/leftint,leftfit->GetParameter(1), leftfit->GetParameter(2),leftfit->GetParameter(3),leftfit->GetParameter(4));
   normleftfit->SetTitle("Normalized fits");

   TF1 * normrightfit = new TF1("normrightfit", "[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4])",xlow,xhigh);
   normrightfit->SetParameters(rightfit->GetParameter(0)/rightint,rightfit->GetParameter(1), rightfit->GetParameter(2),rightfit->GetParameter(3),rightfit->GetParameter(4));

   TF1 * normsigfit =  new TF1("normsigfit", "[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4])",xlow,xhigh);
   normsigfit->SetParameters(signalfit->GetParameter(0)/sigint,signalfit->GetParameter(1), signalfit->GetParameter(2),signalfit->GetParameter(3),signalfit->GetParameter(4));

   //Sum of histograms and fits, also fit the sum independently to compare them all.  sumfit is the fit of hsum, fitsum is the sum rightfit and leftfit
   TH1 * hsum = new TH1F("hsum","hsum",h1->GetNbinsX(),low,high);
   hsum->Add(h1);
   hsum->Add(h2);
   //Make new function to fit hsum and then fit it
   //TF1 * sumfit = new TF1("sumfit", "[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])",xlow,xhigh);
   //Guess
   //sumfit->SetParameters(250,0,2,4,0);
   //Fix exponential shift
   // sumfit->FixParameter(4,leftfit->GetParameter(4));
   //Fix the x intercept to be the same as rightfit
   //sumfit->FixParameter(1,rightfit->GetParameter(1));
   //hsum->Fit(sumfit,"","",low,high);

   //Renormalization constants
   double Nleft = h1->Integral(h1->FindBin(xlow),h1->FindBin(xhigh));
   double Nright = h2->Integral(h2->FindBin(xlow),h2->FindBin(xhigh));

   //Sum leftfit and rightfit
   TF1 * fitsum = new TF1("fitsum","[0]*(((-x+[2])/[3])>-[1])*exp(-pow((-x+[2]),2)/(2*pow([3],2))) + [0]*(((-x+[2])/[3])<=-[1])*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*pow(([4]/abs([1])-abs([1]) - (-x+[2])/[3]),-[4]) + [5]*(((-x+[7])/[8])>-[6])*exp(-pow((-x+[7]),2)/(2*pow([8],2))) + [5]*(((-x+[7])/[8])<=-[6])*pow(([9]/[6]),[9])*exp(-pow([6],2)/2)*pow(([9]/abs([6])-abs([6]) - (-x+[7])/[8]),-[9])",xlow,xhigh);

   fitsum->SetParameter(0,normleftfit->GetParameter(0)*Nleft*h1->GetBinWidth(3));
   fitsum->SetParameter(1,leftfit->GetParameter(1));
   fitsum->SetParameter(2,leftfit->GetParameter(2));
   fitsum->SetParameter(3,leftfit->GetParameter(3));
   fitsum->SetParameter(4,leftfit->GetParameter(4));
   fitsum->SetParameter(5,normrightfit->GetParameter(0)*Nright*h2->GetBinWidth(3));
   fitsum->SetParameter(6,rightfit->GetParameter(1));
   fitsum->SetParameter(7,rightfit->GetParameter(2));
   fitsum->SetParameter(8,rightfit->GetParameter(3));
   fitsum->SetParameter(9,rightfit->GetParameter(4));

   //Renormalize the signal fit to the data
   cout << "Nleft+Nright=" << Nleft+Nright <<endl;
   TF1 * datanormsigfit = new TF1("datanormsigfit", "[0]*( (((x-[2])/[3])<-[1])*exp(-pow((x-[2]),2)/(2*pow([3],2))) + (((x-[2])/[3])>-[1]))*pow(([4]/[1]),[4])*exp(-pow([1],2)/2)*(pow([4]/abs([1])-abs([1]) + (x-[2])/[3],-[4]))",xlow,xhigh);
   datanormsigfit->SetParameters(normsigfit->GetParameter(0)*(Nleft+Nright)*h1->GetBinWidth(3),normsigfit->GetParameter(1),normsigfit->GetParameter(2),normsigfit->GetParameter(3),normsigfit->GetParameter(4));


   //Aesthetics- set all the colors
   signalfit->SetLineColor(1);
   normsigfit->SetLineColor(1);
   h1->SetMarkerColor(3);
   h1->SetMarkerSize(1.5);
   h2->SetMarkerSize(1.5);
   h2->SetMarkerColor(2);
   leftfit->SetLineColor(3);
   normleftfit->SetLineColor(3);
   rightfit->SetLineColor(2);
   normrightfit->SetLineColor(2);
   fitsum->SetLineColor(kGreen+3);
   //sumfit->SetLineColor(4);
   h1->SetMarkerStyle(23);
   h2->SetMarkerStyle(8);
   datanormsigfit->SetLineColor(1);
   hsum->SetMarkerColor(4);
   hsum->SetMarkerStyle(8);

   //Make legend for pads 1 and 2
   TLegend * legend = new TLegend(0.6,0.4,1,0.7);
   legend->AddEntry(h1, "Left Sideband", "lp");
   legend->AddEntry(h2,"Right Sideband", "lp");
   legend->AddEntry(signalfit, "Signal Region Fit","l");

   //Make a legend for pad3
   TLegend * legend3 =new  TLegend(0.6,0.4,1,0.7);
   legend3->AddEntry(fitsum,"Sum of left and right fits", "l");
   //legend3->AddEntry(sumfit,"Fit of the summed histogram", "lp");
   legend3->AddEntry(datanormsigfit,"Average parameter fit", "l");

   gStyle->SetOptStat(0);

   //Draw and save things
   TCanvas * c1 = new TCanvas("c1","c1",1500,650);
   //c1->Divide(2,2);
   c1->cd();
   h1->SetTitle("Unnormalized fits and data- Crystal Ball");
   h1->Draw("P");
   h2->Draw("SAMEP");
   signalfit->GetYaxis()->SetRangeUser(0,signalfit->GetMaximum()*1.1);
   signalfit->Draw("SAME");
   legend->Draw("SAME");
   leftfit->Draw("SAME");
   rightfit->Draw("SAME");

   TCanvas * c2 = new TCanvas("c2","c2",1500,650);
   c2->cd();
   normleftfit->GetYaxis()->SetRangeUser(0,normsigfit->GetMaximum()*1.3);
   normleftfit->Draw();
   normrightfit->Draw("SAME");
   normsigfit->Draw("SAME");
   legend->Draw("SAME");

   TCanvas * c3 = new TCanvas("c3","c3",1500,650);
   c3->cd();
   hsum->SetTitle("Comparison between new and old method- Crystal Ball");
   fitsum->GetYaxis()->SetRangeUser(0,fitsum->GetMaximum()*1.1);
   hsum->Draw("P");
   //sumfit->Draw("SAME");
   fitsum->Draw("SAME");
   datanormsigfit->Draw("SAME");
   legend3->Draw("SAME");

   //Create the titles we want for each plot
   std::string title1 = title + "CB_Unnormalized.pdf";
   const char * titlepointer1=title1.c_str();

   std::string title2 = title + "CB_Normalized.pdf";
   const char * titlepointer2=title2.c_str();

   std::string title3 = title + "CB_Comparison.pdf";
   const char * titlepointer3=title3.c_str();

   c1->Modified();
   c1->Update();
   c1->SaveAs(titlepointer1);

   c2->Modified();
   c2->Update();
   c2->SaveAs(titlepointer2);

   c3->Modified();
   c3->Update();
   c3->SaveAs(titlepointer3);

   cout << "Left fit chisqr/NDF " << leftfit->GetChisquare()/leftfit->GetNDF()<< endl <<"Probability of left sideband fit " << leftfit->GetProb() << endl << "right fit chisqr/NDF " << rightfit->GetChisquare()/rightfit->GetNDF()<< endl << "Probability of right sideband fit " << rightfit->GetProb() << endl;

   return normsigfit;


   }//end else

}
