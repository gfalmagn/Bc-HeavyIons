#include <TF1.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TLatex.h>

// Yen-Jie: systematics table for D meson
// Unit: In percentage

// =============================================================================================================
const int nPtBins = 14;
double PtBins[nPtBins+1] = { 2, 3, 4, 5, 6, 8, 10, 12.5,15.0, 20, 25, 30, 40, 60, 100};

// =============================================================================================================
// D meson decay
// =============================================================================================================
//double DKpiBRUncertainty	= 1.29;			// from PDG 2012
double DKpiBRUncertainty	= 1.02;			// from PDG 2016

// =============================================================================================================
// pp uncertainty
// =============================================================================================================

// Normalization uncertainty
double ppLumiUncertainty 	= 2.3;			// 5% for the moment, to be updated (4/7/2016)

// Point-to-point
double ppTrackingEfficiency 	        = 4;   		// single track systematics from D* studies //4
double PbPbTrackingEfficiency0100 	= 6;   		// single track systematics from D* studies //5
double PbPbTrackingEfficiency010 	= 6.5;   	// single track systematics from D* studies //5
TH1D*  ppSignalExtraction;				// 
TH1D*  ppMesonSelection;				// 
TH1D*  ppBFeedDownCorrection; 				// 
TH1D*  ppTriggerEfficiency; 				// 

TF1* fPPPtShape = new TF1("fPPPtShapeSig","[0]+[1]/(x)+[2]/x/x+[3]*x");


// =============================================================================================================
// PbPb uncertainty
// =============================================================================================================

// Normalization uncertainty
double PbPbNMBUncertainty	= 2;			// uncertainty associated with minbias events,
							// used in RAA for pT < 20 GeV (4/7/2016)
double TAAUncertainty0to100HI	= 2.8;			// Updated number (4/7/2016)
double TAAUncertainty0to100LO	= 3.4;			// Updated number (4/7/2016)
double TAAUncertainty0to10HI	= 1.9;			// Updated number (4/7/2016)
double TAAUncertainty0to10LO	= 3.0;			// Updated number (4/7/2016)
//double PbPbTrigger		= 2.0;			// Statistical uncertainty of the zero-coefficient of the linear fit
double PbPbLumiUncertainty	= 2;			// 10% for the moment, to be updated (from Daniel), NOT used

// Point-to-point

TH1D*  PbPbSignalExtraction;				// (4/7/2016)
TH1D*  PbPbMesonSelection;				// (4/7/2016)
TH1D*  PbPbBFeedDownCorrection;
TH1D*  PbPbTriggerEfficiency;

TF1 *fPbPbPtShape = new TF1("fPbPbPtShapeSig","[0]+[1]/(x)+[2]/x/x+[3]*x");

TH1D*  RAABFeedDownCorrection;				// (4/7/2016)


bool initialized = 0;

void initializationPP()
{
  ppBFeedDownCorrection->SetBinContent(1,	10);
  ppBFeedDownCorrection->SetBinContent(2,	10);
  ppBFeedDownCorrection->SetBinContent(3,	10);
  ppBFeedDownCorrection->SetBinContent(4,	10);
  ppBFeedDownCorrection->SetBinContent(5,	10);
  ppBFeedDownCorrection->SetBinContent(6,	10);
  ppBFeedDownCorrection->SetBinContent(7,	10);
  ppBFeedDownCorrection->SetBinContent(8,	10);
  ppBFeedDownCorrection->SetBinContent(9,	10);
  ppBFeedDownCorrection->SetBinContent(10,	10);
  ppBFeedDownCorrection->SetBinContent(11,	10);
  ppBFeedDownCorrection->SetBinContent(12,	10);
  ppBFeedDownCorrection->SetBinContent(13,	10);
  ppBFeedDownCorrection->SetBinContent(14,	10);
  
  ppMesonSelection->SetBinContent(1,		3.6);
  ppMesonSelection->SetBinContent(2,		3.6);
  ppMesonSelection->SetBinContent(3,		3.6);
  ppMesonSelection->SetBinContent(4,		3.6);
  ppMesonSelection->SetBinContent(5,		3.6);
  ppMesonSelection->SetBinContent(6,		3.6);
  ppMesonSelection->SetBinContent(7,		3.6);
  ppMesonSelection->SetBinContent(8,		3.6);
  ppMesonSelection->SetBinContent(9,		3.6);
  ppMesonSelection->SetBinContent(10,		0.5);
  ppMesonSelection->SetBinContent(11,		0.5);
  ppMesonSelection->SetBinContent(12,		0.5);
  ppMesonSelection->SetBinContent(13,		0.5);
  ppMesonSelection->SetBinContent(14,		0.5);
  
  ppSignalExtraction->SetBinContent(1,		8.2);
  ppSignalExtraction->SetBinContent(2,		7.6);
  ppSignalExtraction->SetBinContent(3,		3.5);
  ppSignalExtraction->SetBinContent(4,		3.0);
  ppSignalExtraction->SetBinContent(5,		3.0);
  ppSignalExtraction->SetBinContent(6,		1.7);
  ppSignalExtraction->SetBinContent(7,		2.1);
  ppSignalExtraction->SetBinContent(8,		2.1);
  ppSignalExtraction->SetBinContent(9,		4.0);
  ppSignalExtraction->SetBinContent(10,		2.0);
  ppSignalExtraction->SetBinContent(11,	        1.6);
  ppSignalExtraction->SetBinContent(12,	        2.3);
  ppSignalExtraction->SetBinContent(13,	        2.8);
  ppSignalExtraction->SetBinContent(14,	        5.2); 

  ppTriggerEfficiency->SetBinContent(1,		0.0);
  ppTriggerEfficiency->SetBinContent(2,		0.0);
  ppTriggerEfficiency->SetBinContent(3,		0.0);
  ppTriggerEfficiency->SetBinContent(4,		0.0);
  ppTriggerEfficiency->SetBinContent(5,		0.0);
  ppTriggerEfficiency->SetBinContent(6,		0.0);
  ppTriggerEfficiency->SetBinContent(7,		0.0);
  ppTriggerEfficiency->SetBinContent(8,		0.0);
  ppTriggerEfficiency->SetBinContent(9,		0.0);
  ppTriggerEfficiency->SetBinContent(10,	1.0);
  ppTriggerEfficiency->SetBinContent(11,	1.0);
  ppTriggerEfficiency->SetBinContent(12,	1.0);
  ppTriggerEfficiency->SetBinContent(13,	1.0);
  ppTriggerEfficiency->SetBinContent(14,	1.0);
  
  fPPPtShape->SetParameters(0.999265,-0.0458006,-0.181359,0);
}

void initializationPbPbCent0100()
{
  PbPbBFeedDownCorrection->SetBinContent(1,	10);
  PbPbBFeedDownCorrection->SetBinContent(2,	10);
  PbPbBFeedDownCorrection->SetBinContent(3,	10);
  PbPbBFeedDownCorrection->SetBinContent(4,	10);
  PbPbBFeedDownCorrection->SetBinContent(5,	10);
  PbPbBFeedDownCorrection->SetBinContent(6,	10);
  PbPbBFeedDownCorrection->SetBinContent(7,	10);
  PbPbBFeedDownCorrection->SetBinContent(8,	10);
  PbPbBFeedDownCorrection->SetBinContent(9,	10);
  PbPbBFeedDownCorrection->SetBinContent(10,	10);
  PbPbBFeedDownCorrection->SetBinContent(11,	10);
  PbPbBFeedDownCorrection->SetBinContent(12,	10);
  PbPbBFeedDownCorrection->SetBinContent(13,	10);
  PbPbBFeedDownCorrection->SetBinContent(14,	10);
  
  RAABFeedDownCorrection->SetBinContent(1,	10);
  RAABFeedDownCorrection->SetBinContent(2,	10);
  RAABFeedDownCorrection->SetBinContent(3,	10);
  RAABFeedDownCorrection->SetBinContent(4,	10);
  RAABFeedDownCorrection->SetBinContent(5,	10);
  RAABFeedDownCorrection->SetBinContent(6,	10);
  RAABFeedDownCorrection->SetBinContent(7,	10);
  RAABFeedDownCorrection->SetBinContent(8,	10);
  RAABFeedDownCorrection->SetBinContent(9,	10);
  RAABFeedDownCorrection->SetBinContent(10,	10);
  RAABFeedDownCorrection->SetBinContent(11,	10);
  RAABFeedDownCorrection->SetBinContent(12,	10);
  RAABFeedDownCorrection->SetBinContent(13,	10);
  RAABFeedDownCorrection->SetBinContent(14,	10);
  
  PbPbMesonSelection->SetBinContent(1,		3.5);
  PbPbMesonSelection->SetBinContent(2,		3.5);
  PbPbMesonSelection->SetBinContent(3,		3.5);
  PbPbMesonSelection->SetBinContent(4,		3.5);
  PbPbMesonSelection->SetBinContent(5,		3.5);
  PbPbMesonSelection->SetBinContent(6,		3.5);
  PbPbMesonSelection->SetBinContent(7,		3.5);
  PbPbMesonSelection->SetBinContent(8,		3.5);
  PbPbMesonSelection->SetBinContent(9,		3.5);
  PbPbMesonSelection->SetBinContent(10,		2.7);
  PbPbMesonSelection->SetBinContent(11,		2.7);
  PbPbMesonSelection->SetBinContent(12,		2.7);
  PbPbMesonSelection->SetBinContent(13,		2.7);
  PbPbMesonSelection->SetBinContent(14,		2.7);
  
  PbPbSignalExtraction->SetBinContent(1,	4.8);
  PbPbSignalExtraction->SetBinContent(2,	2.0);
  PbPbSignalExtraction->SetBinContent(3,	3.0);
  PbPbSignalExtraction->SetBinContent(4,	2.3);
  PbPbSignalExtraction->SetBinContent(5,	1.7);
  PbPbSignalExtraction->SetBinContent(6,	1.7);
  PbPbSignalExtraction->SetBinContent(7,	1.3);
  PbPbSignalExtraction->SetBinContent(8,	1.3);
  PbPbSignalExtraction->SetBinContent(9,	6.5);
  PbPbSignalExtraction->SetBinContent(10,	7.1);
  PbPbSignalExtraction->SetBinContent(11,	9.4);
  PbPbSignalExtraction->SetBinContent(12,	7.5);
  PbPbSignalExtraction->SetBinContent(13,	6.0); //was 3.3% before smoothing
  PbPbSignalExtraction->SetBinContent(14,	7.0);

  PbPbTriggerEfficiency->SetBinContent(1,	0.0);
  PbPbTriggerEfficiency->SetBinContent(2,	0.0);
  PbPbTriggerEfficiency->SetBinContent(3,	0.0);
  PbPbTriggerEfficiency->SetBinContent(4,	0.0);
  PbPbTriggerEfficiency->SetBinContent(5,	0.0);
  PbPbTriggerEfficiency->SetBinContent(6,	0.0);
  PbPbTriggerEfficiency->SetBinContent(7,	0.0);
  PbPbTriggerEfficiency->SetBinContent(8,	0.0);
  PbPbTriggerEfficiency->SetBinContent(9,	0.0);
  PbPbTriggerEfficiency->SetBinContent(10,	2.0);
  PbPbTriggerEfficiency->SetBinContent(11,	2.0);
  PbPbTriggerEfficiency->SetBinContent(12,	2.0);
  PbPbTriggerEfficiency->SetBinContent(13,	2.0);
  PbPbTriggerEfficiency->SetBinContent(14,	2.0);
  
  fPbPbPtShape->SetParameters(0.984161,0.0593406,-0.3992,0.000271564);
}

void initializationPbPbCent010()
{
  PbPbBFeedDownCorrection->SetBinContent(1,	10.);
  PbPbBFeedDownCorrection->SetBinContent(2,	10.);
  PbPbBFeedDownCorrection->SetBinContent(3,	10.);
  PbPbBFeedDownCorrection->SetBinContent(4,	10.);
  PbPbBFeedDownCorrection->SetBinContent(5,	10.);
  PbPbBFeedDownCorrection->SetBinContent(6,	10.);
  PbPbBFeedDownCorrection->SetBinContent(7,	10.);
  PbPbBFeedDownCorrection->SetBinContent(8,	10.);
  PbPbBFeedDownCorrection->SetBinContent(9,	10.);
  PbPbBFeedDownCorrection->SetBinContent(10,	10.);
  PbPbBFeedDownCorrection->SetBinContent(11,	10.);
  PbPbBFeedDownCorrection->SetBinContent(12,	10.);
  PbPbBFeedDownCorrection->SetBinContent(13,	10.);
  PbPbBFeedDownCorrection->SetBinContent(14,	10.);

  RAABFeedDownCorrection->SetBinContent(1,	10.);
  RAABFeedDownCorrection->SetBinContent(2,	10.);
  RAABFeedDownCorrection->SetBinContent(3,	10.);
  RAABFeedDownCorrection->SetBinContent(4,	10.);
  RAABFeedDownCorrection->SetBinContent(5,	10.);
  RAABFeedDownCorrection->SetBinContent(6,	10.);
  RAABFeedDownCorrection->SetBinContent(7,	10.);
  RAABFeedDownCorrection->SetBinContent(8,	10.);
  RAABFeedDownCorrection->SetBinContent(9,	10.);
  RAABFeedDownCorrection->SetBinContent(10,	10.);
  RAABFeedDownCorrection->SetBinContent(11,	10.);
  RAABFeedDownCorrection->SetBinContent(12,	10.);
  RAABFeedDownCorrection->SetBinContent(13,	10.);
  RAABFeedDownCorrection->SetBinContent(14,	10.);
  
  PbPbMesonSelection->SetBinContent(1,		8.1);
  PbPbMesonSelection->SetBinContent(2,		8.1);
  PbPbMesonSelection->SetBinContent(3,		8.1);
  PbPbMesonSelection->SetBinContent(4,		8.1);
  PbPbMesonSelection->SetBinContent(5,		8.1);
  PbPbMesonSelection->SetBinContent(6,		8.1);
  PbPbMesonSelection->SetBinContent(7,		8.1);
  PbPbMesonSelection->SetBinContent(8,		8.1);
  PbPbMesonSelection->SetBinContent(9,		8.1);
  PbPbMesonSelection->SetBinContent(10,		1.7);
  PbPbMesonSelection->SetBinContent(11,		1.7);
  PbPbMesonSelection->SetBinContent(12,		1.7);
  PbPbMesonSelection->SetBinContent(13,		1.7);
  PbPbMesonSelection->SetBinContent(14,		1.7);
  
  PbPbSignalExtraction->SetBinContent(1,	12.8);
  PbPbSignalExtraction->SetBinContent(2,	4.3);
  PbPbSignalExtraction->SetBinContent(3,	5.8);
  PbPbSignalExtraction->SetBinContent(4,	5.4);
  PbPbSignalExtraction->SetBinContent(5,	3.7);
  PbPbSignalExtraction->SetBinContent(6,	3,7);
  PbPbSignalExtraction->SetBinContent(7,	3.4);
  PbPbSignalExtraction->SetBinContent(8,	3.4);
  PbPbSignalExtraction->SetBinContent(9,	12.0); // TO BE FIXED, TEMPORARY SET TO THE RESULT IN THE BIN 9
  PbPbSignalExtraction->SetBinContent(10,	12.0);
  PbPbSignalExtraction->SetBinContent(11,	12.0);      //was 8.6
  PbPbSignalExtraction->SetBinContent(12,	12.7);
  PbPbSignalExtraction->SetBinContent(13,	10.1);    //was 6.5
  PbPbSignalExtraction->SetBinContent(14,	17.5); 
  
  PbPbTriggerEfficiency->SetBinContent(1,	0.0);
  PbPbTriggerEfficiency->SetBinContent(2,	0.0);
  PbPbTriggerEfficiency->SetBinContent(3,	0.0);
  PbPbTriggerEfficiency->SetBinContent(4,	0.0);
  PbPbTriggerEfficiency->SetBinContent(5,	0.0);
  PbPbTriggerEfficiency->SetBinContent(6,	0.0);
  PbPbTriggerEfficiency->SetBinContent(7,	0.0);
  PbPbTriggerEfficiency->SetBinContent(8,	0.0);
  PbPbTriggerEfficiency->SetBinContent(9,	0.0);
  PbPbTriggerEfficiency->SetBinContent(10,	2.0);
  PbPbTriggerEfficiency->SetBinContent(11,	2.0);
  PbPbTriggerEfficiency->SetBinContent(12,	2.0);
  PbPbTriggerEfficiency->SetBinContent(13,	2.0);
  PbPbTriggerEfficiency->SetBinContent(14,	2.0);
  
  fPbPbPtShape->SetParameters(1.00862,-0.277991,0.325087,0.);
}

void initialization(double centL=0,double centH=0)
{  
  ppBFeedDownCorrection = new TH1D("ppBFeedDownCorrection","",nPtBins,PtBins);
  ppMesonSelection = new TH1D("ppMesonSelection","",nPtBins,PtBins);
  ppSignalExtraction = new TH1D("ppSignalExtraction","",nPtBins,PtBins);
  ppTriggerEfficiency = new TH1D("ppTriggerEfficiency","",nPtBins,PtBins);

  PbPbBFeedDownCorrection = new TH1D("PbPbBFeedDownCorrection","",nPtBins,PtBins);
  PbPbMesonSelection = new TH1D("PbPbMesonSelection","",nPtBins,PtBins);
  PbPbSignalExtraction = new TH1D("PbPbSignalExtraction","",nPtBins,PtBins);
  PbPbTriggerEfficiency = new TH1D("PbPbTriggerEfficiency","",nPtBins,PtBins);

  RAABFeedDownCorrection = new TH1D("RAABFeedDownCorrection","",nPtBins,PtBins);

  initializationPP();
  if (centL==0&&centH==100) initializationPbPbCent0100();
  if (centL==0&&centH==10) initializationPbPbCent010();
  initialized=1;  
}

void deleteinitial()
{  
  delete ppBFeedDownCorrection;
  delete ppMesonSelection;
  delete ppSignalExtraction;
  delete ppTriggerEfficiency;

  delete PbPbBFeedDownCorrection;
  delete PbPbMesonSelection;
  delete PbPbSignalExtraction;
  delete PbPbTriggerEfficiency;

  delete RAABFeedDownCorrection;
}

// =============================================================================================================
// RAA systematics
// =============================================================================================================
float normalizationUncertaintyForRAA(double centL=0,double centH=100,bool isupper=true)
{
  double sys = 0;
  sys+=ppLumiUncertainty*ppLumiUncertainty;
  sys+=PbPbNMBUncertainty*PbPbNMBUncertainty;
  
  if(centL==0&&centH==10) 
    {
      // 0-10%
      double TAAUncertainty0to10 = isupper?TAAUncertainty0to10HI:TAAUncertainty0to10LO;
      sys+=TAAUncertainty0to10*TAAUncertainty0to10;
    } 
  else
    {
      // 0-100%
      double TAAUncertainty0to100 = isupper?TAAUncertainty0to100HI:TAAUncertainty0to100LO;
      sys+=TAAUncertainty0to100*TAAUncertainty0to100;
    }
  return sqrt(sys);
}

float systematicsForRAA(double pt,double centL=0,double centH=100, double HLT=0, int stage=0)
{
  initialization(centL,centH);
  
  double sys=0;
  
  if(pt<2) return 0;
  
  if(stage==1) return sqrt(sys);
  
  sys+= PbPbSignalExtraction->GetBinContent(PbPbSignalExtraction->FindBin(pt))*
    PbPbSignalExtraction->GetBinContent(PbPbSignalExtraction->FindBin(pt));
  sys+= ppSignalExtraction->GetBinContent(ppSignalExtraction->FindBin(pt))*
    ppSignalExtraction->GetBinContent(ppSignalExtraction->FindBin(pt));
  
  if(stage==2) return sqrt(sys);
  
  sys+=ppMesonSelection->GetBinContent(ppMesonSelection->FindBin(pt))*
    ppMesonSelection->GetBinContent(ppMesonSelection->FindBin(pt));
  sys+=PbPbMesonSelection->GetBinContent(PbPbMesonSelection->FindBin(pt))*
    PbPbMesonSelection->GetBinContent(PbPbMesonSelection->FindBin(pt));
  
  sys+=fPPPtShape->Eval(pt)*fPPPtShape->Eval(pt);
  sys+=fPbPbPtShape->Eval(pt)*fPbPbPtShape->Eval(pt);
  
  //sys+=(ppTrackingEfficiency*2)*(ppTrackingEfficiency*2);
  //sys+=(PbPbTrackingEfficiency*2)*(PbPbTrackingEfficiency*2);
  
  if(centL==0&&centH==10)  
    {
      //sys+=(ppTrackingEfficiency*2)*(ppTrackingEfficiency*2);
      sys+=(PbPbTrackingEfficiency010*2)*(PbPbTrackingEfficiency010*2);      
    }	    
  if(centL==0&&centH==100) 
    {
      //sys+=(ppTrackingEfficiency*2)*(ppTrackingEfficiency*2);
      sys+=(PbPbTrackingEfficiency0100*2)*(PbPbTrackingEfficiency0100*2);      
    }
  if (stage==3) return sqrt(sys);
  
  //   sys+= ppBFeedDownCorrection->GetBinContent(ppBFeedDownCorrection->FindBin(pt))*
  //         ppBFeedDownCorrection->GetBinContent(ppBFeedDownCorrection->FindBin(pt));
  
  sys+=(RAABFeedDownCorrection->GetBinContent(RAABFeedDownCorrection->FindBin(pt)))*(RAABFeedDownCorrection->GetBinContent(RAABFeedDownCorrection->FindBin(pt)));
  //sys+=(PbPbBFeedDownCorrection->GetBinContent(PbPbBFeedDownCorrection->FindBin(pt))/2)*(PbPbBFeedDownCorrection->GetBinContent(PbPbBFeedDownCorrection->FindBin(pt))/2);

  sys+= PbPbTriggerEfficiency->GetBinContent(PbPbTriggerEfficiency->FindBin(pt))*
    PbPbTriggerEfficiency->GetBinContent(PbPbTriggerEfficiency->FindBin(pt));
  sys+= ppTriggerEfficiency->GetBinContent(ppTriggerEfficiency->FindBin(pt))*
    ppTriggerEfficiency->GetBinContent(ppTriggerEfficiency->FindBin(pt));


  deleteinitial();

  return sqrt(sys);
}

// =============================================================================================================
// RCP systematics
// =============================================================================================================
float normalizationUncertaintyForRCP(double centL=0,double centH=100)
{
   return 0;
}

float systematicsForRCP(double pt, double HLT=0,double centL=0,double centH=100)
{
  //if (!initialized && centL==0&&centH==100) initializationPbPbCent0100();
  // if (!initialized && centL==0&&centH==10) initializationPbPbCent010();
   return 0.2;

}


// =============================================================================================================
// cross-section systematics for pp
// =============================================================================================================
float normalizationUncertaintyForPP()
{
   return sqrt((DKpiBRUncertainty*DKpiBRUncertainty)+(ppLumiUncertainty*ppLumiUncertainty));
}

float systematicsPP(double pt, double HLT=0,int stage=0)
{
  initialization();
  double sys=0;
  
  if(stage==1) return sqrt(sys);
  
  sys+=ppSignalExtraction->GetBinContent(ppSignalExtraction->FindBin(pt))* 
    ppSignalExtraction->GetBinContent(ppSignalExtraction->FindBin(pt));
  
  if(stage==2) return sqrt(sys);
  
  sys+=(ppTrackingEfficiency*2)*(ppTrackingEfficiency*2);
  sys+=ppLumiUncertainty*ppLumiUncertainty;
  sys+=ppMesonSelection->GetBinContent(ppMesonSelection->FindBin(pt))* 
    ppMesonSelection->GetBinContent(ppMesonSelection->FindBin(pt));
  sys+=fPPPtShape->Eval(pt)*fPPPtShape->Eval(pt);
  
  if(stage==3) return sqrt(sys);
  
  sys+=ppBFeedDownCorrection->GetBinContent(ppBFeedDownCorrection->FindBin(pt))* 
    ppBFeedDownCorrection->GetBinContent(ppBFeedDownCorrection->FindBin(pt));
  
  sys+=ppTriggerEfficiency->GetBinContent(ppTriggerEfficiency->FindBin(pt))* 
    ppTriggerEfficiency->GetBinContent(ppTriggerEfficiency->FindBin(pt));
  
  deleteinitial();
  
  return sqrt(sys);
}

// =============================================================================================================
// cross-section systematics for PbPb
// =============================================================================================================
float normalizationUncertaintyForPbPb(double centL=0,double centH=100,bool isupper=true)
{
  double sys = ((DKpiBRUncertainty*DKpiBRUncertainty)+(PbPbNMBUncertainty*PbPbNMBUncertainty));
  if(centL==0&&centH==10)
    {
      // 0-10%
      double TAAUncertainty0to10 = isupper?TAAUncertainty0to10HI:TAAUncertainty0to10LO;
      sys+=TAAUncertainty0to10*TAAUncertainty0to10;
    }
  else
    {
      // 0-100%
      double TAAUncertainty0to100 = isupper?TAAUncertainty0to100HI:TAAUncertainty0to100LO;
      sys+=TAAUncertainty0to100*TAAUncertainty0to100;
    }   
  return sqrt(sys);
}


float systematicsPbPb(double pt, double centL=0,double centH=100, double HLT=0, int stage=0)
{
  initialization(centL,centH);
  
  double sys=0;
  
  if(stage==1) return sqrt(sys);
  
  sys+= PbPbSignalExtraction->GetBinContent(PbPbSignalExtraction->FindBin(pt))* 
    PbPbSignalExtraction->GetBinContent(PbPbSignalExtraction->FindBin(pt));
  
  if(stage==2) return sqrt(sys);
  if(pt<20) sys+=PbPbNMBUncertainty*PbPbNMBUncertainty;
  else sys+=PbPbLumiUncertainty*PbPbLumiUncertainty;
  
  if(centL==0&&centH==100) sys+=(PbPbTrackingEfficiency0100*2)*(PbPbTrackingEfficiency0100*2);
  if(centL==0&&centH==10) sys+=(PbPbTrackingEfficiency010*2)*(PbPbTrackingEfficiency010*2);

  sys+=PbPbMesonSelection->GetBinContent(PbPbMesonSelection->FindBin(pt))* 
    PbPbMesonSelection->GetBinContent(PbPbMesonSelection->FindBin(pt));
  
  sys+=fPbPbPtShape->Eval(pt)*fPbPbPtShape->Eval(pt);

  //sys+=TAAUncertainty0to100*TAAUncertainty0to100;
  
  if(stage==2) return sqrt(sys);
  sys+=PbPbBFeedDownCorrection->GetBinContent(PbPbBFeedDownCorrection->FindBin(pt))* 
    PbPbBFeedDownCorrection->GetBinContent(PbPbBFeedDownCorrection->FindBin(pt));
  
  sys+= PbPbTriggerEfficiency->GetBinContent(PbPbTriggerEfficiency->FindBin(pt))* 
    PbPbTriggerEfficiency->GetBinContent(PbPbTriggerEfficiency->FindBin(pt));
  
  deleteinitial();

  return sqrt(sys);
}


// =============================================================================================================
// Drawer
// =============================================================================================================
void drawSys(double x1,double y1, double x2,double y2, int color = 1)
{
  TLine *l1 = new TLine(x1,y1/100.,x2,y2/100.);
  TLine *l2 = new TLine(x1,-y1/100.,x2,-y2/100.);
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l1->SetLineColor(color);
  l2->SetLineColor(color);
  l1->Draw();
  l2->Draw();
  
}

// =============================================================================================================
// Plot systematics for RAA
// =============================================================================================================
void plotSystematicsRAA(double centL=0,double centH=10)
{
  
  TCanvas*canvas=new TCanvas("canvas","canvas",600,600);//550,500
  canvas->cd();
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.185);//0.200
  canvas->SetRightMargin(0.045);
  canvas->SetTopMargin(0.080);
  canvas->SetBottomMargin(0.150);
  canvas->SetFrameBorderMode(0);
  canvas->SetLogx();

  TH2F* hempty=new TH2F("hempty","",50,1,100.,10.,-0.65,0.65);
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetYaxis()->SetTitle("Systematical Uncertainty");
  hempty->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
  hempty->GetXaxis()->SetTitleOffset(1.40);//0.9
  hempty->GetYaxis()->SetTitleOffset(1.45);//1.
  hempty->GetXaxis()->SetTitleSize(0.05);//0.045
  hempty->GetYaxis()->SetTitleSize(0.05);//0.045
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.050);//0.035
  hempty->GetYaxis()->SetLabelSize(0.050);//0.035  
  hempty->GetXaxis()->SetLabelOffset(0.01);
  hempty->Draw();

   drawSys(1,0, 1,normalizationUncertaintyForRAA(centL,centH),2);
   drawSys(1,normalizationUncertaintyForRAA(centL,centH), 1.5,normalizationUncertaintyForRAA(centL,centH),2);
   drawSys(1.5,0, 1.5,normalizationUncertaintyForRAA(centL,centH),2);


   drawSys(2,0, 2,systematicsForRAA(2,centL,centH,0,0),1);


   for (double i=2;i<100;i+=0.1)
   {      
     drawSys(i,systematicsForRAA(i,centL,centH,0,0), i+0.1,systematicsForRAA(i+0.1,centL,centH,0,0),1);
     //      drawSys(i,systematicsForRAA(i,0,1), i+0.1,systematicsForRAA(i+0.1,0,1),2);
     drawSys(i,sqrt((systematicsForRAA(i,centL,centH,0,2)*systematicsForRAA(i,centL,centH,0,2))-(systematicsForRAA(i,centL,centH,0,1)*systematicsForRAA(i,centL,centH,0,1))),
             i+0.1,sqrt((systematicsForRAA(i+0.1,centL,centH,0,2)*systematicsForRAA(i+0.1,centL,centH,0,2))-(systematicsForRAA(i+0.1,centL,centH,0,1)*systematicsForRAA(i+0.1,centL,centH,0,1))),4);
     drawSys(i,sqrt((systematicsForRAA(i,centL,centH,0,3)*systematicsForRAA(i,centL,centH,0,3))-(systematicsForRAA(i,centL,centH,0,2)*systematicsForRAA(i,centL,centH,0,2))),
             i+0.1,sqrt((systematicsForRAA(i+0.1,centL,centH,0,3)*systematicsForRAA(i+0.1,centL,centH,0,3))-(systematicsForRAA(i+0.1,centL,centH,0,2)*systematicsForRAA(i+0.1,centL,centH,0,2))),kGreen+2);
     drawSys(i,sqrt((systematicsForRAA(i,centL,centH,0,0)*systematicsForRAA(i,centL,centH,0,0))-(systematicsForRAA(i,centL,centH,0,3)*systematicsForRAA(i,centL,centH,0,3))),
             i+0.1,sqrt((systematicsForRAA(i+0.1,centL,centH,0,0)*systematicsForRAA(i+0.1,centL,centH,0,0))-(systematicsForRAA(i+0.1,centL,centH,0,3)*systematicsForRAA(i+0.1,centL,centH,0,3))),kMagenta);     
   }
   
   TH1D *h1 = new TH1D("h1","",100,0,1);
   h1->SetLineWidth(2); h1->SetLineColor(1);
   TH1D *h2 = new TH1D("h2","",100,0,1);
   h2->SetLineWidth(2); h2->SetLineColor(2);
   TH1D *h4 = new TH1D("h4","",100,0,1);
   h4->SetLineWidth(2); h4->SetLineColor(4);
   TH1D *h5 = new TH1D("h5","",100,0,1);
   h5->SetLineWidth(2); h5->SetLineColor(kGreen+2);
    TH1D *h6 = new TH1D("h6","",100,0,1);
   h6->SetLineWidth(2); h6->SetLineColor(kMagenta);
    
  TLatex* texlumi = new TLatex(0.19,0.936,"25.8 pb^{-1} (5.02 TeV pp) + 404 #mub^{-1} (5.02 TeV PbPb)");
  texlumi->SetNDC();
  //texlumi->SetTextAlign(31);
  texlumi->SetTextFont(42);
  texlumi->SetTextSize(0.036);
  texlumi->SetLineWidth(2);
  texlumi->Draw();
  TLatex* texcms = new TLatex(0.22,0.90,"CMS");
  texcms->SetNDC();
  texcms->SetTextAlign(13);
  texcms->SetTextFont(62);//61
  texcms->SetTextSize(0.06);
  texcms->SetLineWidth(2);
  texcms->Draw();
  TLatex* texpre = new TLatex(0.22,0.84,"Performance");
  texpre->SetNDC();
  texpre->SetTextAlign(13);
  texpre->SetTextFont(52);
  texpre->SetTextSize(0.04);
  texpre->SetLineWidth(2);
  texpre->Draw();

  TLatex * texY = new TLatex(0.5,0.8324607,"D^{0} R_{AA}, |y| < 1");//0.2612903,0.8425793
  texY->SetNDC();
  texY->SetTextColor(1);
  texY->SetTextFont(42);
  texY->SetTextSize(0.045);
  texY->SetLineWidth(2);
  texY->Draw();

  TString texper="%";
  TLatex * tlatexeff2 = new TLatex(0.5268456,0.7678883,Form("Centrality %.0f-%.0f%s",centL,centH,texper.Data()));//0.2612903,0.8425793
  tlatexeff2->SetNDC();
  tlatexeff2->SetTextColor(1);
  tlatexeff2->SetTextFont(42);
  tlatexeff2->SetTextSize(0.045);
  tlatexeff2->SetLineWidth(2);
  tlatexeff2->Draw();

   TLegend *leg = new TLegend(0.2147651,0.1762653,0.7818792,0.3717277);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.04);
   leg->AddEntry(h2,"Overall Normalization (N_{MB}, Lumi)","l");
   leg->AddEntry(h1,"Total Systematics","l");
   leg->AddEntry(h4,"Signal Extraction","l");
   leg->AddEntry(h5,"D Meson Selection and Correction","l");
   leg->AddEntry(h6,"B feed down subtraction","l");
   leg->Draw();
   canvas->SaveAs(Form("SystematicSummaryPbPb_Cent%d.pdf",(int)centH));
}

// =============================================================================================================
// Plot systematics for cross section
// =============================================================================================================
void plotSystematicsPP()
{
  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);//550,500
  canvas->cd();
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.185);//0.200
  canvas->SetRightMargin(0.045);
  canvas->SetTopMargin(0.080);
  canvas->SetBottomMargin(0.150);
  canvas->SetFrameBorderMode(0);
  canvas->SetLogx();

  TH2F* hempty=new TH2F("hempty","",50,1,100.,10.,-0.65,0.65);
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetYaxis()->SetTitle("Systematical Uncertainty");
  hempty->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
  hempty->GetXaxis()->SetTitleOffset(1.40);//0.9
  hempty->GetYaxis()->SetTitleOffset(1.45);//1.
  hempty->GetXaxis()->SetTitleSize(0.05);//0.045
  hempty->GetYaxis()->SetTitleSize(0.05);//0.045
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.050);//0.035
  hempty->GetYaxis()->SetLabelSize(0.050);//0.035  
  hempty->GetXaxis()->SetLabelOffset(0.01);
  hempty->Draw();

   drawSys(1,0, 1,normalizationUncertaintyForPP(),2);
   drawSys(1,normalizationUncertaintyForPP(), 1.5,normalizationUncertaintyForPP(),2);
   drawSys(1.5,0, 1.5,normalizationUncertaintyForPP(),2);


   drawSys(2,0, 2,systematicsPP(2),1);

   for (double i=2;i<100;i+=0.1)
   {      
      drawSys(i,systematicsPP(i,0,0), i+0.1,systematicsPP(i+0.1,0,0),1);
//      drawSys(i,systematicsPP(i,0,1), i+0.1,systematicsPP(i+0.1,0,1),2);
      drawSys(i,sqrt((systematicsPP(i,0,2)*systematicsPP(i,0,2))-(systematicsPP(i,0,1)*systematicsPP(i,0,1))),
              i+0.1,sqrt((systematicsPP(i+0.1,0,2)*systematicsPP(i+0.1,0,2))-(systematicsPP(i+0.1,0,1)*systematicsPP(i+0.1,0,1))),4);
      drawSys(i,sqrt((systematicsPP(i,0,3)*systematicsPP(i,0,3))-(systematicsPP(i,0,2)*systematicsPP(i,0,2))),
              i+0.1,sqrt((systematicsPP(i+0.1,0,3)*systematicsPP(i+0.1,0,3))-(systematicsPP(i+0.1,0,2)*systematicsPP(i+0.1,0,2))),kGreen+2);
      drawSys(i,sqrt((systematicsPP(i,0,0)*systematicsPP(i,0,0))-(systematicsPP(i,0,3)*systematicsPP(i,0,3))),
              i+0.1,sqrt((systematicsPP(i+0.1,0,0)*systematicsPP(i+0.1,0,0))-(systematicsPP(i+0.1,0,3)*systematicsPP(i+0.1,0,3))),kMagenta);

   }

   TH1D *h1 = new TH1D("h1","",100,0,1);
   h1->SetLineWidth(2); h1->SetLineColor(1);
   TH1D *h2 = new TH1D("h2","",100,0,1);
   h2->SetLineWidth(2); h2->SetLineColor(2);
   TH1D *h4 = new TH1D("h4","",100,0,1);
   h4->SetLineWidth(2); h4->SetLineColor(4);
   TH1D *h5 = new TH1D("h5","",100,0,1);
   h5->SetLineWidth(2); h5->SetLineColor(kGreen+2);
   TH1D *h6 = new TH1D("h6","",100,0,1);
   h6->SetLineWidth(2); h6->SetLineColor(kMagenta);
    
  //TLatex* texlumi = new TLatex(0.19,0.936,"25.8 pb^{-1} (5.02 TeV pp) + 404 #mub^{-1} (5.02 TeV PbPb)");
  TLatex* texlumi = new TLatex(0.35,0.936,"25.8 pb^{-1} (5.02 TeV pp)");
  texlumi->SetNDC();
  //texlumi->SetTextAlign(31);
  texlumi->SetTextFont(42);
  texlumi->SetTextSize(0.045);
  texlumi->SetLineWidth(2);
  texlumi->Draw();
  TLatex* texcms = new TLatex(0.22,0.90,"CMS");
  texcms->SetNDC();
  texcms->SetTextAlign(13);
  texcms->SetTextFont(62);//61
  texcms->SetTextSize(0.06);
  texcms->SetLineWidth(2);
  texcms->Draw();
  TLatex* texpre = new TLatex(0.22,0.84,"Performance");
  texpre->SetNDC();
  texpre->SetTextAlign(13);
  texpre->SetTextFont(52);
  texpre->SetTextSize(0.04);
  texpre->SetLineWidth(2);
  texpre->Draw();

  TLatex * texY = new TLatex(0.5,0.8324607,"D^{0} d#sigma / dp_{T}, |y| < 1");//0.2612903,0.8425793
  texY->SetNDC();
  texY->SetTextColor(1);
  texY->SetTextFont(42);
  texY->SetTextSize(0.045);
  texY->SetLineWidth(2);
  texY->Draw();
  
  TLegend *leg = new TLegend(0.2147651,0.1762653,0.7818792,0.3717277);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h2,"Overall Normalization (Lumi + BR)","l");
  leg->AddEntry(h1,"Total Systematics","l");
  leg->AddEntry(h4,"Signal Extraction","l");
  leg->AddEntry(h5,"D Meson Selection and Correction","l");
  leg->AddEntry(h6,"B feed down subtraction","l");
  leg->Draw();
  
  canvas->SaveAs("SystematicSummaryPP.pdf");
}


void plotNormalisationUnc()
{
  std::cout<<"normalisation uncertainty RAA 0-100="<<normalizationUncertaintyForRAA(0,100)<<std::endl;
  std::cout<<"normalisation uncertainty RAA 0-10="<<normalizationUncertaintyForRAA(0,10)<<std::endl;
  std::cout<<"normalisation uncertainty pp="<<normalizationUncertaintyForPP()<<std::endl;
  std::cout<<"normalisation uncertainty PbPb 0-100="<<normalizationUncertaintyForPbPb(0,100)<<std::endl;
  std::cout<<"normalisation uncertainty PbPb 0-10="<<normalizationUncertaintyForPbPb(0,10)<<std::endl; 
}
