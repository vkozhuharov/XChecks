#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define USEROOT
#define USEMASS

#define DEBUGALL
#define DEBUG

#define CUSTOMFIT
#define SAVEALLPLOTS



#ifdef USEROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"
#include "TFitResult.h"

TRandom myRandom;
#endif

double errCorr = 0.;



const int NPMAX = 47;  //Number of energy points
const int NS = 10;   //Number of points in the signal
const float ErrScale = 1.;

int NP = 0;

double eBins[NPMAX];
double massBins[NPMAX];


const double eMass = 0.511;  // in MeV

typedef struct {
  double val;
  double stat;
  double syst;
  double err;
} expParameter_t;


typedef struct {
  expParameter_t E;
  expParameter_t mass;     //SQRT(s)
  expParameter_t bkg;      // Bkg from MC
  expParameter_t signal;   //n2cl from data
  expParameter_t pot;      //NPoT
  expParameter_t potCorr;  
  expParameter_t accCorr;
  expParameter_t potMC;
} expMeasurement_t;


expMeasurement_t expData[NPMAX];

typedef struct {
  //For the original graph
  double chi2;
  double pullsRMS;
  double meanVal;
  int NDF;
  double prob;
  double fitPars[2];//[2];  //second index - [0]: value of the par; [1] - error

  //For the graph with masked region
  double minChi2;
  int posMask;
  int NDFMask;
  double probMask;
  double massMask;
  double pullsRMSMask;
  double meanValMask;
  double fitParsMask[2];//[2];
  double fitParsMaskCorr[2][2];//[2];

  double pulls[NPMAX-NS];
  double pulls_ov_sigma[NPMAX-NS];

#ifdef   SAVEALLPLOTS
  TGraphErrors grN2ovNPoT;
  TGraphErrors grN2ovNPoTCorr;
  TGraphErrors grN2ovNPoTCorrMask;
#endif
  TGraphErrors grPullsOvSigma;

} graphTestRes_t; 

typedef struct {
  double X17mass;
  double gve;
  graphTestRes_t grRes;
} results_t;

results_t expRes;


 typedef struct {
  expMeasurement_t data[NPMAX];
  results_t res;
  int success;
} experiment_t;

std::vector<experiment_t> expCollection;

experiment_t modExp;


#define NFITPARS 2
const Double_t fitFuncInitalPars[NFITPARS] = {1.,0};
Double_t fitFunc(Double_t *x, Double_t *par){
// return par[0] + par[1]*x[0];
  //return par[0];
  return par[0] + par[1]*(x[0] - 16.92);
  //  return par[0] + par[1]*(x[0]-16.9) + par[2]*(x[0]-16.9)*(x[0] - 16.9);
}


static Double_t SignalShape(double signalpeak,double bes,double lorewidth,double massn, double sqrts){ // return number of signal events produced(not the accepted) / POT at sqrts(s) at gven = 1
  const double me = 0.511;// MeV
  const double Wb = -7E-6;// MeV
  double eRes = (massn*massn - 2*me*me)/(2.*(me+Wb));
  double eBeam = (sqrts*sqrts - 2*me*me)/(2.*me); // s = 2me^2 + 2meEbeam -> eBeam = (s-2me^2) / 2me
  return signalpeak*TMath::Voigt(eBeam-eRes,bes*eRes,lorewidth*2,5); // last input was 4 but with low BES can have rounding problems
}



void init(){
  for(int i = 0;i<NPMAX;i++) {
    memset(expData,0,NPMAX* (sizeof(expMeasurement_t)));
  }
}

void smearParameter(expParameter_t &par){
  par.val += par.err * myRandom.Gaus();
}

void calcParError( expParameter_t &par  ){
  par.err = std::sqrt(par.stat*par.stat + par.syst*par.syst); 
}

void generateTestData(){
  
  for (int i = 0;i<NPMAX;i++) {
    expData[i].E.val = 275. + 0.7*i;
    
    expData[i].mass.val = std::sqrt(2.*  expData[i].E.val * eMass);



    expData[i].signal.val = 40000;
    expData[i].signal.stat = std::sqrt(  expData[i].signal.val);
    calcParError(expData[i].signal);
    smearParameter(expData[i].signal);    

    expData[i].pot.val = 1e10;
    expData[i].pot.stat = std::sqrt(expData[i].pot.val);
    expData[i].pot.syst = 0.01*expData[i].pot.val;
    calcParError(expData[i].pot);
    smearParameter(expData[i].pot);

    expData[i].potCorr.val = 0.878;
    expData[i].potCorr.err = 0.;//0.08*expData[i].potCorr.val;

    expData[i].accCorr.val = 1.;
    expData[i].accCorr.err = 0.01;
    smearParameter(expData[i].accCorr);  
  }
}


//=========================================================================
void writeParameter(std::ofstream &fs, expParameter_t &par){
  fs.write((char *)&(par.val),sizeof (par.val));
  fs.write((char *)&(par.stat),sizeof (par.stat));
  fs.write((char *)&(par.syst),sizeof (par.syst));
  fs.write((char *)&(par.err),sizeof (par.err));

}


void writeDataPoint(std::ofstream &fs,expMeasurement_t &point){
  writeParameter(fs,point.E);
  writeParameter(fs,point.mass);
  writeParameter(fs,point.bkg);
  writeParameter(fs,point.signal);
  writeParameter(fs,point.pot);  
  writeParameter(fs,point.potCorr);
  writeParameter(fs,point.accCorr);
}


void writeAllData(std::string fname){
  std::ofstream fs(fname, std::ios::out | std::ios::binary );
  for(int i = 0;i<NPMAX;i++) {
    writeDataPoint(fs,(expData[i]));
  }
  fs.close();
}


//=========================================================================
void readParameter(std::ifstream &fs, expParameter_t &par){
  fs.read((char *)&(par.val),sizeof (par.val));
  fs.read((char *)&(par.stat),sizeof (par.stat));
  fs.read((char *)&(par.syst),sizeof (par.syst));
  fs.read((char *)&(par.err),sizeof (par.err));

}

void readDataPoint(std::ifstream &fs,expMeasurement_t &point){
  readParameter(fs,point.E);
  readParameter(fs,point.mass);
  readParameter(fs,point.bkg);
  readParameter(fs,point.signal);
  readParameter(fs,point.pot);  
  readParameter(fs,point.potCorr);
  readParameter(fs,point.accCorr);
}

void readAllData(std::string fname){
  std::ifstream fs(fname, std::ios::in | std::ios::binary );
  for(int i = 0;i<NPMAX;i++) {
    readDataPoint(fs,(expData[i]));
  }
  fs.close();
}

//=========================================================================
void printParameter(expParameter_t &par){
  std::cout << std::setprecision(4);
  std::cout <<" "<< par.val << " " << par.stat << " " << par.syst ;

}

void printDataPoint(expMeasurement_t &point){
  printParameter(point.E);
  printParameter(point.mass);
  printParameter(point.bkg);
  printParameter(point.signal);
  printParameter(point.pot);  
  printParameter(point.potCorr);
  printParameter(point.accCorr);
  std::cout << std::endl;
}

void printAllData(){
  for(int i = 0;i<NPMAX;i++) {
    printDataPoint((expData[i]));
  }
}

//=========================================================================
// Divide and multiply graphs with errors


int divideGraphErrors(  TGraphErrors *g1,TGraphErrors *g2, TGraphErrors *gres ) {
  int nPoints = g1->GetN();
  if(nPoints != g2->GetN()) return 0;
  
  for(int i = 0;i<nPoints;i++) {
    if(g2->GetY()[i] != 0) {
      gres->AddPoint(g1->GetX()[i],g1->GetY()[i]/ g2->GetY()[i] );      
      gres->SetPointError( gres->GetN()-1,g1->GetEX()[i], std::sqrt(
						    g1->GetEY()[i]* g1->GetEY()[i] / (g2->GetY()[i] * g2->GetY()[i]  )
						    +  (  g1->GetY()[i]*g2->GetEY()[i]/ (g2->GetY()[i] * g2->GetY()[i])) * (g1->GetY()[i]*g2->GetEY()[i]/ (g2->GetY()[i] * g2->GetY()[i]))
						     )  )  ;
    } else {
      gres->AddPoint(g1->GetX()[i],0);
    }
  }
  return  nPoints;
}


int multiplyGraphErrors(  TGraphErrors *g1,TGraphErrors *g2, TGraphErrors *gres ) {
  int nPoints = g1->GetN();
  if(nPoints != g2->GetN()) return 0;
  
  for(int i = 0;i<nPoints;i++) {
    gres->AddPoint(g1->GetX()[i],  g1->GetY()[i]*g2->GetY()[i] );
    
    gres->SetPointError(gres->GetN()-1,g1->GetEX()[i],std::sqrt( g1->GetY()[i]* g2->GetEY()[i] *g1->GetY()[i]* g2->GetEY()[i] +
						  g2->GetY()[i]* g1->GetEY()[i] *g2->GetY()[i]* g1->GetEY()[i]));
      
  }
  return  nPoints;
}


//=========================================================================

void printGraphErrors(TGraphErrors *g){
  std::cout << "=========================" << std::endl;  
  for(int i = 0;i<g->GetN();i++) {
    std::cout << g->GetX()[i] << "+-" <<  g->GetEX()[i]<< "   " << g->GetY()[i] << "+-" <<  g->GetEY()[i]<< std::endl;
  }
  std::cout << "=========================" << std::endl;
}



//=========================================================================

void checkGraph(TGraphErrors *gr,graphTestRes_t &res) {


  double minChi2=1e6;
  int iMaskStart = 0;
  TF1 *fitF ;

  //Fit function to describe the graph
#ifdef CUSTOMFIT
  //  TF1 *
  fitF = new TF1("fitF",fitFunc,gr->GetX()[0]-0.5,gr->GetX()[gr->GetN()-1]+0.5,NFITPARS);


  fitF->SetParameters(fitFuncInitalPars);
  gr->Fit(fitF,"q");
#else

  //{
  gr->Fit("pol1","q");
  //  TF1*
  fitF = gr->GetFunction("pol1");
#endif

  TH1F *hPulls = new TH1F("hPulsSignal","Pulls with respect to a fit",
                          50,-5*gr->GetRMS(2),5*gr->GetRMS(2));
  for(int i = 0;i< gr->GetN();i++) {
    hPulls->Fill(gr->GetY()[i] - fitF->Eval(gr->GetX()[i]));
  }
  res.pullsRMS = hPulls->GetRMS();
  res.meanVal  = gr->GetMean(2);
  res.chi2 = fitF->GetChisquare();
  res.NDF = fitF->GetNDF();
  res.prob = fitF->GetProb();

  res.fitPars[0] = fitF->GetParameter(0);
 // res.fitPars[0][1] = fitF->GetParError(0);
  if(fitF->GetNpar() > 1)
    res.fitPars[1] = fitF->GetParameter(1);
 // res.fitPars[1][1] = fitF->GetParError(1);
   if(hPulls) delete hPulls;
#ifndef CUSTOMFIT
   // }
#endif
  for(int imask = 0;NS > 0 && imask<=gr->GetN() - NS;imask++){
    //copy the graph
    TGraphErrors grPart(*gr);
    //Remove NS points from the graph
    for(int i=0;i<NS;i++){
      grPart.RemovePoint(imask);
    }
#ifdef CUSTOMFIT
    fitF->SetParameters(fitFuncInitalPars);
    grPart.Fit(fitF,"q");
#else
    grPart.Fit("pol1","q");
    //    TF1*
    fitF = grPart.GetFunction("pol1");
#endif

    double chi2 = fitF->GetChisquare();
    if (chi2<minChi2){
      minChi2 = chi2;
      iMaskStart = imask;
    }
  }



  TGraphErrors grPart(*gr);
  for(int i=0;i<NS;i++){
    grPart.RemovePoint(iMaskStart);
  }


#ifdef CUSTOMFIT
  fitF->SetParameters(fitFuncInitalPars);
  TFitResultPtr fitResults = grPart.Fit(fitF,"Sq");
#else
  TFitResultPtr fitResults = grPart.Fit("pol1","Sq");
  //  TF1*
  fitF = grPart.GetFunction("pol1");
#endif

  minChi2 = fitF->GetChisquare();

#ifdef DEBUG
  std::cout << "==" << gr->GetName() << "==  "<< "Best chi2 achieved: " << minChi2 << "  NDF:  " <<  fitF->GetNDF()
#ifdef DEBUGALL    
	    << "  with mask starting at:  "
	    << iMaskStart
	    << "  Starting from chi2:   "
	    << res.chi2
#endif
	    << std::endl;
#endif


  double averageValue = res.meanValMask = grPart.GetMean(2); //Get mean of the values on the Y axis

  TH1F *hPullsMask = new TH1F("hPulsSignalMask","Pulls with respect to a fit",
                            50,-5*grPart.GetRMS(2),5*grPart.GetRMS(2));
  TH1F *hPullsOvSigmaMask = new TH1F("hPulsOvSigmaSignalMask","Pulls over sigma with respect to a fit",
                            50,-10,10);
  if(grPart.GetN() != NP - NS) {
    std::cout << "Wrong number of remaining points: " << NP << "   " << NS
	      << " Graph points:  " << grPart.GetN() << std::endl;      
  }

  for(int i = 0;i< grPart.GetN();i++) {
    hPullsMask->Fill(grPart.GetY()[i] - fitF->Eval(grPart.GetX()[i]));
    res.pulls[i] = grPart.GetY()[i] - fitF->Eval(grPart.GetX()[i]);
    res.pulls_ov_sigma[i] = (grPart.GetY()[i] - fitF->Eval(grPart.GetX()[i])) /grPart.GetEY()[i] ;
    hPullsOvSigmaMask->Fill(  (grPart.GetY()[i] - fitF->Eval(grPart.GetX()[i])) /grPart.GetEY()[i]  );
#ifdef     SAVEALLPLOTS
    res.grPullsOvSigma.SetPoint( res.grPullsOvSigma.GetN() , grPart.GetX()[i],  (grPart.GetY()[i] - fitF->Eval(grPart.GetX()[i])) /grPart.GetEY()[i]  );
#endif
    
#ifndef  SAVEALLPLOTS
    res.grPullsOvSigma.SetPoint( res.grPullsOvSigma.GetN() , i,  (grPart.GetY()[i] - fitF->Eval(grPart.GetX()[i])) /grPart.GetEY()[i]  );
#endif


  }

#ifdef DEBUGALL
  TGraphErrors * dgDraw = new TGraphErrors (grPart);

  TCanvas *c1 = new TCanvas();c1->cd(); dgDraw->Draw("AP"); //hPulls->Draw();
  // getchar();
#endif

  res.pullsRMSMask = hPullsMask->GetRMS();
  
#ifdef DEBUG
 std::cout << "==" << gr->GetName() << "==  "<< "Mean value on Y:   " << averageValue << "  RMS:  " << res.pullsRMSMask << " Relative spread: " << res.pullsRMSMask/averageValue
   //	   <<  std::endl;
   // std::cout << "==" << gr->GetName() << "=="
	   << "  RMS of the pulls/sigma:   "<< hPullsOvSigmaMask->GetRMS() <<  std::endl;

 std::cout << "==" << gr->GetName()
	   << "==  Const:   "<< fitF->GetParameter(0) << " +- " << fitF->GetParError(0) 
	   << "   Slope:  " << fitF->GetParameter(1) << " +- " << fitF->GetParError(1)
   //<< "   Quadr. term:  " << fitF->GetParameter(2) << " +- " << fitF->GetParError(2)
	   <<  std::endl;
 
 std::cout << "==" << gr->GetName() << "==  "
	   << " Covvariance matrix   " << std::endl
   	   << fitResults->GetCovarianceMatrix()(0,0) <<"   " << fitResults->GetCovarianceMatrix()(0,1)
	   << std::endl
  	   << fitResults->GetCovarianceMatrix()(1,0) <<"   " << fitResults->GetCovarianceMatrix()(1,1)
 	   <<  std::endl;

 std::cout << "==" << gr->GetName() << "==  "
	   << " Correlation matrix   " << std::endl
   	   << fitResults->GetCorrelationMatrix()(0,0) <<"   " << fitResults->GetCorrelationMatrix()(0,1)
	   << std::endl
  	   << fitResults->GetCorrelationMatrix()(1,0) <<"   " << fitResults->GetCorrelationMatrix()(1,1)
 	   <<  std::endl;



 

#endif

  //Save the results for further usage
  res.meanValMask = averageValue; 
  res.minChi2 =  minChi2;
  res.posMask = iMaskStart;

  // if(iMaskStart == NP-NS  ) {
  //   res.massMask = grPart.GetX()[iMaskStart-1] + ((NS+1)/2)*0.02;
  // } else {
  //   res.massMask = grPart.GetX()[iMaskStart] + ((NS+1)/2)*0.02;
  // }
  
  res.massMask = gr->GetX()[iMaskStart+NS/2];

  if(res.massMask < 0. ) {
    std::cout << "res.massMask:   " <<  res.massMask <<  "   iMaskStart:  " << iMaskStart
      //	      << "   averageValue   " << averageValue
	      << "    grPart.GetX()[iMaskStart]:  " << grPart.GetX()[iMaskStart-1]
	      << std::endl;
  }
  
  res.fitParsMask[0] = fitF->GetParameter(0);
  res.NDFMask = fitF->GetNDF();
  res.probMask = fitF->GetProb();
//  res.fitParsMask[0][1] = fitF->GetParError(0);
  if(fitF->GetNpar() > 1)
    res.fitParsMask[1] = fitF->GetParameter(1);
//  res.fitParsMask[1][1] = fitF->GetParError(1);

  
#ifdef SAVEALLPLOTS
  res.grN2ovNPoTCorrMask = grPart;
#endif

#ifdef CUSTOMFIT  
  if(fitF) delete fitF;
#endif
  
  if(hPullsMask) delete hPullsMask;
  if(hPullsOvSigmaMask) delete hPullsOvSigmaMask;

}




void performExpSystChecks(experiment_t &exp){
  expMeasurement_t * expData =  exp.data;

#ifdef DEBUGALL
  std::cout << "================= New experiment ===================== " << std::endl;
  std::cout << "Experiment X17 mass: " << exp.res.X17mass << std::endl;
#endif

  TGraphErrors *grSignal = new TGraphErrors();
  TGraphErrors *grPot = new TGraphErrors();
  TGraphErrors *grPotCorr = new TGraphErrors();
  TGraphErrors *grAccCorr = new TGraphErrors();


  
  //initialize the graphs
#ifndef USEMASS

  for(int i = 0;i<NP;i++) {
    grSignal->AddPoint(expData[i].E.val, expData[i].signal.val); grSignal->SetPointError(grSignal->GetN()-1, expData[i].E.err,expData[i].signal.err);
    grPot->AddPoint(expData[i].E.val, expData[i].pot.val); grPot->SetPointError ( grPot->GetN()-1,  expData[i].E.err,expData[i].pot.err);
    grPotCorr->AddPoint(expData[i].E.val, expData[i].potCorr.val); grPotCorr->SetPointError(grPotCorr->GetN()-1,expData[i].E.err,expData[i].potCorr.err);
    grAccCorr->AddPoint(expData[i].E.val, expData[i].accCorr.val); grAccCorr->SetPointError(grAccCorr->GetN()-1,expData[i].E.err,expData[i].accCorr.err);
  }
#else
  for(int i = 0;i<NP;i++) {
    grSignal->AddPoint(expData[i].mass.val, expData[i].signal.val); grSignal->SetPointError(grSignal->GetN()-1, expData[i].mass.err,expData[i].signal.err);
    grPot->AddPoint(expData[i].mass.val, expData[i].pot.val); grPot->SetPointError ( grPot->GetN()-1,  expData[i].mass.err,expData[i].pot.err);
    grPotCorr->AddPoint(expData[i].mass.val, expData[i].potCorr.val); grPotCorr->SetPointError(grPotCorr->GetN()-1,expData[i].mass.err,expData[i].potCorr.err);
    grAccCorr->AddPoint(expData[i].mass.val, expData[i].accCorr.val); grAccCorr->SetPointError(grAccCorr->GetN()-1,expData[i].mass.err,expData[i].accCorr.err);
  }



#endif

  
  TGraphErrors *grSignalPot = new TGraphErrors();
  TGraphErrors *grSignalPotAccCorr = new TGraphErrors();
  TGraphErrors *grSignalPotPotCorr = new TGraphErrors();
  TGraphErrors *grSignalPotAccCorrPotCorr = new TGraphErrors();
  TGraphErrors *grSignalPotPotCorrAccCorr = new TGraphErrors();

  
  
  divideGraphErrors(grSignal,grPot, grSignalPot);
  multiplyGraphErrors(grSignalPot,grAccCorr,grSignalPotAccCorr);
  multiplyGraphErrors(grSignalPot,grPotCorr,grSignalPotPotCorr);
  multiplyGraphErrors(grSignalPotPotCorr,grAccCorr,grSignalPotAccCorrPotCorr);
  multiplyGraphErrors(grSignalPotAccCorr,grPotCorr,grSignalPotPotCorrAccCorr);



//  printGraphErrors(grSignalPot );
//  printGraphErrors(grSignalPotAccCorr );
//  printGraphErrors(grSignalPotAccCorrPotCorr );
//  printGraphErrors(grSignalPotPotCorrAccCorr );
  

  grSignalPot->SetName("SignalPot");
  grSignalPotAccCorr->SetName("SignalPotAccCorr");
  grSignalPotAccCorrPotCorr->SetName("SignalPotAccCorrPotCorr");
  grSignalPotPotCorrAccCorr->SetName("SignalPotPotCorrAccCorr");

#ifdef SAVEALLPLOTS
  exp.res.grRes.grN2ovNPoT = (* grSignalPot );
  exp.res.grRes.grN2ovNPoTCorr = (* grSignalPotAccCorr );
#endif
  
//  checkGraph(grSignalPot,exp.res.grRes[0]);
  checkGraph(grSignalPotAccCorr,exp.res.grRes);
//  checkGraph(grSignalPotAccCorrPotCorr,exp.res.grRes[2]);
//  checkGraph(grSignalPotPotCorrAccCorr,exp.res.grRes[3]);

#ifdef SAVEALLPLOTS  

  exp.res.grRes.grN2ovNPoT.SetTitle(TString::Format("grN2ovNPoT_%03f_%06f",exp.res.X17mass,exp.res.gve));
  exp.res.grRes.grN2ovNPoT.SetName(TString::Format("grN2ovNPoT_%03f_%06f",exp.res.X17mass,exp.res.gve));

  exp.res.grRes.grN2ovNPoTCorr.SetTitle(TString::Format("grN2ovNPoTCorr_%03f_%06f",exp.res.X17mass,exp.res.gve));
  exp.res.grRes.grN2ovNPoTCorr.SetName(TString::Format("grN2ovNPoTCorr_%03f_%06f",exp.res.X17mass,exp.res.gve));


   exp.res.grRes.grN2ovNPoTCorrMask.SetTitle(TString::Format("grN2ovNPoTCorrMask_%03f_%06f",exp.res.X17mass,exp.res.gve));
  exp.res.grRes.grN2ovNPoTCorrMask.SetName(TString::Format("grN2ovNPoTCorrMask_%03f_%06f",exp.res.X17mass,exp.res.gve));
#endif

  exp.res.grRes.grPullsOvSigma.SetTitle(TString::Format("grPullsOvSigma_%03f_%06f",exp.res.X17mass,exp.res.gve));
  exp.res.grRes.grPullsOvSigma.SetName(TString::Format("grPullsOvSigma_%03f_%06f",exp.res.X17mass,exp.res.gve));

 

  if(0){ 

  grSignalPot->Fit("pol1","q");
  // TCanvas *c1 = new TCanvas();c1->cd(); grSignalPot->Draw();
  // return;
  
  TF1* fSignalPot=grSignalPot->GetFunction("pol1");
  std::cout << "Chi2:    "  << fSignalPot->GetChisquare() << std::endl;
  TH1F *hPulsSignalPot = new TH1F("hPulsSignalPot","Pulls with respect to a line fit",50,-5*grSignalPot->GetRMS(2),5*grSignalPot->GetRMS(2));
  for(int i = 0;i<grSignalPot->GetN();i++) {
    hPulsSignalPot->Fill(grSignalPot->GetY()[i] - fSignalPot->Eval(grSignalPot->GetX()[i]));
  }
  
  TCanvas *c2 = new TCanvas();c2->cd(); hPulsSignalPot->Draw();



  grSignalPotAccCorr->Fit("pol1","q");
  //  TCanvas *c1 = new TCanvas();c1->cd(); grSignalPotAccCorr->Draw();
  TF1* fSignalPotAccCorr=grSignalPotAccCorr->GetFunction("pol1");
  std::cout << "Chi2:    "  << fSignalPotAccCorr->GetChisquare() << std::endl;
  TH1F *hPulsSignalPotAccCorr = new TH1F("hPulsSignalPotAccCorr","Pulls with respect to a line fit",50,-5*grSignalPotAccCorr->GetRMS(2),5*grSignalPotAccCorr->GetRMS(2));
  for(int i = 0;i<grSignalPotAccCorr->GetN();i++) {
    hPulsSignalPotAccCorr->Fill(grSignalPotAccCorr->GetY()[i] - fSignalPotAccCorr->Eval(grSignalPotAccCorr->GetX()[i]));
  }
  
  TCanvas *c3 = new TCanvas();c3->cd(); hPulsSignalPotAccCorr->Draw();



  grSignalPotAccCorrPotCorr->Fit("pol1","q");
  //  TCanvas *c1 = new TCanvas();c1->cd(); grSignalPotAccCorrPotCorr->Draw();
  TF1* fSignalPotAccCorrPotCorr=grSignalPotAccCorrPotCorr->GetFunction("pol1");
  std::cout << "Chi2:    "  << fSignalPotAccCorrPotCorr->GetChisquare() << std::endl;
  TH1F *hPulsSignalPotAccCorrPotCorr = new TH1F("hPulsSignalPotAccCorrPotCorr","Pulls with respect to a line fit",50,-5*grSignalPotAccCorrPotCorr->GetRMS(2),5*grSignalPotAccCorrPotCorr->GetRMS(2));
  for(int i = 0;i<grSignalPotAccCorrPotCorr->GetN();i++) {
    hPulsSignalPotAccCorrPotCorr->Fill(grSignalPotAccCorrPotCorr->GetY()[i] - fSignalPotAccCorrPotCorr->Eval(grSignalPotAccCorrPotCorr->GetX()[i]));
  }

  
  
  
  TCanvas *c4 = new TCanvas();c4->cd(); hPulsSignalPotAccCorrPotCorr->Draw();

  }
  
  
//  if (0) {
  if (grSignal) delete grSignal;
  if (grPot) delete grPot;
  if (grPotCorr) delete grPotCorr;
  if (grAccCorr) delete grAccCorr;
  if (grSignalPot) delete grSignalPot;
  if (grSignalPotAccCorr) delete grSignalPotAccCorr;
  if (grSignalPotPotCorr) delete grSignalPotPotCorr;
  if (grSignalPotAccCorrPotCorr) delete grSignalPotAccCorrPotCorr;
  if (grSignalPotPotCorrAccCorr) delete grSignalPotPotCorrAccCorr;
  //}

}




//void readMCTestData(std::string fData, std::string fInput){
void readExpDataFromGraphs(TGraphErrors *grNEvents, 
			   TGraphErrors *grNPoT, 
			   TGraphErrors *grEffBkg,
			   TGraphErrors *grSigEff,
			   int mod = 0){

  std::cout << "Mode for data input: " << mod << std::endl;
//  TFile *ParsIn = new TFile(fInput.c_str(),"READ");
//  TFile *DataIn = new TFile(fData.c_str(),"READ");

//  TGraphErrors *grNEvents = (TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000790_0");
  // TGraphErrors *grNEvents = (TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000592_0");
  //TGraphErrors *grNEvents = (TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000379_0");

//  TGraphErrors *grNPoT   = (TGraphErrors *) ParsIn->Get("fPotGraphUsed");
//  TGraphErrors *grEffBkg = (TGraphErrors *) ParsIn->Get("fBkgGraphUsed");
//  TGraphErrors *grSigEff = (TGraphErrors *) ParsIn->Get("fEffiGraphUsed");
#ifdef DEBUGALL
  std::cout << "  " << grNPoT << "  " << grEffBkg << "  " << grSigEff << "  " << std::endl;
#endif

  if(! grNPoT  || !grEffBkg || !grSigEff) {
    std::cout << "=== ERROR === Not all input data structures exist!" << std::endl;
 //   std::cout << "  " << grNPoT << "  " << grEffBkg << "  " << grSigEff << "  " << std::endl;
  }

  //Fill in the experiment measurements structure
  int nPoints = NP = grNEvents->GetN();
  //#ifdef DEBUGALL
  std::cout << "Number of experimental points to be used: " << nPoints<< std::endl;
 //#endif

  if (nPoints != NP ) {
    std::cout << "Number of points is mismatched between expected and provided, better stop" << std::endl;
    return;
  }

#ifdef DEBUGALL
  std::cout << grNPoT->GetN() << "  " << grEffBkg->GetN() << "  " << grSigEff->GetN() << std::endl;
#endif

  if (grNPoT->GetN() != nPoints 
    || grEffBkg->GetN() != nPoints
    || grSigEff->GetN() != nPoints) {
      std::cout << "==== ERROR ====" <<  "Mismatched number of points!" << std::endl;
      return;
  }

  //The input structures seem to be self-consistent ... at least

  experiment_t exp;
  expMeasurement_t *expData;
  
  if(mod == 1) {
    expData = modExp.data;
  } else {
    expData = exp.data;
  }

#ifdef DEBUGALL
  std::cout << "Filling in the necessary structures " << std::endl;
#endif

  for(int i = 0;i<nPoints; i++){
    std::cout << "Reading grNEvents X point " << i << std::endl;
    expData[i].mass.val = grNEvents->GetX()[i];
    std::cout << "Reading grNEvents Y point " << i << std::endl;
    
    expData[i].signal.val = grNEvents->GetY()[i];

    expData[i].signal.stat = std::sqrt(  expData[i].signal.val);
    //    expData[i].signal.syst = expData[i].signal.val * errCorr/100;
    expData[i].signal.val += myRandom.Gaus()* expData[i].signal.val* errCorr/100 ;
    calcParError(expData[i].signal);
    if(mod==2) {
      expData[i].signal.err /= ErrScale;
    }

    
    
    std::cout << "Reading grNPoT Y point " << i << std::endl;

    expData[i].pot.val = grNPoT->GetY()[i]*1e10;
    std::cout << "Reading grNPoT EY point " << i << std::endl;
    if(mod != 2) {
      expData[i].pot.err = grNPoT->GetEY()[i]*1e10;
      std::cout << "Reading grNPoT EX point " << i << std::endl;
      expData[i].mass.err = grNPoT->GetEX()[i];
    } else {
      expData[i].pot.err = modExp.data[i].pot.err/ErrScale;
      std::cout << "Reading grNPoT EX point " << i << std::endl;
      expData[i].mass.err = modExp.data[i].mass.err/ErrScale;

    }
    std::cout << "Reading grEffBkg point " << i << std::endl;

    expData[i].accCorr.val = 1./grEffBkg->GetY()[i];
    //    expData[i].accCorr.err = grEffBkg->GetEY()[i] / expData[i].accCorr.val * expData[i].accCorr.val )  ;
    if(mod != 2) {

      expData[i].accCorr.err = (grEffBkg->GetEY()[i] / grEffBkg->GetY()[i])    * expData[i].accCorr.val   ;
    } else {
      expData[i].accCorr.err = modExp.data[i].accCorr.err/ErrScale;
    }
    //    std::cout << "Read grEffBkg point " << i << std::endl;

    expData[i].potCorr.val = 1.;// 0.878;
    expData[i].potCorr.err = 0.;//0.08*expData[i].potCorr.val;
    

  }

  if(mod != 1)
    expCollection.push_back(exp);

#ifdef DEBUGALL

  std::cout << "=== Reading of data points finished for this experiment === " << std::endl;
#endif

//  ParsIn->Close();
//  DataIn->Close();

 // if(ParsIn) delete ParsIn;
 // if(DataIn) delete DataIn;



}



void readMCTestDataFiles(std::string fData, std::string fInput){
  TFile *ParsIn = new TFile(fInput.c_str(),"READ");
  TFile *DataIn = new TFile(fData.c_str(),"READ");


  TGraphErrors *grNPoT   = (TGraphErrors *) ParsIn->Get("fPotGraphUsed");
  TGraphErrors *grEffBkg = (TGraphErrors *) ParsIn->Get("fBkgGraphUsed");
  TGraphErrors *grSigEff = (TGraphErrors *) ParsIn->Get("fEffiGraphUsed");



  TIter nextkey( DataIn->GetListOfKeys() );
  TKey *key;

  while ( (key = (TKey*)nextkey())) {
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( TGraph::Class() ) ) {
      TGraphErrors *gr = (TGraphErrors *) obj;
      const char *grName = gr->GetName();
      char buf[256]; 
      strcpy(buf,grName);
#ifdef DEBUG
      std::cout << "Found graph with name: " << grName << "  " ;//std::endl;
#endif
      char *ptr;
      //Have to parse the mass and the coupling from the name:
      char* name = strtok_r(buf,"_",&ptr); 
      if(strcmp(name,"gNObs")!= 0) continue;
      
      strtok_r(0,"_",&ptr);
      char *mass = strtok_r(0,"_",&ptr); 
      strtok_r(0,"_",&ptr);
      char *coupling = strtok_r(0,"_",&ptr);
      
      double m = atof(mass); 
      double eps = atof(coupling);

#ifdef DEBUG
      std::cout << "Mass: " << m << " Coupling " << eps << std::endl;
#endif

      TGraphErrors *grNEvents =   gr; //(TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000790_0");  

      readExpDataFromGraphs(grNEvents, grNPoT, grEffBkg, grSigEff);

#ifdef DEBUGALL
      std::cout << "ExpCollection size: " << expCollection.size() << std::endl;
#endif

      expCollection[expCollection.size()-1].res.X17mass = m;
      expCollection[expCollection.size()-1].res.gve = eps;
    }

  }



//  TGraphErrors *grNEvents = (TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000790_0");  


//  readExpDataFromGraphs(grNEvents, grNPoT, grEffBkg, grSigEff);




  ParsIn->Close();
  DataIn->Close();

  if(ParsIn) delete ParsIn;
  if(DataIn) delete DataIn;

}

void performSystChecks(){
  std::cout << "Checking the consistency of " << expCollection.size() << " experiments" << std::endl;
  for(int i = 0;i<expCollection.size();i++){
#ifdef DEBUG
    std::cout << "Processing experiment " << i+1 << " with X Mass " << expCollection[i].res. X17mass << "  and coupling " << expCollection[i].res.gve << std::endl;
#endif

    performExpSystChecks(expCollection[i]);
  }
}

void analyzeSystChecks(){
  std::cout << " Number of experiments to analyze:  " << expCollection.size() << std::endl;

  int nExp=0;
  int nSuccess=0;
  

  TFile *outHistoFile = new TFile("SystCheckResults.root","RECREATE");
  outHistoFile->cd();
  
  TDirectory *hDir = outHistoFile->mkdir("Histos");
  outHistoFile->cd();
  TDirectory *expDir = outHistoFile->mkdir("Experiments");

#ifdef SAVEALLPLOTS
  outHistoFile->cd();
  TDirectory *N2ovNPoTDir = outHistoFile->mkdir("N2ovNPoT");

  outHistoFile->cd();
  TDirectory *N2ovNPoTCorrDir = outHistoFile->mkdir("N2ovNPoTCorr");

  outHistoFile->cd();
  TDirectory *N2ovNPoTCorrMaskDir = outHistoFile->mkdir("N2ovNPoTCorrMask");
#endif

  outHistoFile->cd();
  TDirectory *PullsOvSigmaDir = outHistoFile->mkdir("PullsOvSigma");
  
  hDir->cd();
  
  TH1F *hChi2Mask = new TH1F("hChi2Mask","Chi2 distribution for the fit of the points with masked region", 1000,0.0,1000.0);

  TH1F *hMeanValMask = new TH1F("HMeanValMask","Mean value after masking",100,0.95,1.05);

  TH1F *hConstFitMask = new TH1F ("hConstFitMask","The constant parameter of the fit after masking",100,0.95,1.05);
 
  TH1F *hSlopeFitMask = new TH1F("hSlopeFitMask","Slope parameter of the fit after masking",100,-0.2,0.2);

  TH2F *hSlopeVsConstMask = new TH2F("hSlopeVsConstMask","Correlation between the slope and the constant after masking",100,-0.1,0.1,100,0.975,1.025);

  TH1F *hPullsRMSMask = new TH1F("hPullsRMSMask","RMS of the pulls of the data points wrt fit after masking",100,0.0,.1);

  TH1F *hProbMask = new TH1F("hProbMask", "Probability of the fit after masking",1000,0.0,1.);

#ifndef SAVEALLPLOTS
  TH1::AddDirectory(false);
#endif
  
  TH1F *hChi2 = new TH1F("hChi2","Chi2 distribution for the fit of the original data points", 1000,0.0,1000.0);
  TH1F *hMeanVal = new TH1F("HMeanVal","Mean value before masking",100,0.95,1.05);
  TH1F *hConstFit = new TH1F ("hConstFit","The constant parameter of the fit",100,0.95,1.05);
  TH1F *hSlopeFit = new TH1F("hSlopeFit","Slope parameter of the fit",100,-0.2,0.2);
  TH2F *hSlopeVsConst = new TH2F("hSlopeVsConst","Correlation between the slope and the constant",100,-0.1,0.1,100,0.975,1.025);
  TH1F *hPullsRMS = new TH1F("hPullsRMS","RMS of the pulls of the data points wrt fit",100,0.0,.1);
  TH1F *hProb = new TH1F("hProb", "Probability of the fit",1000,0.0,1.);

  TH1F *hMassMask = new TH1F("hMassMask","Center of the masked region",100,16.01,18.01);
  TH2F *hMassMaskVsMass = new TH2F("hMassMaskVsMass","X17 mass vs center of the masked region",100,16.01,18.01,100,16.01,18.01);

  TH1F *hPosMask = new TH1F("hPosMask","Center of the masked region",47,0,47);
  TH2F *hPosMaskVsMass = new TH2F("hPosMaskVsMass","X17 mass vs center of the masked region",100,16.01,18.01,47,0,47);
  TH2F *hPosMaskVsMass2 = new TH2F("hPosMaskVsMass2","X17 mass vs center of the masked region",47,0,47,100,16.01,18.01);


  

  
  TH1F *hTrueRecoMassDiff = new TH1F("hTrueRecoMassDiff","Difference between the true and the center of masked region mass",100,-2.,2.);

  TH2F *hConstMX17Gve = new TH2F("hConstMX17Gve","Constant parameter as a function of MX17 and Gve",100,16.01,18.01,100,0.,1e-3);

  TH2F *hConstMaskMX17Gve = new TH2F("hConstMaskMX17Gve","Constant parameter as a function of MX17 and Gve",100,16.01,18.01,100,0.,1e-3);

  TH2F *hTrueRecoMassDiffMX17Gve = new TH2F("hTrueRecoMassDiffMX17Gve","Obtained mass difference as a function of MX17 and Gve",100,16.01,18.01,100,0.,1e-3);
  for (int i = 0;i<100;i++){
    for(int j=0;j<100;j++){
      hTrueRecoMassDiffMX17Gve->SetBinContent(i+1,j+1,-10);
    }
  }
 

  TH2F *hChi2MX17Gve = new TH2F("hChi2MX17Gve","Obtained Chi2 as a function of MX17 and Gve",100,16.01,18.01,100,0.,1e-3);
  TH2F *hChi2MaskMX17Gve = new TH2F("hChi2MaskMX17Gve","Obtained Chi2 as a function of MX17 and Gve after masking",100,16.01,18.01,100,0.,1e-3);

  TH2F *hProbMX17Gve = new TH2F("hProbMX17Gve","Obtained probability as a function of MX17 and Gve",100,16.01,18.01,100,0.,1e-3);
  TH2F *hProbMaskMX17Gve = new TH2F("hProbMaskMX17Gve","Obtained probability as a function of MX17 and Gve after masking",100,16.01,18.01,100,0.,1e-3);


  //Graphs to be stored
  TGraph2D *grConstMX17Gve = new TGraph2D(); 
  grConstMX17Gve->SetNameTitle("grConstMX17Gve","Constant parameter as a function of MX17 and Gve");

  TGraph2D *grConstMaskMX17Gve = new TGraph2D(); 
  grConstMaskMX17Gve->SetNameTitle("grConstMaskMX17Gve","Constant parameter as a function of MX17 and Gve");
  
  TGraph2D *grSlopeMX17Gve = new TGraph2D(); 
  grSlopeMX17Gve->SetNameTitle("grSlopeMX17Gve","Slope parameter as a function of MX17 and Gve");

  TGraph2D *grSlopeMaskMX17Gve = new TGraph2D(); 
  grSlopeMaskMX17Gve->SetNameTitle("grSlopeMaskMX17Gve","Slope parameter as a function of MX17 and Gve");

  TGraph2D *grTrueRecoMassDiffMX17Gve  = new TGraph2D(); 
  grTrueRecoMassDiffMX17Gve->SetNameTitle("grTrueRecoMassDiffMX17Gve","Obtained mass difference as a function of MX17 and Gve");

  TGraph2D *grChi2MX17Gve   = new TGraph2D(); 
  grChi2MX17Gve->SetNameTitle("grChi2MX17Gve","Obtained Chi2 as a function of MX17 and Gve");

  TGraph2D *grChi2MaskMX17Gve   = new TGraph2D(); 
  grChi2MaskMX17Gve->SetNameTitle("grChi2MaskMX17Gve","Obtained Chi2 as a function of MX17 and Gve after masking");



  for(int i = 0; i < expCollection.size();i++){
    nExp++;
    hChi2->Fill(expCollection[i].res.grRes.chi2);
    hChi2Mask->Fill(expCollection[i].res.grRes.minChi2);

    hMeanVal->Fill(expCollection[i].res.grRes.meanVal);
    hMeanValMask->Fill(expCollection[i].res.grRes.meanValMask);

    if( expCollection[i].res.grRes.meanVal < 0.7 ) {
      std::cout << "Mass " << expCollection[i].res.X17mass << " Coupling: " << expCollection[i].res.gve << std::endl;
    }

    hConstFit->Fill(expCollection[i].res.grRes.fitPars[0]);
    hConstFitMask->Fill(expCollection[i].res.grRes.fitParsMask[0]);

    hSlopeFit->Fill(expCollection[i].res.grRes.fitPars[1]);
    hSlopeFitMask->Fill(expCollection[i].res.grRes.fitParsMask[1]);

    hSlopeVsConst->Fill(expCollection[i].res.grRes.fitPars[1],expCollection[i].res.grRes.fitPars[0]);
    hSlopeVsConstMask->Fill(expCollection[i].res.grRes.fitParsMask[1],expCollection[i].res.grRes.fitParsMask[0]);

    hPullsRMS->Fill(expCollection[i].res.grRes.pullsRMS);
    hPullsRMSMask->Fill(expCollection[i].res.grRes.pullsRMSMask);

    hProb->Fill(expCollection[i].res.grRes.prob);
    hProbMask->Fill(expCollection[i].res.grRes.probMask);
    if(expCollection[i].res.grRes.probMask > 0.2){
      nSuccess++;
    };

    hMassMask->Fill( expCollection[i].res.grRes.massMask);
    hMassMaskVsMass->Fill( expCollection[i].res.X17mass,expCollection[i].res.grRes.massMask);

    hPosMask->Fill( expCollection[i].res.grRes.posMask);
    hPosMaskVsMass->Fill( expCollection[i].res.X17mass,expCollection[i].res.grRes.posMask);
    hPosMaskVsMass2->Fill( expCollection[i].res.grRes.posMask,expCollection[i].res.X17mass);

    
    hTrueRecoMassDiff->Fill( expCollection[i].res.X17mass-expCollection[i].res.grRes.massMask);
    //    if(fabs(expCollection[i].res.X17mass-expCollection[i].res.grRes.massMask) > 5.) {
      //Debugging checks...
#ifdef DEBUGALL
      std::cout << "True Mass " << expCollection[i].res.X17mass << " Coupling: " << expCollection[i].res.gve
		<< "   Found Mass: " << expCollection[i].res.grRes.massMask 
		<< std::endl;
#endif
      
      //    }

    int binGve = (int) ((expCollection[i].res.gve ) / 1e-5 ) + 1 ;
    int binMass = (int) ((expCollection[i].res.X17mass - 16.01) / 0.02) + 1;
    hConstMX17Gve->SetBinContent(binMass,binGve,expCollection[i].res.grRes.fitPars[0]);
    grConstMX17Gve->SetPoint(grConstMX17Gve->GetN(),expCollection[i].res.X17mass,expCollection[i].res.gve,expCollection[i].res.grRes.fitPars[0]);

    hConstMaskMX17Gve->SetBinContent(binMass,binGve,expCollection[i].res.grRes.fitParsMask[0]);
    grConstMaskMX17Gve->SetPoint(grConstMaskMX17Gve->GetN(),expCollection[i].res.X17mass,expCollection[i].res.gve,expCollection[i].res.grRes.fitParsMask[0]);

    grSlopeMX17Gve->SetPoint(grSlopeMX17Gve->GetN(),expCollection[i].res.X17mass,expCollection[i].res.gve,expCollection[i].res.grRes.fitPars[1]);
    grSlopeMaskMX17Gve->SetPoint(grSlopeMaskMX17Gve->GetN(),expCollection[i].res.X17mass,expCollection[i].res.gve,expCollection[i].res.grRes.fitParsMask[1]);

    
    hTrueRecoMassDiffMX17Gve->SetBinContent(binMass,binGve, expCollection[i].res.X17mass-expCollection[i].res.grRes.massMask);
    grTrueRecoMassDiffMX17Gve->SetPoint(grTrueRecoMassDiffMX17Gve->GetN(),expCollection[i].res.X17mass,expCollection[i].res.gve,expCollection[i].res.X17mass-expCollection[i].res.grRes.massMask);

    hChi2MX17Gve->SetBinContent(binMass,binGve,expCollection[i].res.grRes.chi2);
    grChi2MX17Gve->SetPoint(grChi2MX17Gve->GetN(),expCollection[i].res.X17mass,expCollection[i].res.gve,expCollection[i].res.grRes.chi2);

    hChi2MaskMX17Gve->SetBinContent(binMass,binGve,expCollection[i].res.grRes.minChi2);
    grChi2MaskMX17Gve->SetPoint(grChi2MaskMX17Gve->GetN(),expCollection[i].res.X17mass,expCollection[i].res.gve,expCollection[i].res.grRes.minChi2);

    hProbMX17Gve->SetBinContent(binMass,binGve,expCollection[i].res.grRes.prob);
    hProbMaskMX17Gve->SetBinContent(binMass,binGve,expCollection[i].res.grRes.probMask);

    expDir->cd();
    //#ifndef SAVEALLPLOTS
    TH1::AddDirectory(true);
  //#endif
    TH1F *hPuls = new TH1F(TString::Format("Pulls_%03f_%06f",expCollection[i].res.X17mass,expCollection[i].res.gve),"Pulls",100,-0.1,0.1);
    TH1F *hPuls_ov_sigma = new TH1F(TString::Format("Pulls_wrt_line_%03f_%06f",expCollection[i].res.X17mass,expCollection[i].res.gve),"PullsLine",100,-10,10);
    for(int p = 0;p<NP-NS;p++){
      hPuls->Fill(expCollection[i].res.grRes.pulls[p]);
      hPuls_ov_sigma->Fill(expCollection[i].res.grRes.pulls_ov_sigma[p]);
    }

#ifdef SAVEALLPLOTS
    N2ovNPoTDir->cd();
    expCollection[i].res.grRes.grN2ovNPoT.Write();

    N2ovNPoTCorrDir->cd();
    expCollection[i].res.grRes.grN2ovNPoTCorr.Write();
    
    N2ovNPoTCorrMaskDir->cd();
    expCollection[i].res.grRes.grN2ovNPoTCorrMask.Write();
#endif
    
    PullsOvSigmaDir->cd();
    expCollection[i].res.grRes.grPullsOvSigma.Write();


    
    
  }


  std::cout << "=======  All exp:  " << nExp
	    << "   Successful:  " <<  nSuccess
	    << "   Success Probability:    " << (1. * nSuccess )/nExp
	    << "   Fail Probability:    " << (1. * (nExp - nSuccess) )/nExp
	    << "=======  " <<  std::endl;
  

  outHistoFile->cd();

#ifdef SAVEALLPLOTS
  
  grConstMaskMX17Gve->Write();
  grChi2MaskMX17Gve->Write();
  grSlopeMaskMX17Gve->Write();
  
  grSlopeMX17Gve->Write();  
  grTrueRecoMassDiffMX17Gve->Write();
  grChi2MX17Gve->Write();
  grConstMX17Gve->Write();  
#endif
  
  outHistoFile->Write();
  outHistoFile->Close();
}

int readExpDataFile(std::string fData, int mod = 0){
  TFile *DataIn = new TFile(fData.c_str(),"READ");


  
  TGraphErrors *grNPoT   = (TGraphErrors *) DataIn->Get("gPoT");
  TGraphErrors *grEffBkg = (TGraphErrors *) DataIn->Get("gBkg");
  TGraphErrors *grSigEff = (TGraphErrors *) DataIn->Get("gEff");
  TGraphErrors *grNEvents= (TGraphErrors *) DataIn->Get("gNObs");
  if(grNEvents == NULL)
    grNEvents= (TGraphErrors *) DataIn->Get("gN2NoBkgOvN2");

  if(grNEvents == NULL)
    exit(0);
  
  
  grNPoT->Sort();grNPoT->Scale(1e-10);
  
  grEffBkg->Sort();
  grSigEff->Sort();
  grNEvents->Sort();
  
  readExpDataFromGraphs(grNEvents, grNPoT, grEffBkg, grSigEff,mod);
  
#ifdef DEBUGALL
  std::cout << "ExpCollection size: " << expCollection.size() << std::endl;
#endif
  
  // expCollection[expCollection.size()-1].res.X17mass = 0.;
  // expCollection[expCollection.size()-1].res.gve = 0.;

  return 0;
}

int readMCBkgDataFile(std::string fData){
  TFile *DataIn = new TFile(fData.c_str(),"READ");

  TDirectory *DirIn = (TDirectory *) DataIn -> Get("BkgOnly");

  static float mass = 16.5;
  static float gve = 1e-4;
  
  TGraphErrors *grNPoT   ;
  TGraphErrors *grEffBkg ;
  TGraphErrors *grSigEff ;
  TGraphErrors *grNEvents;

  for (int iMC=0;iMC < 150;iMC++) {
    
  
    grNPoT   = (TGraphErrors *) DirIn->Get(TString::Format("gNPOT_%d",iMC));
    grEffBkg = (TGraphErrors *) DirIn->Get(TString::Format("gNBkg_%d",iMC));
    grSigEff = (TGraphErrors *) DirIn->Get(TString::Format("gEffi_%d",iMC));
    grNEvents= (TGraphErrors *) DirIn->Get(TString::Format("gNObs_%d",iMC));

    
  if(grNEvents == NULL)
    grNEvents= (TGraphErrors *) DataIn->Get("gN2NoBkgOvN2");

  if(grNEvents == NULL)
    exit(0);
  
  
  grNPoT->Sort();//grNPoT->Scale(1e-10); 
  grEffBkg->Sort();
  grSigEff->Sort();
  grNEvents->Sort();
  
  readExpDataFromGraphs(grNEvents, grNPoT, grEffBkg, grSigEff,2);
  
#ifdef DEBUGALL
  std::cout << "ExpCollection size: " << expCollection.size() << std::endl;
#endif
  mass+= 0.02;
  gve+=1e-5;
  expCollection[expCollection.size()-1].res.X17mass = mass;
  expCollection[expCollection.size()-1].res.gve = gve;
  }
  return 0;
}


int main(int argc, char **argv) {
  init();

  int data; 
  char dataFile[56];
  char mcFileSigBkg[56];
  char mcFilePotEff[56];

  int mc = 0;
  int dataFilePresent = 0;
  int mcFileSigBkgPresent=0;
  int mcFilePotEffPresent=0;
  int mcFileBkgOnlyPresent = 0;

  
  int bkg = 0;
  char mcFileBkgOnly[56];
    
  int opt;

  int nDataTypes = 0;
  
  while ((opt = getopt(argc, argv, "di:ms:p:bf:e:")) != -1) {
    switch (opt) {
    case 'd':
      data = 1;
      nDataTypes++;
     break;
    case 'e':
      errCorr = atof(optarg);
      printf("Using error correction of: %lf\n",errCorr);
     break;
     
    case 'i':      
      printf("Using data file: %s\n",optarg);
      strcpy(dataFile,optarg);
      dataFilePresent=1;
      //      nDataTypes++;
      break;
    case 'm':
      mc=1;
      nDataTypes++;
      break;
    case 's':
      printf("Using MC generated series of experiments: %s\n",optarg);
      strcpy(mcFileSigBkg,optarg);
      mcFileSigBkgPresent=1;
      break;
    case 'p':
      printf("Using MC generated PoT and determined efficiencies: %s\n",optarg);
      strcpy(mcFilePotEff ,optarg);
      mcFilePotEffPresent=1;
      break;
    case 'b':
      bkg = 1;
      break;

    case 'f':
      printf("Using background only PoT and determined efficiencies: %s\n",optarg);
      strcpy(mcFileBkgOnly ,optarg);
      nDataTypes++;
      mcFileBkgOnlyPresent =1;
      break;
      
      
    default: /* '?' */
      fprintf(stderr, "Usage:\n  \t%s [-d dataFile] \n\t%s [-m -s SignaBkgFile -p PoTFile]\n",
	      argv[0],argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  
  
  // return 0;;  
    
  //generateTestData();
 // readMCTestDataFiles("data/selectedSPlusB.root","data/usedInput.root");
  //   readMCTestDataFiles("data/geneSignalPlusBkg.root","data/inputFilesUsed.root");
  
  //  writeAllData("output.dat");
  //readAllData("output.dat");
  //  printAllData();

  if((data == 0 && mc == 0 && bkg == 0 ) ||
     //     (  dataFilePresent == 1 && mc == 1 )
     (nDataTypes > 1)
     ) {
    std::cout << " You need either MC or data to be processed " << " .... EXITING" << std::endl;
    exit(0);
  }

  
  
  if(data && dataFilePresent) {
    readExpDataFile(dataFile);
  }
  
  if(mc) {
    if(mcFileSigBkgPresent == 1 && mcFilePotEffPresent ==1) {
      readMCTestDataFiles(mcFileSigBkg,mcFilePotEff);
    } else {
      std::cout << "Please provide the two data files with MC experiments and PoT, signam and Bkg efficiency..." << " ... EXITING" << std::endl;
      exit(0);
    }
  }
  
  if(bkg) {
    readExpDataFile(dataFile,1);
    readMCBkgDataFile(mcFileBkgOnly);
    
  }
  
  
  
  performSystChecks();
  
  
  analyzeSystChecks();
  
  
  
  
  return 0;
}
