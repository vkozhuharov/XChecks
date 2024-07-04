#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <iomanip>
#include <cmath>

#define USEROOT
#define USEMASS
#define DEBUG

#ifdef USEROOT
#include "TH1F.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"

TRandom myRandom;
#endif





const int NP = 47;  //Number of energy points
const int NS = 14;   //Number of points in the signal


double eBins[NP];
double massBins[NP];


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


expMeasurement_t expData[NP];


typedef struct {
  double X17mass;
  double gve;
//  double 

} results_t;

results_t expRes;


typedef struct {
  expMeasurement_t *data;
  results_t *res;
} experiment_t;





Double_t fitFunc(Double_t *x, Double_t *par){
  return par[0] + par[1]*x[0];
}




void init(){
  for(int i = 0;i<NP;i++) {
    memset(expData,0,NP* (sizeof(expMeasurement_t)));
  }
}

void smearParameter(expParameter_t &par){
  par.val += par.err * myRandom.Gaus();
}

void calcParError( expParameter_t &par  ){
  par.err = std::sqrt(par.stat*par.stat + par.syst*par.syst); 
}

void generateTestData(){
  
  for (int i = 0;i<NP;i++) {
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
  for(int i = 0;i<NP;i++) {
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
  for(int i = 0;i<NP;i++) {
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
  for(int i = 0;i<NP;i++) {
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

void checkGraph(TGraphErrors *gr) {
  
  //Fit function to describe the graph
  auto *fitF = new TF1("fitF",fitFunc,gr->GetX()[0],gr->GetX()[gr->GetN()-1],2);
  double minChi2=1e6;
  int iMaskStart = 0;

  for(int imask = 0;imask<gr->GetN() - NS;imask++){
    //copy the graph
    TGraphErrors grPart(*gr);
    //Remove NS points from the graph
    for(int i=0;i<NS;i++){
      grPart.RemovePoint(imask);
    }
    fitF->SetParameters(0.,0.);
    grPart.Fit(fitF,"q");
    
    double chi2 = fitF->GetChisquare();
    if (chi2<minChi2){
      minChi2 = chi2;
      iMaskStart = imask;
    }
  }
  std::cout << "==" << gr->GetName() << "==  "<< "Best chi2 achieved: " << minChi2 << "  NDF:  " <<  fitF->GetNDF()   << "  with mask starting at:  "
        << iMaskStart << std::endl;
  TGraphErrors grPart(*gr);
  for(int i=0;i<NS;i++){
    grPart.RemovePoint(iMaskStart);
  }

  double averageValue = grPart.GetMean(2); //Get mean of the values on the Y axis

  TH1F *hPulls = new TH1F("hPulsSignal","Pulls with respect to a fit",
                            50,-5*grPart.GetRMS(2),5*grPart.GetRMS(2));

  for(int i = 0;i< grPart.GetN();i++) {
    hPulls->Fill(grPart.GetY()[i] - fitF->Eval(grPart.GetX()[i]));
  }

  TGraphErrors * dgDraw = new TGraphErrors (grPart);

  TCanvas *c1 = new TCanvas();c1->cd(); dgDraw->Draw(); //hPulls->Draw();
  // getchar();
  double pullsRMS = hPulls->GetRMS();
  
  std::cout << "==" << gr->GetName() << "==  "<< "Mean value on Y:   " << averageValue << "  RMS:  " << pullsRMS << " Relative spread: " << pullsRMS/averageValue <<  std::endl;

  if(fitF) delete fitF;
  if(hPulls) delete hPulls;

}




void performSystChecks(){

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



  printGraphErrors(grSignalPot );
  printGraphErrors(grSignalPotAccCorr );
  printGraphErrors(grSignalPotAccCorrPotCorr );
  printGraphErrors(grSignalPotPotCorrAccCorr );
  

  grSignalPot->SetName("SignalPot");
  grSignalPotAccCorr->SetName("SignalPotAccCorr");
  grSignalPotAccCorrPotCorr->SetName("SignalPotAccCorrPotCorr");
  grSignalPotPotCorrAccCorr->SetName("SignalPotPotCorrAccCorr");

  checkGraph(grSignalPot);
  checkGraph(grSignalPotAccCorr);
  checkGraph(grSignalPotAccCorrPotCorr);
  checkGraph(grSignalPotPotCorrAccCorr);
 
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


void readMCTestData(std::string fData, std::string fInput){
  TFile *ParsIn = new TFile(fInput.c_str(),"READ");
  TFile *DataIn = new TFile(fData.c_str(),"READ");

  TGraphErrors *grNEvents = (TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000790_0");
  // TGraphErrors *grNEvents = (TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000592_0");
  //TGraphErrors *grNEvents = (TGraphErrors *) DataIn->Get("gNObs_Mass_16.800000_gve_0.000379_0");

  TGraphErrors *grNPoT   = (TGraphErrors *) ParsIn->Get("fPotGraphUsed");
  TGraphErrors *grEffBkg = (TGraphErrors *) ParsIn->Get("fBkgGraphUsed");
  TGraphErrors *grSigEff = (TGraphErrors *) ParsIn->Get("fEffiGraphUsed");
  std::cout << "  " << grNPoT << "  " << grEffBkg << "  " << grSigEff << "  " << std::endl;

  if(! grNPoT  || !grEffBkg || !grSigEff) {
    std::cout << "=== ERROR === Not all input data structures exist!" << std::endl;
 //   std::cout << "  " << grNPoT << "  " << grEffBkg << "  " << grSigEff << "  " << std::endl;
  }

  //Fill in the experiment measurements structure
  int nPoints = grNEvents->GetN();
  std::cout << "Number of experimental points to be used: " << nPoints<< std::endl;

  if (nPoints != NP ) {
    std::cout << "Number of points is mismatched between expected and provided, better stop" << std::endl;
    return;
  }

  std::cout << grNPoT->GetN() << "  " << grEffBkg->GetN() << "  " << grSigEff->GetN() << std::endl;

  if (grNPoT->GetN() != nPoints 
    || grEffBkg->GetN() != nPoints
    || grSigEff->GetN() != nPoints) {
      std::cout << "==== ERROR ====" <<  "Mismatched number of points!" << std::endl;
      return;
  }
  std::cout << "Filling in the necessary structures " << std::endl;
  for(int i = 0;i<nPoints; i++){
    expData[i].mass.val = grNEvents->GetX()[i];

    expData[i].signal.val = grNEvents->GetY()[i];
    expData[i].signal.stat = std::sqrt(  expData[i].signal.val);
    calcParError(expData[i].signal);

    expData[i].pot.val = grNPoT->GetY()[i]*1e10;
    expData[i].pot.err = grNPoT->GetEY()[i]*1e10;
    expData[i].mass.err = grNPoT->GetEX()[i];

    expData[i].accCorr.val = 1./grEffBkg->GetY()[i];
    expData[i].accCorr.err = grEffBkg->GetEY()[i] / (expData[i].accCorr.val * expData[i].accCorr.val )  ;

    expData[i].potCorr.val = 0.878;
    expData[i].potCorr.err = 0.;//0.08*expData[i].potCorr.val;

  }


  std::cout << "=== Reading of data points finished === " << std::endl;

  ParsIn->Close();
  DataIn->Close();

  if(ParsIn) delete ParsIn;
  if(DataIn) delete DataIn;

}

int main() {
  init();

  //generateTestData();
  readMCTestData("data/selectedSPlusB.root","data/usedInput.root");

  //readMCTestData("data/geneSignalPlusBkg.root","data/inputFilesUsed.root");

//  writeAllData("output.dat");

  //readAllData("output.dat");

  //  printAllData();
  performSystChecks();
  
  return 0;
}
