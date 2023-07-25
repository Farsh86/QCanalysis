#include <sstream>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <string>
#include <iostream>
#include <vector>
#include <ctime>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TKey.h>
#include <TLine.h>
#include <THnSparse.h>

using namespace std;

//
// Function to compare two hitmaps saved in THnSparse --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2){

  std::array<long int,5> noisypix = {0, 0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int iybin=1; iybin<=512; iybin++){
    href->GetAxis(1)->SetRange(iybin,iybin); //select a row
    h2->GetAxis(1)->SetRange(iybin,iybin);//select a row
    TH1D *hrefproj = (TH1D*)href->Projection(0);
    TH1D *h2proj = (TH1D*)h2->Projection(0);
    for(int ixbin=1; ixbin<=9216; ixbin++){

      if(hrefproj->GetBinContent(ixbin)==1 && h2proj->GetBinContent(ixbin)==1){//dead in both runs
        noisypix[2]++;
      }
      else if(hrefproj->GetBinContent(ixbin)==1 && h2proj->GetBinContent(ixbin)==0){//dead only in ref run
        noisypix[0]++;
      }
      else if(hrefproj->GetBinContent(ixbin)==0 && h2proj->GetBinContent(ixbin)==1){//dead only in second run
        noisypix[1]++;
      }
      else continue;
    }
    delete hrefproj;
    delete h2proj;
  }

  href->GetAxis(1)->SetRange(1,512);
  h2->GetAxis(1)->SetRange(1,512);

  return noisypix;
}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,5> CompareNoisyTwoRuns(THnSparse *href, THnSparse *h2){

  std::array<long int,5> noisypix = {0, 0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int ibinref=0; ibinref<href->GetNbins(); ibinref++){
    int coordref[2];
    double bincref = href->GetBinContent(ibinref, coordref);
    bool isfound = false;
    for(int ibin2=0; ibin2<h2->GetNbins(); ibin2++){
      int coord2[2];
      double binc2 = h2->GetBinContent(ibin2,coord2);

      if((coordref[0]==coord2[0] && coordref[1]==coord2[1]) && (bincref>2 && binc2>2)){//noisy in both runs
        noisypix[2]++;
        isfound=true;
        break;
      }

    }
    if(!isfound && bincref>2)
      noisypix[0]++; //noisy only in ref run
  }

  for(int ibin2=0; ibin2<h2->GetNbins(); ibin2++){
    int coord2[2];
    double binc2 = h2->GetBinContent(ibin2, coord2);
    bool isfound = false;
    for(int ibinref=0; ibinref<href->GetNbins(); ibinref++){
      int coordref[2];
      double bincref = href->GetBinContent(ibinref,coordref);

      if((coordref[0]==coord2[0] && coordref[1]==coord2[1]) && (bincref>2 && binc2>2)){//noisy in both runs
        isfound=true;
        break;
      }

    }
    if(!isfound && binc2>2)
      noisypix[1]++; //noisy only in second run
  }
  /*cout<<"Stave/run: "<<h2->GetName()<<endl;
  cout<<"Ref run:  "<<noisypix[0]<<endl;
  cout<<"Sec run:  "<<noisypix[1]<<endl;
  cout<<"Both run: "<<noisypix[2]<<endl;*/

  return noisypix;
}

//
//Function to return the number of active (i.e. enabled) chips in a stave
//
int GetNchipsActive(THnSparse *hmap, int maxchip){
  int ix1=1;
  int ix2=1024;
  int activechips = maxchip;

  while(ix2<=hmap->GetAxis(0)->GetNbins()){
    hmap->GetAxis(0)->SetRange(ix1,ix2);
    TH2F *hproj = (TH2F*)hmap->Projection(1,0);
    if((hproj->Integral(1,1024,1,512)<1e-15))
      activechips--;
    ix1=ix2+1;
    ix2+=1024;
    delete hproj;
  }
  hmap->GetAxis(0)->SetRange(1,9216);//reset range
  return activechips;
}

//
//Set Style of a TGraph
//
void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
//Set Style of a TH1
//
void SetStyle(TH1 *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Set Style of TGraphErrors
//
void SetStyle(TGraphErrors *ge, Color_t col){
  ge->SetMarkerStyle(0);
  ge->SetMarkerColor(col);
  ge->SetFillColor(col);
  ge->SetLineColor(col);
}


//
//Function to return the number of runs in which a stave has no hits
//
int GetNrunsWOhits(TH2 *hFhrStv){
  int runswohits = 0;
  for(int iy=1; iy<=hFhrStv->GetNbinsY(); iy++){
    if(hFhrStv->GetBinContent(1,iy)<1e-15)
      runswohits++;
  }
  return runswohits;
}

//
// Get current date and time
//
string GetCurrentDateTime(bool useLocalTime) {
    stringstream currentDateTime;
    // current date/time based on current system
    time_t ttNow = time(0);
    tm * ptmNow;

    if (useLocalTime)
        ptmNow = localtime(&ttNow);
    else
        ptmNow = gmtime(&ttNow);

    currentDateTime << 1900 + ptmNow->tm_year;

    //month
    if (ptmNow->tm_mon < 9)
        //Fill in the leading 0 if less than 10
        currentDateTime << "0" << 1 + ptmNow->tm_mon;
    else
        currentDateTime << (1 + ptmNow->tm_mon);

    //day
    if (ptmNow->tm_mday < 10)
        currentDateTime << "0" << ptmNow->tm_mday << "_";
    else
        currentDateTime <<  ptmNow->tm_mday << "_";

    //hour
    if (ptmNow->tm_hour < 10)
        currentDateTime << "0" << ptmNow->tm_hour;
    else
        currentDateTime << ptmNow->tm_hour;

    //min
    if (ptmNow->tm_min < 10)
        currentDateTime << "0" << ptmNow->tm_min;
    else
        currentDateTime << ptmNow->tm_min;

    //sec
    if (ptmNow->tm_sec < 10)
        currentDateTime << "0" << ptmNow->tm_sec;
    else
        currentDateTime << ptmNow->tm_sec;


    return currentDateTime.str();
}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<float,nMasked+1> GetFHRwithMasking(THnSparse *hmap, const int nchips, double ntrig, TH2 *hhotmap){

  array<float,nMasked+1> fhrstave;

  //clone
  THnSparse *hmapclone = (THnSparse*)hmap->Clone(Form("%s_clone",hmap->GetName()));
  cout<<hmap->GetName()<<" --> "<<hmapclone->GetNbins()<<" noisy pixels"<<endl;

  int nRealMasked = nMasked;//hmapclone->GetNbins() < nMasked ? hmapclone->GetNbins() : nMasked;

  for(int iter=0; iter<nRealMasked+1; iter++){

    TH1F *hproj = (TH1F*)hmapclone->Projection(1);
    long int totalhits = hproj->Integral();
    float fhr = nchips==0 ? 0. : (float)totalhits / (512.*1024.*nchips*ntrig);

    if(ntrig<0) fhr=0.;
    fhrstave[iter] = fhr;

    int coord[2];
    double max = -1.;
    int x=0,y=0;
    long int binwithmax = 0;

    for(int ibin=0; ibin<hmapclone->GetNbins(); ibin++){
      double bincontent = hmapclone->GetBinContent(ibin, coord);
      //if(bincontent<1e-4) continue;
      if(bincontent>max){
        max=bincontent;
        binwithmax = ibin;
        x=coord[0];
        y=coord[1];
      }
    }

    if(nchips && totalhits!=0) hmapclone->SetBinContent(binwithmax,0.);
    if(totalhits!=0 && iter<nMasked/50) hhotmap->SetBinContent((x-1)/4, (y-1)/4, 1);// to avoid a marker in 0,0 for empty histos //iter<nMasked/50 to avoid crowded plot
    delete hproj;
  }

  delete hmapclone;

  return fhrstave;


}
