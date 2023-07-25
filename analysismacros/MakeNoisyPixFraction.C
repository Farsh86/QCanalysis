#include "fstream"
#include "iostream"
#include "TH1.h"
#include "TCanvas.h"
#include "string"
#include "TStyle.h"

using namespace std;

void MakeNoisyPixFraction(string filename, string runnum){
  gStyle->SetOptStat(0000);
  int nstaves[] = {12,16,20,24,30,42,48};
  int npix[] = {512*1024*9,512*1024*9,512*1024*9,512*1024*14*8,512*1024*14*8,512*1024*14*14,512*1024*14*14};
  int nstavestot=nstaves[0]+nstaves[1]+nstaves[2]+nstaves[3]+nstaves[4]+nstaves[5]+nstaves[6];
  ifstream infile(filename.c_str());
  string stave;
  int noisypix;
  TH1D *hhits = new TH1D("hhits",Form("Noisy pix fraction per stave run%s; ;noisy pix (per mille)",runnum.c_str()), nstavestot,0,(double)nstavestot);
  while(infile>>stave>>noisypix){
    string layer = stave.substr(1,1);
    string stavenum = stave.substr(3,2);
    int bin = 0;
    if(layer=="0"){
      bin = stoi(stavenum)+1;
    }
    else{
      for(int ilay=0; ilay<=stoi(layer)-1; ilay++)
        bin += nstaves[ilay];
      bin+=stoi(stavenum)+1;
    }
    hhits->Fill(bin, ((double)noisypix/(double)npix[stoi(layer)])*1000.);
  }

  //labels
  for(int ilay=0; ilay<7; ilay++){
    for(int is=0; is<nstaves[ilay]; is++){
      int bin=0;
      if(!ilay){
        bin = is+1;
      }
      else{
        for(int il=0; il<=ilay-1; il++)
          bin += nstaves[il];
        bin+=is+1;
      }
      hhits->GetXaxis()->SetBinLabel(bin, Form("L%d_%02d",ilay,is));
      hhits->SetBinError(bin,1e-20);
    }
  }

  TCanvas *c = new TCanvas("mycnv","mycnv",1850,600);
  c->SetMargin(0.0654,0.0054,0.101576,0.0998);
  c->SetLogy();
  c->SetTickx();
  c->SetTicky();
  c->SetGridy();
  c->cd();
  hhits->SetMarkerStyle(20);
  hhits->SetMarkerSize(1.2);
  hhits->Draw("PE");

  c->SaveAs(Form("../yaml/noise_masks/fractions/noisypix_fraction_%s.pdf",runnum.c_str()));
}
