#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "TString.h"
#include <string>
#include <array>
#include <map>
#include <iostream>
#include "TFile.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"

using namespace std;

void readvresetd(string treepath, string run, int o2chipid){

  //open file
  TFile *infile = TFile::Open(treepath.c_str());
  TTree *thtree = (TTree*)infile->Get("ITS_calib_tree");
  short int chipid[1024], row[1024], vresetd[1024], hits[1024];
  TH2F *h2pix[1024];
  TH2F *hsummary = new TH2F("h2pix",Form("CHIPID %d - row 100;VRESETD (DAC);",o2chipid),29,100,245,1024,-0.5,1023.5);
  hsummary->SetStats(0);
  for(int i=0; i<1024; i++){
    h2pix[i] = new TH2F(Form("h2pix_%d",i),Form("(row,col) = (100,%d);VRESETD (DAC);",i),29,100,245,1,0,1);
    h2pix[i]->SetStats(0);
  }

  thtree->SetBranchAddress("chipid",&chipid[0]);
  thtree->SetBranchAddress("row",&row[0]);
  thtree->SetBranchAddress("n_hits",&hits[0]);
  thtree->SetBranchAddress("vresetd",&vresetd[0]);

  long int nentries = thtree->GetEntries();
  o2::itsmft::ChipMappingITS mp;
  int lay,sta,ssta,mod,chipInMod;

  for(long int ie=0; ie<nentries; ie++){
    thtree->GetEntry(ie);
    short int crow = row[0];
    short int cchipid = chipid[0];
    mp.expandChipInfoHW(cchipid, lay,sta, ssta, mod,  chipInMod);
    for(int icol=0; icol<1024; icol++){
      if(cchipid==o2chipid){
        h2pix[icol]->SetBinContent(((vresetd[icol]-100)/5)+1, 1, hits[icol]);
        hsummary->SetBinContent(((vresetd[icol]-100)/5)+1, icol+1, hits[icol]);
      }
    }
  }
  TCanvas cpix0("cpix0","cpix0");
  hsummary->Draw("colz");
  cpix0.SaveAs(Form("VRESET1D_300_plots_chip_%d_run%s.pdf[",o2chipid,run.c_str()));
  cpix0.SaveAs(Form("VRESET1D_300_plots_chip_%d_run%s.pdf",o2chipid,run.c_str()));
  for(int i=0; i<1024; i++){
    TCanvas cpix("cpix","cpix");
    h2pix[i]->Draw("COLZ");
    cpix.SaveAs(Form("VRESET1D_300_plots_chip_%d_run%s.pdf",o2chipid,run.c_str()));
    if(i==1023) cpix.SaveAs(Form("VRESET1D_300_plots_chip_%d_run%s.pdf]",o2chipid,run.c_str()));
  }

}
