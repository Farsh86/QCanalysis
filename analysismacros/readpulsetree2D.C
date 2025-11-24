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

void readpulsetree(string treepath, string run, int o2chipid){

  //open file
  TFile *infile = TFile::Open(treepath.c_str());
  TTree *thtree = (TTree*)infile->Get("ITS_calib_tree");
  short int chipid[1024], row[1024], sdel[1024], hits[1024];
  unsigned char charge[1024];
  TH2F *h2pix;
  // TH2F *hsummary = new TH2F("h2pix",Form("CHIPID %d - row 100;Strobe delay (#mus);Charge (e^{-});Pixel",o2chipid),100,0,10,1024,-0.5,1023.5);
  // hsummary->SetStats(0);

  h2pix = new TH2F("h2pix",Form(";Strobe delay (x25ns);Charge (DAC);"),201,0,2010,171,-0.5,170.5);
  h2pix->SetStats(0);

  TH1D *hchipid = new TH1D("hchipid","Chip IDs; ChipID; Counts",24120, -0.5, 24119.5);

  // for(int i=0; i<1024; i++){
  //   h2pix[i] = new TH2F("h2pix",Form("(row,col) = (100,%d); Strobe delay (#mus);Charge (e^{-});",i),100,0,10,1,0,1);
  //   h2pix[i]->SetStats(0);
  // }

  thtree->SetBranchAddress("chipid",&chipid[0]);
  thtree->SetBranchAddress("row",&row[0]);
  thtree->SetBranchAddress("n_hits",&hits[0]);
  thtree->SetBranchAddress("strobedel",&sdel[0]);
  thtree->SetBranchAddress("charge",&charge[0]);

  long int nentries = thtree->GetEntries();
  o2::itsmft::ChipMappingITS mp;
  int lay,sta,ssta,mod,chipInMod;

  for(long int ie=0; ie<nentries; ie++){
    thtree->GetEntry(ie);
    short int crow = row[0];
    short int cchipid = chipid[0];
    //mp.expandChipInfoHW(cchipid, lay,sta, ssta, mod,  chipInMod);
    if(hchipid->GetBinContent(cchipid+1) < 0.01)
      hchipid->Fill(cchipid);
    if(cchipid==o2chipid && crow == 3){
      cout<<crow<<" "<<sdel[1]<<" "<<charge[1]<<" "<<hits[1]<<endl;
      h2pix->Fill(sdel[1], ((int)charge[1]), hits[1]);

    // for(int icol=0; icol<1024; icol++){
      // if(cchipid==o2chipid){
      //   h2pix[icol]->SetBinContent((sdel[icol]/4)+1, 1, hits[icol]);
      //   hsummary->SetBinContent((sdel[icol]/4)+1, icol+1, hits[icol]);
      // }
    }
  }
  TCanvas cpix0("cpix0","cpix0");
  // hsummary->Draw("colz");
  // cpix0.SaveAs(Form("PLplots_run%s.pdf[",run.c_str()));
  // cpix0.SaveAs(Form("PLplots_run%s.pdf",run.c_str()));


  // for(int i=0; i<1024; i++){
  TCanvas *cpix = new TCanvas("cpix","cpix");
  h2pix->Draw("COLZ");
  cpix->SaveAs(Form("PL2Dplots_run%s.pdf",run.c_str()));
    // if(i==1023) cpix.SaveAs(Form("PLplots_run%s.pdf]",run.c_str()));
  // }
  TCanvas *cchip = new TCanvas("cchip","cchip");
  hchipid->Draw("HIST");
  cchip->SaveAs(Form("Chips_run%s.pdf",run.c_str()));

}
