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
  TH2F *h2pix[1024];
  TH2F *hsummary = new TH2F("h2pix",Form("CHIPID %d - row 100;Strobe delay (#mus);Pixel",o2chipid),401,-0.5,400.5,1024,-0.5,1023.5);
  hsummary->SetStats(0);
  for(int i=0; i<1024; i++){
    h2pix[i] = new TH2F(Form("h2pix_row100_col%d",i),Form("(row,col) = (100,%d); Strobe delay (#mus);",i),401,-0.5,400.5,1,0.5,1.5);
    h2pix[i]->SetStats(0);
  }

  thtree->SetBranchAddress("chipid",&chipid[0]);
  thtree->SetBranchAddress("row",&row[0]);
  thtree->SetBranchAddress("n_hits",&hits[0]);
  thtree->SetBranchAddress("strobedel",&sdel[0]);

  long int nentries = thtree->GetEntries();
  o2::itsmft::ChipMappingITS mp;
  int lay,sta,ssta,mod,chipInMod;

  for(long int ie=0; ie<nentries; ie++){
    thtree->GetEntry(ie);
    short int crow = row[0];
    short int cchipid = chipid[0];
    //mp.expandChipInfoHW(cchipid, lay,sta, ssta, mod,  chipInMod);
    for(int icol=0; icol<1024; icol++){
      if(cchipid==o2chipid){
        h2pix[icol]->Fill(sdel[icol], 1, hits[icol]);
        hsummary->Fill(sdel[icol], icol+1, hits[icol]);
	if (icol == 0)
	  cout<<sdel[icol]<<endl;
      }
    }
  }
  TFile *outfl = new TFile(Form("PLplot_file_run%s_chip%d.root",run.c_str(), o2chipid), "RECREATE");
  TCanvas cpix0("cpix0","cpix0");
  hsummary->Draw("colz");
  cpix0.SaveAs(Form("PLplots_run%s_chip%d.pdf[",run.c_str(), o2chipid));
  cpix0.SaveAs(Form("PLplots_run%s_chip%d.pdf",run.c_str(), o2chipid));
  for(int i=0; i<1024; i++){
    TCanvas cpix("cpix","cpix");
    h2pix[i]->Draw("COLZ");
    h2pix[i]->Write();
    cpix.SaveAs(Form("PLplots_run%s_chip%d.pdf",run.c_str(), o2chipid));
    if(i==1023) cpix.SaveAs(Form("PLplots_run%s_chip%d.pdf]",run.c_str(),o2chipid));
  }
  outfl->Close();

}
