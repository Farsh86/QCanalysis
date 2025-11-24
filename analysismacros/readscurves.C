#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include <string>
#include <array>
#include <map>
#include <unordered_map>
#include <iostream>
#include "TFile.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"

//const int chipIDmin[7] = {0, 108, 252, 432,3120, 6480, 14712};
//const int chipIDmax[7] = {107, 251, 431,3119, 6479, 14711, 24119};

//MAIN
void readscurves(std::string treepath = "../Data/tree.root", int runnum = 111111, int chipsel = 369){

  //gStyle->SetOptStat(0000);
  //open file
  TFile *infile = TFile::Open(treepath.c_str());
  TTree *thtree = (TTree*)infile->Get("s-curve-points");
  short int chipid[1024], row[1024];
  unsigned char chg[1024], hits[1024];
  const short int nStaves[7] = {12,16,20,24,30,42,48};
  thtree->SetBranchAddress("chipid",&chipid[0]);
  thtree->SetBranchAddress("row",&row[0]);
  thtree->SetBranchAddress("hits",&hits[0]);
  thtree->SetBranchAddress("chg",&chg[0]);

  // histos 
  TFile *outfl = new TFile(Form("s-curves_run%d.root",runnum), "RECREATE");
  TH1F *hs[11][128]; // we store by default 1 pix every 8 for each row
  for(int irow = 0; irow< 11; irow++)
    for(int ipix=0; ipix<128; ipix++){
      hs[irow][ipix] = new TH1F(Form("hs_row%d_col%d",5 + irow*50, ipix*8), Form("scurve_row%d_col%d;chg(DACs);hits",5 + irow*50, ipix*8),51,-0.5,50.5);
      hs[irow][ipix]->SetStats(0);
    }

  long int nentries = thtree->GetEntries();
  o2::itsmft::ChipMappingITS mp;
  int lay,sta,ssta,mod,chipInMod;
  //loop on tree entries
  short int thischip, thisrow;
  int nrows = 0;
  for(long int i = 0; i<nentries; i++){
    thtree->GetEntry(i);
   
    thischip = chipid[0];
    thisrow = row[0];
    
    if(thischip == chipsel){
      for(int icol=0; icol<1024; icol++)
        hs[thisrow/50][icol/8]->Fill((int)chg[icol],(int)hits[icol]);
    }
  }

  // Save
  TCanvas c1("mycnv1","mycnv1",1000,800);
  hs[0][0]->Draw("HIST");
  for(int irow = 0; irow< 11; irow++){
    for(int ipix=0; ipix<128; ipix++){
      if (irow == 0 && ipix ==0) continue;
      hs[irow][ipix]->Draw("HIST same");
    }
  }
  c1.SaveAs(Form("scurves_run%d_chip%d.pdf[",runnum,chipsel));
  c1.SaveAs(Form("scurves_run%d_chip%d.pdf",runnum,chipsel));

  for(int irow = 0; irow< 11; irow++){
    for(int ipix=0; ipix<128; ipix++){ 
      TCanvas c("mycnv","mycnv",1000,800);
      hs[irow][ipix]->Draw("HIST");
      c.SaveAs(Form("scurves_run%d_chip%d.pdf",runnum,chipsel));
      if(irow==10 && ipix==127) c.SaveAs(Form("scurves_run%d_chip%d.pdf]",runnum,chipsel));
    }
  }
}
