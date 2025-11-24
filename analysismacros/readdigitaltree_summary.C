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
#include <iostream>
#include <fstream>
#include <set>
#include "TFile.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TTree.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"

using namespace std;

const int chipIDmin[7] = {0, 108, 252, 432,3120, 6480, 14712};
const int chipIDmax[7] = {107, 251, 431,3119, 6479, 14711, 24119};

short int chipToLayer(int chipid);

//MAIN
void readdigitaltree_summary(std::string treepaths = "file.dat", std::string runnum = "XXXXXX", std::string title = "DIG .... runnum ....", std::string stavemap = "stave for which you want thr map", int id = 0){

  //Grid connection open 
  TGrid::Connect("alien://");

  //open file with grid paths 
  ifstream infile(treepaths.c_str());
  
  const short int nStaves[7] = {12,16,20,24,30,42,48};
  std::map<short int,std::array<int,512>> noisypix;
  std::map<short int,std::array<int,512>> goodpix; //good pix (50 hits) per dcol
  std::map<short int,std::array<int,512>> deadpix; //dead pix (0 hits) per dcol
  std::map<short int,std::array<int,512>> ineffpix; //ineff pix (0 hits) per dcol
  std::map<short int,std::vector<short int>> trackrow;
  std::map<std::string, std::vector<int>> dcolsP2;

  TH1F *hRowsCountIB = new TH1F("hRowsCountIB","Rows collected per chip in IB; Chip ID; Number of rows",432,-0.5,431.5);
  TH1F *hRowsCountML = new TH1F("hRowsCountML","Rows collected per chip in ML; Chip ID; Number of rows",6048,431.5,6479.5);
  TH1F *hRowsCountOL = new TH1F("hRowsCountOL","Rows collected per chip in OL; Chip ID; Number of rows",17640,6479.5,24119.5);
  TH1F *hRowsCountDistIB = new TH1F("hRowsCountDistIB","Rows collected per chip in IB; Number of rows; Counts",520,-0.5,519.5);
  TH1F *hRowsCountDistML = new TH1F("hRowsCountDistML","Rows collected per chip in ML; Number of rows; Counts",520,-0.5,519.5);
  TH1F *hRowsCountDistOL = new TH1F("hRowsCountDistOL","Rows collected per chip in OL; Number of rows; Counts",520,-0.5,519.5);
  hRowsCountIB->SetStats(0); hRowsCountML->SetStats(0); hRowsCountOL->SetStats(0);

  //stave map
  TH2F *hHitMapStave = 0;
  if(stavemap.find("L0")!=string::npos || stavemap.find("L1")!=string::npos || stavemap.find("L2")!=string::npos){
    hHitMapStave = new TH2F("hHitMapStave",Form("Hit map %s; Col; Row",stavemap.c_str()), 1024*9, -0.5, 1024*9-0.5, 512,-0.5,512-0.5);
  }
  else if(stavemap.find("L3")!=string::npos || stavemap.find("L4")!=string::npos){
    hHitMapStave = new TH2F("hHitMapStave",Form("Hit map %s; Col; Row",stavemap.c_str()), 1024*28, 0, 1024*28, 512*4,0,512*4);
  }
  else{
    hHitMapStave = new TH2F("hHitMapStave",Form("Hit map %s; Col; Row",stavemap.c_str()), 1024*49, 0, 1024*49, 512*4,0,512*4);
  }
  hHitMapStave->SetStats(0);

  o2::itsmft::ChipMappingITS mp;
  int lay,sta,ssta,mod,chipInMod;
  int hitcut = 50;
  string treefile; 
  //start loop on file paths
  while(infile>>treefile){
    TFile *infiletree = TFile::Open(treefile.c_str());
    TTree *thtree = (TTree*)infiletree->Get("ITS_calib_tree");
    short int chipid[1024], row[1024], hits[1024];
    long int nentries = thtree->GetEntries();
    thtree->SetBranchAddress("chipid",&chipid[0]);
    thtree->SetBranchAddress("row",&row[0]);
    thtree->SetBranchAddress("n_hits",&hits[0]); 
    //loop on tree entries
    for(long int i = 0; i<nentries; i++){
      thtree->GetEntry(i);
      short int crow = row[0];
      short int cchipid = chipid[0];
      mp.expandChipInfoHW(cchipid, lay,sta, ssta, mod,  chipInMod);
      short int layer = lay;
      std::string staveid = Form("L%d_%02d",lay,sta);

      //count rows
      int countthisrow=0;
      if(trackrow.count(cchipid)){//check for duplicates
        if(trackrow[cchipid].size()>0){
          for(int ir=0; ir<(int)trackrow[cchipid].size(); ir++){
            if(crow==trackrow[cchipid][ir]) {
              countthisrow++;
              break;
            }
          }
        }
      }

      if(cchipid<432){
        if(!countthisrow)
          hRowsCountIB->Fill((int)cchipid);
      }
      else if(cchipid<6480){
        if(!countthisrow)
          hRowsCountML->Fill(cchipid);
      }
      else{
        if(!countthisrow)
          hRowsCountOL->Fill(cchipid);
      }

      trackrow[cchipid].push_back(crow);//save rows

      if(layer<0) {
        cout<<"Layer is negative! ChipID: "<<cchipid<<endl;
        continue;
      }
      for(int icol=0; icol<1024; icol++){//loop on columns
        //store noisy pix
        if(hits[icol] > hitcut) { 
          noisypix[cchipid][icol/2]++;
        }
        else if(hits[icol] == 50) {
          goodpix[cchipid][icol/2]++;
        }
        else if(hits[icol] == 0) {
          deadpix[cchipid][icol/2]++;
        } else {
	  ineffpix[cchipid][icol/2]++;
	}
        //Fill hit map of selected stave (NOTE: no cuts applied on the tree!!)
        if(staveid==stavemap){
          if(layer<3){
            //std::cout<<crow+1<<" "<<icol<<" "<<thr[icol]<<std::endl;
            if(!countthisrow)
              hHitMapStave->SetBinContent(icol + 1 + 1024*chipInMod, crow+1, hits[icol]);
          }
          else {
            int colbin = chipInMod < 7 ? icol + 1 + chipInMod*1024 + (mod-1)*1024*7 : icol + 1 + (14-chipInMod)*1024 + (mod-1)*1024*7;
            int rowbin = 1;
            if(!ssta) rowbin = chipInMod < 7 ? crow+1 : crow+1+512;
            else rowbin = chipInMod < 7 ? crow+1+1024 : crow+1+1536;
            if(!countthisrow)
              hHitMapStave->SetBinContent(colbin, rowbin, hits[icol]);
          }
        }  
      }//end loop on cols
    }//end loop on the tree entries
  }//end while loop on file paths 

  //extract the bad dcols looking into noisy pix 
  //TH2 for IB,ML,OL
  TH2F *hBadDcolsIB = new TH2F("hbaddcolsIB",Form("%s - bad dcols; Chip; Stave",title.c_str()),9,-0.5,8.5,48,-0.5,47.5);
  TH2F *hBadDcolsML = new TH2F("hbaddcolsML",Form("%s - bad dcols; Chip; Stave",title.c_str()),112,-0.5,111.5,54,-0.5,53.5);
  TH2F *hBadDcolsOL = new TH2F("hbaddcolsOL",Form("%s - bad dcols; Chip; Stave",title.c_str()),196,-0.5,195.5,90,-0.5,89.5);
  TH2F *hBadPixIB = new TH2F("hbadpixIB",Form("%s - bad pix; Chip; Stave",title.c_str()),9,-0.5,8.5,48,-0.5,47.5);
  TH2F *hBadPixML = new TH2F("hbadpixML",Form("%s - bad pix; Chip; Stave",title.c_str()),112,-0.5,111.5,54,-0.5,53.5);
  TH2F *hBadPixOL = new TH2F("hbadpixOL",Form("%s - bad pix; Chip; Stave",title.c_str()),196,-0.5,195.5,90,-0.5,89.5);
  hBadDcolsIB->SetStats(0);hBadDcolsML->SetStats(0);hBadDcolsOL->SetStats(0);
  hBadPixIB->SetStats(0);hBadPixML->SetStats(0);hBadPixOL->SetStats(0);
  
  TFile *outTree = new TFile(Form("../Data/tree_dcols_P2_%s_%d.root",runnum.c_str(),id),"RECREATE");
  TTree mTreeP2("dcols_P2","dcols_P2");
  char longchipid[17];
  int dcolid, n_noisypix, n_goodpix, n_deadpix, n_ineffpix,o2chipid;
  //mTreeP2.Branch("chipid",longchipid,"longchipid[17]/C");
  mTreeP2.Branch("chipid",&o2chipid,"o2chipid/I");
  mTreeP2.Branch("dcol",&dcolid,"dcolid/I");
  mTreeP2.Branch("noisypix", &n_noisypix, "n_noisypix/I");
  mTreeP2.Branch("goodpix", &n_goodpix, "n_goodpix/I");
  mTreeP2.Branch("deadpix", &n_deadpix, "n_deadpix/I");
  mTreeP2.Branch("ineffpix", &n_ineffpix, "n_ineffpix/I");
  
  // format of long chipid: L0_00_{hs}_MOD{mod}_CXX
  cout<<"Filling the tree"<<endl;
  for(int ichip = 0; ichip<24120; ichip++){
    for(int i=0; i<512; i++){
      if(noisypix.count(ichip) || goodpix.count(ichip) || deadpix.count(ichip) || ineffpix.count(ichip)){
        mp.expandChipInfoHW(ichip, lay,sta, ssta, mod,  chipInMod);
        strncpy(longchipid, Form("L%d_%02d_%d_MOD%d_C%02d",lay,sta,ssta,mod,chipInMod), 17);
        dcolid = i;
	o2chipid=ichip;
        n_noisypix = (noisypix.count(ichip)) ? noisypix[ichip][i] : 0;
        n_goodpix = (goodpix.count(ichip)) ? goodpix[ichip][i] : 0;
        n_deadpix = (deadpix.count(ichip)) ? deadpix[ichip][i] : 0;
        n_ineffpix = (ineffpix.count(ichip)) ? ineffpix[ichip][i] : 0;
        mTreeP2.Fill();
      }
    }
  }

  /* TO BE REVISED
  for( const auto& pair : noisypix ) {
    int dcolshits[512] = { 0 };
    for(int i=0; i<(int)pair.second.size(); i++) {
      int dcol = ((pair.second[i] - pair.second[i]%1000) / 1000) / 2;
      dcolshits[dcol]++;
    }
    //count dcols 
    int dcolcounter = 0;
    for(int i=0; i<512; i++) {
      if(dcolshits[i] > 0) {
        //Fill the tree with the chipid and dcolID
        //mp.expandChipInfoHW(pair.first, lay,sta, ssta, mod,  chipInMod);
        //strncpy(longchipid, Form("L%d_%02d_%d_MOD%d_C%02d",lay,sta,ssta,mod,chipInMod), 17);
        //dcolid = i;
        //n_noisypix = dcolshits[i];
        //n_goodpix = goodpix[pair.first][i];
        //n_deadpix = deadpix[pair.first][i];
	//n_ineffpix = ineffpix[pair.first][i];
        //mTreeP2.Fill();
      }
      if(dcolshits[i] > 50) {
        dcolcounter++;
      }
    }
    //count single pixels not belonging to bad dcols
    int noisypixcounter = 0; 
    for(int i=0; i<(int)pair.second.size(); i++) {
      int dcol = ((pair.second[i] - pair.second[i]%1000) / 1000) / 2;
      if(dcolshits[i] <= 50) noisypixcounter++;
    }
    //fill plot with bad dcols per chip 
    mp.expandChipInfoHW(pair.first, lay,sta, ssta, mod,  chipInMod);
    if(lay<3){
      int stabin = !lay ? sta : lay==1 ? sta+nStaves[lay-1] : sta+nStaves[lay-2]+nStaves[lay-1];
      hBadDcolsIB->SetBinContent(chipInMod+1, stabin+1, dcolcounter);
      hBadPixIB->SetBinContent(chipInMod+1, stabin+1, noisypixcounter);
    }
    else{
      if(lay<5){
        int chipbin = chipInMod<7 ? (chipInMod+1)+14*(mod-1)+ssta*56 : chipInMod+14*(mod-1)+ssta*56;
        int stabin = lay==3 ? sta : sta+nStaves[lay-1];
        hBadDcolsML->SetBinContent(chipbin, stabin+1, dcolcounter);
        hBadPixML->SetBinContent(chipbin, stabin+1, noisypixcounter);
      }
      else {
        int chipbin = chipInMod<7 ? (chipInMod+1)+14*(mod-1)+ssta*98 : chipInMod+14*(mod-1)+ssta*98;
        int stabin = lay==5 ? sta : sta+nStaves[lay-1];
        hBadDcolsOL->SetBinContent(chipbin, stabin+1, dcolcounter);
        hBadPixOL->SetBinContent(chipbin, stabin+1, noisypixcounter);
      }
    }
  } // end loop on noisy pix map
  */
  mTreeP2.Write();
  outTree->Close();
  delete outTree;

  //calculate distributions of number of rows
  for(int ibin=1; ibin<=hRowsCountIB->GetNbinsX(); ibin++)
    if(hRowsCountIB->GetBinContent(ibin)>0) hRowsCountDistIB->Fill(hRowsCountIB->GetBinContent(ibin)); //to do the merging part properly
  for(int ibin=1; ibin<=hRowsCountML->GetNbinsX(); ibin++)
    if(hRowsCountML->GetBinContent(ibin)>0) hRowsCountDistML->Fill(hRowsCountML->GetBinContent(ibin));
  for(int ibin=1; ibin<=hRowsCountOL->GetNbinsX(); ibin++)
    if(hRowsCountOL->GetBinContent(ibin)>0) hRowsCountDistOL->Fill(hRowsCountOL->GetBinContent(ibin));

  TFile *outROOT = new TFile(Form("../Plots/out_digscan_%d.root",id),"RECREATE");
  hHitMapStave->Write();
  hBadPixIB->Write();
  hBadPixML->Write();
  hBadPixOL->Write();
  hBadDcolsIB->Write();
  hBadDcolsML->Write();
  hBadDcolsOL->Write();
  hRowsCountIB->Write();
  hRowsCountML->Write();
  hRowsCountOL->Write();
  hRowsCountDistIB->Write();
  hRowsCountDistML->Write();
  hRowsCountDistOL->Write();
  outROOT->Close();
  delete outROOT;
 
}

short int chipToLayer(int chipid){
  for(short int il=0; il<7; il++){
    if(chipid>=chipIDmin[il] && chipid<=chipIDmax[il]){
      return il;
    }
  }
  return -1;
}
