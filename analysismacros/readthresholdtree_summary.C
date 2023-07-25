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
#include "ITSMFTReconstruction/ChipMappingITS.h"

const int chipIDmin[7] = {0, 108, 252, 432,3120, 6480, 14712};
const int chipIDmax[7] = {107, 251, 431,3119, 6479, 14711, 24119};

short int chipToLayer(int chipid);

//MAIN
void readthresholdtree_summary(std::string treepath = "../Data/tree.root", std::string title = "THR .... runnum ....", std::string stavemap = "stave for which you want thr map"){

  //gStyle->SetOptStat(0000);
  //open file
  TFile *infile = TFile::Open(treepath.c_str());
  TTree *thtree = (TTree*)infile->Get("ITS_calib_tree");
  short int chipid[1024], row[1024], thr[1024];
  char noise[1024];
  bool success[1024];
  const short int nStaves[7] = {12,16,20,24,30,42,48};
  int runnum = std::stoi(treepath.substr(treepath.find_last_of('/')+1,6));
  std::unordered_map<int,std::array<int,5>> countsum;
  std::unordered_map<int,std::array<int,2>> countsuccess_nolimits;
  std::unordered_map<int,std::array<int,2>> countsuccess_highnoise;
  std::unordered_map<short int,std::vector<short int>> trackrow;

  TH2F *hChipOutRange = new TH2F("hChipOutRange","Hits out of injected range (all superimposed); Col; Row",1024,-0.5,1023.5,512,-0.5,511.5);
  TH1F *hThrLayer[7];
  for(int il=0; il<7; il++){
    hThrLayer[il] = new TH1F(Form("hThrLayer%d",il), Form("Thr dist per pixel - L%d; Threshold (e-); Counts",il), 120, 30, 400);
  }
  TH1F *hUnSuccess_nolimits = new TH1F("hUnSuccess_nolimits","Un-success rate per chip - S-curve with no limits; Chip ID; unsuccess rate (%)",24120,-0.5,24119.5);
  hUnSuccess_nolimits->SetStats(0);
  TH1F *hUnSuccess_highnoise = new TH1F("hUnSuccess_highnoise","Un-success rate per chip - S-curve with noise>= 15 e-; Chip ID; unsuccess rate (%)",24120,-0.5,24119.5);
  hUnSuccess_highnoise->SetStats(0);
  TH2F *hThrChip = new TH2F("hThrChip","Threshold distribution per chip; ChipID; Threshold (e-)", 24120, -0.5, 24119.5, 100,30,400);
  hThrChip->SetStats(0);
  TH2F *hNoiseChip = new TH2F("hNoiseChip","Noise distribution per chip; ChipID; Noise (e-)", 24120, -0.5, 24119.5, 30,0,40);
  hNoiseChip->SetStats(0);

  TH1F *hRowsCountIB = new TH1F("hRowsCountIB","Rows collected per chip in IB; Chip ID; Number of rows",432,-0.5,431.5);
  TH1F *hRowsCountML = new TH1F("hRowsCountML","Rows collected per chip in ML; Chip ID; Number of rows",6048,431.5,6479.5);
  TH1F *hRowsCountOL = new TH1F("hRowsCountOL","Rows collected per chip in OL; Chip ID; Number of rows",17640,6479.5,24119.5);
  TH1F *hRowsCountDistIB = new TH1F("hRowsCountDistIB","Rows collected per chip in IB; Number of rows; Counts",20,-0.5,19.5);
  TH1F *hRowsCountDistML = new TH1F("hRowsCountDistML","Rows collected per chip in ML; Number of rows; Counts",20,-0.5,19.5);
  TH1F *hRowsCountDistOL = new TH1F("hRowsCountDistOL","Rows collected per chip in OL; Number of rows; Counts",20,-0.5,19.5);
  hRowsCountIB->SetStats(0); hRowsCountML->SetStats(0); hRowsCountOL->SetStats(0);

  //stave map
  TH2F *hThrMapStave = 0;
  if(stavemap.find("L0")!=string::npos || stavemap.find("L1")!=string::npos || stavemap.find("L2")!=string::npos){
    hThrMapStave = new TH2F("hThrMapStave",Form("Thr map %s; Col; Row",stavemap.c_str()), 1024*9, -0.5, 1024*9-0.5, 512,-0.5,512-0.5);
  }
  else if(stavemap.find("L3")!=string::npos || stavemap.find("L4")!=string::npos){
    hThrMapStave = new TH2F("hThrMapStave",Form("Thr map %s; Col; Row",stavemap.c_str()), 1024*28, 0, 1024*28, 512*2,0,512*2);
  }
  else{
    hThrMapStave = new TH2F("hThrMapStave",Form("Thr map %s; Col; Row",stavemap.c_str()), 1024*49, 0, 1024*49, 512*2,0,512*2);
  }
  hThrMapStave->SetStats(0);

  thtree->SetBranchAddress("chipid",&chipid[0]);
  thtree->SetBranchAddress("row",&row[0]);
  thtree->SetBranchAddress("thr",&thr[0]);
  thtree->SetBranchAddress("noise",&noise[0]);
  thtree->SetBranchAddress("success",&success[0]);

  long int nentries = thtree->GetEntries();
  o2::itsmft::ChipMappingITS mp;
  int lay,sta,ssta,mod,chipInMod;
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
          if(crow==trackrow[cchipid][ir]) countthisrow++;
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

      //Fill threshold map of selected stave (NOTE: no cuts applied on the tree!!)
      if(staveid==stavemap){
        if(layer<3){
          //std::cout<<crow+1<<" "<<icol<<" "<<thr[icol]<<std::endl;
          hThrMapStave->SetBinContent(icol + 1 + 1024*chipInMod, crow+1, thr[icol]);
        }
        else {
          int colbin = chipInMod < 7 ? icol + 1 + chipInMod*1024 + (mod-1)*1024*7 : icol + 1 + (14-chipInMod)*1024 + (mod-1)*1024*7;
          int rowbin = 1;
          if(!ssta) rowbin = chipInMod < 7 ? crow+1 : crow+1+512;
          else rowbin = chipInMod < 7 ? crow+1+1024 : crow+1+1536;
          hThrMapStave->SetBinContent(colbin, rowbin, thr[icol]);
        }
      }

      //count unsuccess per chipid
      if(!success[icol] && thr[icol]<1e-9 && noise[icol]<1e-9){
        if(!countsuccess_nolimits.count(cchipid)){
          countsuccess_nolimits[cchipid][0] = 1; //un-success
          countsuccess_nolimits[cchipid][1] = 1; //count entries
        }
        else {
          countsuccess_nolimits[cchipid][0] ++; //un-success
          countsuccess_nolimits[cchipid][1] ++; //count entries
        }
      }
      else if(/*!success[icol] &&*/ (int)noise[icol]>=15){
        if(!countsuccess_highnoise.count(cchipid)){
          countsuccess_highnoise[cchipid][0] = 1; //un-success
          countsuccess_highnoise[cchipid][1] = 1; //count entries
        }
        else {
          countsuccess_highnoise[cchipid][0] ++; //un-success
          countsuccess_highnoise[cchipid][1] ++; //count entries
        }
      }
      else {

        if(!countsuccess_nolimits.count(cchipid)){
          countsuccess_nolimits[cchipid][1] = 1; //count entries
        }
        else {
          countsuccess_nolimits[cchipid][1] ++; //count entries
        }

        if(!countsuccess_highnoise.count(cchipid)){
          countsuccess_highnoise[cchipid][1] = 1; //count entries
        }
        else {
          countsuccess_highnoise[cchipid][1] ++; //count entries
        }

      }

      //Calculate sums for averages
      if(thr[icol]<1e-9 && (int)noise[icol]<1e-9){
        continue; //skip if didn't find upper and lower limits in the fit
      }
      //if(noise[icol]>=15){
        //cout<<"ChipID: "<<cchipid<<" - (row,col)="<<"("<<crow<<","<<icol<<") has noise = "<<(int)noise[icol]<<" electrons"<<endl;
      //}
      hThrLayer[layer]->Fill(thr[icol]);
      hThrChip->Fill(cchipid+1, thr[icol]);
      hNoiseChip->Fill(cchipid+1, (int)noise[icol]);
      if(!countsum.count(cchipid)){
        countsum[cchipid][0] = thr[icol];
        countsum[cchipid][1] = noise[icol];
        countsum[cchipid][4] = 1;
      }
      else {
        countsum[cchipid][0] += thr[icol];
        countsum[cchipid][1] += noise[icol];
        countsum[cchipid][4] ++;
      }
    }
  }

  //TH2 for IB
  TH2F *hAvgThr = new TH2F("hthrIB",Form("%s; Chip; Stave",title.c_str()),9,-0.5,8.5,48,-0.5,47.5);
  TH2F *hAvgNoi = new TH2F("hnoiIB",Form("%s; Chip; Stave)",title.c_str()),9,-0.5,8.5,48,-0.5,47.5);
  TH2F *hAvgThrML = new TH2F("hthrML",Form("%s; Chip; Stave",title.c_str()),112,-0.5,111.5,54,-0.5,53.5);
  TH2F *hAvgNoiML = new TH2F("hnoiML",Form("%s; Chip; Stave)",title.c_str()),112,-0.5,111.5,54,-0.5,53.5);
  TH2F *hAvgThrOL = new TH2F("hthrOL",Form("%s; Chip; Stave",title.c_str()),196,-0.5,195.5,90,-0.5,89.5);
  TH2F *hAvgNoiOL = new TH2F("hnoiOL",Form("%s; Chip; Stave)",title.c_str()),196,-0.5,195.5,90,-0.5,89.5);
  hAvgThr->SetStats(0);hAvgNoi->SetStats(0);hAvgThrML->SetStats(0);hAvgNoiML->SetStats(0);hAvgThrOL->SetStats(0);hAvgNoiOL->SetStats(0);

  //loop over map, calculate averages, convert chipid, fill plot
  for(auto const& [chipID, t_arr] : countsum){
    if(!t_arr[4]) continue;
    mp.expandChipInfoHW(chipID, lay,sta, ssta, mod,  chipInMod);
    if(lay<3){//IB
      int stabin = !lay ? sta : lay==1 ? sta+nStaves[lay-1] : sta+nStaves[lay-2]+nStaves[lay-1];
      hAvgThr->SetBinContent(chipInMod+1, stabin+1, t_arr[0] / t_arr[4]);
      hAvgNoi->SetBinContent(chipInMod+1, stabin+1, t_arr[1] / t_arr[4]);
    }
    else{
      if(lay<5){//ML
        int chipbin = chipInMod<7 ? (chipInMod+1)+14*(mod-1)+ssta*56 : chipInMod+14*(mod-1)+ssta*56;
        int stabin = lay==3 ? sta : sta+nStaves[lay-1];
        hAvgThrML->SetBinContent(chipbin, stabin+1, t_arr[0] / t_arr[4]);
        hAvgNoiML->SetBinContent(chipbin, stabin+1, t_arr[1] / t_arr[4]);
      }
      else {
        int chipbin = chipInMod<7 ? (chipInMod+1)+14*(mod-1)+ssta*98 : chipInMod+14*(mod-1)+ssta*98;
        int stabin = lay==5 ? sta : sta+nStaves[lay-1];
        hAvgThrOL->SetBinContent(chipbin, stabin+1, t_arr[0] / t_arr[4]);
        hAvgNoiOL->SetBinContent(chipbin, stabin+1, t_arr[1] / t_arr[4]);
      }
    }
  }

  /*TFile *outROOT = new TFile(Form("../Data/out_threshold_run%d.root",runnum),"RECREATE");
  hAvgThr->SetName(Form("hAvgThrIB_run%d_1600000000000",runnum));
  hAvgThrML->SetName(Form("hAvgThrML_run%d_1600000000001",runnum));
  hAvgThrOL->SetName(Form("hAvgThrOL_run%d_1600000000001",runnum));
  hAvgThr->Write();
  hAvgThrOL->Write();
  hAvgThrML->Write();
  outROOT->Close();
  delete outROOT;*/

  //Fill histograms showing averages per HIC (and not per chip)
  TH2F *hAvgThrML_hic = new TH2F("hthrML_hic",Form("%s; Module; Stave",title.c_str()),8,-0.5,7.5,54,-0.5,53.5);
  TH2F *hAvgNoiML_hic = new TH2F("hnoiML_hic",Form("%s; Module; Stave)",title.c_str()),8,-0.5,7.5,54,-0.5,53.5);
  TH2F *hAvgThrOL_hic = new TH2F("hthrOL_hic",Form("%s; Module; Stave",title.c_str()),14,-0.5,13.5,90,-0.5,89.5);
  TH2F *hAvgNoiOL_hic = new TH2F("hnoiOL_hic",Form("%s; Module; Stave)",title.c_str()),14,-0.5,13.5,90,-0.5,89.5);
  hAvgThrML_hic->SetStats(0);hAvgNoiML_hic->SetStats(0);hAvgThrOL_hic->SetStats(0);hAvgNoiOL_hic->SetStats(0);
  ///// ML
  double sumhicthr = 0., sumhicnoi = 0.;
  double countchips = 0.;
  int div = 0;
  for(int istave=0; istave<hAvgThrML->GetNbinsY(); istave++){
    for(int ichip=0; ichip<=hAvgThrML->GetNbinsX(); ichip++){
      if(ichip/14 != div){
        hAvgThrML_hic->SetBinContent(ichip/14,istave+1,!countchips ? 0. : sumhicthr/countchips);
        hAvgNoiML_hic->SetBinContent(ichip/14,istave+1,!countchips ? 0. : sumhicnoi/countchips);
        sumhicthr=0;
        sumhicnoi=0;
        countchips=0;
        div = ichip/14;
        if(div==8) break;
      }
      sumhicthr+=hAvgThrML->GetBinContent(ichip+1,istave+1);
      sumhicnoi+=hAvgNoiML->GetBinContent(ichip+1,istave+1);
      countchips+=hAvgThrML->GetBinContent(ichip+1,istave+1)<1e-10 ? 0. : 1.;
    }
  }
  sumhicthr=0.;
  sumhicnoi = 0.;
  countchips = 0.;
  div = 0;
  /////OL
  for(int istave=0; istave<hAvgThrOL->GetNbinsY(); istave++){
    for(int ichip=0; ichip<=hAvgThrOL->GetNbinsX(); ichip++){
      if(ichip/14 != div){
        hAvgThrOL_hic->SetBinContent(ichip/14,istave+1,!countchips ? 0. : sumhicthr/countchips);
        hAvgNoiOL_hic->SetBinContent(ichip/14,istave+1,!countchips ? 0. : sumhicnoi/countchips);
        sumhicthr=0;
        sumhicnoi=0;
        countchips=0;
        div = ichip/14;
        if(div==14) break;
      }
      sumhicthr+=hAvgThrOL->GetBinContent(ichip+1,istave+1);
      sumhicnoi+=hAvgNoiOL->GetBinContent(ichip+1,istave+1);
      countchips+=hAvgThrOL->GetBinContent(ichip+1,istave+1)<1e-10 ? 0. : 1.;
    }
  }

  //Fill histo showing unsuccess rate in the case: no limits
  for(auto const& [chipID, t_arr] : countsuccess_nolimits){
    hUnSuccess_nolimits->SetBinContent(chipID+1, ((float)t_arr[0] / (float)t_arr[1]) * 100.);
  }
  //Fill histo showing unsuccess rate in the case: high noise
  for(auto const& [chipID, t_arr] : countsuccess_highnoise){
    hUnSuccess_highnoise->SetBinContent(chipID+1, ((float)t_arr[0] / (float)t_arr[1]) * 100.);
  }

  //calculate distributions of number of rows
  for(int ibin=1; ibin<=hRowsCountIB->GetNbinsX(); ibin++)
    hRowsCountDistIB->Fill(hRowsCountIB->GetBinContent(ibin));
  for(int ibin=1; ibin<=hRowsCountML->GetNbinsX(); ibin++)
    hRowsCountDistML->Fill(hRowsCountML->GetBinContent(ibin));
  for(int ibin=1; ibin<=hRowsCountOL->GetNbinsX(); ibin++)
    hRowsCountDistOL->Fill(hRowsCountOL->GetBinContent(ibin));

  //Draw row counter plots
  for(int i=0; i<6; i++){
    TCanvas cRowCount("cRowCount","cRowCount", 1000, 700);
    switch(i){
      case 0: hRowsCountIB->Draw("HIST"); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf[",runnum)); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf",runnum)); break;
      case 1: hRowsCountML->Draw("HIST"); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf",runnum)); break;
      case 2: hRowsCountOL->Draw("HIST"); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf",runnum)); break;
      case 3: hRowsCountDistIB->Draw("HIST"); cRowCount.SetLogy(); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf",runnum)); break;
      case 4: hRowsCountDistML->Draw("HIST"); cRowCount.SetLogy(); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf",runnum)); break;
      case 5: hRowsCountDistOL->Draw("HIST"); cRowCount.SetLogy(); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf",runnum)); cRowCount.SaveAs(Form("../Plots/RowCount_run%d.pdf]",runnum)); break;
    }
  }

  //Draw threshold map of selected stave
  TCanvas cMapStave("cMapStave","cMapStave",1800, 400);
  hThrMapStave->Draw("colz");
  //if(stavemap.find("L0")!=string::npos || stavemap.find("L1")!=string::npos || stavemap.find("L2")!=string::npos){
    //hThrMapStave->GetYaxis()->SetRangeUser(4,26);
  //}
  cMapStave.SaveAs(Form("../Plots/%s_thr_map_run%d.png",stavemap.c_str(),runnum));

  //Draw unsuccess rate
  for(int i=0; i<2; i++){
    TCanvas cUnSucc("cUnSucc","cUnSucc",1800, 700);
    cUnSucc.SetLogy();
    i==0 ? hUnSuccess_nolimits->Draw("HIST") : hUnSuccess_highnoise->Draw("HIST");
    if(!i) cUnSucc.SaveAs(Form("../Plots/UnsuccessRate_run%d.pdf[",runnum));
    cUnSucc.SaveAs(Form("../Plots/UnsuccessRate_run%d.pdf",runnum));
    if(i==1) cUnSucc.SaveAs(Form("../Plots/UnsuccessRate_run%d.pdf]",runnum));
  }


  //Draw plots with threshold and noise distributions per chip id
  for(int i=0; i<2; i++){
    TCanvas cchips("cDistChipid","cDistChipid",1800, 700);
    cchips.SetLogz();
    i==0 ? hThrChip->Draw("COLZ") : hNoiseChip->Draw("COLZ");
    if(!i) cchips.SaveAs(Form("../Plots/Distributions_per_chip_run%d.pdf[",runnum));
    cchips.SaveAs(Form("../Plots/Distributions_per_chip_run%d.pdf",runnum));
    if(i==1) cchips.SaveAs(Form("../Plots/Distributions_per_chip_run%d.pdf]",runnum));
  }

  //lines for layers
  TLine *l01 = new TLine(-0.5,11.5,8.5,11.5);
  TLine *l12 = new TLine(-0.5,27.5,8.5,27.5);
  TLine *l34 = new TLine(-0.5,23.5,111.5,23.5);
  TLine *l56 = new TLine(-0.5,41.5,197.5,41.5);
  TLatex *lat = new TLatex();
  TLatex *lattitle = new TLatex();
  lat->SetTextFont(42);
  lattitle->SetTextFont(42);
  lattitle->SetTextAngle(90);
  lattitle->SetTextSize(0.04);
  lattitle->SetNDC();

  TCanvas *cTHRIB = new TCanvas("cTHRIB","cTHRIB",1100,1800);
  cTHRIB->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgThr->Draw("colz");
  hAvgThr->SetMaximum(170);
  hAvgThr->SetMinimum(50);
  hAvgThr->GetXaxis()->SetTickLength(0.01);
  hAvgThr->GetYaxis()->SetTickLength(0.02);
  l01->Draw();
  l12->Draw();
  lat->DrawLatex(0,9,"L0");
  lat->DrawLatex(0,25,"L1");
  lat->DrawLatex(0,45,"L2");
  lattitle->DrawLatex(0.98,0.76,"Threshold (#it{e^{-}})");
  cTHRIB->SaveAs(Form("../Plots/Threshold_IB_run%d.pdf", runnum));


  TCanvas *cNOIIB = new TCanvas("cNOIIB","cNOIIB",1100,1800);
  cNOIIB->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgNoi->Draw("colz");
  hAvgNoi->SetMaximum(7);
  hAvgNoi->SetMinimum(1);
  hAvgNoi->GetXaxis()->SetTickLength(0.01);
  hAvgNoi->GetYaxis()->SetTickLength(0.02);
  l01->Draw();
  l12->Draw();
  lat->DrawLatex(0,9,"L0");
  lat->DrawLatex(0,25,"L1");
  lat->DrawLatex(0,45,"L2");
  lattitle->DrawLatex(0.97,0.8,"Noise (#it{e^{-}})");
  cNOIIB->SaveAs(Form("../Plots/Noise_IB_run%d.pdf", runnum));

  ///
  /// ML
  ///
  TCanvas *cTHRML = new TCanvas("cTHRML","cTHRML",1100,1800);
  cTHRML->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgThrML->Draw("colz");
  hAvgThrML->SetMaximum(280);
  hAvgThrML->SetMinimum(50);
  hAvgThrML->GetXaxis()->SetTickLength(0.01);
  hAvgThrML->GetYaxis()->SetTickLength(0.02);
  l34->Draw();
  lat->DrawLatex(1,20,"L3");
  lat->DrawLatex(1,50,"L4");
  lattitle->DrawLatex(0.98,0.76,"Threshold (#it{e^{-}})");
  cTHRML->SaveAs(Form("../Plots/Threshold_ML_run%d.pdf", runnum));


  TCanvas *cNOIML = new TCanvas("cNOIML","cNOIML",1100,1800);
  cNOIML->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgNoiML->Draw("colz");
  hAvgNoiML->SetMaximum(7);
  hAvgNoiML->SetMinimum(1);
  hAvgNoiML->GetXaxis()->SetTickLength(0.01);
  hAvgNoiML->GetYaxis()->SetTickLength(0.02);
  l34->Draw();
  lat->DrawLatex(1,20,"L3");
  lat->DrawLatex(1,50,"L4");
  lattitle->DrawLatex(0.97,0.8,"Noise (#it{e^{-}})");
  cNOIML->SaveAs(Form("../Plots/Noise_ML_run%d.pdf", runnum));

  TCanvas *cTHRML_hic = new TCanvas("cTHRML_hic","cTHRML_hic",1100,1800);
  cTHRML_hic->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgThrML_hic->Draw("colz");
  hAvgThrML_hic->SetMaximum(280);
  hAvgThrML_hic->SetMinimum(50);
  hAvgThrML_hic->GetXaxis()->SetTickLength(0.01);
  hAvgThrML_hic->GetYaxis()->SetTickLength(0.02);
  l34->Draw();
  lat->DrawLatex(1,20,"L3");
  lat->DrawLatex(1,50,"L4");
  lattitle->DrawLatex(0.98,0.76,"Threshold (#it{e^{-}})");
  cTHRML_hic->SaveAs(Form("../Plots/Threshold_ML_HIC_run%d.pdf", runnum));


  TCanvas *cNOIML_hic = new TCanvas("cNOIML_hic","cNOIML_hic",1100,1800);
  cNOIML_hic->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgNoiML_hic->Draw("colz");
  hAvgNoiML_hic->SetMaximum(7);
  hAvgNoiML_hic->SetMinimum(1);
  hAvgNoiML_hic->GetXaxis()->SetTickLength(0.01);
  hAvgNoiML_hic->GetYaxis()->SetTickLength(0.02);
  l34->Draw();
  lat->DrawLatex(1,20,"L3");
  lat->DrawLatex(1,50,"L4");
  lattitle->DrawLatex(0.97,0.8,"Noise (#it{e^{-}})");
  cNOIML_hic->SaveAs(Form("../Plots/Noise_ML_HIC_run%d.pdf", runnum));

  ///
  /// OL
  ///
  TCanvas *cTHROL = new TCanvas("cTHROL","cTHROL",1100,1800);
  cTHROL->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgThrOL->Draw("colz");
  hAvgThrOL->SetMaximum(280);
  hAvgThrOL->SetMinimum(50);
  hAvgThrOL->GetXaxis()->SetTickLength(0.01);
  hAvgThrOL->GetYaxis()->SetTickLength(0.02);
  l56->Draw();
  lat->DrawLatex(1,35,"L5");
  lat->DrawLatex(1,84,"L6");
  lattitle->DrawLatex(0.98,0.76,"Threshold (#it{e^{-}})");
  cTHROL->SaveAs(Form("../Plots/Threshold_OL_run%d.pdf", runnum));


  TCanvas *cNOIOL = new TCanvas("cNOIOL","cNOIOL",1100,1800);
  cNOIOL->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgNoiOL->Draw("colz");
  hAvgNoiOL->SetMaximum(7);
  hAvgNoiOL->SetMinimum(1);
  hAvgNoiOL->GetXaxis()->SetTickLength(0.01);
  hAvgNoiOL->GetYaxis()->SetTickLength(0.02);
  l56->Draw();
  lat->DrawLatex(1,35,"L5");
  lat->DrawLatex(1,84,"L6");
  lattitle->DrawLatex(0.97,0.8,"Noise (#it{e^{-}})");
  cNOIOL->SaveAs(Form("../Plots/Noise_OL_run%d.pdf", runnum));

  TCanvas *cTHROL_hic = new TCanvas("cTHROL_hic","cTHROL_hic",1100,1800);
  cTHROL_hic->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgThrOL_hic->Draw("colz");
  hAvgThrOL_hic->SetMaximum(280);
  hAvgThrOL_hic->SetMinimum(50);
  hAvgThrOL_hic->GetXaxis()->SetTickLength(0.01);
  hAvgThrOL_hic->GetYaxis()->SetTickLength(0.02);
  l56->Draw();
  lat->DrawLatex(1,35,"L5");
  lat->DrawLatex(1,84,"L6");
  lattitle->DrawLatex(0.98,0.76,"Threshold (#it{e^{-}})");
  cTHROL_hic->SaveAs(Form("../Plots/Threshold_OL_HIC_run%d.pdf", runnum));


  TCanvas *cNOIOL_hic = new TCanvas("cNOIOL_hic","cNOIOL_hic",1100,1800);
  cNOIOL_hic->SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgNoiOL_hic->Draw("colz");
  hAvgNoiOL_hic->SetMaximum(7);
  hAvgNoiOL_hic->SetMinimum(1);
  hAvgNoiOL_hic->GetXaxis()->SetTickLength(0.01);
  hAvgNoiOL_hic->GetYaxis()->SetTickLength(0.02);
  l56->Draw();
  lat->DrawLatex(1,35,"L5");
  lat->DrawLatex(1,84,"L6");
  lattitle->DrawLatex(0.97,0.8,"Noise (#it{e^{-}})");
  cNOIOL_hic->SaveAs(Form("../Plots/Noise_OL_HIC_run%d.pdf", runnum));

  //Draw thrshold distributions
  for(int il=0; il<7; il++){
    TCanvas cThrLay;
    hThrLayer[il]->Draw("HIST");
    //cThrLay.SetLogy();
    if(!il) cThrLay.SaveAs(Form("../Plots/ThresholdDist_run%d.pdf[",runnum));
    cThrLay.SaveAs(Form("../Plots/ThresholdDist_run%d.pdf",runnum));
    if(il==6) cThrLay.SaveAs(Form("../Plots/ThresholdDist_run%d.pdf]",runnum));
  }

}

short int chipToLayer(int chipid){
  for(short int il=0; il<7; il++){
    if(chipid>=chipIDmin[il] && chipid<=chipIDmax[il]){
      return il;
    }
  }
  return -1;
}
