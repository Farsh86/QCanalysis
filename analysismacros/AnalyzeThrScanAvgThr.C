#include <string>
#include <iostream>
#include <vector>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>
#include <TGraph.h>

using namespace std;

void SetStyle(TH1 *h, Int_t col);
void DoAnalysis(string filepath, const int nChips, bool isIB, int layernum);

//
// MAIN
//
void AnalyzeThrScanAvgThr(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  string layernum = fpath.substr(fpath.find("Layer")+5, 1);

  if(stoi(layernum)>=0 && stoi(layernum)<=2) nchips = 9;
  else if (stoi(layernum)==3 && stoi(layernum)==4) nchips = 54*2;
  else nchips = 98*2;

  bool isIB;
  if(nchips==9) isIB=kTRUE;
  else isIB=kFALSE;

  //Call
  DoAnalysis("../Data/"+fpath, nchips, isIB, stoi(layernum));

}

//
//Set Style
//
void SetStyle(TGraph *gr, Int_t col){
  gr->SetLineColor(col);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.4);
  gr->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, bool isIB, int layernum){

  gStyle->SetOptStat(0000);

  //std::vector<TH2S*> hmaps;
  //std::vector<string> timestamps, runnumbers;
  int nTimes=0;
  Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  //Read the tree containing the average thresholds for each chip
  TFile *infl = new TFile(filepath.c_str());

  TTree *tr = (TTree*)infl->Get("thrscan_thr");
  int run, stave, chip;
  float avgthrchip;
  int numofstaves = 6;
  switch(layernum){
    case 0:
      numofstaves = 6;
      break;
    case 1:
      numofstaves = 8;
      break;
    case 2:
      numofstaves = 10;
      break;
    default:
      break;
  }
  TGraph *gr_stave[numofstaves];
  for(int istave=0; istave<numofstaves; istave++)
    gr_stave[istave] = new TGraph();

  tr->SetBranchAddress("runnum", &run);
  tr->SetBranchAddress("stavenum", &stave);
  tr->SetBranchAddress("chipnum", &chip);
  tr->SetBranchAddress("avgchipthr", &avgthrchip);
  Long64_t nentries = tr->GetEntries();
  int prevstave, prevrun, istave=0, irun=0, ix=0;
  float avgthrstave = 0;
  float minthr=1e5, maxthr=-1e5;// minimum and maximum threshold to set the correct scale of the plot
  vector<int> allruns, allstavenum;
  for(Long64_t i=0; i<nentries; i++){
    tr->GetEntry(i);
    //int pos = TMath::Floor((float)i/((float)nChips*(float)numofstaves));
    if(!i){
      prevstave = stave;
      prevrun = run;
      allruns.push_back(run);
      allstavenum.push_back(stave);
    }
    if(stave!=prevstave && run==prevrun){//change of stave but not of run!
      prevstave=stave;
      avgthrstave/=(float)nChips;
      if(avgthrstave>maxthr) maxthr=avgthrstave;
      else if(avgthrstave<minthr) minthr=avgthrstave;
      gr_stave[istave] -> SetPoint(irun, gr_stave[istave]->GetN(), avgthrstave);
      avgthrstave=0.;
      istave++;
      if(i < nChips*numofstaves) allstavenum.push_back(stave);
    }
    else if(stave!=prevstave && run!=prevrun){//change of stave and of run!
      prevstave=stave;
      prevrun=run;
      avgthrstave/=(float)nChips;
      if(avgthrstave>maxthr) maxthr=avgthrstave;
      else if(avgthrstave<minthr) minthr=avgthrstave;
      gr_stave[istave] -> SetPoint(irun, gr_stave[istave]->GetN(), avgthrstave);
      avgthrstave=0.;
      irun++;
      istave=0;
      allruns.push_back(run);
    }
    else if(i==nentries-1){//for the last entry
      avgthrstave+=avgthrchip; //sum the last number in the tree
      avgthrstave/=(float)nChips;
      if(avgthrstave>maxthr) maxthr=avgthrstave;
      else if(avgthrstave<minthr) minthr=avgthrstave;
      gr_stave[istave] -> SetPoint(irun, gr_stave[istave]->GetN(), avgthrstave);
    }
    avgthrstave+=avgthrchip;
    //cout<<run<<"  "<<stave<<"  "<<chip<<"  "<<avgthrchip<<endl;
  }

  //Set style (for OB needs to be modified)
  for(int istave=0; istave<numofstaves; istave++)
    SetStyle(gr_stave[istave], col[istave]);

  //Set labels creating a fake histogram
  int npoints = gr_stave[0]->GetN();
  TH1F *hfake = new TH1F("hfake", "; Run; Avg. Stave Threshold (DAC)", npoints, -0.5, (double)npoints-0.5);
  for(int ir=0; ir<(int)allruns.size(); ir++)
      hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", allruns[ir]));


  //Draw
  TCanvas *canvas = new TCanvas();
  canvas->cd();
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetMargin(0.0988,0.1,0.194,0.0993);
  TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
  for(int istave=0; istave<numofstaves; istave++)
    leg->AddEntry(gr_stave[istave], Form("Stv%d",allstavenum[istave]), "lp");
  hfake->Draw();
  hfake->GetYaxis()->SetRangeUser(minthr-0.05*minthr, maxthr+0.05*maxthr);
  hfake->GetXaxis()->SetTitleOffset(2.8);
  hfake->SetTitle(Form("%s Layer-%d, %s",isIB?"IB":"OB",layernum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  for(int istave=0; istave<numofstaves; istave++)
    gr_stave[istave]->Draw("P L same");
  leg->Draw();
  canvas->SaveAs(Form("../Plots/%s_layer%d_thresholds_vs_run_%s.pdf", isIB?"IB":"OB", layernum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas->SaveAs(Form("../Plots/%s_layer%d_thresholds_vs_run_%s.root", isIB?"IB":"OB", layernum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
}
