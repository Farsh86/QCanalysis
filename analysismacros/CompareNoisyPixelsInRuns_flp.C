#include <string>
#include <iostream>
#include <vector>
#include <array>
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
#include <TGraphErrors.h>
#include <TLine.h>
#include <TText.h>
#include <TSystem.h>
#include <TTree.h>

using namespace std;

//Functions
std::array<long int,4> CompareTwoRuns(vector<array<int,2>> run1/*ref run*/, vector<array<int,2>> run2);
void SetStyle(TGraphErrors *ge, Color_t col);
void DoAnalysis(string filepath, const int nChips, bool isIB, long int refrun, int layernum);
array<int,6> Preparestavearray(const int numofstaves, vector<array<int,2>> noisypix);
void SetStyle2(TGraph *gr, Int_t col);

void CompareNoisyPixelsInRuns_flp(){
  string fpath;
  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
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

  cout<<"Available (non-empty) runs in your file:\n"<<endl;
  int run;
  int prevrun;
  TFile *infile=new TFile(Form("../Data/%s",fpath.c_str()));
  TTree *tr = (TTree*)infile->Get("fhitscan");
  tr->SetBranchAddress("runnum", &run);
  Long64_t nentries = tr->GetEntries();

  for(Long64_t i=0; i<nentries; i++){
    tr->GetEntry(i);
    if(!i){
      prevrun=run;
    }
    if(prevrun!=run){
      cout<<prevrun<<endl;
      prevrun=run;
    }
    if(i==nentries-1)
      cout<<run<<endl;
  }


  long int runnum;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>runnum;

  DoAnalysis("../Data/"+fpath, nchips, isIB, runnum, stoi(layernum));
}

//
// Analysis
//
void DoAnalysis(string filepath, const int nChips, bool isIB, long int refrun, int layernum){

  gStyle->SetOptStat(0000);

  std::vector<string> timestamps, runnumbers;
  std::vector<array<long int,4>> noisypixinfo;
  Int_t color[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  //Read the file and the list of plots with entries
  TFile *infl = new TFile(filepath.c_str());
  TTree *tr = (TTree*)infl->Get("fhitscan");
  int run, stave, chip, col, row;
  float hits;
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

  tr->SetBranchAddress("runnum", &run);
  tr->SetBranchAddress("stavenum", &stave);
  tr->SetBranchAddress("chipnum", &chip);
  tr->SetBranchAddress("row", &row);
  tr->SetBranchAddress("col", &col);
  tr->SetBranchAddress("hits", &hits);
  Long64_t nentries = tr->GetEntries();
  vector<vector<array<int, 2>>> noisypixperrun;
  int prevrun;
  vector<array<int,2>> noisypix;
  array<int,2> colrow;
  vector<int> allruns;
  for(Long64_t i=0; i<nentries; i++){//loop on tree entries
    tr->GetEntry(i);
    //int pos = TMath::Floor((float)i/((float)nChips*(float)numofstaves));
    if(!i){
      prevrun = run;
    }
    else if(run!=prevrun){//change of stave and of run!
      allruns.push_back(prevrun);
      prevrun=run;
      noisypixperrun.push_back(noisypix);
      noisypix.clear();
    }
    else if(i==nentries-1){//for the last entry
      prevrun=run;
      noisypixperrun.push_back(noisypix);
      noisypix.clear();
      allruns.push_back(run);
    }
    if(hits/(300*50000)<1e-5) continue; //temporary cut
    colrow[0]=col;
    colrow[1]=row;
    if(hits>1) noisypix.push_back(colrow); // >1 is to reduce contribution from cosmics
  }

  //find position of the ref run
  int posrefrun=-1;
  for(int irun=0; irun<(int)allruns.size(); irun++){
    if(refrun==allruns[irun]){
      posrefrun=irun;
      break;
    }
  }

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  vector<long int> runlabel;
  for(int irun=0; irun<(int)allruns.size(); irun++){
    if(allruns[irun]==refrun) continue; // do not compare refrun with itself
    noisypixinfo.push_back(CompareTwoRuns(noisypixperrun[posrefrun], noisypixperrun[irun]));
    noisypixinfo[noisypixinfo.size()-1][0] = allruns[irun];//save timestamp on run2 (the ref one is known)
    runlabel.push_back(allruns[irun]);
    //cout<<allruns[irun]<<endl;
    //cout<<noisypixinfo[noisypixinfo.size()-1][0]<<"  "<<noisypixinfo[noisypixinfo.size()-1][1]<<"  "<<noisypixinfo[noisypixinfo.size()-1][2]<<"  "<<noisypixinfo[noisypixinfo.size()-1][3]<<endl;
  }

  //////////////////////////////////////////////////
  ///////////Make plot with #dead pix vs run////////
  //////////////////////////////////////////////////
  const int NSTAVES = 6; //FIXME
  TGraph *grnoisyperrun[NSTAVES];
  array<int,NSTAVES> noisypixperstave;
  for(int is=0; is<numofstaves; is++)
    grnoisyperrun[is]=new TGraph();

  double minnoi=1e20, maxnoi=-1;
  for(int irun=0; irun<(int)allruns.size(); irun++){
    noisypixperstave=Preparestavearray(NSTAVES, noisypixperrun[irun]);
    for(int is=0; is<numofstaves; is++){
      if(noisypixperstave[is]<minnoi) minnoi=noisypixperstave[is];
      if(noisypixperstave[is]>maxnoi) maxnoi=noisypixperstave[is];
      grnoisyperrun[is]->SetPoint(irun, irun, noisypixperstave[is]);
    }
  }
  //Style
  for(int is=0; is<numofstaves; is++)
    SetStyle2(grnoisyperrun[is], color[is]);

  //Set labels creating a fake histogram
  int npoints1 = grnoisyperrun[0]->GetN();
  TH1F *hfake1 = new TH1F("hfake1", "; Run; #Noisy pixels", npoints1, -0.5, (double)npoints1-0.5);
  for(int ir=0; ir<(int)allruns.size(); ir++)
      hfake1->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", allruns[ir]));


  //Draw
  TCanvas *canvas1 = new TCanvas();
  canvas1->cd();
  canvas1->SetTickx();
  canvas1->SetTicky();
  //canvas1->SetLogy();
  canvas1->SetMargin(0.0988,0.1,0.194,0.0993);
  TLegend *leg1 = new TLegend(0.904, 0.197,0.997,0.898);
  for(int istave=0; istave<numofstaves; istave++)
    leg1->AddEntry(grnoisyperrun[istave], Form("Stv%d",istave+numofstaves), "lp");
  hfake1->Draw();
  hfake1->GetYaxis()->SetRangeUser(minnoi-0.2*minnoi, maxnoi+0.2*maxnoi);
  hfake1->GetXaxis()->SetTitleOffset(2.8);
  hfake1->GetYaxis()->SetTitleOffset(1.4);
  hfake1->SetTitle(Form("%s Layer-%d - Noisy pixels per run, %s",isIB?"IB":"OB",layernum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  for(int istave=0; istave<numofstaves; istave++)
    grnoisyperrun[istave]->Draw("P L same");
  leg1->Draw();
  canvas1->SaveAs(Form("../Plots/%s_layer%d_noisypix_vs_run_%s.pdf", isIB?"IB":"OB", layernum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas1->SaveAs(Form("../Plots/%s_layer%d_noisypix_vs_run_%s.root", isIB?"IB":"OB", layernum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));


  ////////////////////////////////////////////
  ////////////////Make plot with bars/////////
  ////////////////////////////////////////////
  TGraphErrors *ge_nref = new TGraphErrors();
  TGraphErrors *ge_n2 = new TGraphErrors();
  TGraphErrors *ge_ncom1 = new TGraphErrors();
  TGraphErrors *ge_ncom2 = new TGraphErrors();
  double xshift = 3.;
  double max = -1;
  double min = 1e35;

  for(int icomp=0; icomp<(int)noisypixinfo.size(); icomp++){//first the older data and last the most recent
    //first couple of bar on the left
    int ipoint = icomp;
    if(!ipoint) xshift=1.;
    else xshift = 3.;
    ge_nref->SetPoint(ipoint, ipoint*xshift, (double)noisypixinfo[icomp][3]/2.+(double)noisypixinfo[icomp][1]/2.);
    ge_nref->SetPointError(ipoint, 0.5, (double)noisypixinfo[icomp][1]/2.);
    ge_ncom1->SetPoint(ipoint, ipoint*xshift, 0.);
    ge_ncom1->SetPointError(ipoint, 0.5, (double)noisypixinfo[icomp][3]/2.);
    if((double)noisypixinfo[icomp][3]/2.+(double)noisypixinfo[icomp][1]/2.+(double)noisypixinfo[icomp][1]/2. > max) max = (double)noisypixinfo[icomp][3]/2.+(double)noisypixinfo[icomp][1]/2.+(double)noisypixinfo[icomp][1]/2.;

    //second couple of bar on the right
    ge_n2->SetPoint(ipoint, ipoint*xshift+1, -(double)noisypixinfo[icomp][3]/2.-(double)noisypixinfo[icomp][2]/2.);
    ge_n2->SetPointError(ipoint, 0.5, (double)noisypixinfo[icomp][2]/2.);
    ge_ncom2->SetPoint(ipoint, ipoint*xshift+1, 0.);
    ge_ncom2->SetPointError(ipoint, 0.5, (double)noisypixinfo[icomp][3]/2.);
    if(-(double)noisypixinfo[icomp][3]/2.-(double)noisypixinfo[icomp][2]/2.-(double)noisypixinfo[icomp][2]/2. < min) min = -(double)noisypixinfo[icomp][3]/2.-(double)noisypixinfo[icomp][2]/2.-(double)noisypixinfo[icomp][2]/2.;
  }

  //Style
  SetStyle(ge_nref, kBlue);
  SetStyle(ge_ncom1, kBlack);
  SetStyle(ge_ncom2, kBlack);
  SetStyle(ge_n2, kRed+2);

  //Draw
  TCanvas *canvas = new TCanvas("mycanvas", "mycanvas", 1300, 800);
  canvas->SetMargin(0.0092, 0.1271, 0.1759, 0.0996);
  canvas->cd();

  //fake histo (just for the axes)
  double x2,y2;
  ge_ncom2->GetPoint(ge_ncom2->GetN()-1, x2,y2);
  TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
  hfake->Draw();
  //canvas->SetLogy();
  hfake->SetTitle(Form("run%06ld compared to all",refrun));
  ge_nref->Draw("P E2 same");
  ge_ncom1->Draw("E2 same");
  ge_ncom2->Draw("E2 same");
  ge_n2->Draw("E2 same");
  hfake->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
  hfake->GetYaxis()->SetLabelColor(kWhite);
  hfake->GetYaxis()->SetTickLength(0.005);
  TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
  lineref->SetLineColor(kGray-1);
  lineref->SetLineStyle(2);
  lineref->Draw("same");

  //draw labels on x axis
  //TAxis *ax = ge_nref->GetHistogram()->GetXaxis();
  int counter = 0;
  for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
    hfake->GetXaxis()->SetBinLabel(k, Form("run%06ld", runlabel[counter]));
    counter++;
  }

  //draw legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref, "#splitline{#noisy pix}{ref. run only}", "f");
  leg->AddEntry(ge_n2, "#splitline{#noisy pix}{2nd run only}", "f");
  leg->AddEntry(ge_ncom1, "#splitline{#noisy pix}{both}");
  leg->Draw("same");

  canvas->SaveAs(Form("../Plots/%s-Layer%d_NoisyPixComparison_run%ld_compared_to_run_%s.pdf", isIB?"IB":"OB",layernum, refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas->SaveAs(Form("../Plots/%s-Layer%d_NoisyPixComparison_run%ld_compared_to_run_%s.root", isIB?"IB":"OB",layernum, refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

}

//
// Evaluate dead pix per stave on layer
//
array<int,6> Preparestavearray(const int numofstaves, vector<array<int,2>> noisypix){

  array<int,6> numnoi;
  for(int i=0; i<numofstaves; i++)
    numnoi[i]=0;

  for(int id=0; id<(int)noisypix.size(); id++){
    //check the rows
    int realrow=1000;
    int rown=noisypix[id][1];
    int count=0;
    while(realrow>512){
      realrow=rown-512*count;
      count++;
    }
    numnoi[count-1]++;
  }
  return numnoi;
}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,4> CompareTwoRuns(vector<array<int,2>> run1/*ref run*/, vector<array<int,2>> run2){

  std::array<long int,4> noisypixinfo = {0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int ipix=0; ipix<(int)run1.size(); ipix++){
    int col1 = run1[ipix][0];
    int row1 = run1[ipix][1];
    bool isthere = false;

    for(int ipixs=0; ipixs<(int)run2.size(); ipixs++){
      if(run2[ipixs][0]==col1 && run2[ipixs][1]==row1){//noisy in both runs
        noisypixinfo[3]++;
        isthere=true;
        break;
      }
    }
    if(!isthere){//noisy only in ref run (run1)
      noisypixinfo[1]++;
    }
  }

  for(int ipix=0; ipix<(int)run2.size(); ipix++){
    int col2 = run2[ipix][0];
    int row2 = run2[ipix][1];
    bool isthere = false;

    for(int ipixs=0; ipixs<(int)run1.size(); ipixs++){
      if(run1[ipixs][0]==col2 && run1[ipixs][1]==row2){//noisy in both runs
        isthere=true;
        break;
      }
    }
    if(!isthere){//noisy only in 2nd run (run2)
      noisypixinfo[2]++;
    }
  }

  return noisypixinfo;
}

//
// Style
//
void SetStyle(TGraphErrors *ge, Color_t col){
  ge->SetMarkerStyle(0);
  ge->SetMarkerColor(col);
  ge->SetFillColor(col);
  ge->SetLineColor(col);
}

//
//Set Style
//
void SetStyle2(TGraph *gr, Int_t col){
  gr->SetLineColor(col);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.4);
  gr->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}
