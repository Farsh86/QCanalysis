#include <fstream>
#include <iostream>
#include <map>
#include "TH1.h"
#include "TCanvas.h"
#include "string"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

using namespace std;

void MakeNoisyPixFraction(string filename, string runnum){
  gStyle->SetOptStat(0000);
  int nstaves[] = {12,16,20,24,30,42,48};
  int npix[] = {512*1024*9,512*1024*9,512*1024*9,512*1024*14*8,512*1024*14*8,512*1024*14*14,512*1024*14*14};
  int nstavestot=nstaves[0]+nstaves[1]+nstaves[2]+nstaves[3]+nstaves[4]+nstaves[5]+nstaves[6];
  ifstream infile(filename.c_str());
  string hs;
  int layer, stave, mod, chipinmod, row, col;
  std::map<std::string, int> totalnoisypix; // stave id --> noisy pix
  TH1D *hhits = new TH1D("hhits",Form("Percentage of noisy pixels per stave in ITS2 - Cosmic run %s - ITS2 framing 67 kHz - Recorded readout frames (ROF): 27.5 #times 10^{6} - Stave average thresholds: 100 #it{e}^{-};Stave ID;Noisy pixels (%%)",runnum.c_str()), nstavestot,0,(double)nstavestot);
  string intro;
  TLine *lineIBOB = new TLine(48,3e-4,48, 4e-2);
  lineIBOB->SetLineStyle(2);
  TLatex *lat = new TLatex();
  lat->SetTextFont(42);

  //arrows
  TArrow *ar50pix = new TArrow(0,0.00953,9,0.00953,0.01,"<|");
  //ar50pix->SetAngle(180);
  ar50pix->SetLineColor(kRed);
  ar50pix->SetFillColor(kRed);

  //box
  TBox *b1 = new TBox(57.59,0.013,92.6,0.033);
  b1->SetFillColor(kWhite);
  TBox *b2 = new TBox(4.15,0.0128,39.2,0.033);
  b2->SetFillColor(kWhite);
  TBox *b3 = new TBox(9.63,0.0083,23.27,0.011);
  b3->SetFillColor(kWhite);

  std::getline(infile,intro); //first line
  cout<<"Reading "<<intro<<endl;
  while(infile>>stave>>layer>>hs>>mod>>chipinmod>>row>>col){
    string staveid = std::to_string(stave)+"_"+std::to_string(layer);
    totalnoisypix[staveid]++;
  }
  infile.close();

  for (auto const& stv : totalnoisypix){
    int bin = 0;
    int l = std::stoi(stv.first.substr(stv.first.find("L")+1,1));
    int s = std::stoi(stv.first.substr(stv.first.find("_")+1, 2));
    if(!l){ // if it's L0
      bin = std::stoi(stv.first.substr(stv.first.find("_")+1, 2));
    }
    else{ // for all the other layers
      for(int ilay=0; ilay<l; ilay++)
        bin += nstaves[ilay];
      bin+=s;
    }
    hhits->Fill(bin, ((double)stv.second/(double)npix[l])*100.);
    //hhits->GetXaxis()->SetBinLabel(bin+1, Form("L%d_%02d",l,s));
    hhits->SetBinError(bin+1,1e-20);
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
  hhits->SetMarkerColor(kRed);
  hhits->Draw("PE");
  hhits->GetYaxis()->SetLabelSize(0.045);
  hhits->GetXaxis()->SetLabelSize(0.045);
  hhits->GetYaxis()->SetRangeUser(3e-4,4e-2);
  hhits->GetYaxis()->SetTitleSize(0.045);
  hhits->GetXaxis()->SetTitleSize(0.045);
  hhits->GetYaxis()->SetTitleOffset(0.8);
  hhits->GetXaxis()->SetTitleOffset(0.85);
  hhits->GetYaxis()->SetMoreLogLabels();
  b1->Draw("same");
  b2->Draw("same");
  lat->DrawLatex(4.37,0.022,"#splitline{Inner Barrel}{pixels with >10^{-2} hits/ROF}");
  lat->DrawLatex(57.7,0.022,"#splitline{Outer Barrel}{pixels with >10^{-6} hits/ROF}");
  lat->SetTextSize(0.03);
  //lat->DrawLatex(4.37,0.13,"#it{L0,L1,L2: 4.7#times10^{6} pixels/stave}");
  //lat->DrawLatex(57.7,0.13,"#it{L3,L4: 58.7#times10^{6} pixels/stave}");
  //lat->DrawLatex(57.7,0.103,"#it{L5,L6: 102.8#times10^{6} pixels/stave}");
  lineIBOB->Draw();
  ar50pix->Draw();
  lat->SetTextSize(0.05);
  //lat->DrawLatex(159,0.0146,"#bf{ALICE Performance}");
  lat->SetTextSize(0.035);
  lat->SetTextColor(kRed);

  lat->DrawLatex(10,0.009,"50 pixels/chip");

  c->SaveAs(Form("../Plots/noisypix_fraction_%s.pdf",runnum.c_str()));
}
