#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void MakePlotITHRVCASN(string sourcefile = "file.txt", bool isVCASN = kTRUE) {
  const short int nStaves[7] = {12,16,20,24,30,42,48};
  TH2F *hAvgIB = new TH2F("hAvgIB",Form("%s; Chip; Stave",isVCASN ? "VCASN":"ITHR"),9,-0.5,8.5,48,-0.5,47.5);
  TH2F *hAvgML = new TH2F("hAvgML",Form("%s; Chip; Stave",isVCASN ? "VCASN":"ITHR"),112,-0.5,111.5,54,-0.5,53.5);
  TH2F *hAvgOL = new TH2F("hAvgOL",Form("%s; Chip; Stave",isVCASN ? "VCASN":"ITHR"),196,-0.5,195.5,90,-0.5,89.5);
  hAvgIB->SetStats(0); hAvgML->SetStats(0); hAvgOL->SetStats(0);

  TH1F *hDistLayer[7];
  for(int i=0; i<7; i++){
    hDistLayer[i] = new TH1F(Form("hDistLayer%s_L%d",isVCASN ? "VCASN":"ITHR", i),Form("Layer %d - %s distribution chip by chip; %s (DAC); Counts",i,isVCASN ? "VCASN":"ITHR", isVCASN ? "VCASN":"ITHR"), 60,0,120);
  }

  //read file and fill plots
  ifstream infl(sourcefile.c_str());
  string intro;
  getline(infl, intro, '\n');
  int layer, stave, mod, chipInMod, value; //value can be ITHR or VCASN
  string hs;
  while(infl>>layer>>stave>>hs>>mod>>chipInMod>>value){
    hDistLayer[layer]->Fill(value);
    int ssta = hs=="L" ? 0:1;
    if(layer<3){//IB
      int stabin = !layer ? stave : layer==1 ? stave+nStaves[layer-1] : stave+nStaves[layer-2]+nStaves[layer-1];
      hAvgIB->SetBinContent(chipInMod+1, stabin+1, value);
    }
    else{
      if(layer<5){//ML
        int chipbin = chipInMod<7 ? (chipInMod+1)+14*(mod-1)+ssta*56 : chipInMod+14*(mod-1)+ssta*56;
        int stabin = layer==3 ? stave : stave+nStaves[layer-1];
        hAvgML->SetBinContent(chipbin, stabin+1, value);
      }
      else {
        int chipbin = chipInMod<7 ? (chipInMod+1)+14*(mod-1)+ssta*98 : chipInMod+14*(mod-1)+ssta*98;
        int stabin = layer==5 ? stave : stave+nStaves[layer-1];
        hAvgOL->SetBinContent(chipbin, stabin+1, value);
      }
    }
  }

  //tlines
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


  TCanvas cnv1("cnv1","cnv1",1100,1800);
  cnv1.SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgIB->Draw("colz");
  hAvgIB->SetMaximum(70);
  hAvgIB->SetMinimum(40);
  l01->Draw();
  l12->Draw();
  lat->DrawLatex(0,9,"L0");
  lat->DrawLatex(0,25,"L1");
  lat->DrawLatex(0,45,"L2");
  lattitle->DrawLatex(0.98,0.76,Form("%s (DAC)",isVCASN ? "VCASN":"ITHR"));
  cnv1.SaveAs(Form("%s_plots.pdf[", isVCASN ? "VCASN":"ITHR"));
  cnv1.SaveAs(Form("%s_plots.pdf", isVCASN ? "VCASN":"ITHR"));

  TCanvas cnv2("cnv2","cnv2",1100,1800);
  cnv2.SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgML->Draw("colz");
  hAvgML->SetMaximum(70);
  hAvgML->SetMinimum(40);
  hAvgML->GetXaxis()->SetTickLength(0.01);
  hAvgML->GetYaxis()->SetTickLength(0.02);
  l34->Draw();
  lat->DrawLatex(1,20,"L3");
  lat->DrawLatex(1,50,"L4");
  lattitle->DrawLatex(0.98,0.76,Form("%s (DAC)",isVCASN ? "VCASN":"ITHR"));
  cnv2.SaveAs(Form("%s_plots.pdf", isVCASN ? "VCASN":"ITHR"));

  TCanvas cnv3("cnv3","cnv3",1100,1800);
  cnv3.SetMargin(0.099, 0.1632, 0.1007, 0.0999);
  hAvgOL->Draw("colz");
  hAvgOL->SetMaximum(70);
  hAvgOL->SetMinimum(40);
  hAvgOL->GetXaxis()->SetTickLength(0.01);
  hAvgOL->GetYaxis()->SetTickLength(0.02);
  l56->Draw();
  lat->DrawLatex(1,35,"L5");
  lat->DrawLatex(1,84,"L6");
  lattitle->DrawLatex(0.98,0.76,Form("%s (DAC)",isVCASN ? "VCASN":"ITHR"));
  cnv3.SaveAs(Form("%s_plots.pdf", isVCASN ? "VCASN":"ITHR"));

  for(int i=0; i<7; i++){
    TCanvas cnv4("cnv4","cnv4");
    cnv4.SetLogy();
    hDistLayer[i]->Draw("HIST");
    cnv4.SaveAs(Form("%s_plots.pdf", isVCASN ? "VCASN":"ITHR"));
    if(i==6) cnv4.SaveAs(Form("%s_plots.pdf]", isVCASN ? "VCASN":"ITHR"));
  }

}
