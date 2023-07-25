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
#include <TGraph.h>
#include <TLine.h>

using namespace std;

void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns, long int refrun);

//
// MAIN
//
void CompareLayerDeadPixels(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*THRMAPS_DEADPIXMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  bool isIB;
  if(fpath.find("IB")!=string::npos){
    isIB = kTRUE;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    if(stoi(layernum)>=0 && stoi(layernum)<=2) nchips = 9;
    else if (stoi(layernum)==3 && stoi(layernum)==4) nchips = 54*2;
    else nchips = 98*2;
    if(nchips==9) isIB=kTRUE;
    else isIB=kFALSE;
  }

  //Choose whether to skip runs
  string skipans, skipruns;
  cout<<endl;
  cout<<"Would you like to skip some run(s)? [y/n] ";
  cin>>skipans;
  if(skipans=="y" || skipans=="Y"){
    cout<<endl;
    cout<<"Specify run number(s) separated by comma (no white spaces!):";
    cin>>skipruns;
    cout<<endl;
  }
  else
    skipruns=" ";

  cout<<"Available runs in your file:\n"<<endl;
  TFile *infile=new TFile(fpath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TKey *key;
  TObject *obj;
  vector<string> laynums;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    string objtitle = (string)obj->GetTitle();
    if(objname.find("Stv")!=string::npos) break;
    if(objtitle.find("Threshold")!=string::npos) continue;
    string runnum =  objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    laynums.push_back(laynum);
    //cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
    if((int)laynums.size()>1 && laynum!=laynums[laynums.size()-2])
      break;

    cout<<runnum<<endl;
  }


  long int refrun;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>refrun;


  //Call
  DoAnalysis(fpath, nchips, isIB, skipruns, refrun);

}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns, long int refrun){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers, laynums;
  int nTimes=0, nRuns=1;

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  TH2 *h2;
  vector<int> posrefrun;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    string objtitle = (string)obj->GetTitle();
    if(objname.find("Stv")!=string::npos) break;
    if(objtitle.find("Threshold")!=string::npos) continue;

    h2 = (TH2*)obj;
    //if(!h2->GetEntries()) continue;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(h2);

    //cout<<runnum<<"  "<<timestamp<<endl;
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);

    laynums.push_back(laynum);

    //position of ref run
    if(stol(runnum)==refrun)
      posrefrun.push_back(nTimes);

    nTimes++;
    if(nTimes>1 && laynum==laynums[laynums.size()-2])
      nRuns++;
    else nRuns=1;
  }

  const int nLayers = (int)hmaps.size()==nRuns ? 1 : stoi(laynums[laynums.size()-1])+1;


  //find max for plots (same for all staves)
  double max = -1.;
  for(int ihist=0; ihist<(int)hmaps.size(); ihist++)
    for(int ix=1; ix<=hmaps[ihist]->GetNbinsX();ix++)
      for(int iy=1; iy<=hmaps[ihist]->GetNbinsY();iy++)
        if(hmaps[ihist]->GetBinContent(ix,iy)>max)
          max=hmaps[ihist]->GetBinContent(ix,iy);

  double maxlimit = max+0.5*max;
  int nbins = (int)(max+0.5*max-0.9)*1.0;

  TH2D *hCorr[nLayers];

  for(int ilay=0; ilay<nLayers; ilay++)
    hCorr[ilay] = new TH2D(Form("hCorr_L%s",laynums[ilay*nRuns].c_str()), Form("Layer-%s - DeadPix corr. %s - Ref. run: %ld; # Dead Pixel per Chip (run%ld); # Dead Pixel per Chip (runs)",laynums[ilay*nRuns].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str(),refrun,refrun), nbins, 0.9, maxlimit, nbins, 0.9, maxlimit);

  int ilayer=nLayers-1;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){
    if(stol(runnumbers[ihist])==refrun){
      if(ihist>0)
        if(laynums[ihist-1]!=laynums[ihist]){// in case the ref run is the first in the list of all layers
          ilayer--;
        }
      continue; //skip ref run
    }
    for(int ibinx=1; ibinx<=hmaps[ihist]->GetNbinsX(); ibinx++){
      for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){
        double data_refrun = hmaps[posrefrun[ilayer]]->GetBinContent(ibinx,ibiny);
        double data_run = hmaps[ihist]->GetBinContent(ibinx,ibiny);
        hCorr[ilayer]->Fill(data_refrun, data_run);
        //cout<<fhr_refrun<<"   "<<fhr_run<<endl;
      }
    }
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
  }

  gStyle->SetPalette(1);
  //Draw
  TLine *line = new TLine(1,1,maxlimit,maxlimit);
  line->SetLineStyle(2);
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogx();
    canvas->SetLogy();
    canvas->SetTickx();
    canvas->SetTicky();

    hCorr[ilay]->Draw("COLZ");
    line->Draw("same");
    hCorr[ilay]->GetXaxis()->SetTitleOffset(1.2);
    canvas->SaveAs(Form("../Plots/Layer%s_deadpixcorr_refrun%ld_%s.pdf", laynums[ilay*nRuns].c_str(), refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_deadpixcorr_refrun%ld_%s.root", laynums[ilay*nRuns].c_str(), refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    delete canvas;
  }

}
