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
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/CcdbDatabase.h"
#include "inc/ccdb.h"

using namespace std;
using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;

void DoAnalysis(string filepath, int IBorOB, string skipruns, long int refrun, bool ccdb_upload);

//
// MAIN
//
void CompareLayerOccupancy(){
  string fpath;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  int IBorOB;
  bool ccdb_upload;
  //IBorOB = 0 if I want to check all IB layers                                                                   
  //IBorOB = 1 if I want to check all OB layers                                                                    
  //IBorOB = 2 if I want to check all IB + OB layers or if I want to check a single layer                         

  if(fpath.find("IB")!=string::npos){
    IBorOB = 0;
  }
  else if (fpath.find("OB")!=string::npos){
    IBorOB = 1;
  }
  else if (fpath.find("all")!=string::npos){
    IBorOB = 2;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    IBorOB = 2;
  }

  //Choose whether to skip runs
  string skipans, skipruns, CCDB_up;
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

  cout<<"Would you like to upload the output to ccdb? [y/n] ";
  cin>>CCDB_up;
  cout<<endl;
  if(CCDB_up =="y"||CCDB_up =="Y") ccdb_upload= true;
  else ccdb_upload= false;

if(ccdb_upload)SetTaskName(__func__);




  cout<<"Available runs in your file:\n"<<endl;
  TFile *infile=new TFile(fpath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TKey *key;
  TObject *obj;
  TH2 *h2;
  vector<string> laynums;
  vector<int> RunList[7];
  vector<int> CommonRunList;
  int LastLayerInList=-1;

  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")!=string::npos) break;
    string runnum =  objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    h2 = (TH2*)obj;
    if(!h2->GetEntries()) {
      //      cout << "Layer: " << laynum << " Run number: " << runnum << " --> No entries " <<endl;
      continue;
    }

    LastLayerInList = stoi(laynum);
    laynums.push_back(laynum);
    RunList[stoi(laynum)].push_back(stoi(runnum));
    if ((int)laynums.size()==1 || ((int)laynums.size()>1 && laynum!=laynums[laynums.size()-2]))  {
      //cout << "\nLayer: " << laynum << "\nRuns: " <<endl;
    }
    //    cout<<runnum<<endl;
  }

  bool isCommon=1;
  for (Int_t i = 0; i< (int)RunList[LastLayerInList].size(); i++){
    isCommon=1; //let's suppose the run i is common to all layers
    for (Int_t lay = 0; lay<7; lay++){ 
      if (isCommon==0) break; // if a run is missing for one layer, than it's not a common run
      if (lay == LastLayerInList) continue;
      if ((int)RunList[lay].size() == 0) continue; //the layer has no runs
      for (Int_t l = 0; l< (int)RunList[lay].size(); l++){
	if (RunList[LastLayerInList][i] == RunList[lay][l])  break;
	else {
	  if (l== (int)RunList[lay].size()-1) {
	    isCommon=0; //the run is not common to one layer, therefore it's not a common run
	    break;
	  }
	}
      }
    }
    if (isCommon) CommonRunList.push_back(RunList[LastLayerInList][i]);
  }

  cout << "Common runs: " << endl;
  for (int i=0; i<(int)CommonRunList.size(); i++){
    cout << CommonRunList[i] << endl;
  }

  long int refrun;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>refrun;


  //Call
  DoAnalysis(fpath, IBorOB, skipruns, refrun,ccdb_upload);

}

//
// Analyse data
//
void DoAnalysis(string filepath, int IBorOB, string skipruns, long int refrun, bool ccdb_upload){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers, laynums;
  int nTimes=0, nRuns=1;
  int nRunsB[7]={-1};
  for (int ilay=0; ilay < 7; ilay++){
    nRunsB[ilay] =-1;
  }
  int nRunsTot=0;
  int nRunsTotFixed=0;
  int nLayersInput =1;

//Setting up the connection to the ccdb database

//      CcdbDatabase* ccdb;
//      if(ccdb_upload) ccdb = SetupConnection();       ~To-Do- Currently not working      
  std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");

  auto* ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());

  ccdb->connect(ccdbport.c_str(), "", "", "");


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
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")!=string::npos) break;
    h2 = (TH2*)obj;
    if(!h2->GetEntries()) continue;
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
    if(nRunsB[stoi(laynum)]==-1) nRunsB[stoi(laynum)]=0;
    if(nTimes>1){
      if (laynum==laynums[laynums.size()-2]) {
        nRunsB[stoi(laynum)]++;
      }
      else nLayersInput++;
    }
  }

  int nLayers;
  if (nLayersInput==1) nLayers =1;
  else {
    if (IBorOB==0) nLayers = 3;
    else if (IBorOB==1) nLayers = 4;
    else nLayers = 7;
  }
  cout << "nLayers= " << nLayers << endl;
  for (int ilay=0; ilay < 7; ilay++){
    //    cout <<     nRunsB[ilay]<< endl;
  }

  TH2D *hCorr[nLayers];
  //bins
  double exponent = -14.;
  double base = 10;
  double bins[100];
  bins[0] = 1e-14;
  for(int i=1; i<100; i++){
    bins[i] = i%9==0 ?TMath::Power(base, exponent) : bins[i-1] + TMath::Power(base, exponent);
    if((i+1)%9==0) {
      exponent++;
    }
  }
  int ilayEff=0;
  for(int ilay=0; ilay<nLayers; ilay++) {
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    if (ilay>0 && nRunsB[ilayEff-1]!=-1) nRunsTot += (nRunsB[ilayEff-1]+1);
    if (nRunsB[ilayEff] ==-1) continue;
    //    cout << "ilay " << ilay << " ilayEff: " << ilayEff << " nRunsTot " << nRunsTot << endl;
    hCorr[ilay] = new TH2D(Form("hCorr_L%s",laynums[nRunsTot].c_str()), Form("Layer-%s - FHR corr. %s - Ref. run: %ld; FHR (run%ld); FHR (runs)",laynums[nRunsTot].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str(),refrun,refrun), 99, bins, 99, bins);
  }

  int ilayer=0;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){
    if (IBorOB==1)  ilayer = stoi(laynums[ihist])-3;
    else ilayer = stoi(laynums[ihist]);
    if (nLayersInput==1) ilayer = 0;
    //    cout << "ilayer " << ilayer << endl;
    //    cout <<"laynum " << laynums[ihist] << " runnumber " << runnumbers[ihist]<< endl;
    if(stol(runnumbers[ihist])==refrun){
      continue; //skip ref run
    }
    //cout << "\n run number" << runnumbers[ihist] << endl;
    for(int ibinx=1; ibinx<=hmaps[ihist]->GetNbinsX(); ibinx++){
      for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){
        double fhr_refrun = hmaps[posrefrun[ilayer]]->GetBinContent(ibinx,ibiny);
        double fhr_run = hmaps[ihist]->GetBinContent(ibinx,ibiny);
        hCorr[ilayer]->Fill(fhr_refrun, fhr_run);
	//cout<<"stave " << ibiny-1 << " chip " << ibinx-1 << "-> fhr ref run: " << fhr_refrun<<", fhr run: "<<fhr_run<<endl;
      }
    }
  }

  gStyle->SetPalette(1);
  //Draw
  nRunsTot=0;
  TLine *line = new TLine(1e-14,1e-14,1e-3,1e-3);
  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    if (ilay>0 && nRunsB[ilayEff-1]!=-1) nRunsTot += (nRunsB[ilayEff-1]+1);
    if (nRunsB[ilayEff] ==-1) continue;
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogy();
    canvas->SetLogx();
    canvas->SetTickx();
    canvas->SetTicky();

    hCorr[ilay]->Draw("COLZ");
    line->Draw("same");
    hCorr[ilay]->GetXaxis()->SetTitleOffset(1.2);
if(ccdb_upload){
    string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str());
    string canvas_name = Form("Layer%s_fakehitratecorr",laynums[nRunsTot].c_str());
    string Reference_run = to_string(refrun);

    canvas->SetName(canvas_name.c_str());
        auto mo= std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%s",laynums[nRunsTot].c_str()),TaskClass, DetectorName,std::stoi(runnumbers.front()),Runperiod);
        mo->addMetadata("ReferenceRunNumber",Reference_run.c_str());
	mo->setIsOwner(false);
        ccdb->storeMO(mo);
	}
canvas->SaveAs(Form("../Plots/Layer%s_fakehitratecorr_refrun%ld_%s.pdf", laynums[nRunsTot].c_str(), refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_fakehitratecorr_refrun%ld_%s.root", laynums[nRunsTot].c_str(), refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    delete canvas;
  }

//Disconnencting the interface
//if(ccdb_upload)
 ccdb->disconnect();
}
