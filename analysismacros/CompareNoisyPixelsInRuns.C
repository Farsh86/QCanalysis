#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TText.h>
#include <TSystem.h>
#include <TKey.h>
#include <THnSparse.h>
#include <TClass.h> 
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/CcdbDatabase.h"
#include "inc/ccdb.h"
#include "inc/constants.h"
#include <TAttAxis.h>
#include <TAxis.h>

using namespace std;
vector<int> RunList[7];
using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;

//Functions
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2);
void SetStyle(TGraphErrors *ge, Color_t col);
void DoAnalysis(string filepath, const int nChips, string skipruns, long int refrun, int IBorOB, bool ccdb_upload);

void CompareNoisyPixelsInRuns(){
  string fpath;
  bool ccdb_upload;
  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  //IBorOB = 0 if I want to check all IB layers                                                              
  //IBorOB = 1 if I want to check all OB layers                                                              
  //IBorOB = 2 if I want to check all IB + OB layers or if I want to check a single layer                    

  int IBorOB=0;
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
  TObject *obj;
  TKey *key;
  vector<string> stavenums;
  vector<string> laynums;  
  vector<int> CommonRunList;
  int LastLayerInList=-1;

  int counter=0;
  bool isNewLayer=0;
  bool isNewStave=0;
  while((key = ((TKey*)next()))){
    obj=key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
         && (!obj->InheritsFrom("THnSparse"))) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    string runnum =  objname.substr(objname.find("run")+3, 6);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    string laynum = objname.substr(objname.find("L")+1,1);

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user
    if(stvnum.find("_")!=string::npos)
      stvnum = objname.substr(objname.find("Stv")+3,1);
    stavenums.push_back(stvnum);
    //    if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-2])  break;

    laynums.push_back(laynum);
    LastLayerInList = stoi(laynum);

    //CREATE THE LIST OF AVAILABLE RUNS FOR EACH LAYER (the input plots are ordered in tis way: for each Layer and Stave the full list of runnumber is given)
    isNewLayer=0;
    if ((int)laynums.size()==1 || ((int)laynums.size()>1 && laynum!=laynums[laynums.size()-2]))  { //new layer
      //      cout << "\nLayer: " << laynum << "\nRuns: " <<endl;
      isNewLayer=1;
      isNewStave=0;
    }
    if (isNewLayer) {
      //      cout<<"run: "<<runnum<<" laynum: "<<laynum<< " stavenum " << stvnum<< endl; //first stave of new layer
      RunList[stoi(laynum)].push_back(stoi(runnum));
    }
    else {
      if ((int)stavenums.size()>1 && stvnum==stavenums[stavenums.size()-2]){
	if (!isNewStave) {
	  RunList[stoi(laynum)].push_back(stoi(runnum));
	  //	  cout<<"run: "<<runnum<<" laynum: "<<laynum<< " stavenum " << stvnum<< endl;
	}
      }
      else{
	isNewStave=1;
      }
    }
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

  cout << "\n\nCommon runs: " << endl;
  for (int i=0; i<(int)CommonRunList.size(); i++){
    cout << CommonRunList[i] << endl;
  }

  long int refrun;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>refrun;

  DoAnalysis(fpath, nchips, skipruns, refrun, IBorOB, ccdb_upload);
}

//
// Analysis
//
void DoAnalysis(string filepath, const int nChips, string skipruns, long int refrun, int IBorOB, bool ccdb_upload){

  gStyle->SetOptStat(0000);

  std::vector<THnSparse*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums;
  vector<int> posrefrun;
  int nTimes=0;
  int nRunsB[7]={-1};
  for (int ilay=0; ilay < 7; ilay++){
    nRunsB[ilay] =-1;
  }
  int nRunsTot=0;
  int nRunsTotFixed=0;
  int nLayersInput =0;

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
  THnSparse *hsparse;
  bool isNewStave=0;
  bool isNewLayer=0;
  while((key=((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D or sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    hsparse = (THnSparse*)obj;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = objname.substr(objname.find("Stv")+3,1);
    }

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(hsparse);

    if(stol(runnum)==refrun) //position of refence run
      posrefrun.push_back((int)hmaps.size()-1);

    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    laynums.push_back(laynum);
    stavenums.push_back(stvnum);

    nTimes++;
    if(nRunsB[stoi(laynum)]==-1) {
      nRunsB[stoi(laynum)]=0;
      nLayersInput++;
    }

    isNewLayer=0;
    if ((int)laynums.size()==1 || ((int)laynums.size()>1 && laynum!=laynums[laynums.size()-2]))  { //new layer
      //      cout << "\nLayer: " << laynum << "\nRuns: " <<endl;
      isNewLayer=1;
      isNewStave=0;
    }
    if (isNewLayer) {
      //      cout << "isNewLayer " << endl;
      //      cout<<"run: "<<runnum<<" laynum: "<<laynum<< " stavenum " << stvnum<< endl; //first stave of new layer
      nRunsB[stoi(laynum)]++;
    }
    else {
      //      cout << "isNewStave " << isNewStave << endl;
      if ((int)stavenums.size()>1 && stvnum==stavenums[stavenums.size()-2]){
	if (!isNewStave) {
	  nRunsB[stoi(laynum)]++;
	  //  cout<<"run: "<<runnum<<" laynum: "<<laynum<< " stavenum " << stvnum<< endl;
	}
      }
      else{
	isNewStave=1;
      }
    }
  }

  //  cout << "nLayersInput " << nLayersInput << endl;
  int nLayers;
  if (nLayersInput==1) nLayers =1;
  else {
    if (IBorOB==0) nLayers = 3;
    else if (IBorOB==1) nLayers = 4;
    else nLayers = 7;
  }

  //  cout << "nLayers= " << nLayers << endl;
  int nRunsMax =0; //maximum number of runs among different layers
  for (int ilay=0; ilay < 7; ilay++){
    //    cout <<  "nRunsB " <<    nRunsB[ilay]<< endl;
    if (nRunsB[ilay] > nRunsMax) nRunsMax = nRunsB[ilay];
  }

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  vector<int> runlabel[7];
  int istave = 0;
  int ilayer = 0;
  long int first[7][nRunsMax], second[7][nRunsMax], both[7][nRunsMax];
  int irun=0;
  vector<array<long int,5>> noisypix;
  int ilayEff = 0;
  for(int ilay=0; ilay<7; ilay++) {
    for(int i=0; i<nRunsMax; i++){
      first[ilay][i]=0; second[ilay][i]=0; both[ilay][i]=0;
    }
  }

  Bool_t isFirstLayer =1;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run   
    if(runnumbers[ihist].find(std::to_string(refrun))!=string::npos){
      if(ihist>0){// in case ref run is the first into the list of runs
        if(stavenums[ihist-1]!=stavenums[ihist]){
          irun=0;
        }
      }
      continue;
    }

    if (IBorOB==1)  ilayer = stoi(laynums[ihist])-3;
    else ilayer = stoi(laynums[ihist]);
    if (nLayersInput==1) ilayer = 0;
    istave = stoi(stavenums[ihist]);

    noisypix.push_back(CompareTwoRuns(hmaps[posrefrun[istave]], hmaps[ihist]));
    noisypix[noisypix.size()-1][3] = stol(stavenums[ihist]);
    noisypix[noisypix.size()-1][4] = stol(laynums[ihist]);
    first[ilayer][irun]+=noisypix[noisypix.size()-1][0];
    second[ilayer][irun]+=noisypix[noisypix.size()-1][1];
    both[ilayer][irun]+=noisypix[noisypix.size()-1][2];
    irun++;

    if(stavenums[ihist]=="0"){
      runlabel[stoi(laynums[ihist])].push_back(stoi(runnumbers[ihist]));
    }

    if(ihist>0){
      if(stavenums[ihist-1]!=stavenums[ihist]){
	irun=0;
      }
    }
  }//end loop on histograms

  //Make plot for each layer and for each stave in the root file
  TGraphErrors *ge_nref[nLayers];
  TGraphErrors *ge_n2[nLayers];
  TGraphErrors *ge_ncom1[nLayers];
  TGraphErrors *ge_ncom2[nLayers];
  TGraphErrors *ge_nref_stave[nLayers][100];//100 is an abitrary large number of staves
  TGraphErrors *ge_n2_stave[nLayers][100];
  TGraphErrors *ge_ncom1_stave[nLayers][100];
  TGraphErrors *ge_ncom2_stave[nLayers][100];
  double xshift = 3.;
  double max[nLayers];
  double min[nLayers];
  double scalingFactor = 1e2;

  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;

    ge_nref[ilay] = new TGraphErrors();
    ge_n2[ilay] = new TGraphErrors();
    ge_ncom1[ilay] = new TGraphErrors();
    ge_ncom2[ilay] = new TGraphErrors();
    max[ilay] = -1.;
    min[ilay] = 1e35;

    for(int ir=0; ir<nRunsB[ilayEff]-1; ir++){//first the older data and last the most recent
      //first couple of bar on the left
      if(!ir) xshift=1.;
      else xshift = 3.;
      ge_nref[ilay]->SetPoint(ir, ir*xshift, (double)both[ilay][ir]/2.+((double)first[ilay][ir]/2.) / scalingFactor);
      ge_nref[ilay]->SetPointError(ir, 0.5, ((double)first[ilay][ir]/2.) / scalingFactor);
      ge_ncom1[ilay]->SetPoint(ir, ir*xshift, 0.);
      ge_ncom1[ilay]->SetPointError(ir, 0.5, (double)both[ilay][ir]/2.);
      if((double)both[ilay][ir]/2.+((double)first[ilay][ir]) / scalingFactor > max[ilay]) max[ilay] = (double)both[ilay][ir]/2.+((double)first[ilay][ir]) / scalingFactor;

      //second couple of bar on the right
      ge_n2[ilay]->SetPoint(ir, ir*xshift+1, -(double)both[ilay][ir]/2.-((double)second[ilay][ir]/2.) / scalingFactor);
      ge_n2[ilay]->SetPointError(ir, 0.5, ((double)second[ilay][ir]/2.) / scalingFactor);
      ge_ncom2[ilay]->SetPoint(ir, ir*xshift+1, 0.);
      ge_ncom2[ilay]->SetPointError(ir, 0.5, (double)both[ilay][ir]/2.);
      if(-(double)both[ilay][ir]/2.-((double)second[ilay][ir]) / scalingFactor < min[ilay]) min[ilay] = -(double)both[ilay][ir]/2.-((double)second[ilay][ir]) / scalingFactor;
    }//end first loop on runs

    //Style
    SetStyle(ge_nref[ilay], kBlue);
    SetStyle(ge_ncom1[ilay], kBlack);
    SetStyle(ge_ncom2[ilay], kBlack);
    SetStyle(ge_n2[ilay], kRed+2);
  }//end loop on layers

  //FILL PLOTS FOR EACH STAVE
  //fill the plots for each stave
  int cnt = 0;
  double maxs[nLayers][100];
  double mins[nLayers][100];
  for(int ilay=nLayers-1; ilay>=0; ilay--){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    for(int is=nStavesInLayAll[ilayEff]-1; is>=0; is--){
      maxs[ilay][is] = -1.;
      mins[ilay][is] = 1e35;
      ge_nref_stave[ilay][is] = new TGraphErrors();
      ge_ncom1_stave[ilay][is] = new TGraphErrors();
      ge_n2_stave[ilay][is] = new TGraphErrors();
      ge_ncom2_stave[ilay][is] = new TGraphErrors();
      for(int ir=0; ir<nRunsB[ilayEff]-1; ir++){//first the older data and last the most recent
        if(!ir) xshift=1.;
        else xshift = 3.;
        ge_nref_stave[ilay][is]->SetPoint(ir, ir*xshift, (double)noisypix[cnt][2]/2.+((double)noisypix[cnt][0]/2.) / scalingFactor);
        ge_nref_stave[ilay][is]->SetPointError(ir, 0.5, ((double)noisypix[cnt][0]/2.) / scalingFactor);
        ge_ncom1_stave[ilay][is]->SetPoint(ir, ir*xshift, 0.);
        ge_ncom1_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
        if((double)noisypix[cnt][2]/2.+((double)noisypix[cnt][0]) / scalingFactor > maxs[ilay][is]) maxs[ilay][is] = (double)noisypix[cnt][2]/2.+((double)noisypix[cnt][0]) / scalingFactor;

        //second couple of bar on the right
        ge_n2_stave[ilay][is]->SetPoint(ir, ir*xshift+1, -(double)noisypix[cnt][2]/2.-((double)noisypix[cnt][1]/2.) / scalingFactor);
        ge_n2_stave[ilay][is]->SetPointError(ir, 0.5, ((double)noisypix[cnt][1]/2.) / scalingFactor);
        ge_ncom2_stave[ilay][is]->SetPoint(ir, ir*xshift+1, 0.);
        ge_ncom2_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
        if(-(double)noisypix[cnt][2]/2.-((double)noisypix[cnt][1]) / scalingFactor < mins[ilay][is]) mins[ilay][is] = -(double)noisypix[cnt][2]/2.-((double)noisypix[cnt][1]) / scalingFactor;
        cnt++;
      }
      SetStyle(ge_nref_stave[ilay][is], kBlue);
      SetStyle(ge_ncom1_stave[ilay][is], kBlack);
      SetStyle(ge_ncom2_stave[ilay][is], kBlack);
      SetStyle(ge_n2_stave[ilay][is], kRed+2);
    }
  }

  //Legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref[0], Form("#splitline{#noisy pix}{#splitline{ref. run only}{#divide %.2f}}", scalingFactor), "f");
  leg->AddEntry(ge_n2[0], Form("#splitline{#noisy pix}{#splitline{2nd run only}{#divide %.2f}}",scalingFactor), "f");
  leg->AddEntry(ge_ncom1[0], "#splitline{#noisy pix}{both}");

  //Draw plot for each layer
  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    TCanvas *canvas = new TCanvas(Form("mycanvas_%d",ilayEff), Form("mycanvas_%d",ilayEff), 1300, 800);
    canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
    canvas->cd();

    //fake histo (just for the axes)
    double x2,y2;
    ge_ncom2[ilay]->GetPoint(ge_ncom2[ilay]->GetN()-1, x2,y2);
    TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
    //draw labels on x axis
    int counter = 0;
    for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
      hfake->GetXaxis()->SetBinLabel(k, Form("run%i", runlabel[ilayEff][counter]));
      counter++;
    }
    hfake->Draw();
    //canvas->SetLogy();
    hfake->SetTitle(Form("Layer-%i - %s%06ld compared to all",ilayEff, filepath.find("run")==string::npos? "":"run",refrun));
   
    ge_nref[ilay]->Draw("E2 same");
    ge_ncom1[ilay]->Draw("E2 same");
    ge_ncom2[ilay]->Draw("E2 same");
    ge_n2[ilay]->Draw("E2 same");
    hfake->GetYaxis()->SetRangeUser(min[ilay]+0.1*min[ilay], max[ilay]+0.1*max[ilay]);
    hfake->GetYaxis()->SetTickLength(0.005);
    TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
    lineref->SetLineColor(kGray-1);
    lineref->SetLineStyle(2);
    lineref->Draw("same");

    //draw legend
    leg->Draw("same");
    string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str());
    string Reference_run = to_string(refrun);
    if(ccdb_upload){
      string canvas_name = Form("Layer%d_NoisyPixComparison_Allstaves",ilay);
      canvas->SetName(canvas_name.c_str());
      auto mo= std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%d",ilay),TaskClass, DetectorName,std::stoi(runnumbers.front()),Runperiod);
      mo->addMetadata("ReferenceRunNumber",Reference_run.c_str());
      mo->setIsOwner(false);
      ccdb->storeMO(mo);	
    }
    canvas->SaveAs(Form("../Plots/Layer%i_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf", ilayEff,filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%i_NoisyPixComparison_%s%ld_compared_to_run_%s.root", ilayEff,filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

    delete canvas;
    delete hfake;
    delete lineref;
  }//end loop on layers

  //Draw plot for each stave
  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str());
    string Reference_run = to_string(refrun);

    for(int is=0; is<nStavesInLayAll[ilayEff]; is++){
      TCanvas *canvas = new TCanvas(Form("mycanvas_%d_%d",ilayEff,is), Form("mycanvas_%d_%d",ilayEff,is), 1300, 800);
      canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
      canvas->cd();

      //fake histo (just for the axes)
      double x2,y2;
      ge_ncom2_stave[ilay][is]->GetPoint(ge_ncom2_stave[ilay][is]->GetN()-1, x2,y2);
      TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
      //draw labels on x axis
      int counter = 0;
      for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
        hfake->GetXaxis()->SetBinLabel(k, Form("run%i", runlabel[ilayEff][counter]));
        counter++;
      }
      hfake->Draw();
      hfake->SetTitle(Form("Layer-%i - Stave-%d - %s%06ld compared to all", ilayEff, is, filepath.find("run")==string::npos? "":"run",refrun));
      ge_nref_stave[ilay][is]->Draw("E2 same");
      ge_ncom1_stave[ilay][is]->Draw("E2 same");
      ge_ncom2_stave[ilay][is]->Draw("E2 same");
      ge_n2_stave[ilay][is]->Draw("E2 same");
      hfake->GetYaxis()->SetRangeUser(mins[ilay][is]+0.1*mins[ilay][is], maxs[ilay][is]+0.1*maxs[ilay][is]);
      hfake->GetYaxis()->SetTickLength(0.005);
      TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
      lineref->SetLineColor(kGray-1);
      lineref->SetLineStyle(2);
      lineref->Draw("same");
      leg->Draw("same");

      TString LayerTitle="";
      if (IBorOB==0) LayerTitle = "AllLayersIB_";
      else if (IBorOB==1) LayerTitle = "AllLayersOB_";
      else {
	if (nLayers==1) LayerTitle = Form("Layer%i_",ilayEff);
	else LayerTitle = "AllLayers_";
      }
      
      if(ccdb_upload){
	canvas->SetName(Form("Layer%d-Stave%d_NoisyPixComparison",ilay,is));
	auto mo2= std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%d",ilay),TaskClass, DetectorName,std::stoi(runnumbers.front()),Runperiod);
        mo2->addMetadata("ReferenceRunNumber",Reference_run.c_str());
        mo2->setIsOwner(false);
        ccdb->storeMO(mo2);

      }
      if(!ilay && !is) canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf[", LayerTitle.Data(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf", LayerTitle.Data() ,filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      if(ilay==nLayers-1 && is==nStavesInLayAll[ilayEff]-1) canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf]", LayerTitle.Data(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

      delete canvas;
      delete hfake;
      delete lineref;
    }
  }
//Disconnencting the interface
//if(ccdb_upload)
 ccdb->disconnect();
}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2){

  std::array<long int,5> noisypix = {0, 0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int ibinref=0; ibinref<href->GetNbins(); ibinref++){
    int coordref[2];
    double bincref = href->GetBinContent(ibinref, coordref);
    bool isfound = false;
    for(int ibin2=0; ibin2<h2->GetNbins(); ibin2++){
      int coord2[2];
      double binc2 = h2->GetBinContent(ibin2,coord2);

      if((coordref[0]==coord2[0] && coordref[1]==coord2[1]) && (bincref>2 && binc2>2)){//noisy in both runs
        noisypix[2]++;
        isfound=true;
        break;
      }

    }
    if(!isfound && bincref>2)
      noisypix[0]++; //noisy only in ref run
  }

  for(int ibin2=0; ibin2<h2->GetNbins(); ibin2++){
    int coord2[2];
    double binc2 = h2->GetBinContent(ibin2, coord2);
    bool isfound = false;
    for(int ibinref=0; ibinref<href->GetNbins(); ibinref++){
      int coordref[2];
      double bincref = href->GetBinContent(ibinref,coordref);

      if((coordref[0]==coord2[0] && coordref[1]==coord2[1]) && (bincref>2 && binc2>2)){//noisy in both runs
        isfound=true;
        break;
      }

    }
    if(!isfound && binc2>2)
      noisypix[1]++; //noisy only in second run
  }

  return noisypix;
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
