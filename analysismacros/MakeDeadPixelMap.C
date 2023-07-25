#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "inc/itsAnalysis.hh"
#include "inc/constants.h"
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/CcdbDatabase.h"
#include "inc/ccdb.h"

void MakeDeadPixelMap(){

bool ccdb_upload;
string CCDB_up;

 cout<<"Would you like to upload the output to ccdb? [y/n] ";
  cin>>CCDB_up;
  cout<<endl;
 if(CCDB_up =="y"||CCDB_up =="Y") ccdb_upload= true;
  else ccdb_upload= false;

if(ccdb_upload)SetTaskName(__func__);

std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");

auto ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());

  ccdb->connect(ccdbport.c_str(), "", "", "");

  itsAnalysis myAnalysis("Hits on Layer"); // Change to "Hits on Layer" if using FHR as fake data

  auto nLayers      = myAnalysis.nLayers();     // int of number of layers
  auto laynums      = myAnalysis.Layers();      //vec of layers
  auto runNumbers   = myAnalysis.Runs();        //vec of run numbers
  auto hmaps        = myAnalysis.loadedHists(); // all histograms for layers and runs needed
  auto nRuns        = myAnalysis.nRuns();

  std::unique_ptr<TFile> myFile( TFile::Open(Form("DeadPixMapResults_run%s_to_run%s.root",runNumbers[nRuns-1].c_str(),runNumbers[0].c_str()), "RECREATE") );

  for (string layer : laynums){ // loop over layers
    int ilay=stoi(layer);
    int nStavesInLay = myAnalysis.stavesInLayer(ilay);
    
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay),800,1200);

    cnv.SetTopMargin(0.4);
    if (ilay>=3) cnv.Divide(1,nStavesInLay/2,0,0);
    else cnv.Divide(1,nStavesInLay,0,0);

    int iHists = -1;

    for(int istave=0; istave<nStavesInLay; istave++){
      iHists ++;
      TH2F *hHotMap;
      if (ilay==3 || ilay==4) TH2F *hHotMap = new TH2F(Form("hHotMap_L%s_Stv%i",layer.c_str(),istave), "; ; ", 28672/4,0.5,28672.5,2048/4,0.5,2048.5);
      else if (ilay>=5) TH2F *hHotMap = new TH2F(Form("hHotMap_L%s_Stv%i",layer.c_str(),istave), "; ; ", 50176/4,0.5,50176.5,2048/4,0.5,2048.5);
      else TH2F *hHotMap = new TH2F(Form("hHotMap_L%s_Stv%i",layer.c_str(),istave), "; ; ", 9216,0.5,9216.5,512,0.5,512.5);

      int cnt = 0;

      for(auto hist : myAnalysis.loadLayerSparse(stoi(layer))){ //loop hist for layer
        if(hist->GetEntries()>1e4){
          if(istave==0)cout<<hist->GetName()<<" skipped because has more than 10000 entries (probably a bad run)."<<endl;
          continue;
        }
        int stave = myAnalysis.getStaveNumber(hist);
        if(stave!=istave) continue;
        if(!cnt) {hHotMap =(TH2F*)hist->Projection(1,0);
                  if(ilay>=3) hHotMap->Rebin2D(4,4); 
        }
        else {
          TH2F *htemp = (TH2F*)hist->Projection(1,0);
          if(ilay>=3) htemp->Rebin2D(4,4);
          hHotMap->Add(htemp);
          delete htemp;
        }
        cnt++;
      }

      hHotMap->SetTitle(Form("hHotMap_lay%i_stv%i",ilay,istave));
      hHotMap->SetName(Form("hHotMap_lay%i_stv%i",ilay,istave));
      hHotMap->SetMarkerStyle(20);
      hHotMap->SetMarkerSize(0.6);
      hHotMap->SetMarkerColor(kRed);
      hHotMap->SetLineColor(kRed);
      hHotMap->Write();

      cnv.cd(iHists+1);
      cnv.GetPad(iHists+1)->SetTickx();
      cnv.GetPad(iHists+1)->SetTicky();
      cnv.GetPad(iHists+1)->SetRightMargin(0.01);
      if(!iHists) cnv.GetPad(iHists+1)->SetTopMargin(0.1);
      hHotMap->SetTitle(" ");

      hHotMap->GetXaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetLabelSize(0.13);
      hHotMap->GetXaxis()->SetLabelSize(0.13);

      if(iHists>0){
        hHotMap->GetXaxis()->SetLabelOffset(999);
        hHotMap->GetXaxis()->SetTickLength(0.05);
        hHotMap->GetXaxis()->SetNdivisions(530);
      }
      else{
        hHotMap->GetXaxis()->SetLabelOffset(0.003);
        hHotMap->GetXaxis()->SetTickLength(0.05);
        hHotMap->GetXaxis()->SetNdivisions(530);
      }
      hHotMap->SetStats(0);
      hHotMap->DrawCopy("P X+");
	if(ccdb_upload){
	string Runperiod = Form("from_run%s_to_run%s",runNumbers.back().c_str(),runNumbers[0].c_str());
	hHotMap->SetTitle(Form("Dead pixel map of Layer%i Stave %i",ilay,istave));
	hHotMap->GetXaxis()->SetLabelOffset(0.003);
	hHotMap->SetTitleSize(0.05);		
	hHotMap->GetYaxis()->SetLabelSize(0.11);
	auto mo = std::make_shared<o2::quality_control::core::MonitorObject>(hHotMap, TaskName+Form("/Layer%d",ilay), TaskClass, DetectorName,std::stoi(runNumbers.back()),Runperiod);
       mo->setIsOwner(false);
       ccdb->storeMO(mo);
		}

	hHotMap->SetTitle(" ");
	hHotMap->GetXaxis()->SetLabelOffset(999);
      TLatex lat;
      lat.SetTextAngle(90);
      lat.SetNDC();
      lat.SetTextSize(0.15);
      lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));

      delete hHotMap;

      cout<<"Layer: "<<layer<<" stv: "<<istave<<" / "<<nStavesInLay-1<<endl;

      if(((ilay>=3) && (istave+1)==nStavesInLay/2) || (istave+1)==nStavesInLay){
          cnv.cd();
          TLatex lat;
          lat.SetNDC();
          lat.SetTextSize(0.03);
          lat.DrawLatex(0.01,0.98,Form("L%s",layer.c_str()));

          if ((ilay>=3) && (istave+1)==nStavesInLay/2) cnv.SaveAs(Form("../Plots/Layer%i_Top_Deadpixmap_run%s_to_run%s.pdf",ilay,runNumbers[nRuns-1].c_str(),runNumbers[0].c_str()));
          else if ((ilay>=3) && (istave+1)==nStavesInLay) cnv.SaveAs(Form("../Plots/Layer%i_Bottom_Deadpixmap_run%s_to_run%s.pdf",ilay,runNumbers[nRuns-1].c_str(),runNumbers[0].c_str()));
          else cnv.SaveAs(Form("../Plots/Layer%i_Deadpixmap_run%s_to_run%s.pdf",ilay,runNumbers[nRuns-1].c_str(),runNumbers[0].c_str()));

          cnt=0;
          iHists = -1;
      }

    } // End loop on staves
  } // End loop on layers
 ccdb->disconnect();
} // end of qMakeDeadPixelMap()
