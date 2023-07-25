#include "inc/itsAnalysis.hh"
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/CcdbDatabase.h"
#include "inc/ccdb.h"

void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
}

Int_t col[] = {810, 807, 797, 827, 417, 841, 868, 867, 860, 602, 921, 874};

// Main function
void CompareLayerThresholds(){

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


  itsAnalysis myAnalysis("Threshold");

  auto laynums      = myAnalysis.Layers();      //vec of layers
  auto runNumbers   = myAnalysis.Runs();        //vec<string> of run numbers
  auto hmaps        = myAnalysis.loadedHists(); // all histograms for layers and runs needed
  auto nRuns        = myAnalysis.nRuns();

  cout<<"Runs that will be used for comparisons: "<<endl;
  for (auto x: runNumbers){
    cout<<x<<endl;
  }
  string refrun;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>refrun;

  TH2D *hCorr[7];
  TH2 *hRef[7];
  for (string layer : laynums){ // loop over layers
    int ilay = stoi(layer);
    hCorr[ilay] = new TH2D(Form("hCorr_L%s",layer.c_str()), Form("Layer-%s - THR corr. for run%s to run%s - Ref. run: %s; Chip Threshold (run%s); Chip Threshold (runs)",layer.c_str(),runNumbers[nRuns-1].c_str(),runNumbers[0].c_str(),refrun.c_str(),refrun.c_str()), 150, 7.5, 13.5, 150, 7.5, 13.5);
    hCorr[ilay]->SetStats(0);
    
    for(auto hist : myAnalysis.loadLayer(ilay)){
        if(myAnalysis.getRunNumber(hist)==refrun){
            hRef[ilay] = hist;
        }
      }
    for(auto hist : myAnalysis.loadLayer(ilay)){
      for(int ibinx=1; ibinx<=hist->GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<=hist->GetNbinsY(); ibiny++){
          double data_refrun = hRef[ilay]->GetBinContent(ibinx,ibiny);
          double data_run = hist->GetBinContent(ibinx,ibiny);
          hCorr[ilay]->Fill(data_refrun, data_run);
        }
      }
    }
  }

  gStyle->SetPalette(1);
  //Draw
  TLine *line = new TLine(7.5,7.5,13.5,13.5);
  for (string layer : laynums){ // loop over layers
    int ilay = stoi(layer);    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    hCorr[ilay]->Draw("COLZ");
    line->Draw("same");
    hCorr[ilay]->GetXaxis()->SetTitleOffset(1.2);
    if(ccdb_upload){
	string Runperiod = Form("from_run%s_to_run%s",runNumbers[0].c_str(),runNumbers[nRuns-1].c_str());
 //   string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str()); //This should be used for actual data      
	canvas->SetName(Form("Layer%s_Threshold_correlation",layer.c_str()));			
	auto mo1= std::make_shared<o2::quality_control::core::MonitorObject>(canvas,TaskName+Form("/Layer%s",layer.c_str()), TaskClass, DetectorName,std::stoi(runNumbers.back()),Runperiod);
 	mo1->addMetadata("ReferenceRunNumber",refrun.c_str());
        mo1->setIsOwner(false);
        ccdb->storeMO(mo1);
}
    canvas->SaveAs(Form("../Plots/Layer%s_thresholdcorr_run%s_to_run%s_referencerun%s.pdf", layer.c_str(),runNumbers[0].c_str(),runNumbers[nRuns-1].c_str(),refrun.c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_thresholdcorr_run%s_to_run%s_referencerun%s.root", layer.c_str(),runNumbers[0].c_str(),runNumbers[nRuns-1].c_str(),refrun.c_str()));
    delete canvas;
  }
//Disconnencting from the database
 ccdb->disconnect();
} // end of CompareLayerThresholds()
