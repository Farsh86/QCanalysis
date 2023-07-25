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

void AnalyzeLayerDeadPixels(){

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


itsAnalysis myAnalysis("DeadPixel");

  auto laynums      = myAnalysis.Layers();      //vec of layers
  auto runNumbers   = myAnalysis.Runs();        //vec of run numbers
  auto nRuns        = myAnalysis.nRuns();
  auto hmaps        = myAnalysis.loadedHists(); // all histograms for layers and runs needed

  TGraph *trend[7][100][2]; //list of trends
  for (int ilayer = 0; ilayer < 7; ++ilayer){ //loop over layers
    for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //loop over number of histograms
      for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){ // Loop over y-bins (staves)
        for(int hs=0; hs<=1; hs++){
          trend[ilayer][ibiny-1][hs] = new TGraph();

        }
      }
    }
  }

  int staveswithdead[7][2];
  double maxtot=-1.;

  TH1F *hproj = new TH1F();
  for (string layer : laynums){ // loop over layers
    staveswithdead[stoi(layer)][0] = 0;
    staveswithdead[stoi(layer)][1] = 0;
    int irun = myAnalysis.nRuns()-1; // Count back to have right order in histogram
    for(auto hist : myAnalysis.loadLayer(stoi(layer))){ //loop hist for layer
      for(int ibiny=1; ibiny<=hist->GetNbinsY(); ibiny++){//loop on y bins (staves)
        TH1D *hproj = hist->ProjectionX("proj",ibiny,ibiny); //single stave
        trend[stoi(layer)][ibiny-1][0]->SetName(Form("gr_L%s_HS-upper_stave%d",layer.c_str(),ibiny-1));
        trend[stoi(layer)][ibiny-1][1]->SetName(Form("gr_L%s_HS-lower_stave%d",layer.c_str(),ibiny-1));


        if(stoi(layer)<=2){
          int nChips = myAnalysis.nChips(stoi(layer));
          int chipswithdeadpix = 0;

          for(int ibinx=1; ibinx<=hist->GetNbinsX(); ibinx++){
            if(hist->GetBinContent(ibinx,ibiny)>0) chipswithdeadpix++;
          }
          if(chipswithdeadpix>0){
            cout<<"Layer: "<<layer<<" Stave: "<<ibiny-1<<" Run: "<<myAnalysis.getRunNumber(hist)<<" --> # Chips with dead pix: "<<chipswithdeadpix<<endl;
            staveswithdead[stoi(layer)][0]++;
          }

          trend[stoi(layer)][ibiny-1][0]->SetPoint(irun, irun, hproj->Integral());
          if(hproj->Integral()>maxtot)
            maxtot=hproj->Integral();
        }
        
        if(stoi(layer)>=3){
          int nChips = myAnalysis.nChips(stoi(layer));
          int chipswithdeadpix_upper = 0;
          int chipswithdeadpix_lower = 0;

          for(int ibinx=1; ibinx<=hist->GetNbinsX()/2; ibinx++){//evaluate the number of disabled chips
            if(hist->GetBinContent(ibinx,ibiny)>0) chipswithdeadpix_upper++;
          }
          for(int ibinx=hist->GetNbinsX()/2+1; ibinx<=hist->GetNbinsX(); ibinx++){//evaluate the number of disabled chips
            if(hist->GetBinContent(ibinx,ibiny)>0) chipswithdeadpix_lower++;
          }

          if(chipswithdeadpix_upper>0){
            cout<<"Layer: "<<layer<<" Stave HS-upper: "<<ibiny-1<<" Run: "<<myAnalysis.getRunNumber(hist)<<" --> # Chips with dead pix: "<<chipswithdeadpix_upper<<endl;
            staveswithdead[stoi(layer)][0]++;
          }
          if(chipswithdeadpix_lower>0){
            cout<<"Layer: "<<layer<<" Stave HS-lower: "<<ibiny-1<<" Run: "<<myAnalysis.getRunNumber(hist)<<" --> # Chips with dead pix: "<<chipswithdeadpix_lower<<endl;
            staveswithdead[stoi(layer)][0]++;
          }
          trend[stoi(layer)][ibiny-1][0]->SetPoint(irun, irun, hproj->Integral(0,hist->GetNbinsX()/2));
          trend[stoi(layer)][ibiny-1][1]->SetPoint(irun, irun, hproj->Integral(hist->GetNbinsX()/2+1,hist->GetNbinsX()));
          if(hproj->Integral(0,hist->GetNbinsX()/2)>maxtot)
            maxtot=hproj->Integral(0,hist->GetNbinsX()/2);
          if(hproj->Integral(hist->GetNbinsX()/2+1,hist->GetNbinsX())>maxtot)
            maxtot=hproj->Integral(hist->GetNbinsX()/2+1,hist->GetNbinsX());
        }

        // Set style IB
        if (stoi(layer) <=2){
          if((ibiny-1)<hist->GetNbinsY()/2)
            SetStyle(trend[stoi(layer)][ibiny-1][0], col[ibiny-1], 24);
          else
            SetStyle(trend[stoi(layer)][ibiny-1][0], col[ibiny-1-hist->GetNbinsY()/2], 26);
        }
        //Set Style OB
        for(int hs=0; hs<=1; hs++){
          if(stoi(layer) ==3 || stoi(layer)==4){
            if((ibiny-1)<hist->GetNbinsY()/3)
              SetStyle(trend[stoi(layer)][ibiny-1][hs], col[ibiny-1], 24);
            else if ((ibiny-1)<hist->GetNbinsY()*2/3)
              SetStyle(trend[stoi(layer)][ibiny-1][hs], col[ibiny-1-(hist->GetNbinsY()*1/3)], 26);
            else if ((ibiny-1)<hist->GetNbinsY()*3/3)
              SetStyle(trend[stoi(layer)][ibiny-1][hs], col[ibiny-1-(hist->GetNbinsY()*2/3)], 25);
          }
          else if (stoi(layer) ==5 || stoi(layer)==6){
            if((ibiny-1)<hist->GetNbinsY()/4)
              SetStyle(trend[stoi(layer)][ibiny-1][hs], col[ibiny-1], 24);
            else if ((ibiny-1)<hist->GetNbinsY()*2/4)
              SetStyle(trend[stoi(layer)][ibiny-1][hs], col[ibiny-1-(hist->GetNbinsY()*1/4)], 26);
            else if ((ibiny-1)<hist->GetNbinsY()*3/4)
              SetStyle(trend[stoi(layer)][ibiny-1][hs], col[ibiny-1-(hist->GetNbinsY()*2/4)], 25);
            else if ((ibiny-1)<hist->GetNbinsY()*4/4)
              SetStyle(trend[stoi(layer)][ibiny-1][hs], col[ibiny-1-(hist->GetNbinsY()*3/4)], 30);
          }
        }
      }
      irun--;
    }
  }

  int npoints = myAnalysis.nRuns();
  TH1F *hfake = new TH1F("hfake", "; # Dead Pixels", npoints, -0.5, (double)npoints-0.5);
  for(int ir=0; ir<npoints; ir++) // Set bin labels to run numbers
      hfake->GetXaxis()->SetBinLabel(npoints-(ir), Form("run%06d",stoi(myAnalysis.Runs()[ir])));

  for (string layer : laynums){ // loop over layers
    if(stoi(layer)<=2){
      TCanvas *canvas = new TCanvas();
      canvas->cd();
      canvas->SetLogy();
      canvas->SetTickx();
      canvas->SetTicky();
      canvas->SetMargin(0.0988,0.1,0.194,0.0993);
      TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
      hfake->GetYaxis()->SetRangeUser(8e-1, maxtot+0.5*maxtot);
      hfake->GetXaxis()->SetTitleOffset(2.8);
      hfake->SetStats(0);
      hfake->SetTitle(Form("Layer-%s (%i Staves with DeadPixels), from run%s to run%s; Run number; #Dead pixels",layer.c_str(),(int)ceil(staveswithdead[stoi(layer)][0]/nRuns),runNumbers.back().c_str(),runNumbers[0].c_str()));
      hfake->Draw();
      for(int istave=0; istave<myAnalysis.stavesInLayer(stoi(layer)); istave++){
        leg->AddEntry(trend[stoi(layer)][istave][0], Form("Stv%d",istave), "p");
        trend[stoi(layer)][istave][0]->Draw("P same");   
      }
      leg->Draw("same");
      if(ccdb_upload){
	string Runperiod = Form("from_run%s_to_run%s",runNumbers.back().c_str(),runNumbers[0].c_str());		
 //   string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str()); //This should be used for actual data 
	canvas->SetName(Form("Layer%s_Dead_pixels",layer.c_str()));
	auto mo1= std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%s",layer.c_str()), TaskClass, DetectorName,std::stoi(runNumbers.back()),Runperiod);
        mo1->setIsOwner(false);
        ccdb->storeMO(mo1);
}
	canvas->SaveAs(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.pdf", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      canvas->SaveAs(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.root", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      delete canvas;
      delete leg;
    }

    if(stoi(layer)>=3){
      for(int hs=0; hs<=1; hs++){
        TCanvas *canvas = new TCanvas();
        canvas->cd();
        canvas->SetLogy();
        canvas->SetTickx();
        canvas->SetTicky();
        canvas->SetMargin(0.0988,0.1,0.194,0.0993);
        TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
        leg->SetNColumns(2);
        hfake->GetYaxis()->SetRangeUser(8e-1, maxtot+0.5*maxtot);
        hfake->GetXaxis()->SetTitleOffset(2.8);
        hfake->SetStats(0);
        hfake->Draw();
        for(int istave=0; istave<myAnalysis.stavesInLayer(stoi(layer)); istave++){
          leg->AddEntry(trend[stoi(layer)][istave][hs], Form("Stv%d",istave), "p");
          trend[stoi(layer)][istave][hs]->Draw("P same");      
        }
        leg->Draw("same");
        if(hs==0){
          hfake->SetTitle(Form("Layer-%s HS-upper (%i Staves with DeadPixels), from run%s to run%s;Run number; #Dead pixels",layer.c_str(),(int)ceil(staveswithdead[stoi(layer)][0]/nRuns),runNumbers.back().c_str(),runNumbers[0].c_str()));
if(ccdb_upload){
        string Runperiod = Form("from_run%s_to_run%s",runNumbers.back().c_str(),runNumbers[0].c_str());
//   string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str()); //This should be used for actual data 
        canvas->SetName(Form("Layer%s_HS_Upper_Dead_pixels",layer.c_str()));
        auto mo2= std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%s",layer.c_str()), TaskClass, DetectorName,std::stoi(runNumbers.back()),Runperiod);
        mo2->setIsOwner(false);
        ccdb->storeMO(mo2);
}
         canvas->SaveAs(Form("../Plots/Layer%s_HS-upper_DeadPixels_run%s-run%s.pdf", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      canvas->SaveAs(Form("../Plots/Layer%s_HS-upper_DeadPixels_run%s-run%s.root", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
        }
        if(hs==1){
          hfake->SetTitle(Form("Layer-%s HS-lower (%i Staves with DeadPixels), from run%s to run%s;Run number; #Dead pixels",layer.c_str(),(int)ceil(staveswithdead[stoi(layer)][1]/nRuns),runNumbers.back().c_str(),runNumbers[0].c_str()));
if(ccdb_upload){
        string Runperiod = Form("from_run%s_to_run%s",runNumbers.back().c_str(),runNumbers[0].c_str());  
//   string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str()); //This should be used for actual data 
        canvas->SetName(Form("Layer%s_HS_Lower_Dead_pixels",layer.c_str()));
        auto mo3= std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%s",layer.c_str()), TaskClass, DetectorName,std::stoi(runNumbers.back()),Runperiod);
        mo3->setIsOwner(false);
        ccdb->storeMO(mo3);
}
         canvas->SaveAs(Form("../Plots/Layer%s_HS-lower_DeadPixels_run%s-run%s.pdf", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      canvas->SaveAs(Form("../Plots/Layer%s_HS-lower_DeadPixels_run%s-run%s.root", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
        }
        delete canvas;
        delete leg;
      }
    }
  }

  gStyle->SetPalette(1);
  for (string layer : laynums){ // loop over layers
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    auto histos = myAnalysis.loadLayer(stoi(layer));
    for(int i = histos.size()-1; i >= 0; i--){ //loop over number of histograms
      auto hist = histos[i];
      hist->Draw("colz");
      hist->SetMinimum(8e-1);
      hist->SetMaximum(maxtot+0.5*maxtot);
      hist->SetTitle(Form("Layer%s Run%s (%i/%i);Chip Number; Stave Number",layer.c_str(),myAnalysis.getRunNumber(hist).c_str(),(int)histos.size()-i,(int)histos.size()));
      hist->GetZaxis()->SetTitle("# Dead Pixels");
      // Save frames of GIF
      canvas->Print(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.gif+40", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      // Save last frame as gif++ so gif loops
      if (i==0)canvas->Print(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.gif++", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
    }
  }
//Disconnencting from the database
 ccdb->disconnect();

} // end of analyzeLayerDeadPixels()

