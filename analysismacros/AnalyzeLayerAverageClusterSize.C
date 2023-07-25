#include <string>
#include <iostream>
#include <vector>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TKey.h>
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/CcdbDatabase.h"
#include "inc/ccdb.h"

//using namespace o2::framework;
using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;
using namespace std;

void SetStyle(TGraph *gObject, int color, Style_t marker);
void DoAnalysis(string sFile_Path, string sSkip_Runs, int IBorOB, bool ccdb_upload);

//___________________________________________________
// MAIN
void AnalyzeLayerAverageClusterSize(){

   string sFile_Path;
   // int nchips = 9;
   cout << "\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n" << endl;
   gSystem->Exec("ls ../Data/*ClusterTask* -Art | tail -n 500"); //@ prints list of available files for this task
   cout << "\nCopy file name: ";
   cin >> sFile_Path; //@ put path to file(s)
   cout << endl;

   int IBorOB;
   bool ccdb_upload;

   if(sFile_Path.find("IB") != string::npos){
      IBorOB = 0;
   } else if (sFile_Path.find("OB") != string::npos){
      IBorOB = 1;
   } else if (sFile_Path.find("all") != string::npos){
      IBorOB = 2;
   } else {
      // string layernum = sFile_Path.substr(sFile_Path.find("Layer") + 5, 1);
      IBorOB = 2;
   }

   string sSkip_Answer, sSkip_Runs, CCDB_up;
   cout << endl;
   cout << "Would you like to skip some run(s)? [y/n]";
   cin  >> sSkip_Answer;

   if(sSkip_Answer == "y" || sSkip_Answer == "Y"){
      cout << endl;
      cout << "Specify run number(s) separated by comma (no white spaces!):";
      cin  >> sSkip_Runs;
      cout << endl;
   } else {
      sSkip_Runs = " ";
   }

   cout << "Would you like to upload the output to ccdb? [y/n]";
   cin  >> CCDB_up;
   cout << endl;

   if(CCDB_up =="y" || CCDB_up =="Y") ccdb_upload = true;
   else ccdb_upload = false;

   if(ccdb_upload) SetTaskName(__func__);

   //Call
   DoAnalysis(sFile_Path, sSkip_Runs, IBorOB, ccdb_upload);
}

//___________________________________________________
//Set Style
void SetStyle(TGraph *gObject, int color, Style_t marker){
   gObject->SetLineColor(color);
   gObject->SetMarkerStyle(marker);
   gObject->SetMarkerSize(1.4);
   gObject->SetMarkerColor(color);
}

//___________________________________________________
// Analyse data
void DoAnalysis(string sFile_Path, string sSkip_Runs, int IBorOB, bool ccdb_upload){

   //___________________________________________________
   //Auxiliary variables
   gStyle->SetOptStat(0000);
   int nChips {0};
   int nTimes {0};
   int nRunsTot {0};
   int nInputLayers {1};
   int numberOfRuns_PerLayer[7] { -1, -1, -1, -1, -1, -1, -1 };

   vector<TH2*> vInput_Histo;
   vector<string> /*timestamps, */ vNumberOfRuns, vLayerNumber; // timestamps is not used
   int color[] {810, 807, 797, 827, 417, 841, 868, 867, 860, 602, 921, 874};

   //___________________________________________________
   //Setting up the connection to the ccdb database

   // CcdbDatabase* ccdb;
   // if(ccdb_upload) ccdb = SetupConnection();	~To-Do- Currently not working
   unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");
   auto ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());
   ccdb->connect(ccdbport.c_str(), "", "", "");

   //___________________________________________________
   //Read the file and the list of plots with entries
   TFile *Input_File = new TFile(sFile_Path.c_str(), "READ");

   TList *Input_List = Input_File->GetListOfKeys();
   TKey *key;
   TObject *Input_Object;
   TIter next(Input_List);

   while((key = ((TKey*)next()))){
      Input_Object = key->ReadObj();

      if((strcmp(Input_Object->IsA()->GetName(), "TProfile") != 0) && (!Input_Object->InheritsFrom("TH2")) && (!Input_Object->InheritsFrom("TH1"))){
         cout << "<W> Object " << Input_Object->GetName() << " is not 1D or 2D histogram : will not be converted" << endl;
      }

      string sDummy = (string) Input_Object->GetTitle();
      if(sDummy.find("Average Cluster Size") == string::npos) continue;

      string objname = (string) Input_Object->GetName();
      cout << "... Reading " << Input_Object->GetName() << " with title \"" << Input_Object->GetTitle() << "\"" << endl;

      //Extraction of "run number", "time stamp" and "layer number"
      // string timestamp = objname.find("run") == string::npos ? objname.substr(objname.find("_", 2) + 1, 13) : objname.substr(objname.find("_", 6) + 1, 13);
      string sRunNumber   = objname.find("run") == string::npos ? "norun" : objname.substr(objname.find("run") + 3, 6);
      string sLayerNumber = objname.substr(objname.find("L") + 1, 1);
      if(sSkip_Runs.find(sRunNumber) != string::npos) continue; //eventually skip runs if specified

      // timestamps.push_back(timestamp);
      vInput_Histo.push_back((TH2*) Input_Object);
      vNumberOfRuns.push_back(sRunNumber);
      vLayerNumber.push_back(sLayerNumber);

      // cout<<"run: "<<sRunNumber<<"   timestamp: "<<timestamp<<"    sLayerNumber: "<<sLayerNumber<<endl;

      nTimes++; //@ suggestion: needed to separate contribution for a given layer from different runs
      if(numberOfRuns_PerLayer[stoi(sLayerNumber)] == -1) numberOfRuns_PerLayer[stoi(sLayerNumber)] = 0; //@ by default "-1"
      if(nTimes > 1){
         if(sLayerNumber == vLayerNumber[vLayerNumber.size() - 2]){ //@ if previous layer number was the same or not. If not, it means we have at least 2 layers
            numberOfRuns_PerLayer[stoi(sLayerNumber)]++;
         } else {
            nInputLayers++; //@ difinition comes from name
         }
      }
   }

   if (nTimes == 0) {cout << "\nInput files contains no plots\n" << endl; return;}

   //___________________________________________________
   int nLayers;
   if(nInputLayers == 1){
      nLayers = 1;
   } else {
      switch(IBorOB){
         case 0:  nLayers = 3; break;
         case 1:  nLayers = 4; break;
         default: nLayers = 7;
      }
   }

   //___________________________________________________
   TGraph *gAverageClusterSize[nLayers][50][2]; //@ [second index] -> number of bins on Y-axis (or number of staves), [third index] -> number of halfs (for OB 2 halfs)

   int ilayer = 0;
   for(int ihist = (int) vInput_Histo.size() - 1; ihist >= 0; ihist--){ //@ from old to new
      if(IBorOB == 1){
         ilayer = stoi(vLayerNumber[ihist]) - 3;
      } else {
         ilayer = stoi(vLayerNumber[ihist]);
      }

      if (nInputLayers == 1) ilayer = 0;
      // cout << "ilayer " <<ilayer << endl;

      for(int ibiny = 1; ibiny <= vInput_Histo[ihist]->GetNbinsY(); ibiny++){ //@ number of bins correspond to # of staves
         for (int b = 0; b <= 1; b++){
            gAverageClusterSize[ilayer][ibiny-1][b] = new TGraph();
         }
      }
   } // @ end of loop over histograms

   //___________________________________________________
   int irun = 0;
   int ChipMin {1};
   int ChipMax {1};
   int nStaveParts {0};

   for(int ihist = (int) vInput_Histo.size() - 1; ihist >= 0; ihist--){ //start from the last in order to have the runs from the oldest to the newest

      //@ can be discarded since the same code is written at line 247
      if(IBorOB == 1){
         ilayer = stoi(vLayerNumber[ihist]) - 3;
      } else {
         ilayer = stoi(vLayerNumber[ihist]);
      }
      if (nInputLayers == 1) ilayer = 0;

      if(stoi(vLayerNumber[ihist]) < 3){
         nChips = vInput_Histo[ihist]->GetNbinsX(); //IB layers
      } else {
         nChips = vInput_Histo[ihist]->GetNbinsX()/2.0; //@ division of OB layers into 2 parts
      }

      if(stoi(vLayerNumber[ihist]) >= 3){
         nStaveParts = 1;
      } else {
         nStaveParts = 0;
      }

      for(int ibiny = 1; ibiny <= vInput_Histo[ihist]->GetNbinsY(); ibiny++){ //loop on y bins (stave)s
         //cout << "\n stave number: " << ibiny-1 << endl;
         int nDeadchips = 0;

         for(int StavePart = 0; StavePart <= nStaveParts; StavePart++){ //loop over the two Half Staves for OB

            gAverageClusterSize[ilayer][ibiny-1][StavePart]->SetName(Form("gr_L%s_stave%d_HS%i", vLayerNumber[ihist].c_str(), ibiny - 1, StavePart));

            ChipMin = 1; // @ For IB chips start from 0 to 8
            ChipMax = vInput_Histo[ihist]->GetNbinsX();

            // @ division is performed by number of chips, not staves. What is a reason?
            if(stoi(vLayerNumber[ihist]) >= 3){
               if(StavePart == 0) { //HS Lower
                  ChipMin = 1;
                  ChipMax = vInput_Histo[ihist]->GetNbinsX()/2.0;
               } else { //HS Upper
                  ChipMin = vInput_Histo[ihist]->GetNbinsX()/2.0 + 1;
                  ChipMax = vInput_Histo[ihist]->GetNbinsX();
               }
            }

            nDeadchips = 0;
            for(int ibinx = ChipMin; ibinx <= ChipMax; ibinx++){//evaluate the number of disabled chips
               if(vInput_Histo[ihist]->GetBinContent(ibinx, ibiny) < 1e-15) nDeadchips++; // @ here can be used TH1 *hproj
            }

            if(nDeadchips > 0) cout << "Layer " << vLayerNumber[ihist] << " Stave " << ibiny - 1 << " Run: " << vNumberOfRuns[ihist] << " Number of lanes: " << nChips << " --> Active lanes:" << nChips - nDeadchips << endl;

            double StaveOccupancy = 0.0;

            if(nDeadchips != nChips){
               for(int ibinx = ChipMin; ibinx <= ChipMax; ibinx++){
                  StaveOccupancy += vInput_Histo[ihist]->GetBinContent(ibinx, ibiny);
               }
               gAverageClusterSize[ilayer][ibiny-1][StavePart]->SetPoint(irun, irun, StaveOccupancy/(nChips-nDeadchips)); // @ final output
               // cout << "Stave " << ibiny - 1 << " average cluster occupation " << StaveOccupancy/(nChips-nDeadchips) << endl;
            } else {
               gAverageClusterSize[ilayer][ibiny-1][StavePart]->SetPoint(irun, irun, 0.0);
            }

            double x = 0.0;
            double y = 0.0;

            //___________________________________________________
            if(stoi(vLayerNumber[ihist]) < 3){

               if((ibiny - 1) < vInput_Histo[ihist]->GetNbinsY()/2) {
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1], 24);
               } else {
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1-vInput_Histo[ihist]->GetNbinsY()/2], 26);
               }

            } else if (stoi(vLayerNumber[ihist]) == 3 || stoi(vLayerNumber[ihist]) == 4){

               if((ibiny - 1) < vInput_Histo[ihist]->GetNbinsY()/3){
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1], 24);
               } else if ((ibiny - 1) < vInput_Histo[ihist]->GetNbinsY()*2/3){
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1-vInput_Histo[ihist]->GetNbinsY()/3], 26);
               } else {
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1-vInput_Histo[ihist]->GetNbinsY()*2/3], 25);
               }

            } else if (stoi(vLayerNumber[ihist]) == 5 || stoi(vLayerNumber[ihist]) == 6){
               if((ibiny - 1)<int(vInput_Histo[ihist]->GetNbinsY()/4)){
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1], 24);
               } else if((ibiny - 1) < 2*int(vInput_Histo[ihist]->GetNbinsY()/4)){
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1-int(vInput_Histo[ihist]->GetNbinsY()/4)], 26);
               } else if((ibiny - 1) < 3*vInput_Histo[ihist]->GetNbinsY()/4){
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1-2*int(vInput_Histo[ihist]->GetNbinsY()/4)], 25);
               } else {
                  SetStyle(gAverageClusterSize[ilayer][ibiny-1][StavePart], color[ibiny-1-int(vInput_Histo[ihist]->GetNbinsY()*3/4)], 30);
               }
            }
            // cout << "gAverageClusterSize values for run "<< irun<< " "  << gAverageClusterSize[ilayer][ibiny-1][StavePart]->GetPointY(irun) << endl;
         } //Stavepart end
      } // biny part end
      irun++;
      if(ihist > 0 && vLayerNumber[ihist - 1] != vLayerNumber[ihist]) irun = 0;
   } // end of loop over histograms

   //___________________________________________________
   int ilayEff = 0;
   TH1F *hDummy[7];
   nRunsTot = 0;

   for(int ilay = 0; ilay < nLayers; ilay++){
      if(nLayers == 1){ // @ # of input layers
         ilayEff = stoi(vLayerNumber[0]);
      } else if (IBorOB == 1){ // @ OB part
         ilayEff = ilay + 3;
      } else {
         ilayEff = ilay; //
      }

      if(numberOfRuns_PerLayer[ilayEff] != -1) nRunsTot += numberOfRuns_PerLayer[ilayEff] + 1;

      int npoints = gAverageClusterSize[ilay][0][0]->GetN(); // no sense, since we always fill graph

      // cout << "ilay " << ilay << " ilayEff " << ilayEff << " nRunsTot " << nRunsTot<< endl;

      if(npoints == 0) continue;
      hDummy[ilay] = new TH1F(Form("hClusterOccup_L%i", ilay), "; Run; Average cluster size (pixel/nChips)", npoints, -0.5, (double) npoints - 0.5);

      for(int ir = 0; ir <= numberOfRuns_PerLayer[ilayEff]; ir++){
         hDummy[ilay]->GetXaxis()->SetBinLabel(ir+1, Form("%06d", stoi(vNumberOfRuns[nRunsTot-ir-1]))); // @ run%06d
      }
   }

   //___________________________________________________
   TString sPathOut = Form("../Plots/AverageClusterSize_%s.root", sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root")-sFile_Path.find("from")).c_str());
   TFile *Output_File = new TFile(sPathOut, "RECREATE");

   //___________________________________________________
   //Draw
   nRunsTot = 0;

   for(int ilay = 0; ilay < nLayers; ilay++){
      if(nLayers==1){
         ilayEff = stoi(vLayerNumber[0]);
      } else if(IBorOB==1){
         ilayEff = ilay + 3;
      } else {
         ilayEff = ilay;
      }

      if(ilay > 0 && numberOfRuns_PerLayer[ilayEff-1] != -1) nRunsTot += (numberOfRuns_PerLayer[ilayEff-1] + 1);

      TCanvas *canvas = new TCanvas("canvas", "canvas");
      TCanvas *Secondcanvas = new TCanvas("Secondcanvas", "Secondcanvas");

      if(ilayEff >= 3){
         Secondcanvas->SetName("canvas_HSUpper");
         canvas->SetName("canvas_HSLower");
      }

      //canvas->SetLogy();
      canvas->SetTickx();
      canvas->SetTicky();
      canvas->SetMargin(0.0988,0.15,0.194,0.0993);
      //Secondcanvas->SetLogy();
      Secondcanvas->SetTickx();
      Secondcanvas->SetTicky();
      Secondcanvas->SetMargin(0.0988,0.15,0.194,0.0993);

      Bool_t IsTwoCanvas = 0;
      if(ilayEff >= 3) IsTwoCanvas = 1;

      TLegend *leg = new TLegend(0.857, 0.197,0.997,0.898);
      if(ilayEff >= 3) leg->SetNColumns(2);

      //@@ Stop here
      for(int istave = 0; istave < vInput_Histo[nRunsTot]->GetNbinsY(); istave++){
         leg->AddEntry(gAverageClusterSize[ilay][istave][0], Form("Stv%d",istave), "p");
      }

      hDummy[ilay]->GetYaxis()->SetRangeUser(1, 13);
      hDummy[ilay]->GetXaxis()->SetTitleOffset(2);
      hDummy[ilay]->SetTitle(Form("Layer-%s, %s", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root")-sFile_Path.find("from")).c_str()));
      canvas->cd();

      if(IsTwoCanvas) hDummy[ilay]->SetTitle(Form("Layer-%s, HS Lower, %s", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root")-sFile_Path.find("from")).c_str()));
      hDummy[ilay]->Draw();
      for(int istave = 0; istave < vInput_Histo[nRunsTot]->GetNbinsY(); istave++){
         gAverageClusterSize[ilay][istave][0]->Draw("P same");
         Output_File->WriteTObject(gAverageClusterSize[ilay][istave][0]);
      }

      leg->Draw("same");

      if(!IsTwoCanvas){
         if(ccdb_upload){
            // The number 27 is the sum of the 2*6 digit run numbers+ len("_to_run")+len("from_run")
            string Runperiod = Form("%s",sFile_Path.substr(sFile_Path.find("from"), 27).c_str());
            int RunNumber = std::stoi(sFile_Path.substr(sFile_Path.find("n")+1,6).c_str());
	    string canvas_name = Form("Layer%s_average_cluster_size", vLayerNumber[nRunsTot].c_str());
            canvas->SetName(canvas_name.c_str());
            auto mo1 = std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%s", vLayerNumber[nRunsTot].c_str()), TaskClass, DetectorName, RunNumber, Runperiod); // @ CCDB upload
            mo1->setIsOwner(false);
            ccdb->storeMO(mo1);
         }
         canvas->SaveAs(Form("../Plots/Layer%s_average_cluster_size_%s.pdf", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root") - sFile_Path.find("from")).c_str()));
         canvas->SaveAs(Form("../Plots/Layer%s_average_cluster_size_%s.root", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root") - sFile_Path.find("from")).c_str()));
      } else {
         // The number 27 is the sum of the 2*6 digit run numbers+ len("_to_run")+len("from_run")
         if(ccdb_upload){
            string Runperiod = Form("%s",sFile_Path.substr(sFile_Path.find("from"),27).c_str());
        int RunNumber = std::stoi(sFile_Path.substr(sFile_Path.find("n")+1,6).c_str());
	  string canvas_name2 = Form("Layer%s_average_cluster_size_HSLower", vLayerNumber[nRunsTot].c_str());
            canvas->SetName(canvas_name2.c_str());
            auto mo2 = std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%s", vLayerNumber[nRunsTot].c_str()), TaskClass, DetectorName,RunNumber,Runperiod); // @ CCDB upload
            mo2->setIsOwner(false);
            ccdb->storeMO(mo2);
         }
         canvas->SaveAs(Form("../Plots/Layer%s_average_cluster_size_%s_HSLower.root", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root") - sFile_Path.find("from")).c_str()));
         canvas->SaveAs(Form("../Plots/Layer%s_average_cluster_size_%s_HSLower.pdf", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root") - sFile_Path.find("from")).c_str()));
      }

      if(IsTwoCanvas){
         Secondcanvas->cd();
         hDummy[ilay]->SetTitle(Form("Layer-%s, HS Upper, %s", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root") - sFile_Path.find("from")).c_str()));
         hDummy[ilay]->Draw();

         for(int istave=0; istave<vInput_Histo[nRunsTot]->GetNbinsY(); istave++){
            gAverageClusterSize[ilay][istave][1]->Draw("P same");
            Output_File->WriteTObject(gAverageClusterSize[ilay][istave][1]);
         }
         leg->Draw("same");

         if(ccdb_upload){
            string Secondcanvas_name = Form("Layer%s_average_cluster_size_HSUpper", vLayerNumber[nRunsTot].c_str());

            // The number 27 is the sum of the 2*6 digit run numbers+ len("_to_run")+len("from_run")
            string Runperiod = Form("%s",sFile_Path.substr(sFile_Path.find("from"), 27).c_str());
       int RunNumber = std::stoi(sFile_Path.substr(sFile_Path.find("n")+1,6).c_str());
	     Secondcanvas->SetName(Secondcanvas_name.c_str());
            auto mo3 = std::make_shared<o2::quality_control::core::MonitorObject>(Secondcanvas, TaskName+Form("/Layer%s", vLayerNumber[nRunsTot].c_str()), TaskClass, DetectorName, RunNumber , Runperiod); //  @ CCDB  upload
            mo3->setIsOwner(false);
            ccdb->storeMO(mo3);
         }
         Secondcanvas->SaveAs(Form("../Plots/Layer%s_average_cluster_size_%s_HSUpper.pdf", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root")-sFile_Path.find("from")).c_str()));
         Secondcanvas->SaveAs(Form("../Plots/Layer%s_average_cluster_size_%s_HSUpper.root", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root")-sFile_Path.find("from")).c_str()));
      }

      delete canvas;
      delete Secondcanvas;
      delete leg;
   }
   Output_File->Close();

   //___________________________________________________
   //Make GIF with TH2 for each run and for each layer
   nRunsTot = 0;
   irun = 1;
   gStyle->SetPalette(1);

   for(int ilay = 0; ilay < nLayers; ilay++) {//remove images if they exist already
      if (nLayers==1){
         ilayEff = stoi(vLayerNumber[0]);
      } else if (IBorOB==1){
         ilayEff = ilay + 3;
      } else {
         ilayEff = ilay;
      }

      if(ilay > 0 && numberOfRuns_PerLayer[ilayEff-1] != -1) nRunsTot += (numberOfRuns_PerLayer[ilayEff-1] + 1);
      gSystem->Unlink(Form("../Plots/Layer%s_average_cluster_size_map_%s.gif", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root") - sFile_Path.find("from")).c_str()));
   }

   nRunsTot = vNumberOfRuns.size() - 1; // @

   for(int ihist = (int) vInput_Histo.size()-1; ihist >= 0; ihist--){// start from the last in order to have the runs from the oldest to the newest
      if(ihist != ((int) vInput_Histo.size()-1) && (stoi(vLayerNumber[ihist]) != stoi(vLayerNumber[ihist+1]))){
         if(numberOfRuns_PerLayer[stoi(vLayerNumber[ihist])] != -1) nRunsTot -= (numberOfRuns_PerLayer[stoi(vLayerNumber[ihist])] + 1);
      }

      TCanvas *canvas = new TCanvas();
      canvas->cd();
      canvas->SetLogz();
      canvas->SetTickx();
      canvas->SetTicky();
      canvas->SetRightMargin(0.15);
      vInput_Histo[ihist]->Draw("colz");
      vInput_Histo[ihist]->SetMinimum(1e-1);
      vInput_Histo[ihist]->SetMaximum(1e+2);
      vInput_Histo[ihist]->GetZaxis()->SetTitle("Average cluster size (pixel/nChips)");
      vInput_Histo[ihist]->SetTitle(Form("Layer-%s, Run %06d (%d/%d)", vLayerNumber[nRunsTot].c_str(), stoi(vNumberOfRuns[ihist]), irun, numberOfRuns_PerLayer[stoi(vLayerNumber[ihist])]+1));
      canvas->Print(Form("../Plots/Layer%s_average_cluster_size_map_%s.gif+40", vLayerNumber[nRunsTot].c_str(), sFile_Path.substr(sFile_Path.find("from"), sFile_Path.find(".root")-sFile_Path.find("from")).c_str()));
      irun++;

      if(!ihist && nLayers == 1){
         break;
      }

      if(nLayers > 1){
         if(ihist > 0){
            if(vLayerNumber[ihist-1]!=vLayerNumber[ihist]){
               irun = 1;
            }
         }
      }
      delete canvas;
   }

   cout << "\nROOT file " << sPathOut << " has been created" << endl;

   // Disconnencting the interface
   // if(ccdb_upload)
   ccdb->disconnect();
}
