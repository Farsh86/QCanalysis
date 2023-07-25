#include <string>
#include <iostream>
#include <fstream>
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

//using namespace o2::framework;
using namespace std;

void SetStyle(TGraph *gObject, int color, Style_t marker);
void DoAnalysis(string sFile_Path, string output_file_name);
std::vector<int> GetChipInModList(int layer, int stave, int lane);

//___________________________________________________
// MAIN
void DumpEmptyChipsClusterOccupancy(string output_file_name = "file.txt"){

   string sFile_Path;
   // int nchips = 9;
   cout << "\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n" << endl;
   gSystem->Exec("ls ../Data/*ClusterTask* -Art | tail -n 500"); //@ prints list of available files for this task
   cout << "\nCopy file name: ";
   cin >> sFile_Path; //@ put path to file(s)
   cout << endl;

   //Call
   DoAnalysis(sFile_Path, output_file_name);
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
void DoAnalysis(string sFile_Path, string output_file_name){

   //___________________________________________________
   //Auxiliary variables
   gStyle->SetOptStat(0000);
   int nChips {0};
   int nTimes {0};
   int nRunsTot {0};
   int nInputLayers {1};
   int numberOfRuns_PerLayer[7] { -1, -1, -1, -1, -1, -1, -1 }; //@ name could be wrong

   vector<TH2*> vInput_Histo;
   vector<string> /*timestamps, */ vNumberOfRuns, vLayerNumber; // @ std:: , timestamps is not used

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
      if(sDummy.find("Cluster Occupancy") == string::npos) continue;

      string objname = (string) Input_Object->GetName();
      cout << "... Reading " << Input_Object->GetName() << " with title \"" << Input_Object->GetTitle() << "\"" << endl;

      //Extraction of "run number", "time stamp" and "layer number"
      // string timestamp = objname.find("run") == string::npos ? objname.substr(objname.find("_", 2) + 1, 13) : objname.substr(objname.find("_", 6) + 1, 13);
      string sRunNumber   = objname.find("run") == string::npos ? "norun" : objname.substr(objname.find("run") + 3, 6);
      string sLayerNumber = objname.substr(objname.find("L") + 1, 1);

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

   // Check empty chips
   ofstream outfl (output_file_name.c_str());
   for(int ihist = 0; ihist < (int)vInput_Histo.size(); ihist++){
     for(int istave = 1; istave<=vInput_Histo[ihist]->GetNbinsY(); istave++) {
       for(int ilane = 1; ilane<=vInput_Histo[ihist]->GetNbinsX(); ilane++){
         if(vInput_Histo[ihist]->GetBinContent(ilane, istave) < 1e-8) {//select empty lane or chip
           std::vector<int> o2chipid = GetChipInModList(std::stoi(vLayerNumber[ihist]), istave-1, ilane-1);
           string stavestring = Form("L%s_%02d",vLayerNumber[ihist].c_str(),istave-1);
           string hsstring;
           if(std::stoi(vLayerNumber[ihist])<5) {
             if(ilane-1 <= 8) hsstring = "L";
             else hsstring = "U";
           } else {
             if(ilane-1 < 14) hsstring = "L";
             else hsstring = "U";
           }
           string modulestring;
           if(std::stoi(vLayerNumber[ihist])<3) {
             modulestring = "1";
           } else if (std::stoi(vLayerNumber[ihist])<5){
             if(ilane-1 < 8) modulestring = std::to_string((ilane-1)/2 + 1);
             else modulestring = std::to_string((ilane-1-8)/2 + 1);
           } else {
             if(ilane-1 < 14) modulestring = std::to_string((ilane-1)/2 + 1);
             else modulestring = std::to_string((ilane-1-14)/2 + 1);
           }
           for(int i=0; i<(int)o2chipid.size(); i++){
             outfl<<vNumberOfRuns[ihist]<<" "<<stavestring<<" "<<hsstring<<" "<<modulestring<<" "<<Form("%02d",o2chipid[i])<<"\n";
           }
         } // end if
       } // end loop on lanes
     }//end loop on staves
   }// end loop on hist

   outfl.close();
}

std::vector<int> GetChipInModList(int layer, int stave, int lane) {

  std::vector<int> outlist;

  if(layer<3) {
    outlist.push_back(lane);
  } else {

    if(lane % 2 == 0){
      for(int i=0; i<=6; i++)
        outlist.push_back(i);
    } else {
      for(int i=8; i<=14; i++){
        outlist.push_back(i);
      }
    }
  }

  return outlist;
}
