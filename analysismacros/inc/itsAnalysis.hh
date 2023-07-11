#ifndef itsAnalysis_hh
#define itsAnalysis_hh

#include <string>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TList.h>
#include "TH2.h"
#include "THnSparse.h"
#include "TSystem.h"
#include "TKey.h"

//using namespace std;

//---------------------------------------------------------------------------
// Description:
// Class for loading ROOT files to analyse in ALICE ITS2 Threshold QC
// Questions: bas.hofman@cern.ch
//---------------------------------------------------------------------------

class itsAnalysis {

private :
	vector<string> runNumbers_;
	vector<string> layerNumbers_;
	vector<string> runNumbers_skips;
	vector<string> layerNumbers_skips;

	string skip_ans,skip_runs,skip_layers;

	std::vector<TH2*> hmaps_skips;
	std::vector<THnSparse*> hmaps_skips_sparse;

	std::vector<TH2*> hmaps;
	std::vector<THnSparse*> hmaps_sparse;

public :
	itsAnalysis(string histType) {  // Function to load in a ROOT file
		// We will ask for the file to analyse and load that
		cout<<"\n\n=> Available file(s) for the analysis: \n"<<endl;
		gSystem->Exec("ls ../Data/*THRMAPS_DEADPIXMAPS* -Art | tail -n 500");
		cout<<"\nCopy file name: ";
		
		string filepath;
		cin>>filepath;
		cout<<endl;

		TFile *file = TFile::Open(filepath.c_str());
		if (!file) {std::cout<<"Could not open root file"<<std::endl;}

		TObject *obj; TKey *key;
		TIter next( file->GetListOfKeys());
		while ((key = (TKey *) next())) {
		    obj = file->Get(key->GetName()); // copy object to memory
		    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
		         && (!obj->InheritsFrom("TH2"))
			       && (!obj->InheritsFrom("TH1"))
		         && (!obj->InheritsFrom("THnSparse"))
		       ) { cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;}   
			string objname = (string)obj->GetName();
			string objtitle = (string)obj->GetTitle();
			//if(objname.find("Stv")!=string::npos) break; //hsparse of pixel map has "stave" instead of Stv
			if(objtitle.find(histType)==string::npos) continue;

			if(histType == "THR")	{
				if(objname.find("AvgThrIB")!=string::npos){
					if (std::find(layerNumbers_.begin(), layerNumbers_.end(), to_string(0)) == layerNumbers_.end()) {
						layerNumbers_.push_back(to_string(0)); 
						layerNumbers_.push_back(to_string(1)); 
						layerNumbers_.push_back(to_string(2)); 
					}
					TH2F *h_L0 = new TH2F(("h2_L0_run"+itsAnalysis::getRunNumber((TH2*)obj)).c_str(),"ITS Layer0, Threshold Vs Chip and Stave",9,0,9,  12,0,12);
					TH2F *h_L1 = new TH2F(("h2_L1_run"+itsAnalysis::getRunNumber((TH2*)obj)).c_str(),"ITS Layer1, Threshold Vs Chip and Stave",9,0,9,  16,0,16);
					TH2F *h_L2 = new TH2F(("h2_L2_run"+itsAnalysis::getRunNumber((TH2*)obj)).c_str(),"ITS Layer2, Threshold Vs Chip and Stave",9,0,9,  20,0,20);
					for(int xbin = 0; xbin < 9; xbin++){
						for(int ybin = 0; ybin < 12; ybin++){
							h_L0->SetBinContent(xbin+1,ybin+1,((TH2*)obj)->GetBinContent(xbin+1,ybin+1));
						}
						for(int ybin = 0; ybin < 16; ybin++){
							h_L1->SetBinContent(xbin+1,ybin+1,((TH2*)obj)->GetBinContent(xbin+1,ybin+1+12));
						}
						for(int ybin = 0; ybin < 20; ybin++){
							h_L2->SetBinContent(xbin+1,ybin+1,((TH2*)obj)->GetBinContent(xbin+1,ybin+1+12+16));
						}
					}
					h_L0->SetStats(0); h_L1->SetStats(0);h_L2->SetStats(0);
					hmaps.push_back(h_L0);
					hmaps.push_back(h_L1);
					hmaps.push_back(h_L2);
				}


				if(objname.find("AvgThrML")!=string::npos){
					if (std::find(layerNumbers_.begin(), layerNumbers_.end(), to_string(3)) == layerNumbers_.end()) {
						layerNumbers_.push_back(to_string(3)); 
						layerNumbers_.push_back(to_string(4)); 
					}
					TH2F *h_L3 = new TH2F(("h2_L3_run"+itsAnalysis::getRunNumber((TH2*)obj)).c_str(),"ITS Layer3, Threshold Vs Chip and Stave",8*14,0,8*14,  24,0,24);
					TH2F *h_L4 = new TH2F(("h2_L4_run"+itsAnalysis::getRunNumber((TH2*)obj)).c_str(),"ITS Layer4, Threshold Vs Chip and Stave",8*14,0,8*14,  30,0,30);
					for(int xbin = 0; xbin < 8*14; xbin++){
						for(int ybin = 0; ybin < 24; ybin++){
							h_L3->SetBinContent(xbin+1,ybin+1,((TH2*)obj)->GetBinContent(xbin+1,ybin+1));
						}
						for(int ybin = 0; ybin < 30; ybin++){
							h_L4->SetBinContent(xbin+1,ybin+1,((TH2*)obj)->GetBinContent(xbin+1,ybin+1+24));
						}
					}
					h_L3->SetStats(0); h_L4->SetStats(0);
					hmaps.push_back(h_L3);
					hmaps.push_back(h_L4);
				}


				if(objname.find("AvgThrOL")!=string::npos){
					if (std::find(layerNumbers_.begin(), layerNumbers_.end(), to_string(5)) == layerNumbers_.end()) {
						layerNumbers_.push_back(to_string(5)); 
						layerNumbers_.push_back(to_string(6)); 
					}
					TH2F *h_L5 = new TH2F(("h2_L5_run"+itsAnalysis::getRunNumber((TH2*)obj)).c_str(),"ITS Layer5, Threshold Vs Chip and Stave",14*14,0,14*14,  42,0,42);
					TH2F *h_L6 = new TH2F(("h2_L6_run"+itsAnalysis::getRunNumber((TH2*)obj)).c_str(),"ITS Layer6, Threshold Vs Chip and Stave",14*14,0,14*14,  48,0,48);
					for(int xbin = 0; xbin < 14*14; xbin++){
						for(int ybin = 0; ybin < 42; ybin++){
							h_L5->SetBinContent(xbin+1,ybin+1,((TH2*)obj)->GetBinContent(xbin+1,ybin+1));
						}
						for(int ybin = 0; ybin < 48; ybin++){
							h_L6->SetBinContent(xbin+1,ybin+1,((TH2*)obj)->GetBinContent(xbin+1,ybin+1+42));
						}
					}
					h_L5->SetStats(0); h_L6->SetStats(0);
					hmaps.push_back(h_L5);
					hmaps.push_back(h_L6);
				}
			}

			else{
				//Load in histograms into vector
				if(obj->InheritsFrom("THnSparse"))hmaps_sparse.push_back((THnSparse*)obj);
				else hmaps.push_back((TH2*)obj);

				//Find layerNumbers_ if unique
				string laynum =  objname.find("L")==string::npos ? "":objname.substr(objname.find("L")+1, 1);
				if (std::find(layerNumbers_.begin(), layerNumbers_.end(), laynum) == layerNumbers_.end()) {
				  layerNumbers_.push_back(laynum);
				}
			}

			//Find runNumber if unique
			string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
			if (std::find(runNumbers_.begin(), runNumbers_.end(), runnum) == runNumbers_.end()) {
			  runNumbers_.push_back(runnum);
			}

		} // end of loop over objects

		//Sort Layer numbers:
		sort(layerNumbers_.begin(),layerNumbers_.end());

		cout << "Found run numbers: ";
		for (auto i : runNumbers_){ cout << i <<","; }
		cout << endl << "Found layers: ";
		for (auto i : layerNumbers_){ cout << i <<","; }
		cout << endl<<endl;

		cout<<"Would you like to skip some run(s)? [y/n]";
		cin>>skip_ans;
		if(skip_ans=="y" || skip_ans=="Y"){
		  cout<<"Specify run number(s) separated by comma (no white spaces!):";
		  cin>>skip_runs;
		  cout<<endl;
		}

		for (auto hist: hmaps)	{
			string objname = (string)hist->GetName();
			string run     = objname.substr(objname.find("run")+3, 6);
			if (skip_runs.find(run) == std::string::npos)
				hmaps_skips.push_back((TH2*)hist);
		}

		for (auto hist: hmaps_sparse)	{
			string objname = (string)hist->GetName();
			string run     = objname.substr(objname.find("run")+3, 6);
			if (skip_runs.find(run) == std::string::npos){
				hmaps_skips_sparse.push_back((THnSparse*)hist);
			}
		}

		for (auto run : runNumbers_) {
			if (skip_runs.find(run) == std::string::npos)
				runNumbers_skips.push_back(run);
		}

		cout<<"Would you like to skip some layer(s)? [y/n]";
		cin>>skip_ans;
		if(skip_ans=="y" || skip_ans=="Y"){
		  cout<<"Specify layer number(s) separated by comma (no white spaces!):";
		  cin>>skip_layers;
		  cout<<endl;
		}

		for (auto run : layerNumbers_) {
			if (skip_layers.find(run) == std::string::npos)
				layerNumbers_skips.push_back(run);
		}


	} // end of loadFile()s

	int nChips(int layer) { // Function to return number of chips in a its Layer
		int nChips_ = -999;
		if (layer == 0 || layer == 1 || layer == 2) nChips_ = 9;
		if (layer == 3 || layer == 4) nChips_ = 8;	// Modules instead of chips
		if (layer == 5 || layer == 6) nChips_ = 14; // Modules instead of chips
		return nChips_;
	}

	int stavesInLayer(int layer)	{
		int staves;
		if (layer == 0) staves = 12;
		if (layer == 1) staves = 16;
		if (layer == 2) staves = 20;
		if (layer == 3)	staves = 24;
		if (layer == 4) staves = 30;
		if (layer == 5)	staves = 42;
		if (layer == 6) staves = 48;
		return staves;
	}

	vector<string> Layers() {
		return layerNumbers_skips;
	}

	vector<string> Runs() {
		return runNumbers_skips;
	}

	int nLayers() {
		return layerNumbers_skips.size();
	}

	int nRuns() {
		return runNumbers_skips.size();
	}

	std::vector<TH2*> loadedHists() {
		return hmaps_skips;
	}

	std::vector<TH2*> loadLayer(int layer) {
		std::vector<TH2*> hmaps_per_layer; // empty list, histograms of layer
		for (auto hist: hmaps_skips)	{
			string objname = (string)hist->GetName();
			if (stoi(objname.substr(objname.find("L")+1,1))==layer){ // Check title for layer
				hmaps_per_layer.push_back(hist);
			}
		}
		return hmaps_per_layer;
	}

	std::vector<THnSparse*> loadLayerSparse(int layer) {
		std::vector<THnSparse*> hmaps_per_layer_sparse; // empty list, histograms of layer
		for (auto hist: hmaps_skips_sparse)	{
			string objname = (string)hist->GetName();
			if (stoi(objname.substr(objname.find("L")+1,1))==layer){ // Check title for layer
				hmaps_per_layer_sparse.push_back(hist);
			}
		}
		return hmaps_per_layer_sparse;
	}

	string getRunNumber(TH2* hist) {
		string objname = hist->GetName();
		string runN = objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
		return runN;
	}

	string getRunNumber(THnSparse* hist) {
		string objname = hist->GetName();
		string runN = objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
		return runN;
	}

	int getStaveNumber (THnSparse* hist){
		string hname = hist->GetName();
		string stvnum = hname.substr(hname.find("Stv")+3,2);
		if(stvnum.find("_")!=string::npos){
		  stvnum = hname.substr(hname.find("Stv")+3,1);
		}
		return stoi(stvnum);
	}

}; // end of class def

#endif