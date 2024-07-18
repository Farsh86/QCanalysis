//#include "inc/itsAnalysis.hh"
//#include "QualityControl/PostProcessingInterface.h"
//#include "QualityControl/Reductor.h"
//#include "QualityControl/DatabaseFactory.h"
//#include "QualityControl/RootClassFactory.h"
//#include "QualityControl/DatabaseInterface.h"
//#include "QualityControl/MonitorObject.h"
//#include "QualityControl/QcInfoLogger.h"
//#include "QualityControl/CcdbDatabase.h"
//#include "inc/ccdb.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TDirectory.h"
#include "vector"
#include "string"
#include "set"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include <TSystem.h>
#include <iostream>
#include <filesystem>
#include <list>
//include "stdio.h"



//using namespace o2::quality_control::repository;
//using namespace o2::quality_control::core;
using namespace std;



struct LayerParameters
{
	string name;
	int nFirstStave;
	int nLastStave;
	int nChipsPerStave;
	int nStavesForLayer;
};




std::pair <Double_t, Double_t> GetAvgThresholdsForStaveAndNumberOfFailedNChip(TH2F* h, Int_t stvNum, Int_t numChip) {
	std::pair <Double_t, Double_t> avgFailed; 
	TH1D* hproj = h->ProjectionX("proj", stvNum, stvNum);
	Double_t sum = 0;
	Int_t cNumChip = numChip;
	Int_t failedNumChip = 0;
	for (int i = 1; i <= numChip; i++) {
		Double_t forChip = hproj->GetBinContent(i);
		if (forChip == 0) {
			cNumChip--;
			failedNumChip++;
		}
		sum += forChip;
	}
	Double_t avg = 0;
	if (cNumChip > 0) {
		avg = sum / cNumChip;
	}
	avgFailed.first = avg;
	avgFailed.second = failedNumChip;
	return avgFailed;
}



set <string> RemoveEmptyRuns(set <string> allRunNumbers, set <string> toDelete) {
	for (auto run : toDelete) {
		allRunNumbers.erase(run);
	}
	return allRunNumbers;
}



TGraph* AddStyle(TGraph* result, Int_t color, Int_t style) {
	result->SetMarkerColor(color);
	result->SetMarkerStyle(style);
	result->SetMarkerSize(1.5);
	return result;
}




set <string> GetRunNumbers(TFile* analyzefile)
{
	set <string> result;
	TList* listKeys = analyzefile->GetListOfKeys();
	for (TObject* key : *listKeys) {
		string fileName = key->GetName();
		string sRunNumber = fileName.substr(fileName.find("run") + 3, 6);
		result.insert(sRunNumber);
	}
	return result;
}


void PrintGraphToPDF(int numLayer, TCanvas* canvas, TMultiGraph* graph, TLegend* legend, string nameFile) {
	canvas->SetBottomMargin(0.19);
	canvas->SetTopMargin(0.075);
	string plotsDir = "../Plots/";
	string outputDirectoryPath = plotsDir + nameFile;
	if (!filesystem::exists(plotsDir)) {
		filesystem::create_directory(plotsDir);
	}
	if (!filesystem::exists(outputDirectoryPath)) {
		filesystem::create_directory(outputDirectoryPath);
	}
	if (numLayer == 0) {
		canvas->Print((outputDirectoryPath + Form("/%s.pdf[", nameFile.c_str())).c_str());
	}
	graph->Draw("AP");
	legend->Draw();
	canvas->Draw();
	canvas->Print((outputDirectoryPath + Form("/%s.pdf", nameFile.c_str())).c_str());
	canvas->Print((outputDirectoryPath + Form("/%s_Layer%i.png", nameFile.c_str(), numLayer)).c_str());
	if (numLayer == 10) {
		canvas->Print((outputDirectoryPath + Form("/%s.pdf]", nameFile.c_str())).c_str());
	}
}

vector <LayerParameters> init() {
	LayerParameters l0 = { "IB", 1, 12, 9, 1 };
	LayerParameters l1 = { "IB", 13, 28, 9, 13 };
	LayerParameters l2 = { "IB", 29, 48, 9, 29 };
	LayerParameters l3Top = { "ML", 1, 12, 112, 1 };
	LayerParameters l3Bottom = { "ML", 13, 24, 112, 1 };
	LayerParameters l4Top = { "ML", 25, 39, 112, 25 };
	LayerParameters l4Bottom = { "ML", 40, 54, 112, 25 };
	LayerParameters l5Top = { "OL", 1, 21, 196, 1 };
	LayerParameters l5Bottom = { "OL", 22, 42, 196, 1 };
	LayerParameters l6Top = { "OL", 43, 66, 196, 43 };
	LayerParameters l6Bottom = { "OL", 67, 90, 196, 43 };

	vector <LayerParameters> config = { l0, l1, l2, l3Top, l3Bottom, l4Top, l4Bottom, l5Top, l5Bottom, l6Top, l6Bottom };

	return config;
}



void AnalyzeAllLayerNoiseThresholds() {
	
	vector <LayerParameters> config = init();

	vector <Int_t> colorConfig = { kRed, kGreen, kBlue, kMagenta, kOrange };
	vector <Int_t> styleConfig = { 24,25,26,30,28,32 };
	string filepath;
	cout << endl << "Please enter the full path to the file to analyze" << endl;
	gSystem->Exec("ls ../Data/*all-layers* -Art | tail -n 500");
	cin >> filepath;
	string runNumberSuffix = filepath.substr(filepath.find("_run"), filepath.size() - filepath.find("_run") - 5);
	
	TFile* analyzefile = new TFile(filepath.c_str());
	set <string> allRunNumbers = GetRunNumbers(analyzefile);
	for (int numLayer = 0; numLayer <= 10; numLayer++) {
		LayerParameters configToUse = config[numLayer];
		TMultiGraph* avgLayerThreshold = new TMultiGraph();
		
		TLegend* legend = new TLegend(1, 0.19, 0.9, 0.925);
		
		int colorPointer = 0;
		int stylePointer = 0;
		set <string> emptyRunNumbers = {};
		//emptyRunNumbers.insert("551199");
		for (int currentStv = configToUse.nFirstStave; currentStv <= configToUse.nLastStave; currentStv++) {
			TGraph* stvGraph = new TGraph();
			
			Double_t runSequenceNumber = 1;
			for (auto runNumber : allRunNumbers) {
			//	if (runNumber != "551199") {
					TList* listKeys = analyzefile->GetListOfKeys();
					for (TObject* key : *listKeys) {
						string fileName = key->GetName();
						if (fileName.find(configToUse.name) != string::npos && fileName.find(runNumber) != string::npos) {
							TH2F* h = (TH2F*)analyzefile->Get(key->GetName());
							if (h->Integral() == 0) {
								emptyRunNumbers.insert(runNumber);
								continue;
							}
							std::pair <Double_t, Double_t> avgFailed = GetAvgThresholdsForStaveAndNumberOfFailedNChip(h, currentStv, configToUse.nChipsPerStave);
							Double_t avgThrForStv = avgFailed.first;
							if (avgThrForStv != 0) {
								stvGraph->AddPoint(runSequenceNumber, avgThrForStv);
							}
							runSequenceNumber++;
						}
					}
				}
			//}
			stvGraph = AddStyle(stvGraph, colorConfig[colorPointer], styleConfig[stylePointer]);
		
			int colorConfigSize = colorConfig.size() - 1;
			if (colorPointer == colorConfigSize) {
				colorPointer = 0;
				stylePointer = stylePointer + 1;
			}
			else {
				colorPointer = colorPointer + 1;
			}
			legend->AddEntry(stvGraph, Form("Stv%d", currentStv - configToUse.nStavesForLayer), "p");
			avgLayerThreshold->Add(stvGraph);
			
		}
		string layerName = (numLayer < 3) ? to_string(numLayer) : to_string((numLayer + 1) / 2 + 1) + ((numLayer % 2 == 0) ? "Bottom" : "Top");

		avgLayerThreshold->SetTitle(Form("Average noise of thresholds for stave at layer-%s", layerName.c_str()));
		TAxis* X = avgLayerThreshold->GetXaxis();

		
		set <string> notEmptyRuns = RemoveEmptyRuns(allRunNumbers, emptyRunNumbers);
		vector <string> forBigData;
		int notEmptyRunSize = notEmptyRuns.size(); 
		int k = 0;
		if (notEmptyRunSize > 40) {
			int axisRegisterCoef = (notEmptyRunSize / 40) - 1;
			for (auto i : notEmptyRuns) {
				if (k % (2 + axisRegisterCoef) != 0) {
					forBigData.push_back(i);
				}
				k++;
			}
			
			X->SetNdivisions(notEmptyRunSize, 2 + axisRegisterCoef, 0, true);
			Double_t runSequenceNumber = 2;
			X->ChangeLabel(-1, 315, 0, 10, -1, -1, "");
			X->ChangeLabel(1, 315, 0, 10, -1, -1, "");
			//X->ChangeLabel(2, 315, 0, 10, -1, -1, "");
			
			for (auto runNum : forBigData) {
					int runNumber = stoi(runNum);
					X->ChangeLabel(runSequenceNumber, 300, -1, 10, -1, -1, Form("run%i", runNumber));
				
					runSequenceNumber++;
			}
		}
		else {
			
			X->SetNdivisions((notEmptyRunSize * 2));
			Double_t runSequenceNumber = 2;
			X->ChangeLabel(1, 315, 0, 10, -1, -1, "");
			X->ChangeLabel(notEmptyRunSize + 2, 315, 0, 10, -1, -1, "");
			
			for (auto runNum : notEmptyRuns) {
				int runNumber = stoi(runNum);
				X->ChangeLabel(runSequenceNumber, 300, -1, 10, -1, -1, Form("run%i", runNumber));
			
				runSequenceNumber++;
			}
		}
		avgLayerThreshold->GetYaxis()->SetTitle("Average noise of thresholds for stave");
		
		//string nameFileAvg = Form("AvgNoiseThr_without_551199_%s", runNumberSuffix.c_str());
		
		string nameFileAvg = Form("AvgNoiseThr_%s", runNumberSuffix.c_str());
	
		TCanvas* canvas = new TCanvas(Form("c_L%s", layerName.c_str()), Form("c_L%s", layerName.c_str()), 0, 0, 640, 480);
		PrintGraphToPDF(numLayer, canvas, avgLayerThreshold, legend, nameFileAvg);
		canvas->Close();
	}
}