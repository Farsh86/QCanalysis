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

using namespace std;
int ns[] = {12,16,20,24,30,42,48}; // n staves for each layer

void DoAnalysis(string filepath);
int getlayer(int ifee){
  return (ifee) < 36 ? 0 : (ifee) < 84 ? 1 : (ifee) < 144 ? 2 : (ifee) < 192 ? 3 : (ifee) < 252 ? 4 : (ifee) < 336 ? 5 : 6;
}

int getstave(int lay, int ifee){
  if(!lay) return (ifee % (3*ns[lay])) / 3;
  int ntot = ifee;
  for(int ilay=0; ilay<lay; ilay++)
    ntot-=ilay<3 ? 3*ns[ilay] : 2*ns[ilay];
 return lay<3 ? ntot / 3 :  ntot / 2;
}


//
// MAIN
//
void DumpLaneStatusFlag(){

  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*LaneStatus* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;
  
  //Call
  DoAnalysis(fpath);

}

//
// Analyse data
//
void DoAnalysis(string filepath){

	gStyle->SetOptStat(0000);

	//std::vector<TH2*> herr, hfault, hok, hwarning;
	std::vector<TH2*> herr;
	std::vector<TH2*> hfault;
	std::vector<TH2*> hok;
	std::vector<TH2*> hwarning;
	//std::vector<string> timestamps, runnumbers;
	std::vector<string> timestamps1, runnumbers1, timestamps2, runnumbers2, timestamps3, runnumbers3, timestamps4, runnumbers4;

	//Read the file and the list of plots with entries
	TFile *infile=new TFile(filepath.c_str());
	TList *list = infile->GetListOfKeys();
	TKey *key;
	TObject *obj;
	TIter next(list);
	TH2 *h2err;
	TH2 *h2fault;
	TH2 *h2warning;
	while((key = ((TKey*)next()))){
		obj = key->ReadObj();
		if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
				&& (!obj->InheritsFrom("TH2"))
				&& (!obj->InheritsFrom("TH1"))
			 ) {
			cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
		}
		string objname = (string)obj->GetName();
		cout<<"objname = "<<objname<<endl;

		string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
		string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);

		if(objname.find("LSerror")!=string::npos){
			h2err = (TH2*)obj;
			cout<<"... Reading "<<obj->GetName()<<endl;
			h2err->ClearUnderflowAndOverflow();
			herr.push_back(h2err);
			timestamps1.push_back(timestamp);
			runnumbers1.push_back(runnum);
		}
		if(objname.find("LSfault")!=string::npos){
			h2fault = (TH2*)obj;
			cout<<"... Reading "<<obj->GetName()<<endl;
			h2fault->ClearUnderflowAndOverflow();
			hfault.push_back(h2fault);
			timestamps2.push_back(timestamp);
			runnumbers2.push_back(runnum);
		}
		if(objname.find("LSwarning")!=string::npos){
			h2warning = (TH2*)obj;
                        cout<<"... Reading "<<obj->GetName()<<endl;
			h2warning->ClearUnderflowAndOverflow();
			hwarning.push_back(h2warning);
			timestamps4.push_back(timestamp);
			runnumbers4.push_back(runnum);
		}
	}
	//dump lanes in error 
	std::ofstream flerr(Form("../Plots/lane_dump_error_run%s_to_run%s.txt", runnumbers1[runnumbers1.size()-1].c_str(), runnumbers1[0].c_str()));
	for(int i=0; i<(int)herr.size(); i++){
		for(int ifee = 0; ifee < herr[i]->GetNbinsX(); ifee++) {
			for(int ilane = 0; ilane < herr[i]->GetNbinsY(); ilane++) {
				double binc = herr[i]->GetBinContent(ifee+1, ilane+1);
     			        if(binc>0){
        				int layer = getlayer(ifee);
        				int stave = getstave(layer,ifee);
        				flerr<<"RUN"<<runnumbers1[i]<<" "<<"FEEID"<<ifee<<" --> "<< "L"<<layer<<"_"<<stave<<" lane "<<ilane<<"\n";
      				}
			}
		}
	}
	flerr.close();

	//dump lanes in fault 
	std::ofstream flfault(Form("../Plots/lane_dump_fault_run%s_to_run%s.txt", runnumbers2[runnumbers1.size()-1].c_str(), runnumbers2[0].c_str()));
	for(int i=0; i<(int)hfault.size(); i++){
		for(int ifee = 0; ifee < hfault[i]->GetNbinsX(); ifee++) {
			for(int ilane = 0; ilane < hfault[i]->GetNbinsY(); ilane++) {
				double binc = hfault[i]->GetBinContent(ifee+1, ilane+1);
     			        if(binc>0){
        				int layer = getlayer(ifee);
        				int stave = getstave(layer,ifee);
        				flfault<<"RUN"<<runnumbers2[i]<<" "<<"FEEID"<<ifee<<" --> "<< "L"<<layer<<"_"<<stave<<" lane "<<ilane<<"\n";
      				}
			}
		}
	}
	flfault.close();

	//dump lanes in warning
	std::ofstream flwarn(Form("../Plots/lane_dump_warning_run%s_to_run%s.txt", runnumbers4[runnumbers1.size()-1].c_str(), runnumbers4[0].c_str()));
	for(int i=0; i<(int)hwarning.size(); i++){
		for(int ifee = 0; ifee < hwarning[i]->GetNbinsX(); ifee++) {
			for(int ilane = 0; ilane < hwarning[i]->GetNbinsY(); ilane++) {
				double binc = hwarning[i]->GetBinContent(ifee+1, ilane+1);
     			        if(binc>0){
        				int layer = getlayer(ifee);
        				int stave = getstave(layer,ifee);
        				flwarn<<"RUN"<<runnumbers4[i]<<" "<<"FEEID"<<ifee<<" --> "<< "L"<<layer<<"_"<<stave<<" lane "<<ilane<<"\n";
      				}
			}
		}
	}
        flwarn.close();
} 
