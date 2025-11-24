#include <TFileMerger.h>
#include <TGrid.h>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

void DownloadAndMergeFromGrid(string fname, string run){

  TFileMerger *mgr = new TFileMerger();
  TGrid::Connect("alien://");
  cout<<"Processing file "<<fname.c_str()<<endl;
  ifstream infl(fname.c_str());
  string s;
  int count = 0;
  while(infl>>s){
    mgr->AddFile(s.c_str());
    count++;
  }
  mgr->OutputFile(Form("../Data/%s_merged.root", run.c_str()));
  cout<<endl;
  cout<<"... Merging "<<count<<" files"<<endl;
  cout<<endl;
  mgr->Merge();
}

void DownloadAndMergeFromGrid_step(string fname, string run, int nsteps) {
  TGrid::Connect("alien://");
  ifstream in(fname.c_str());
  string s;

  for(int istep=1; istep<=nsteps; istep++){
    TFileMerger* mg= new TFileMerger();
    int c=0;
    while(in>>s){
      mg->AddFile(s.c_str());
      c++;
      if(c == 60){
        cout<<endl;
        Printf("......Step: %d - Merging %d files", istep, c);
        cout<<endl;
        mg->OutputFile(Form("../Data/%s_merged_%d.root", run.c_str(), istep));
        mg->Merge();
        break;
      }
    }
    delete mg;
  }
}
