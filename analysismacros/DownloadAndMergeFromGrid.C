#include <TFileMerger.h>
#include <TGrid.h>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

void DownloadAndMergeFromGrid(){

  TFileMerger *mgr = new TFileMerger();
  TGrid::Connect("alien://");

  ifstream infl("filemerge.dat");
  string s;
  int count = 0;
  while(infl>>s){
    mgr->AddFile(s.c_str());
    count++;
  }
  mgr->OutputFile("../Data/509243_merged.root");
  cout<<endl;
  cout<<"... Merging "<<count<<" files"<<endl;
  cout<<endl;
  mgr->Merge();
}
