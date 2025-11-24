#include <iostream>
#include <fstream> 
#include <string>
#include <stdio.h>
#include <map>
#include <vector>
#include "TString.h"

using namespace std;

/*Split download list creating files for each epn*/

void SplitDownloadList(string filepath){
  
  ifstream infile(filepath.c_str());
  ofstream outfile;
  string fl;
  map<string,vector<string>> filelistepn;

  int counter = 0, filecounter=0; 
  while(infile>>fl){
    string epn = fl.substr(fl.find("epn"),6);
    filelistepn[epn].push_back(fl);
  }
  
  for( const auto& pair : filelistepn ) {
    outfile.open(Form("filemerge_%d.dat",filecounter));
    for(string s : pair.second){
      outfile<<s<<"\n";
    }
    filecounter++;
    outfile.close();
  }   
}
