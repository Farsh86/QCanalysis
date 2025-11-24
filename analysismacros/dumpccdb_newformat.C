#include "iostream"
#include "fstream"
#include "vector"
#include "sstream"
#include "TFile.h"
#include "string"
#include "ITSMFTReconstruction/ChipMappingITS.h"

using namespace std;

struct CalibrationResStruct {
    int Layer;
    int Stave;
    int Hs;
    int HIC;
    int O2ChipID;
    int ChipID;
    int VCASN;
    int Unsuccess;
  };

vector<std::string> splitString(string s, string delimiter);
CalibrationResStruct CalibrationParser(string input);

void dumpccdb() {
  TFile *infl = new TFile("../Data/calib_scan_VCASN_561406.root");
  vector<char> *v;
  infl->GetObject("ccdb_object",v);
  string fullstring="";
  for(std::vector<char>::iterator it = v->begin(); it != v->end(); ++it) fullstring+=*it;

  //split and dump in txt file
  ofstream outfile ("../Data/vcasn_P2_test_561406.dat");
  outfile<<"Layer "<<" Stave"<<" Hs"<<" HIC"<<" Chip"<<" VCASN"<<endl; //intro
  cout<<"splitting string"<<endl;
  auto splitRes = splitString(fullstring, "O2ChipID:");
  cout<<"split string done"<<endl;
  int counter = 0;
  for (auto StaveStr : splitRes) {
    //cout<<"Srave str: "<<StaveStr<<endl;
    if(StaveStr.empty()) continue;
    CalibrationResStruct result = CalibrationParser(StaveStr);
    if(result.VCASN < 30) {
      cout<<result.Layer<<" "<<result.Stave<<" "<<result.Hs<<" "<<result.HIC<<" "<<result.ChipID<<" --> VCASN: "<<result.VCASN<<", Unsuccess: "<<result.Unsuccess<<endl;
    }
    string hs = result.Hs == 0 ? "L":"U";
    outfile<<result.Layer<<" "<<result.Stave<<" "<<hs<<" "<<result.HIC<<" "<<result.ChipID<<" "<<result.VCASN<<endl; //print last dcol after all pipes
    counter++;
  }
  cout<<"Processed chips: "<<counter<<endl;
  outfile.close();
}

vector<std::string> splitString(string s, string delimiter)
{
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;

  while ((pos_end = s.find(delimiter, pos_start)) != string::npos){
    token = ((delimiter.find("O2ChipID") != string::npos && (pos_start!=0 || pos_end!=0))) ? delimiter : "";
    token += s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
    token = "";
  }

  res.push_back(delimiter.find("O2ChipID") != string::npos ? delimiter+=s.substr(pos_start) : s.substr(pos_start));
  return res;
}

CalibrationResStruct CalibrationParser(string input)
{
  CalibrationResStruct result;
  auto StaveINFO = splitString(input, ",");

  for (string info : StaveINFO) {
    if (info.size() == 0)
      continue;

    std::string name = splitString(info, ":")[0];
    string data = splitString(info, ":")[1];
    if (name == "O2ChipID") {
      result.O2ChipID = std::stoi(data);
      //convert it to stave/hs/hic
      o2::itsmft::ChipMappingITS mp;
      int lay,sta,ssta,mod,chipInMod;
      mp.expandChipInfoHW(result.O2ChipID, lay,sta, ssta, mod,  chipInMod);
      result.Hs = ssta;
      result.Layer = lay;
      result.Stave = sta;
      result.HIC = mod;
      result.ChipID = chipInMod;
    }
    if (name == "VCASN") {
      result.VCASN = std::stoi(data);
    }
    if (name == "Status") {
      result.Unsuccess = std::stoi(data);
    }
  }
  return result;
}
