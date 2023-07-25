#include <unordered_map>
#include <string>
#include <TFile.h>
#include <vector>
#include <fstream>
#include <sstream>
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "DataFormatsITSMFT/NoiseMap.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/Logger.h"
// chips fully disabled from confDB (dead in all runs --> to be added on top of what found from QC)
std::vector<std::string> deadchips_db = {"L3_03L_M3_C03","L3_19U_M4_C04","L4_05U_M1_C03","L4_09L_M2_C02","L4_11U_M2_C06","L4_13U_M4_C11","L4_27U_M2_C06","L5_02L_M5_C03","L5_04L_M4_C06","L5_04U_M1_C13","L5_07U_M5_C04","L5_11L_M6_C06","L5_14L_M3_C04","L5_15L_M6_C11","L5_16L_M3_C06","L5_18U_M7_C06","L5_22L_M5_C11","L5_25L_M7_C09","L5_28U_M1_C09","L5_28U_M1_C13","L5_32U_M1_C05","L5_35U_M1_C10","L5_35U_M5_C05","L5_36L_M2_C04","L5_41L_M1_C14","L6_04L_M5_C13","L6_04L_M6_C13","L6_05U_M2_C05","L6_05U_M7_C06","L6_08L_M7_C11","L6_22L_M4_C02","L6_29L_M4_C12","L6_35L_M5_C02","L6_35L_M5_C03","L6_37U_M4_C12","L6_38U_M3_C14","L6_39L_M6_C13","L6_39U_M4_C05"};

void makeBadMap(int run, std::vector<std::string> chips, bool storetoccdb)
{
  std::unordered_map<std::string, int> chmap;
  o2::itsmft::ChipMappingITS mp;
  o2::itsmft::NoiseMap noiseMap(mp.getNChips());
  std::string fnm = fmt::format("o2-itsmft-NoiseMap_r{}",run);

  auto disableChips = [&](std::string partName)
  {
    //fnm += fmt::format("-{}", partName);
    for (size_t i=0;i<mp.getNChips();i++) {
      auto chname = mp.getChipNameHW(i);
      if (chname.find(partName)==0) {
	      printf("Masking %d %s\n", int(i), chname.c_str());
	      noiseMap.maskFullChip(i);
      }
    }
  };
  auto storeMap = [&](std::string name)
  {
    TFile fl(name.c_str(),"recreate");
    fl.WriteObjectAny(&noiseMap, o2::itsmft::NoiseMap::Class(), "ccdb_object");
  };

  for (auto& ch : chips) {
    disableChips(ch);
  }
  // add fully disabled chips on top
  for(auto& ch : deadchips_db){
    disableChips(ch);
  }
  auto& cm = o2::ccdb::BasicCCDBManager::instance();
  cm.setURL("http://alice-ccdb.cern.ch");
  auto lims = cm.getRunDuration(run);
  if (lims.first == 0) {
    LOG(error) << "failed to extract run limits for run " << run;
  }
  lims.first -= o2::ccdb::CcdbObjectInfo::MINUTE;
  lims.second += 10*o2::ccdb::CcdbObjectInfo::MINUTE;
  fnm += fmt::format("-{}-{}.root", lims.first, lims.second);
  if(storetoccdb){
    std::ofstream scriptfile("run_upload.sh", std::ios_base::app); // open in append mode
    std::string command = fmt::format("o2-ccdb-upload --host http://alice-ccdb.cern.ch -p ITS/Calib/DeadMap -k ccdb_object -f {} --starttimestamp {} --endtimestamp {} -m \"JIRA={};runNumber={}\"", fnm, lims.first, lims.second, "O2-3676", run);
    scriptfile<<command<<"\n";
    scriptfile.close();
    LOGP(info, "upload as: o2-ccdb-upload --host http://alice-ccdb.cern.ch -p ITS/Calib/DeadMap -k ccdb_object -f {} --starttimestamp {} --endtimestamp {} -m \"JIRA={};runNumber={}\"", fnm, lims.first, lims.second, "O2-3676", run);
    storeMap(fnm);
  }
}


void DoDeadChipMaps(std::string inputfilename = "file.txt", bool storetoccdb = false)
{
  std::unordered_map<int,std::vector<string>> mymap;
  std::ifstream infl(inputfilename.c_str());
  int run, mod;
  std::string stave, chipInMod, hs, line;
  for(int line_no = 1; std::getline(infl, line); ++line_no) {
    std::stringstream ss(line);
    ss>>run;
    ss>>stave;
    ss>>hs;
    ss>>mod;
    ss>>chipInMod;
    if(false)
      cout<<run<<" "<<stave<<" "<<hs<<" "<<mod<<" "<<chipInMod<<endl;
    int layer = std::stoi(stave.substr(stave.find("L")+1,1));
    mymap[run].push_back(layer<3 ? stave+"_C"+chipInMod : stave+hs+"_M"+std::to_string(mod)+"_C"+chipInMod);
  }

  //loop on the map
  for(auto const& x : mymap){
    makeBadMap(x.first, x.second, storetoccdb);
  }
}
