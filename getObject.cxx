#include <string>
#include <iostream>
#include <time.h>
#include <algorithm>
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/CcdbDatabase.h"
#include "QualityControl/MonitorObject.h"
#include "CCDB/CcdbApi.h"
#include <TBufferJSON.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

using namespace std;
using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;

//functions to download data
string GetCorrectTS(string selrun, vector<string> runs, vector<string> timestamps);
<<<<<<< Updated upstream
array<string,2> GetLastRunWithTS(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname);
array<string,2> GetRunWithTS24hAgo(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string timestamp);
vector<string> GetGoodRunList(o2::ccdb::CcdbApi ccdbApi, string run1, string run2, string runtype, string qcpathstart);
bool RunShifter(CcdbDatabase *ccdb, int opt, string syncasync);
bool RunExpert(CcdbDatabase *ccdb, int opt, string syncasync);
void DownloadTimestamps(CcdbDatabase *ccdb, o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, long int ts_start, long int ts_end, int lnum);
void DownloadRuns(CcdbDatabase *ccdb, o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, int lnum, vector<string> runlistfromfile);
bool GetListOfHisto(CcdbDatabase* ccdb, string taskname, string objname, vector<long int> timestamps, vector<long int> timestamps2, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers, vector<int>runnumbers2);
bool Download(int choice, CcdbDatabase* ccdb, o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, long int ts_start, long int ts_end, int lnum, vector<string> runlistfromfile);
=======
array<string,2> GetLastRunWithTS(o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname);
array<string,2> GetRunWithTS24hAgo(o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, string timestamp);
vector<string> GetGoodRunList(o2::ccdb::CcdbApi &ccdbApi, string run1, string run2, string runtype, string qcpathstart);
bool RunShifter(CcdbDatabase *ccdb, int opt, string syncasync);
bool RunExpert(CcdbDatabase *ccdb, int opt, string syncasync);
void DownloadTimestamps(CcdbDatabase *ccdb, o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, long int ts_start, long int ts_end, int lnum);
void DownloadRuns(CcdbDatabase *ccdb, o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, int lnum, vector<string> runlistfromfile);
bool GetListOfHisto(o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, vector<long int> timestamps, vector<long int> timestamps2, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers, vector<int>runnumbers2);
bool Download(int choice, CcdbDatabase* ccdb, o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, long int ts_start, long int ts_end, int lnum, vector<string> runlistfromfile);
>>>>>>> Stashed changes
string GetOptName(int opt);
string GetListName(int opt, int ilist);

const int nStavesInLay[7] = {12, 16, 20, 24, 30, 42, 48};
TFile *outputfile;
std::string selpassname;

//to which CCDB we have to connect
// For P2 operations put: ali-qcdb.cern.ch:8083
string ccdbport = "ali-qcdb-gpn.cern.ch:8083";


void getObject(string expOrShift, int whatToDownload)
{

  int dbnum;
  std::cout<<"Select the CCDB from where to download data"<<std::endl;
  std::cout<<"  1. localhost:8083 (make sure to setup lxtunnel properly)"<<std::endl;
  std::cout<<"  2. ali-qcdb-gpn.cern.ch:8083 (default)" << std::endl;
  std::cout<<"  3. ccdb-test.cern.ch:8080" << std::endl;
  std::cout<<"Enter number: ";
  std::cin>>dbnum;
  std::cout<<std::endl;
  switch(dbnum){
    case 1:
      {
        ccdbport = "localhost:8083";
        break;
      }
    case 2:
      {
        ccdbport = "ali-qcdb-gpn.cern.ch:8083";
        break;
      }
    case 3:
      {
        ccdbport = "ccdb-test.cern.ch:8080";
        break;
      }
    default:
      {
        std::cout<<"Unknown number. Exiting..."<<std::endl;
        exit(1);
      }
  }

  int sa;
  string syncasync;
  std::cout<<std::endl;
  std::cout<<"You would like to check synchronous or asynchronous QC plots? "<<std::endl;
  std::cout<<"  1. Sync"<<std::endl;
  std::cout<<"  2. Async"<<std::endl;
  std::cout<<"Enter number: ";
  cin>>sa;
  std::cout<<std::endl;
  switch(sa){
    case 1:
    {
      syncasync = "sync";
      break;
    }
    case 2:
    {
      syncasync = "async";
      break;
    }
    default:
    {
      std::cout<<"Uknown number. Exiting ..."<<std::endl;
      exit(1);
    }
  }

  if(sa==2) {
    std::cout<<"Type the pass name to download as reported in PassName field in qcdb (ex: apass3): ";
    cin>>selpassname;
  }

  std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");
  if(mydb == nullptr) {
    cerr << "Failed to create DatabaseInterface" << endl;
    return;
  }

  auto* ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());
  if(ccdb == nullptr) {
    cerr << "ccdb pointer is null" << endl;
    return;
  }

  ccdb->connect(ccdbport.c_str(), "", "", "");

  if(expOrShift == "expert") {
    RunExpert(ccdb, whatToDownload, syncasync);
  }
  else if (expOrShift == "shifter") {
    RunShifter(ccdb, whatToDownload, syncasync);
  }
  else {
    cerr<<"The expOfShift string "<<expOrShift<<" is not recognized, please check!"<<endl;
  }

  ccdb->disconnect();

  return;
}//end main

//
//Get Last (most recent) Run
//
<<<<<<< Updated upstream
array<string,2> GetLastRunWithTS(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname){
=======
array<string,2> GetLastRunWithTS(o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname){
>>>>>>> Stashed changes

  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  stringstream ss(objectlist);
  string word;
  array<string,2> runts;

  while(ss>>word){
<<<<<<< Updated upstream
    if(word=="Created:"){// take the one related to file creation
=======
    if(word=="Validity:"){// take the one related to file creation
>>>>>>> Stashed changes
      ss>>word;
      runts[0] = word;
    }
    if(word=="RunNumber"){
      ss>>word;
      ss>>word;
      runts[1] = word;
      break;
    }
  }

  return runts;
}

//
// Get Run 24h depending on the actual time stamp
//
<<<<<<< Updated upstream
array<string,2> GetRunWithTS24hAgo(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string timestamp){
=======
array<string,2> GetRunWithTS24hAgo(o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, string timestamp){
>>>>>>> Stashed changes

  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  stringstream ss(objectlist);
  string word;
  array<string,2> runts;
  bool islast = false;

  long int stamp_int_actual = stol(timestamp);
  long int stamp_int_24hago = stamp_int_actual - 86400000; //remove number of milliseconds in 1 day

  while(ss>>word){
<<<<<<< Updated upstream
    if(word=="Created:"){// take the one related to file creation
=======
    if(word=="Validity:"){// take the one related to file creation
>>>>>>> Stashed changes
      ss>>word;
      if(stol(word)>=stamp_int_24hago){
        runts[0] = word;
      }
      else islast=true;
    }
    if(word=="RunNumber"){
      ss>>word;
      ss>>word;
      if(islast) break;
      if(stol(runts[0])>=stamp_int_24hago) runts[1] = word;
    }
  }

  return runts;
}

//
<<<<<<< Updated upstream
// Expert mode
=======
// Shift mode
>>>>>>> Stashed changes
//
bool RunShifter(CcdbDatabase *ccdb, int opt, std::string syncasync){

  std::string qcpathstart;
  if(syncasync == "sync") {
    qcpathstart = "qc/";
  } else if (syncasync == "async") {
    qcpathstart = "qc_async/";
  }

  //Choose the layer number
  int layernum;
  cout<<endl;
  cout<<endl;
  cout<<"Enter the layer number [put -1 for all IB layers]"<<endl;
  cin>>layernum;

  //add error and trigger-flags plot to data
  bool adderrordata = true;
  if(opt==2) adderrordata = false;

  //taskname
  string taskname[4] = {qcpathstart+"ITS/MO/ITSFHR", qcpathstart+"ITS/MO/ITSFHR", qcpathstart+"ITS/MO/ITSFHR", qcpathstart+"ITS/MO/ITSFHR"};

  //set the task name
  switch(opt){
    case 1: {// fake-hit
      taskname[0] = qcpathstart+"ITS/MO/ITSFHR";//L0T, L0B
      taskname[1] = qcpathstart+"ITS/MO/ITSFHR";//L1T, L1B
      taskname[2] = qcpathstart+"ITS/MO/ITSFHR"; //L2T
      taskname[3] = qcpathstart+"ITS/MO/ITSFHR"; //L2B
      break;
    }

    case 2: { //thr scan
<<<<<<< Updated upstream
      taskname[0] = qcpathstart+"ITS/MO/ITSTHRTask0";
      taskname[1] = qcpathstart+"ITS/MO/ITSTHRTask1";
      taskname[2] = qcpathstart+"ITS/MO/ITSTHRTask2T";
      taskname[3] = qcpathstart+"ITS/MO/ITSTHRTask2B";
=======
      taskname[0] = qcpathstart+"ITS/MO/ITSThresholdCalibrationTask";
      taskname[1] = qcpathstart+"ITS/MO/ITSThresholdCalibrationTask";
      taskname[2] = qcpathstart+"ITS/MO/ITSThresholdCalibrationTask";
      taskname[3] = qcpathstart+"ITS/MO/ITSThresholdCalibrationTask";
>>>>>>> Stashed changes
      break;
    }

    default: break;
  }//end of switch



  //CCDB api initialization
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init(ccdbport.c_str());

  //set variables (run interval of timestamp interval)
  array<string,2> runts1;
  array<string,2> runts2;

  //Decide how many different elements are needed
  int nListElements = 2;
  switch(opt){
    case 1: nListElements = 2; break;
    case 2: nListElements = 3; break;
    default: nListElements = 2;
  }
  if(adderrordata && opt==1)
    nListElements+=2;//2 because we add trigger and error plots

  //Output file
  string layername;
  if(layernum==-1)
    layername = "all-IB-layers";
  else
    layername = Form("Layer%d",layernum);

  string suffix = "run";
  string optname = GetOptName(opt);

  //Download depending on the option (opt)
  switch(opt){
    case 1: {

      //run interval definition
      cout<<"Finding runs in: "<<taskname[0]<<"/Occupancy/Layer0/Layer0ChipStave"<<endl;
      runts2 = GetLastRunWithTS(ccdbApi, taskname[0], "Occupancy/Layer0/Layer0ChipStave"); //take a random object name since run-list is the same.
      runts1 = GetRunWithTS24hAgo(ccdbApi, taskname[0], "Occupancy/Layer0/Layer0ChipStave", runts2[0]);
      cout<<"Run interval selected:       "<<runts1[1]<<"-"<<runts2[1]<<endl;
      cout<<"Timestamp interval selected: "<<runts1[0]<<"-"<<runts2[0]<<endl;

      //output file
      outputfile = new TFile(Form("Data/Output_%s_%s_from_%s%s_to_%s%s%s.root", layername.c_str(), optname.c_str(), suffix.c_str(), runts1[1].c_str(), suffix.c_str(), runts2[1].c_str(), adderrordata? "_w_error_and_trig_data":""), "RECREATE");
      outputfile->cd();

      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Occupancy/Layer%d/Layer%dChipStave",layernum,layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[layernum],  objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum, vector<string>());
              break;
            }

            case 1: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9){
                    Download(1, ccdb, ccdbApi, taskname[layernum+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum, vector<string>());
                  }
                  else{
                    Download(1, ccdb, ccdbApi, taskname[layernum], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum, vector<string>());
                  }
                }
                else {
                  Download(1, ccdb, ccdbApi, taskname[layernum], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum, vector<string>());
                }

              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum, vector<string>());
              break;
            }

            case 3: {//error files
              string objname = "General/TriggerVsFeeid";
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum, vector<string>());
              break;
            }
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){

        vector<string> goodrunlist = GetGoodRunList(ccdbApi, runts1[1], runts2[1], "Fhr", qcpathstart);

        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
                Download(1, ccdb, ccdbApi, taskname[ilay], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay,vector<string>());
              }
              break;
            }

            case 1: {
              for(int ilay=0; ilay<=2; ilay++){
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9){
                      Download(1, ccdb, ccdbApi, taskname[ilay+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay,vector<string>());
                    }
                    else{
                      Download(1, ccdb, ccdbApi, taskname[ilay], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay,vector<string>());
                    }
                  }
                  else{
                    Download(1, ccdb, ccdbApi, taskname[ilay], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay, vector<string>());
                  }
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[0], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), 0, vector<string>());
              break;
            }

            case 3: {//trigger files
              string objname = "General/TriggerVsFeeid";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[0], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), 0, vector<string>());
              break;
            }
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }//end case 1
<<<<<<< Updated upstream

    case 2: {//thresholds
      //run interval definition
      cout<<"Finding runs in: "<<taskname[0]<<"/Threshold/Layer0/Threshold_Vs_Chip_and_Stave"<<endl;
      runts2 = GetLastRunWithTS(ccdbApi, taskname[0], "Threshold/Layer0/Threshold_Vs_Chip_and_Stave"); //take a random object name since run-list is the same.
      runts1 = GetRunWithTS24hAgo(ccdbApi, taskname[0], "Threshold/Layer0/Threshold_Vs_Chip_and_Stave", runts2[0]);
      vector<string> goodrunlist = GetGoodRunList(ccdbApi, runts1[1], runts2[1], "Thr", qcpathstart);
      cout<<"Run interval selected:       "<<runts1[1]<<"-"<<runts2[1]<<endl;
      cout<<"Timestamp interval selected: "<<runts1[0]<<"-"<<runts2[0]<<endl;

      outputfile = new TFile(Form("Data/Output%s_%s_from_%s%s_to_%s%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),runts1[1].c_str(), suffix.c_str(), runts2[1].c_str(), adderrordata? "_w_error_and_trig_data":""), "RECREATE");
      outputfile->cd();

      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum,vector<string>());
              break;
            }

            case 1: {
              string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum, vector<string>());
              break;
            }

            case 2: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9)
                    Download(1, ccdb, ccdbApi, taskname[layernum+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum, vector<string>());
                  else
                    Download(1, ccdb, ccdbApi, taskname[layernum], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum, vector<string>());
                }
                else{
                  Download(1, ccdb, ccdbApi, taskname[layernum],  Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum, vector<string>());
                }
              }
              break;
            }

            default: break;
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
                Download(1, ccdb, ccdbApi, taskname[ilay], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay, vector<string>());
              }
              break;
            }

            case 1: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
                Download(1, ccdb, ccdbApi, taskname[ilay], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay, vector<string>());
              }
              break;
            }

            case 2: {
              for(int ilay=0; ilay<=2; ilay++){
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9)
                      Download(1, ccdb, ccdbApi, taskname[ilay+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay, vector<string>());
                    else
                      Download(1, ccdb, ccdbApi, taskname[ilay], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay, vector<string>());
                  }
                  else{
                    Download(1, ccdb, ccdbApi, taskname[ilay], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay, vector<string>());
                  }
                }
              }
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }// end of case 2

  }//end switch

=======
 
    case 2: {//thresholds
        cout << "THR task does not support shifter mode ";
        /*
        //run interval definition
        cout << "Finding runs in: " << taskname[0] << "ThrNoiseChipAverageIB" << endl;
        runts2 = GetLastRunWithTS(ccdbApi, taskname[0], "ThrNoiseChipAverageIB"); //take a random object name since run-list is the same.
        runts1 = GetRunWithTS24hAgo(ccdbApi, taskname[0], "ThrNoiseChipAverageIB", runts2[0]);
        vector<string> goodrunlist = GetGoodRunList(ccdbApi, runts1[1], runts2[1], "Thr", qcpathstart);
        cout << "Run interval selected:       " << runts1[1] << "-" << runts2[1] << endl;
        cout << "Timestamp interval selected: " << runts1[0] << "-" << runts2[0] << endl;

        outputfile = new TFile(Form("Data/Output%s_%s_from_%s%s_to_%s%s.root", layername.c_str(), optname.c_str(), suffix.c_str(), runts1[1].c_str(), suffix.c_str(), runts2[1].c_str()), "RECREATE");
        outputfile->cd();

        if (layernum >= 0) {
                    string objname = "ThrNoiseChipAverageIB";
                    cout << "\nAll data in " << taskname[layernum] + "/" + objname << " between run" << runts1[1] << " and run" << runts2[1] << " are going to be downloaded." << endl;
                    Download(1, ccdb, ccdbApi, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum, vector<string>());
         
        }//end if layernum>=0

        else if (layernum == -1) {
                    for (int ilay = 0; ilay <= 2; ilay++) {
                        string objname = "ThrNoiseChipAverageIB";
                        cout << "\nAll data in " << taskname[ilay] + "/" + objname << " between run" << runts1[1] << " and run" << runts2[1] << " are going to be downloaded." << endl;
                        Download(1, ccdb, ccdbApi, taskname[ilay], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay, vector<string>());
                    }
              
        }//end if layernum==-1
        break;
        */
    }// end of case 2
    
  }//end switch
 
>>>>>>> Stashed changes
  outputfile->Close();
  delete outputfile;

  return 1;
}

//
// Expert mode
//
bool RunExpert(CcdbDatabase *ccdb, int opt, std::string syncasync){
  std::string qcpathstart;
  if(syncasync == "sync") {
    qcpathstart = "qc/";
  } else if (syncasync == "async") {
    qcpathstart = "qc_async/";
  }

  int IBorOB;
  int layernum;
  if(opt==3) layernum=-2;
  else{

  //choose IB, OB or both
  cout << endl;
  cout << endl;
  cout << "Choose beteen IB (0), OB (1) or both (2, = all layers of IB and OB)" << endl;
  cin>> IBorOB;

  //Choose the layer number
  cout<<endl;

  switch(IBorOB){
  case 0: {
    cout<<"Enter the layer number [put -1 for all IB layers]"<<endl;
    cin >>layernum;
    break;
  }
  case 1: {
    cout<<"Enter the layer number [put -1 for all OB layers]"<<endl;
    cin >>layernum;
    break;
  }
  case 2: {
    layernum=-1;
    break;
  }
  default: break;
  }

  }

<<<<<<< Updated upstream
  //  cin>>layernum;

=======
>>>>>>> Stashed changes
  //taskname
  string taskname[8]    = {qcpathstart+"ITS/MO/FHRTask", qcpathstart+"ITS/MO/FHRTask", qcpathstart+"ITS/MO/FHRTask", qcpathstart+"ITS/MO/FHRTask",  qcpathstart+"ITS/MO/FHRTask",  qcpathstart+"ITS/MO/FHRTask",  qcpathstart+"ITS/MO/FHRTask",  qcpathstart+"ITS/MO/FHRTask"};

  //set the task name
  switch(IBorOB){
  case 0: {//Inner Barrel
    switch(opt){
    case 1: {// fake-hit
      taskname[0] = qcpathstart+"ITS/MO/FHRTask";//L0T, L0B
      taskname[1] = qcpathstart+"ITS/MO/FHRTask";//L1T, L1B
      taskname[2] = qcpathstart+"ITS/MO/FHRTask"; //L2T
      taskname[3] = qcpathstart+"ITS/MO/FHRTask"; //L2B
      break;
    }

    case 2: { //thr scan
      taskname[0] = qcpathstart+"ITS/MO/ITSTHRTask0";
      taskname[1] = qcpathstart+"ITS/MO/ITSTHRTask1";
      taskname[2] = qcpathstart+"ITS/MO/ITSTHRTask2T";
      taskname[3] = qcpathstart+"ITS/MO/ITSTHRTask2B";

      break;
    }
    default: break;
    }//end of switch
    break;
  }
  case 1: { //Outer Barrel
    switch(opt){
    case 1: {// fake-hit
      taskname[4] = qcpathstart+"ITS/MO/FHRTask";//L3
      taskname[5] = qcpathstart+"ITS/MO/FHRTask";//L4
      taskname[6] = qcpathstart+"ITS/MO/FHRTask";//L5
      taskname[7] = qcpathstart+"ITS/MO/FHRTask";//L6

      break;
    }

    case 2: { //thr scan
      //path name to be changed for OB
      taskname[4] = qcpathstart+"ITS/MO/ITSTHRTask0";
      taskname[5] = qcpathstart+"ITS/MO/ITSTHRTask1";
      taskname[6] = qcpathstart+"ITS/MO/ITSTHRTask2T";
      taskname[7] = qcpathstart+"ITS/MO/ITSTHRTask2B";

      break;
    }
    default: break;
    }//end of switch
  }
  case 2: { //Inner and Outer Barrel
    switch(opt){
    case 1: {// fake-hit
      taskname[0] = qcpathstart+"ITS/MO/FHRTask";//L0
      taskname[1] = qcpathstart+"ITS/MO/FHRTask";//L1
      taskname[2] = qcpathstart+"ITS/MO/FHRTask";//L2
      taskname[3] = qcpathstart+"ITS/MO/FHRTask";//not used for FHR
      taskname[4] = qcpathstart+"ITS/MO/FHRTask";//L3
      taskname[5] = qcpathstart+"ITS/MO/FHRTask";//L4
      taskname[6] = qcpathstart+"ITS/MO/FHRTask";//L5
      taskname[7] = qcpathstart+"ITS/MO/FHRTask";//L6

      break;
    }

    case 2: { //thr scan
      // path name to be changed
<<<<<<< Updated upstream
      taskname[0] = qcpathstart+"ITS/MO/ITSTHRTask0";
      taskname[1] = qcpathstart+"ITS/MO/ITSTHRTask1";
      taskname[2] = qcpathstart+"ITS/MO/ITSTHRTask2T";
      taskname[3] = qcpathstart+"ITS/MO/ITSTHRTask2B";
      taskname[4] = qcpathstart+"ITS/MO/ITSTHRTask0";
      taskname[5] = qcpathstart+"ITS/MO/ITSTHRTask1";
      taskname[6] = qcpathstart+"ITS/MO/ITSTHRTask2T";
      taskname[7] = qcpathstart+"ITS/MO/ITSTHRTask2B";
=======
      taskname[0] = qcpathstart+"ITS/MO/ITSThresholdCalibrationTask";
>>>>>>> Stashed changes

      break;
    }
    default: break;
    }//end of switch
  }
  default: break;
  }//end of switch

  //Ask whether attach a error report to the results
  bool adderrordata = false;
  string erranswer = "n";
  if(opt==1){//only is user selects to download FHR
    cout<<endl;
    cout<<"Error-flag plots needed? [y/n]"<<endl;
    cin>>erranswer;
    if(erranswer=="y")
      adderrordata = true;
  }

  //Chose about run number or time interval
  int choice;
  cout<<endl;
  cout<<"Download data - Choose an option:\n"<<endl;
  cout<<"1. Enter run numbers (not possible in qc_async)"<<endl;
  cout<<"2. Enter timestamps (not possible in qc_async)"<<endl;
  cout<<"3. Run list as input txt file (must be in Data folder)"<<endl;
  cout<<"[Type option number]: ";
  cin>>choice;

  //set variables (run interval of timestamp interval)
  string run1="0", run2="0";
  time_t ts_start=0, ts_end=0;
  vector<string> nums; //can contain selected run interval or timestamp interval
  vector<string> runlistfromfile;
  switch(choice){
    case 1: {
      cout<<"\n"<<"Run interval [just run number, NO ZEROs in front]"<<endl;
      cout<<"Enter first  run: ";
      cin>> run1;
      cout<<"Enter second run: ";
      cin>> run2;
      nums.push_back(run1);
      nums.push_back(run2);
      break;
    }

    case 2:{
      //Enter time period of the data to be downloaded
      struct tm datetime_start;
      struct tm datetime_end;
      string day, month, hour, min, sec;
      int year;
      cout<<"\n"<<"Time period [format (hour 24H): DD MM YYYY HH MM SS]"<<endl;
      cout<<"Enter start date and time: ";
      cin>> day >> month >> year >> hour >> min >> sec;
      string datetime1 = Form("%d%s%s_%s%s%s",year,month.c_str(),day.c_str(),hour.c_str(),min.c_str(),sec.c_str());
      datetime_start.tm_mday = stoi(day);
      datetime_start.tm_mon = stoi(month)-1;//-1 because january = 0;
      datetime_start.tm_year = year-1900;
      datetime_start.tm_hour = stoi(hour);
      datetime_start.tm_min = stoi(min);
      datetime_start.tm_sec = stoi(sec);
      datetime_start.tm_isdst = -1; //-1 means DST unknown
      ts_start = mktime(&datetime_start);//get timestamp start
      ts_start*=1000; //convert in milliseconds

      cout<<"Enter end date and time  : ";
      cin>> day >> month >> year >> hour >> min >> sec;
      string datetime2 = Form("%d%s%s_%s%s%s",year,month.c_str(),day.c_str(),hour.c_str(),min.c_str(),sec.c_str());
      datetime_end.tm_mday = stoi(day);
      datetime_end.tm_mon = stoi(month)-1;//-1 because january = 0;
      datetime_end.tm_year = year-1900;
      datetime_end.tm_hour = stoi(hour);
      datetime_end.tm_min = stoi(min);
      datetime_end.tm_sec = stoi(sec);
      datetime_end.tm_isdst = -1; //-1 means DST unknown
      ts_end = mktime(&datetime_end);//get timestamp start
      ts_end*=1000; //convert in milliseconds
      printf("Timestamp start (ms): %ld\n", (long) ts_start);
      printf("Timestamp end   (ms): %ld\n", (long) ts_end);
      nums.push_back(datetime1);
      nums.push_back(datetime2);
      break;
    }

    case 3: {
      cout<<"\n"<<"Files potentially containing a run list:"<<endl;
      cout<<endl;
      gSystem->Exec("ls Data/*.txt -Art | tail -n 500");
      string filename;
      cout<<"\n"<<"Enter file name (example: my_run_list.txt): ";
      cin>>filename;
      ifstream inflruns(Form("Data/%s",filename.c_str()));
      string runi;
      while(inflruns>>runi){
        runlistfromfile.push_back(runi); // run list from user
      }
      inflruns.close();
      nums.push_back(runlistfromfile[runlistfromfile.size()-1]);
      nums.push_back(runlistfromfile[0]);
      run1 = runlistfromfile[runlistfromfile.size()-1];
      run2 = runlistfromfile[0];
      break;
    }

    default: break;
  }

  //Decide how many different elements are needed
  int nListElements = 2;
  switch(opt){
    case 1: nListElements = 2; break;
    case 2: nListElements = 3; break;
    default: nListElements = 2;
  }
  if(adderrordata)
    nListElements+=1;

  //CCDB api initialization
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init(ccdbport.c_str());

  //Output file
  string layername;
  if(layernum==-1){
    if(IBorOB==0)    layername = "_all-IB-layers";
    else if (IBorOB==1) layername = "_all-OB-layers";
    else layername = "_all-layers";
  }
  else if(layernum==-2) layername="";
  else
    layername = Form("_Layer%d",layernum);

  int layernumEff=layernum;
  if(IBorOB==1) layernumEff=layernum+1;

  int ilayMin = 0;
  int ilayMax = 2;
  if (IBorOB==1) {
    ilayMin = 3;
    ilayMax =6;
  }
  else if (IBorOB==2) {
    ilayMin = 0;
    ilayMax =6;
  }

  string suffix;
  switch(choice){
    case 1: suffix = "run"; break;
    case 2: suffix = ""; break;
    case 3: suffix = "run"; break;
    default: suffix="";
  }
  string optname = GetOptName(opt);
<<<<<<< Updated upstream
  outputfile = new TFile(Form("Data/Output%s_%s_from_%s%s_to_%s%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),nums[0].c_str(), suffix.c_str(), nums[1].c_str(), adderrordata? "_w_error_data":""), "RECREATE");
=======
  outputfile = new TFile(Form("Data/Output%s_%s_from_%s%s_to_%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),nums[0].c_str(), suffix.c_str(), nums[1].c_str()), "RECREATE");
>>>>>>> Stashed changes
  outputfile->cd();

  //Download depending on the option (opt)
  switch(opt){
    case 1: {
      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Occupancy/Layer%d/Layer%dChipStave",layernum,layernum);
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, taskname[layernumEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end,layernum, runlistfromfile);
              break;
            }

            case 1: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9){//L2B
                    Download(choice, ccdb, ccdbApi, taskname[layernumEff+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                  }
                  else{
                    Download(choice, ccdb, ccdbApi, taskname[layernumEff], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                  }
                }
                else{
                  Download(choice, ccdb, ccdbApi, taskname[layernumEff], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, taskname[layernumEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
              break;
            }

            default: break;
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
<<<<<<< Updated upstream
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
=======
		        int ilayEff=ilay;
		        if (IBorOB==1) ilayEff=ilay+1;
	        	else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
>>>>>>> Stashed changes
                string objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
                cout<<"\nAll data in "<<taskname[ilayEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, taskname[ilayEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay, runlistfromfile);
              }
              break;
            }

            case 1: {
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
<<<<<<< Updated upstream
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
=======
		        int ilayEff=ilay;
		        if (IBorOB==1) ilayEff=ilay+1;
		        else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
>>>>>>> Stashed changes
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9){//L2B
                      Download(choice, ccdb, ccdbApi, taskname[ilayEff+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay, runlistfromfile);
                    }
                    else{
                      Download(choice, ccdb, ccdbApi, taskname[ilayEff], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay, runlistfromfile);
                    }
                  }
                  else{
                    Download(choice, ccdb, ccdbApi, taskname[ilayEff],  Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay, runlistfromfile);
                  }
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, taskname[0], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }//end case 1
<<<<<<< Updated upstream

    case 2: {// thresholds
      vector<string> goodrunlist = GetGoodRunList(ccdbApi, run1, run2, "Thr", qcpathstart); //good run list
      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, taskname[layernumEff],  objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end,layernum,runlistfromfile);
              break;
            }

            case 1: {
              string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, taskname[layernumEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
              break;
            }

            case 2: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9)
                    Download(choice, ccdb, ccdbApi, taskname[layernumEff+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                  else
                    Download(choice, ccdb, ccdbApi, taskname[layernumEff], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                }
                else{
                  Download(choice, ccdb, ccdbApi, taskname[layernumEff], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                }

              }
              break;
            }

            default: break;
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
	      //for(int ilay=0; ilay<=2; ilay++){
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilayEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, taskname[ilayEff], objname, run1, run2, goodrunlist, (long)ts_start, (long)ts_end,ilay, runlistfromfile);
              }
              break;
            }

            case 1: {
	      // for(int ilay=0; ilay<=2; ilay++){
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilayEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, taskname[ilayEff],  objname, run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay, runlistfromfile);
              }
              break;
            }

            case 2: {
	      //              for(int ilay=0; ilay<=2; ilay++){
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9)
                      Download(choice, ccdb, ccdbApi, taskname[ilayEff+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay, runlistfromfile);
                    else
                      Download(choice, ccdb, ccdbApi, taskname[ilayEff], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay, runlistfromfile);
                  }
                  else{
                    Download(choice, ccdb, ccdbApi, taskname[ilayEff],  Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay, runlistfromfile);
                  }
                }
              }
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
=======
    case 2: {// thresholds
        string tskn = "qc/ITS/MO/ITSThresholdCalibrationTask";

                    string objname = "THRChipAverageIB";
                    cout << "\nAll data in " << tskn + "/" + objname << " between run" << run1 << " and run" << run2 << " are going to be downloaded." << endl;
                    Download(choice, ccdb, ccdbApi, tskn, objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                      
                    objname = "THRChipAverageML";
                    cout << "\nAll data in " << tskn + "/" + objname << " between run" << run1 << " and run" << run2 << " are going to be downloaded." << endl;
                    Download(choice, ccdb, ccdbApi, tskn, objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);
                        
                    objname = "THRChipAverageOL";
                    cout << "\nAll data in " << tskn + "/" + objname << " between run" << run1 << " and run" << run2 << " are going to be downloaded." << endl;
                    Download(choice, ccdb, ccdbApi, tskn, objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);

                    string objname = "ThrNoiseChipAverageIB";
                    cout << "\nAll data in " << tskn + "/" + objname << " between run" << run1 << " and run" << run2 << " are going to be downloaded." << endl;
                    Download(choice, ccdb, ccdbApi, tskn, objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);

                    objname = "ThrNoiseChipAverageML";
                    cout << "\nAll data in " << tskn + "/" + objname << " between run" << run1 << " and run" << run2 << " are going to be downloaded." << endl;
                    Download(choice, ccdb, ccdbApi, tskn, objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);

                    objname = "ThrNoiseChipAverageOL";
                    cout << "\nAll data in " << tskn + "/" + objname << " between run" << run1 << " and run" << run2 << " are going to be downloaded." << endl;
                    Download(choice, ccdb, ccdbApi, tskn, objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);

        break;
>>>>>>> Stashed changes
    }// end of case 2

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

<<<<<<< Updated upstream
case 3: { //tracks
=======
    case 3: { //tracks
>>>>>>> Stashed changes
      string trackname = (syncasync == "async") ? "Tracks" : "ITSTrackTask";
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "AngularDistribution", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "EtaDistribution", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "PhiDistribution", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "NClusters", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "VertexZ", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "VertexRvsZ", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "VertexCoordinates", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "NVertexContributors", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "Ntracks", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "AssociatedClusterFraction", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+trackname, "NClustersPerTrackEta", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);

      break;
    }//end of case 3
    //////////////////////////////////////////////////////////
<<<<<<< Updated upstream
case 4: { //FEE
=======
    case 4: { //FEE
>>>>>>> Stashed changes
      cout<<"\nAll data in "<<qcpathstart+"ITS/MO/ITSFEE/LaneStatus/LaneStatusFlagError(ERROR,FAULT,WARNING)"<<" between run"<<run1<<" and run"<<run2<<" downloading..."<<endl;
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE","LaneStatus/laneStatusFlagERROR", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE","LaneStatus/laneStatusFlagFAULT", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE","LaneStatus/laneStatusFlagWARNING", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE","LaneStatus/laneStatusFlagCumulativeERROR", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE","LaneStatus/laneStatusFlagCumulativeFAULT", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE","LaneStatus/laneStatusFlagCumulativeWARNING", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      string objname = "TriggerVsFeeid";
<<<<<<< Updated upstream
      cout<<"\nAll data in "<<qcpathstart<< "ITS/MO/ITSFEE/"<<objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
=======
      cout<<"\nAll data in "<<qcpathstart<< "ITS/MO/ITSFEE/"<<objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;  
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE", objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      objname = "RDHSummaryCumulative";
      cout<<"\nAll data in "<<qcpathstart<< "ITS/MO/ITSFEE/"<<objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;  
>>>>>>> Stashed changes
      Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/ITSFEE", objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0, runlistfromfile);
      break;
    }//end of case 4
    //////////////////////////////////////////////////////////
    case 5:{ //Clusters
      string cluname = (syncasync=="async") ? "Clusters":"ITSClusterTask";
      if(layernum>=0){

         //Cluster occupation
         cout<<"\nAll data in "<< Form("%sITS/MO/%s/Layer%d/ClusterOccupation",qcpathstart.c_str(), cluname.c_str() ,layernum) << " between run" << run1 << " and run" << run2 << " are going to be downloaded." <<endl;
         Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+cluname, Form("Layer%d/ClusterOccupation", layernum), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);

         //Average cluster size
         cout<<"\nAll data in "<< Form("%sITS/MO/%s/Layer%d/AverageClusterSize", qcpathstart.c_str(), cluname.c_str(), layernum) << " between run" << run1 << " and run" << run2 << " are going to be downloaded." <<endl;
         Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+cluname, Form("Layer%d/AverageClusterSize", layernum), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum, runlistfromfile);

      } else if (layernum == -1){

         for(int ilay = ilayMin; ilay <= ilayMax; ilay++){ // "-1" option corresponds to ALL LAYERS in OB or IB, or even for whole ITS

            //Cluster occupation
            cout<<"\nAll data in "<< Form("%sITS/MO/%s/Layer%d/ClusterOccupation", qcpathstart.c_str(), cluname.c_str(),ilay) << " between run" << run1 << " and run" << run2 << " are going to be downloaded." <<endl;
            Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+cluname, Form("Layer%d/ClusterOccupation", ilay), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay, runlistfromfile);

            //Average cluster size
            cout<<"\nAll data in "<< Form("%sITS/MO/%s/Layer%d/AverageClusterSize",qcpathstart.c_str(), cluname.c_str(), ilay) << " between run" << run1 << " and run" << run2 << " are going to be downloaded." <<endl;
            Download(choice, ccdb, ccdbApi, qcpathstart+"ITS/MO/"+cluname, Form("Layer%d/AverageClusterSize", ilay), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay, runlistfromfile);
         }
      }
      break;
   }//end of case 5
 }//end switch

  outputfile->Close();
  delete outputfile;

  return 1;
}

//
// Return run list of GOOD runs = runs in all CCDB paths for TH2
//
<<<<<<< Updated upstream
vector<string> GetGoodRunList(o2::ccdb::CcdbApi ccdbApi, string run1, string run2, string runtype, string qcpathstart){//run type can be "Thr" or "Fhr"

  string objnames[4] = {"Threshold/Layer0/Threshold_Vs_Chip_and_Stave", "Threshold/Layer1/Threshold_Vs_Chip_and_Stave", "Threshold/Layer2/Threshold_Vs_Chip_and_Stave", "Threshold/Layer2/Threshold_Vs_Chip_and_Stave"};
  string tasknames[4] = {qcpathstart+"ITS/MO/ITSTHRTask0", qcpathstart+"ITS/MO/ITSTHRTask1", qcpathstart+"ITS/MO/ITSTHRTask2T", qcpathstart+"ITS/MO/ITSTHRTask2B"};
=======
vector<string> GetGoodRunList(o2::ccdb::CcdbApi &ccdbApi, string run1, string run2, string runtype, string qcpathstart){//run type can be "Thr" or "Fhr"

  string objnames[3] = {"ThrNoiseChipAverageIB", "ThrNoiseChipAverageIB", "ThrNoiseChipAverageIB"};
  string tasknames[3] = {qcpathstart+"ITS/MO/ITSThresholdCalibrationTask", qcpathstart+"ITS/MO/ITSThresholdCalibrationTask", qcpathstart+"ITS/MO/ITSThresholdCalibrationTask"};
>>>>>>> Stashed changes

  //in case of FHR
  if(runtype=="Fhr"){
    objnames[0] = "Occupancy/Layer0/Layer0ChipStave";
    objnames[1] = "Occupancy/Layer1/Layer1ChipStave";
    objnames[2] = "Occupancy/Layer2/Layer2ChipStave";
<<<<<<< Updated upstream
    objnames[3] = "Occupancy/Layer2/Layer2ChipStave";
=======
>>>>>>> Stashed changes

    tasknames[0] = qcpathstart+"ITS/MO/ITSFHR";
    tasknames[1] = qcpathstart+"ITS/MO/ITSFHR";
    tasknames[2] = qcpathstart+"ITS/MO/ITSFHR";
<<<<<<< Updated upstream
    tasknames[3] = qcpathstart+"ITS/MO/ITSFHR";
  }

  vector<vector<string>> runlists;
  for(int i=0; i<4; i++){//get run list between run1 and run2 for all the paths defined above
=======
  }

  vector<vector<string>> runlists;
  for(int i=0; i<3; i++){//get run list between run1 and run2 for all the paths defined above
>>>>>>> Stashed changes
    string objectlist = ccdbApi.list(tasknames[i] + "/" + objnames[i],false,"text/plain");
    stringstream ss(objectlist);
    string word;
    vector<string>runlist;
    while(ss>>word){
<<<<<<< Updated upstream
      if(word=="Created:"){// take the one related to file creation
=======
      if(word=="Validity:"){// take the one related to file creation
>>>>>>> Stashed changes
        ss>>word;
        //alltimestampsALT.push_back(word);
      }
      if(word=="RunNumber"){
        ss>>word;
        ss>>word;
        if(stoi(word)>=stoi(run1) && stoi(word)<=stoi(run2)){
	        if(std::find(runlist.begin(), runlist.end(), word) == runlist.end())//if element doesn't exist already
		        runlist.push_back(word);
	      }
        if(stoi(word)==stoi(run1)) break;
      }
    }
    runlists.push_back(runlist);
  }

  //Extract runs present in all runlists
  vector<string> runsfiltered;
  for(int il=0; il<(int)runlists[0].size(); il++){
    string run = runlists[0][il];
    int cntrun = 1;
    if((int)runlists[2].size()==0) cntrun++;// in case half L2 is missing
    if((int)runlists[3].size()==0) cntrun++;// in case half L2 is missing
    //look for this run in the other lists
    for(int cntlist=1; cntlist<4; cntlist++){
      for(int il2=0; il2<(int)runlists[cntlist].size(); il2++){
        if(run==runlists[cntlist][il2]) {
          cntrun++;
        }
      }
    }
    if(cntrun==4){
      runsfiltered.push_back(run);
      cout<<"Run: "<<run<<endl;
    }
  }

  return runsfiltered;
}


//
// Download data based on timestamps
//
<<<<<<< Updated upstream
void DownloadTimestamps(CcdbDatabase* ccdb, o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, long int ts_start, long int ts_end, int lnum){
=======
void DownloadTimestamps(CcdbDatabase* ccdb, o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, long int ts_start, long int ts_end, int lnum){
>>>>>>> Stashed changes

  //Extract all the time stamps of the object
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  stringstream ss(objectlist);
  string word;
  vector <string> timestamps;
  while(ss>>word){
<<<<<<< Updated upstream
    if(word=="Created:"){// take the one related to file creation
=======
    if(word=="Validity:"){// take the one related to file creation
>>>>>>> Stashed changes
      ss>>word;
      timestamps.push_back(word);
    }
  }

  //filter all the timestamps to get the ones within start and end dates (get most recent!!)
  vector <long int> timestamps_selperiod;
  int counter=0;
  for(int its=0; its<(int)timestamps.size()-1; its++){
    if(stol(timestamps[its])>=ts_start && stol(timestamps[its])<=ts_end){
      counter++;
      if(counter==1) {
        if(its==0)
          timestamps_selperiod.push_back(stol(timestamps[its]));//the first (most recent) is taken by definition if it is really the first in the list
        else if(stol(timestamps[its-1])-stol(timestamps[its])>65000) timestamps_selperiod.push_back(stol(timestamps[its]));//if not the real first, check if it is the most recent (in this way the user can insert whatever date interval)
      }
      if(stol(timestamps[its])-stol(timestamps[its+1])>65000) timestamps_selperiod.push_back(stol(timestamps[its+1]));
    }

    else if(stol(timestamps[its])>ts_end) continue;
    else break;
  }
  cout<<"\n"<<"Selected timestamps: "<<endl;
  for(int i=0; i<(int)timestamps_selperiod.size();i++){
    cout<<timestamps_selperiod[i]<<endl;
  }

  bool isperstave = 0;
  if(objname.find("HITMAP")!=string::npos) isperstave = 1;
<<<<<<< Updated upstream
  GetListOfHisto(ccdb, taskname, objname, timestamps_selperiod, vector<long int>(), lnum, 0, isperstave, vector<int>(), vector<int>());
=======
  GetListOfHisto(ccdbApi, taskname, objname, timestamps_selperiod, vector<long int>(), lnum, 0, isperstave, vector<int>(), vector<int>());
>>>>>>> Stashed changes

  timestamps_selperiod.clear();
}


//
// Download data based on run numbers - available in metadata from 03/07/2019 21.49 (run 582 --> fake hit scan)
//
<<<<<<< Updated upstream
void DownloadRuns(CcdbDatabase* ccdb, o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, int lnum, vector<string> runlistfromfile){
=======
void DownloadRuns(CcdbDatabase* ccdb, o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, int lnum, vector<string> runlistfromfile){
>>>>>>> Stashed changes
  
  bool isQCAsync = taskname.find("qc_async")!=string::npos;
  
  //Extract all the time stamps and run numbers of the object
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  string objectlist2 = " ";
  cout<<endl;
  cout<<endl;
  cout<<"Ready to get files from "<<taskname<<"/"<<objname<<endl;

  stringstream ss(objectlist);
  string word;
  vector<string> alltimestamps, timestamps, runs, passnames;

  //filter normal path but correct timestamps with the one from alternative path (if needed)
  while(ss>>word){
<<<<<<< Updated upstream
    if(word=="Created:"){// take the one related to file creation
=======
    if(word=="Validity:"){// take the one related to file creation
>>>>>>> Stashed changes
      ss>>word;
      alltimestamps.push_back(word);
    }
    if(word=="RunNumber"){
      ss>>word;
      ss>>word;
      if(word.size()==3) continue; //protection for fee task in particular: skip runs with 3 digits.
      runs.push_back(word);
      timestamps.push_back(alltimestamps[alltimestamps.size()-1]);//this keep only the timestamps connected to a run number
      if(!isQCAsync && stoi(word)==stoi(run1)) break; // skip in QCasync since runs are in random order in QCDB
    }
    if(isQCAsync && word=="PassName") { //relevant only for qc_async 
      ss>>word;
      ss>>word;
      passnames.push_back(word); 
    }
  }

  if(!isQCAsync) {
    passnames=runs; // just to assign to passnames some random values 
  }

  //filter all runs to get the ones within run1 and run2 (get most recent for each run!!) --> take from good run list if defined or from user-defined list
  vector <long int> timestamps_selperiod;
  vector <int> runs_selperiod;
  int counter=0;

  for(int irun=0; irun<(int)runs.size(); irun++){
    if(runlistfromfile.size()==0 && (stoi(runs[irun])>=stoi(run1) && stoi(runs[irun])<=stoi(run2)) && (goodrunlist.size()==0 || std::find(goodrunlist.begin(), goodrunlist.end(), runs[irun]) != goodrunlist.end()) && (!isQCAsync)){ // NOT DONE IN QC ASYNC
      counter++;
      if(counter==1){
        runs_selperiod.push_back(std::stoi(runs[irun]));
        timestamps_selperiod.push_back(std::stol(timestamps[irun]));
      }
      else{
        if(runs[irun]==runs[irun-1]) continue;
        else{
          runs_selperiod.push_back(std::stoi(runs[irun]));
          timestamps_selperiod.push_back(std::stol(timestamps[irun]));
        }
      }
    }

    else if (runlistfromfile.size()>0 && (stoi(runs[irun])>=stoi(runlistfromfile[runlistfromfile.size()-1]) && stoi(runs[irun])<=stoi(runlistfromfile[0]))) { //user-defined list
      bool isfound = false;
      for(int imatch = 0; imatch<(int)runlistfromfile.size(); imatch++) {
        if(runs[irun] == runlistfromfile[imatch] && (!isQCAsync || selpassname == passnames[irun])) {
          isfound = true;
          break;
        }
        if(!isQCAsync && stoi(runs[irun]) > stoi(runlistfromfile[0])) {
          break;
        }
        if(!isQCAsync && stoi(runs[irun]) < stoi(runlistfromfile[runlistfromfile.size()-1])){
          break;
        }
      }
      if(isfound) {
        counter++;
        if(counter==1){
          runs_selperiod.push_back(std::stoi(runs[irun]));
          timestamps_selperiod.push_back(std::stol(timestamps[irun]));
        }
        else{
          if(!isQCAsync && runs[irun]==runs[irun-1]) continue;
          if(std::find(runs_selperiod.begin(), runs_selperiod.end(), std::stoi(runs[irun])) != runs_selperiod.end()) continue;
          else{
            runs_selperiod.push_back(std::stoi(runs[irun]));
            timestamps_selperiod.push_back(std::stol(timestamps[irun]));
          }
        }
      } // end of if on isfound
    } // end of else if

  } // end of for loop on all runs

  cout<<"\n"<<"Selected runs and corresponding timestamps: "<<endl;
  for(int i=0; i<(int)runs_selperiod.size();i++){
    cout<<"run"<<runs_selperiod[i]<<" - "<<timestamps_selperiod[i]<<endl;
  }

  bool isperstave = 0;
  if(objname.find("HITMAP")!=string::npos) isperstave = 1;
<<<<<<< Updated upstream
  GetListOfHisto(ccdb, taskname, objname, timestamps_selperiod, vector<long int>(), lnum, 1, isperstave, runs_selperiod, vector<int>());
=======
  GetListOfHisto(ccdbApi, taskname, objname, timestamps_selperiod, vector<long int>(), lnum, 1, isperstave, runs_selperiod, vector<int>());
>>>>>>> Stashed changes

  timestamps_selperiod.clear();
  runs_selperiod.clear();
  timestamps.clear();
  runs.clear();
  alltimestamps.clear();

}

//
// Get correct timestamp
//
string GetCorrectTS(string selrun, vector<string> runs, vector<string> timestamps){

  for(int it=0; it<(int)timestamps.size(); it++){
    if(selrun==runs[it]){
      return timestamps[it];
    }
  }
  return "nots";
}

//
//Get list of histogram inside an object
//
<<<<<<< Updated upstream
bool GetListOfHisto(CcdbDatabase* ccdb, string taskname, string objname, vector<long int> timestamps, vector<long int> timestamps2, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers, vector<int>runnumbers2){
=======
bool GetListOfHisto(o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, vector<long int> timestamps, vector<long int> timestamps2, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers, vector<int>runnumbers2){
>>>>>>> Stashed changes

  string qcpathstart;
  if(taskname.find("qc_async")!=string::npos) {
    qcpathstart = "qc_async";
  } else {
    qcpathstart = "qc";
  }
  //Getting root files from the database and write them to file
  cout<<"\n"<<"... Getting files from the database"<<endl;

  string stvnum = "0";
  if(isperstave){
    stvnum = objname.substr(objname.find("Stave")+5,2);
    if(stvnum.find("/")!=string::npos)
      stvnum = objname.substr(objname.find("Stave")+5,1);
  }

<<<<<<< Updated upstream
 // cout<<"timestamps  size: "<<timestamps.size()<<endl;
 // cout<<"timestamps2 size: "<<timestamps2.size()<<endl;

  //in case L2T is missing (excluded)
  if((int)timestamps.size()==0 && (objname.find("Layer2/Threshold_Vs_Chip_and_Stave")!=string::npos || objname.find("Layer2/DeadPixel_Vs_Chip_and_Stave")!=string::npos)){
    for(int i=0; i<(int)timestamps2.size();i++){

      //MonitorObject *monitor = ccdb->retrieve(taskname, objname, timestamps[i]);

      auto monitor = ccdb->retrieveMO(Form("ITS/MO/ITS%sTask2B",objname.find("Chip_and_Stave")!=string::npos ? "THR":"Raw"), objname, timestamps2.size()>0 ? timestamps2[i]:timestamps[i], { 0, 0, "", "", qcpathstart });

      if (monitor == nullptr) {
        cerr << "Failed to get MonitorObject for timestamp: " << timestamps[i]<< endl;
        return 0;
      }

      TObject *obj = nullptr;
      obj = monitor->getObject();
      monitor->setIsOwner(false);
      string c = obj->ClassName();
      TH2 *h2s = 0x0;

      //if(strstr(c,"TH2")!=nullptr){
      if(c.find("TH2")!=string::npos){
        string histname = "";
        histname = Form("h2_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers2[i]) : "", timestamps2[i]);

        h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
        outputfile->cd();
        h2s->Write();
      }
    }
  }

  //General case
  for(int i=0; i<(int)timestamps.size();i++){

    //MonitorObject *monitor = ccdb->retrieve(taskname, objname, timestamps[i]);

    auto monitor = ccdb->retrieveMO(qcpathstart=="qc_async" ? taskname.substr(9) : taskname.substr(3), objname, timestamps[i], { 0, 0, "", "", qcpathstart });

    if (monitor == nullptr) {
      cerr << "Failed to get MonitorObject for timestamp: " << timestamps[i]<< endl;
      return 0;
    }

    TObject *obj = nullptr;
    obj = monitor->getObject();
    monitor->setIsOwner(false);
    //for L2B only
    TObject *obj2 = nullptr;
    if(/*objname.find("Layer2ChipStave")!=string::npos || objname.find("ErrorFile")!=string::npos || objname.find("TriggerFile")!=string::npos
       || */objname.find("Layer2/Threshold_Vs_Chip_and_Stave")!=string::npos || objname.find("Layer2/DeadPixel_Vs_Chip_and_Stave")!=string::npos){
      auto monitor2 = ccdb->retrieveMO(Form("ITS/MO/ITS%sTask2B",objname.find("Chip_and_Stave")!=string::npos ? "THR":"Raw"), objname, timestamps2.size()>0 ? timestamps2[i]:timestamps[i], { 0, 0, "", "", qcpathstart });
      obj2 = monitor2->getObject();
      monitor2->setIsOwner(false);
    }
    string c = obj->ClassName();
    TH2 *h2s = 0x0;
    TH2 *h2sbis = 0x0; //for L2B only
=======

  //General case
  for(int i=0; i<(int)timestamps.size();i++){	  
    std::map<std::string, std::string> metadata;
    metadata["RunNumber"] = std::to_string(runnumbers[i]);
    auto *obj = ccdbApi.retrieveFromTFileAny<TObject>(taskname+"/"+objname, metadata, timestamps[i]);
    if (obj == nullptr) {
      cerr << "Failed to get TObject for timestamp: " << timestamps[i]<< endl;
      return 0;
    }

    string c = obj->ClassName();
    TH2 *h2s = 0x0; 
>>>>>>> Stashed changes
    TH1 *h1s = 0x0;
    THnSparse *hSp = 0x0;
    TTree *tree;
    //////////////////////////////////
    //Cluster per tracks vs eta
    if(objname.find("NClustersPerTrackEta")!=string::npos){
      string histname = "";
      histname = Form("NClustersPerTrackEta_h2_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h2s->Write();
    }
    //////////////////////////////////
    //Cluster fraction
    if(objname.find("AssociatedClusterFraction")!=string::npos){
      string histname = "";
      histname = Form("AssociatedClusterFraction_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }
    //////////////////////////////////
    //Ntracks per event
    if(objname.find("Ntracks")!=string::npos){
      string histname = "";
      histname = Form("Ntracks_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }
    //////////////////////////////////
    //Vertex Contributors
    if(objname.find("NVertexContributors")!=string::npos){
      string histname = "";
      histname = Form("NVertexContributors_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }
    //////////////////////////////////
    //VertexCoordinates (XvsY)
    if(objname.find("VertexCoordinates")!=string::npos){
      string histname = "";
      histname = Form("VertexCoordinates_h2_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h2s->Write();
    }
    //////////////////////////////////
    //VertexRvsZ
    if(objname.find("VertexRvsZ")!=string::npos){
      string histname = "";
      histname = Form("VertexRvsZ_h2_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h2s->Write();
    }
  //////////////////////////////////
  //VertexZ
    if(objname.find("VertexZ")!=string::npos){
      string histname = "";
      histname = Form("VertexZ_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //AngularDistribution
    if(objname.find("AngularDistribution")!=string::npos){
      string histname = "";
      histname = Form("AngularDistribution_h2_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h2s->Write();
    }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //ClusterUsage
    if(objname.find("ClusterUsage")!=string::npos){
      string histname = "";
      histname = Form("ClusterUsage_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //EtaDistribution
   if(objname.find("EtaDistribution")!=string::npos){
      string histname = "";
      histname = Form("EtaDistribution_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //PhiDistribution
    if(objname.find("PhiDistribution")!=string::npos){
      string histname = "";
      histname = Form("PhiDistribution_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //NClusters
    if(objname == "NClusters"){
      string histname = "";
      histname = Form("NClusters_h1_%s%s_%ld", isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//LaneStatusFlag
if(objname.find("LaneStatus") != string::npos){
    string histname = "";
    if(objname.find("LaneStatus/laneStatusFlagERROR")!=string::npos) histname = Form("h2_LSerror%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
    else if(objname.find("LaneStatus/laneStatusFlagFAULT")!=string::npos) histname = Form("h2_LSfault%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
    else if(objname.find("LaneStatus/laneStatusFlagOK")!=string::npos) histname = Form("h2_LSok%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
    else if(objname.find("LaneStatus/laneStatusFlagWARNING")!=string::npos) histname = Form("h2_LSwarning%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
    else if(objname.find("LaneStatus/laneStatusFlagCumulativeERROR")!=string::npos) histname = Form("h2_LCSerror%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
    else if(objname.find("LaneStatus/laneStatusFlagCumulativeFAULT")!=string::npos) histname = Form("h2_LCSfault%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
    else if(objname.find("LaneStatus/laneStatusFlagCumulativeWARNING")!=string::npos) histname = Form("h2_LCSwarning%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
    h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
    outputfile->cd();
    h2s->Write();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    //if(strstr(c,"TH1")!=nullptr){
    if(c.find("TH1")!=string::npos && taskname.find("Track")==string::npos){
      string histname = "";
      histname = Form("h1_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }
    //if(strstr(c,"TH2")!=nullptr){
    if(c.find("TH2")!=string::npos && taskname.find("Track")==string::npos && objname.find("LaneStatus")==string::npos){
      string histname = "";
      if(objname.find("Error")!=string::npos)
        histname = Form("h2_err%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
      else if(objname.find("Trigger")!=string::npos)
        histname = Form("h2_trg%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
<<<<<<< Updated upstream
      else
        histname = Form("h2_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);

      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
      if(/*objname.find("Layer2ChipStave")!=string::npos || objname.find("ErrorFile")!=string::npos || objname.find("TriggerFile")!=string::npos
        || */objname.find("Layer2/Threshold_Vs_Chip_and_Stave")!=string::npos || objname.find("Layer2/DeadPixel_Vs_Chip_and_Stave")!=string::npos){
        h2sbis = dynamic_cast<TH2*>(obj2->Clone(Form("%s_L2B",histname.c_str())));
        h2s->Add(h2sbis);
      }
=======
      else if(objname.find("RDHSummaryCumulative")!=string::npos)
	    histname = Form("h2_RDH%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
      else if(objname.find("Thr")!=string::npos){
          histname = Form("h2%s%s_%s_%ld", isperstave ? Form("_Stv%s", stvnum.c_str()) : "", isrunknown ? Form("_run%d", runnumbers[i]) : "", objname.substr(objname.size() - 2, 2).c_str(), timestamps[i]);
      }
      else {
	    histname = Form("h2_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
      }
      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str())); 
>>>>>>> Stashed changes
      outputfile->cd();
      h2s->Write();
    }

    if(c.find("THnSparse")!=string::npos){
      string histname = "";
      histname = Form("hsparse_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);

      hSp = dynamic_cast<THnSparse*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      hSp->Write();
    }

    if(c.find("TTree")!=string::npos){
      string treename = Form("ttree_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);;
      //obj->Print();
      tree = dynamic_cast<TTree*>(obj->Clone(treename.c_str()));
      //tree->Print();
      outputfile->cd();
      tree->Write();
    }

    delete h1s;
    delete h2s;
    delete hSp;
    delete obj;
  }

  return 1;
}

//
// Download depending on the choice
//
<<<<<<< Updated upstream
bool Download(int choice, CcdbDatabase* ccdb, o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, long int ts_start, long int ts_end, int lnum, vector<string> runlistfromfile){
=======
bool Download(int choice, CcdbDatabase* ccdb, o2::ccdb::CcdbApi &ccdbApi, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, long int ts_start, long int ts_end, int lnum, vector<string> runlistfromfile){
>>>>>>> Stashed changes

  switch(choice){
    case 1: {
      DownloadRuns(ccdb, ccdbApi, taskname, objname, run1, run2, goodrunlist, lnum, vector<string>());//download with runs
      break;//download runs data
    }

    case 2: {
      DownloadTimestamps(ccdb, ccdbApi, taskname, objname, ts_start, ts_end, lnum); //download with timestamps
      break;
    }

    case 3: {
      DownloadRuns(ccdb, ccdbApi, taskname, objname, run1, run2, goodrunlist, lnum, runlistfromfile);//download with runs from input run list
      break;
    }

    default: return 0;
  }

  return 1;

}

//
// Get name depending on the option
//
string GetOptName(int opt){
  switch(opt){
     case 1: return "FHRMAPS_HITMAPS";
     case 2: return "THRMAPS_DEADPIXMAPS";
     // case 3: return "NOISYPIX_TREE";
     case 3: return "TrackTask";
     case 4: return "LaneStatusFlag";
     case 5: return "ClusterTask";
     default: return "0";
  }
}

//
// Get name to assign to the output list
//
string GetListName(int opt, int ilist){
  switch(opt){
    case 1:{
      if(!ilist)
        return "fhrmaps";
      else
        return "hitmaps";
    }
    default: return "0";
  }
}
