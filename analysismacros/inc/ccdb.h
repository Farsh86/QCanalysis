//This header file can serve as configuration for the databse where the analysis macros push their output to facilitate future updates.

#include <string>
#include "QualityControl/CcdbDatabase.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"

using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;
using namespace std;

const string ccdbport = "ccdb-test.cern.ch:8080";
const string DetectorName = "ITS";
string TaskName = "QC_Offline/";
string TaskClass = "OfflineQC";
/*CcdbDatabase* SetupConnection(){
	std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");
	auto ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());
	ccdb->connect(ccdbport.c_str(), "", "", "");
	return ccdb;}
*/
void SetTaskName(string funcname){TaskName.append(funcname);}
