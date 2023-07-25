#include "QualityControl/CcdbDatabase.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/RootClassFactory.h"
#include "inc/ccdb.h"
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TKey.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
void SetCommonAxis(TH2D *histo);
void SetStyle(TGraph *h, Int_t col, Style_t mkr);
void DoAnalysis(string filepath, const int nChips, string skipruns,
                bool ccdb_upload);
TString SIBorOB[2] = {"IB", "OB"};

int ns[] = {12, 16, 20, 24, 30, 42, 48}; // n staves for each layer

//
// MAIN
//
void AnalyzeLaneStatusFlagExtended() {
  // gStyle->SetPalette(55);
  string fpath;
  int nchips = 9;
  cout << "\n\n=> Available file(s) for the analysis (the last should be the "
          "file you want!): \n"
       << endl;
  gSystem->Exec("ls ../Data/*LaneStatus* -Art | tail -n 500");
  cout << "\nCopy file name: ";
  cin >> fpath;
  cout << endl;
  bool ccdb_upload;
  // Choose whether to skip runs
  string skipans, skipruns, CCDB_up;
  cout << endl;
  cout << "Would you like to skip some run(s)? [y/n] ";
  cin >> skipans;
  if (skipans == "y" || skipans == "Y") {
    cout << endl;
    cout << "Specify run number(s) separated by comma (no white spaces!):";
    cin >> skipruns;
    cout << endl;
  } else
    skipruns = " ";
  cout << "Would you like to upload the output to ccdb? [y/n] ";
  cin >> CCDB_up;
  cout << endl;
  if (CCDB_up == "y" || CCDB_up == "Y")
    ccdb_upload = true;
  else
    ccdb_upload = false;

  if (ccdb_upload)
    SetTaskName(__func__);

  // Call
  DoAnalysis(fpath, nchips, skipruns, ccdb_upload);
}

double GetDispersionOfTriggers(TH2 *H, int biny, int firstx, int lastx){
  // TODO: can be improved playing with projections and TH1 methods
  double xm=0, x2m=0;
  for (int i=firstx; i<lastx+1; i++){
    double z = H->GetBinContent(i,biny)*1.e-6;
    xm += z/(lastx-firstx+1);
    x2m += z*z/(lastx-firstx+1);
  }
  if (xm>0) return TMath::Sqrt(x2m-xm*xm+1.e-8) / xm;  // 1.e-something is a protection needed.
  else return 9.99;
}

//
// Set Style
//
void SetStyle(TGraph *h, Int_t col, Style_t mkr) {
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  // h->SetFillStyle(0);
  // h->SetFillColorAlpha(col,0.8);
}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, string skipruns,
                bool ccdb_upload) {

  gStyle->SetOptStat(0000);

  // std::vector<TH2*> herr, hfault, hok, hwarning;
  std::vector<TH2 *> herr, herrcumulative;
  std::vector<TH2 *> hfault, hfaultcumulative;
  std::vector<TH2 *> hok;
  std::vector<TH2 *> hwarning;
  std::vector<TH2 *> htrg;
  std::vector<string> runnumbersErr, runnumbersFaul, runnumbersWarn, runnumbersCumul, runnumbersTrg; // only runnumbersErr is used later in the macro. TODO: protect in case some of the plots are missing?
  std::vector<string> timestamps1, timestamps2, timestamps3, timestamps4; // not in use
  int col[] = {810, 807, 797, 827, 417, 841, 868, 867, 860, 602};
  //	Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"),
  //TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"),
  //TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"),
  //TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"),
  //TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");

  auto ccdb = dynamic_cast<CcdbDatabase *>(mydb.get());

  ccdb->connect(ccdbport.c_str(), "", "", "");

  // Read the file and the list of plots with entries
  TFile *infile = new TFile(filepath.c_str());
  TList *list = infile->GetListOfKeys();
  TKey *key;
  TObject *obj;
  TIter next(list);
  TH2 *h2err;
  TH2 *h2fault;
  TH2 *h2errcumulative;
  TH2 *h2faultcumulative;
  TH2 *h2warning;
  TH2 *h2trg;
  while ((key = ((TKey *)next()))) {
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(), "TProfile") != 0) &&
        (!obj->InheritsFrom("TH2")) && (!obj->InheritsFrom("TH1"))) {
      cout << "<W> Object " << obj->GetName()
           << " is not 1D or 2D histogram : will not be converted" << endl;
    }

    string objname = (string)obj->GetName();
    // if(objname.find("trg")!=string::npos) continue; //skip trigger flags plot
    // here
    cout << "objname = " << objname << endl;
    // if(objname.find("trg")==string::npos) continue;
    string timestamp = objname.find("run") == string::npos
                           ? objname.substr(objname.find("_", 2) + 1, 13)
                           : objname.substr(objname.find("_", 6) + 1, 13);
    string runnum = objname.find("run") == string::npos
                        ? "norun"
                        : objname.substr(objname.find("run") + 3, 6);

    if (skipruns.find(runnum) != string::npos)
      continue; // eventually skip runs specified by the user
    if (objname.find("LSerror") != string::npos) {
      h2err = (TH2 *)obj;
      cout << "... Reading " << obj->GetName() << endl;
      h2err->ClearUnderflowAndOverflow();
      herr.push_back(h2err);
      timestamps1.push_back(timestamp);
      runnumbersErr.push_back(runnum);
    }
    if (objname.find("LSfault") != string::npos) {
      h2fault = (TH2 *)obj;
      cout << "... Reading " << obj->GetName() << endl;
      h2fault->ClearUnderflowAndOverflow();
      hfault.push_back(h2fault);
      timestamps2.push_back(timestamp);
      runnumbersFaul.push_back(runnum);
    }
    if (objname.find("LSwarning") != string::npos) {
      h2warning = (TH2 *)obj;
      cout << "... Reading " << obj->GetName() << endl;
      h2warning->ClearUnderflowAndOverflow();
      hwarning.push_back(h2warning);
      timestamps4.push_back(timestamp);
      runnumbersWarn.push_back(runnum);
    }
    if (objname.find("trg") != string::npos) {
      h2trg = (TH2 *)obj;
      cout << "... Reading " << obj->GetName() << endl;
      h2trg->ClearUnderflowAndOverflow();
      htrg.push_back(h2trg);
      timestamps3.push_back(timestamp);
      runnumbersTrg.push_back(runnum);
    }
    if (objname.find("LCSerror") != string::npos) {
      h2errcumulative = (TH2 *)obj;
      cout << "... Reading " << obj->GetName() << endl;
      h2errcumulative->ClearUnderflowAndOverflow();
      herrcumulative.push_back(h2errcumulative);
    }
    if (objname.find("LCSfault") != string::npos) {
      h2faultcumulative = (TH2 *)obj;
      cout << "... Reading " << obj->GetName() << endl;
      h2faultcumulative->ClearUnderflowAndOverflow();
      hfaultcumulative.push_back(h2faultcumulative);
    }
  }

  int nRuns = (int)runnumbersErr.size();
  for (int i = 0; i < (int)runnumbersErr.size(); i++)
    cout << "RUNNUM=" << runnumbersErr[i] << endl;

  // sum all the histos in a single histogram (for summary plot) for each layer
  // TODO: does it still make sense to use the reset plots?
  TH2D *hSummary1 = (TH2D *)herr[0]->Clone("hSummary1");
  for (int iplot = 1; iplot < (int)herr.size(); iplot++) {
    hSummary1->Add(herr[iplot]);
  }
  TH2D *hSummary2 = (TH2D *)hfault[0]->Clone("hSummary2");
  for (int iplot = 1; iplot < (int)hfault.size(); iplot++) {
    hSummary2->Add(hfault[iplot]);
  }
  TH2D *hSummary3 = (TH2D *)hwarning[0]->Clone("hSummary3");
  for (int iplot = 1; iplot < (int)hwarning.size(); iplot++) {
    hSummary3->Add(hwarning[iplot]);
  }

  // Draw summary plot
  TCanvas canvas;
  canvas.cd();
  canvas.SetTickx();
  canvas.SetTicky();
  canvas.SetLogz();
  canvas.SetMargin(0.18, 0.2, 0.194, 0.0993);
  canvas.SetRightMargin(0.15);
  hSummary1->SetTitle(
      Form("Lane Status Flag ERROR at the end of the run, %s",
           filepath
               .substr(filepath.find("from"),
                       filepath.find("_w_") - filepath.find("from"))
               .c_str()));
  hSummary1->Draw("colz");
  SetCommonAxis(hSummary1);

  TCanvas canvas12 = (TCanvas)canvas.Clone("canvas12");
  canvas12.SetMargin(0.18, 0.2, 0.194, 0.0993);
  canvas12.SetRightMargin(0.15);
  hSummary2->SetTitle(
      Form("Lane Status Flag FAULT at the end of the run, %s",
           filepath
               .substr(filepath.find("from"),
                       filepath.find("_w_") - filepath.find("from"))
               .c_str()));
  hSummary2->Draw("colz");
  SetCommonAxis(hSummary2);

  TCanvas canvas13 = (TCanvas)canvas.Clone("canvas13");
  canvas13.SetMargin(0.18, 0.2, 0.194, 0.0993);
  canvas13.SetRightMargin(0.15);
  hSummary3->SetTitle(
      Form("Lane Status Flag WARNING at the end of the run, %s",
           filepath
               .substr(filepath.find("from"),
                       filepath.find("_w_") - filepath.find("from"))
               .c_str()));
  hSummary3->Draw("colz");
  SetCommonAxis(hSummary3);



  if (ccdb_upload) {
    string Runperiod =
        Form("%s", filepath.substr(filepath.find("from"), 27).c_str());
    int RunNumber =
        std::stoi(filepath.substr(filepath.find("run") + 3, 6).c_str());
    canvas.SetName("Summary_Lane_Status_Flag_ERROR");
    auto mo_err = std::make_shared<o2::quality_control::core::MonitorObject>(
        &canvas, TaskName, TaskClass, DetectorName, RunNumber, Runperiod);
    mo_err->setIsOwner(false);
    ccdb->storeMO(mo_err);

    canvas12.SetName("Summary_Lane_Status_Flag_FAULT");
    auto mo_fault = std::make_shared<o2::quality_control::core::MonitorObject>(
        &canvas12, TaskName, TaskClass, DetectorName, RunNumber, Runperiod);
    mo_fault->setIsOwner(false);
    ccdb->storeMO(mo_fault);

    canvas13.SetName("Summary_Lane_Status_Flag_WARNING");
    auto mo_warn = std::make_shared<o2::quality_control::core::MonitorObject>(
        &canvas13, TaskName, TaskClass, DetectorName, RunNumber, Runperiod);
    mo_warn->setIsOwner(false);
    ccdb->storeMO(mo_warn);
  }

  canvas.SaveAs(Form("../Plots/LaneStatusFlag_%s.pdf[",
                     filepath
                         .substr(filepath.find("from"),
                                 filepath.find(".root") - filepath.find("from"))
                         .c_str()));
  canvas.SaveAs(Form("../Plots/LaneStatusFlag_%s.pdf",
                     filepath
                         .substr(filepath.find("from"),
                                 filepath.find(".root") - filepath.find("from"))
                         .c_str()));
  canvas12.SaveAs(
      Form("../Plots/LaneStatusFlag_%s.pdf",
           filepath
               .substr(filepath.find("from"),
                       filepath.find(".root") - filepath.find("from"))
               .c_str()));
  canvas13.SaveAs(
      Form("../Plots/LaneStatusFlag_%s.pdf",
           filepath
               .substr(filepath.find("from"),
                       filepath.find(".root") - filepath.find("from"))
               .c_str()));


  // Draw trend plots 3status x 7Layer
  const int NRun = (int)herr.size();
  const int NStatus = 3;
  const int NLayer = 7;
  TString StatusKind[NStatus] = {"ERROR", "FAULT", "WARNING"};
  double binedge_L0[13] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36};
  double binedge_L1[17] = {36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84};
  double binedge_L2[21] = {84,  87,  90,  93,  96,  99,  102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144};
  double binedge_L3[25] = {144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192};
  double binedge_L4[31] = {192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252};
  double binedge_L5[43] = {252, 254, 256, 258, 260, 262, 264, 266, 268, 270, 272, 274, 276, 278, 280, 282, 284, 286, 288, 290, 292, 294, 296, 298, 300, 302, 304, 306, 308, 310, 312, 314, 316, 318, 320, 322, 324, 326, 328, 330, 332, 334, 336};
  double binedge_L6[49] = {336, 338, 340, 342, 344, 346, 348, 350, 352, 354, 356, 358, 360, 362, 364, 366, 368, 370, 372, 374, 376, 378, 380, 382, 384, 386, 388, 390, 392, 394, 396, 398, 400, 402, 404, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432};
  const int nstep[NLayer] = {3, 3, 3, 2, 2, 2, 2}; // 3 FEEID's for a stave in IB, 2 FEEID's for OB
  const int nstave[NLayer] = {12, 16, 20, 24, 30, 42, 48}; // Number of Staves in a Layer
  TH1D *hfeex[NStatus][NRun];
  TH1D *hrebin[NStatus][NRun][NLayer];
  TH1D *htimeratio[NRun];
  double nOrbitDispersion[NRun][NLayer];
  double tmp_nonzero = 0;
  double tmp_nol = 0;
  double hbfNOK = 0;
  for (int iRun = 0; iRun < NRun; iRun++) {
    htimeratio[iRun] = new TH1D(Form("timeratio%d", iRun),
                                Form("timeratio%d", iRun), 432, 0, 432);
    ///////////////the order of runnum reversed in the root file, so reverting
    ///the order for figure.///////////////////////
    int iRev = NRun - iRun - 1;
    hfeex[0][iRun] = (TH1D *)herr[iRev]->ProjectionX(); // total number of error
    hfeex[1][iRun] =
        (TH1D *)hfault[iRev]->ProjectionX(); // total number of fault
    hfeex[2][iRun] =
        (TH1D *)hwarning[iRev]->ProjectionX(); // total number of warning
    for (int iStatus = 0; iStatus < NStatus; iStatus++) {
      hfeex[iStatus][iRun]->Reset(); // reset
      for (int iFid = 0; iFid < 432; iFid++) {
        tmp_nol = 0;
        for (int iLane = 0; iLane < 28; iLane++) {
          if (iStatus == 0)
            tmp_nonzero = herr[iRev]->GetBinContent(iFid + 1, iLane + 1);
          if (iStatus == 1)
            tmp_nonzero = hfault[iRev]->GetBinContent(iFid + 1, iLane + 1);
          if (iStatus == 2)
            tmp_nonzero = hwarning[iRev]->GetBinContent(iFid + 1, iLane + 1);
          if (tmp_nonzero > 0)
            tmp_nol++; // check number of nonzero lanes
        }
        hfeex[iStatus][iRun]->SetBinContent(iFid + 1, tmp_nol);
      }
    }

    for (int iFid = 0; iFid < 432; iFid++) {
      
      hbfNOK = 0;

      for (int iLane = 0; iLane < 28; iLane++) {
        hbfNOK += herrcumulative[iRev]->GetBinContent(iFid + 1, iLane + 1) + hfaultcumulative[iRev]->GetBinContent(iFid + 1, iLane + 1);
      }

      double Norbit = htrg[iRev]->GetBinContent(iFid + 1, 1); // Number of ORBIT triggers seen by this FEEID

      int laneinfeeid = iFid < 145 ? 3 : iFid < 252 ? 8 : 14;
      if (Norbit > 0)
        hbfNOK /= (Norbit * laneinfeeid);
      else
        hbfNOK = -1;

      htimeratio[iRun]->SetBinContent(iFid + 1, hbfNOK);
    }

    nOrbitDispersion[iRun][0] = GetDispersionOfTriggers( htrg[iRev], 1, binedge_L0[0] +1 , binedge_L0[nstave[0]]);
    nOrbitDispersion[iRun][1] = GetDispersionOfTriggers( htrg[iRev], 1, binedge_L1[0] +1 , binedge_L1[nstave[1]]);
    nOrbitDispersion[iRun][2] = GetDispersionOfTriggers( htrg[iRev], 1, binedge_L2[0] +1 , binedge_L2[nstave[2]]);
    nOrbitDispersion[iRun][3] = GetDispersionOfTriggers( htrg[iRev], 1, binedge_L3[0] +1 , binedge_L3[nstave[3]]);
    nOrbitDispersion[iRun][4] = GetDispersionOfTriggers( htrg[iRev], 1, binedge_L4[0] +1 , binedge_L4[nstave[4]]);
    nOrbitDispersion[iRun][5] = GetDispersionOfTriggers( htrg[iRev], 1, binedge_L5[0] +1 , binedge_L5[nstave[5]]);
    nOrbitDispersion[iRun][6] = GetDispersionOfTriggers( htrg[iRev], 1, binedge_L6[0] +1 , binedge_L6[nstave[6]]);

    

    for (int iStatus = 0; iStatus < NStatus; iStatus++) {
      hrebin[iStatus][iRun][0] = (TH1D *)hfeex[iStatus][iRun]->Rebin(
          12, Form("hrebin%d_%d_L0", iStatus, iRun), binedge_L0);
      hrebin[iStatus][iRun][1] = (TH1D *)hfeex[iStatus][iRun]->Rebin(
          16, Form("hrebin%d_%d_L1", iStatus, iRun), binedge_L1);
      hrebin[iStatus][iRun][2] = (TH1D *)hfeex[iStatus][iRun]->Rebin(
          20, Form("hrebin%d_%d_L2", iStatus, iRun), binedge_L2);
      hrebin[iStatus][iRun][3] = (TH1D *)hfeex[iStatus][iRun]->Rebin(
          24, Form("hrebin%d_%d_L3", iStatus, iRun), binedge_L3);
      hrebin[iStatus][iRun][4] = (TH1D *)hfeex[iStatus][iRun]->Rebin(
          30, Form("hrebin%d_%d_L4", iStatus, iRun), binedge_L4);
      hrebin[iStatus][iRun][5] = (TH1D *)hfeex[iStatus][iRun]->Rebin(
          42, Form("hrebin%d_%d_L5", iStatus, iRun), binedge_L5);
      hrebin[iStatus][iRun][6] = (TH1D *)hfeex[iStatus][iRun]->Rebin(
          48, Form("hrebin%d_%d_L6", iStatus, iRun), binedge_L6);
    }
  }

  // PLOTS WITH NUMBER OF LANES IN E/F/W AT THE END OF RUN
  TGraph *grt_L0[NStatus][nstave[0]];
  TGraph *grt_L1[NStatus][nstave[1]];
  TGraph *grt_L2[NStatus][nstave[2]];
  TGraph *grt_L3[NStatus][nstave[3]];
  TGraph *grt_L4[NStatus][nstave[4]];
  TGraph *grt_L5[NStatus][nstave[5]];
  TGraph *grt_L6[NStatus][nstave[6]];

  for (int iStatus = 0; iStatus < NStatus; iStatus++) {
    for (int iLayer = 0; iLayer < NLayer; iLayer++) {
      for (int iStave = 0; iStave < nstave[iLayer]; iStave++) {
        if (iLayer == 0) {
          grt_L0[iStatus][iStave] = new TGraph();
          SetStyle(grt_L0[iStatus][iStave], col[iStave % 10], 24 + iStave / 10);
        }
        if (iLayer == 1) {
          grt_L1[iStatus][iStave] = new TGraph();
          SetStyle(grt_L1[iStatus][iStave], col[iStave % 10], 24 + iStave / 10);
        }
        if (iLayer == 2) {
          grt_L2[iStatus][iStave] = new TGraph();
          SetStyle(grt_L2[iStatus][iStave], col[iStave % 10], 24 + iStave / 10);
        }
        if (iLayer == 3) {
          grt_L3[iStatus][iStave] = new TGraph();
          SetStyle(grt_L3[iStatus][iStave], col[iStave % 10], 24 + iStave / 10);
        }
        if (iLayer == 4) {
          grt_L4[iStatus][iStave] = new TGraph();
          SetStyle(grt_L4[iStatus][iStave], col[iStave % 10], 24 + iStave / 10);
        }
        if (iLayer == 5) {
          grt_L5[iStatus][iStave] = new TGraph();
          SetStyle(grt_L5[iStatus][iStave], col[iStave % 10], 24 + iStave / 10);
        }
        if (iLayer == 6) {
          grt_L6[iStatus][iStave] = new TGraph();
          SetStyle(grt_L6[iStatus][iStave], col[iStave % 10], 24 + iStave / 10);
        }

        for (int iRun = 0; iRun < NRun; iRun++) {
          double binc =
              hrebin[iStatus][iRun][iLayer]->GetBinContent(iStave + 1);
          if (binc < 1e-15)
            binc = -20.;
          if (iLayer == 0)
            grt_L0[iStatus][iStave]->SetPoint(iRun, iRun, binc > 9 ? 9 : binc);
          if (iLayer == 1)
            grt_L1[iStatus][iStave]->SetPoint(iRun, iRun, binc > 9 ? 9 : binc);
          if (iLayer == 2)
            grt_L2[iStatus][iStave]->SetPoint(iRun, iRun, binc > 9 ? 9 : binc);
          if (iLayer == 3)
            grt_L3[iStatus][iStave]->SetPoint(iRun, iRun,
                                              binc > 16 ? 16 : binc);
          if (iLayer == 4)
            grt_L4[iStatus][iStave]->SetPoint(iRun, iRun,
                                              binc > 16 ? 16 : binc);
          if (iLayer == 5)
            grt_L5[iStatus][iStave]->SetPoint(iRun, iRun,
                                              binc > 28 ? 28 : binc);
          if (iLayer == 6)
            grt_L6[iStatus][iStave]->SetPoint(iRun, iRun,
                                              binc > 28 ? 28 : binc);
          cout << "istatus" << iStatus << " irun" << iRun << " ilayer" << iLayer
               << " istave" << iStave << "  val=" << binc << endl;
        }
      }
    }
  }
  ///// PLOT WITH E/F AVERAGED OVER TIME
  TGraph *time_L0[nstave[0]];
  TGraph *time_L1[nstave[1]];
  TGraph *time_L2[nstave[2]];
  TGraph *time_L3[nstave[3]];
  TGraph *time_L4[nstave[4]];
  TGraph *time_L5[nstave[5]];
  TGraph *time_L6[nstave[6]];

  for (int iLayer = 0; iLayer < NLayer; iLayer++) {
    for (int iStave = 0; iStave < nstave[iLayer]; iStave++) {
      if (iLayer == 0) {
        time_L0[iStave] = new TGraph();
        SetStyle(time_L0[iStave], col[iStave % 10], 24 + iStave / 10);
      }
      if (iLayer == 1) {
        time_L1[iStave] = new TGraph();
        SetStyle(time_L1[iStave], col[iStave % 10], 24 + iStave / 10);
      }
      if (iLayer == 2) {
        time_L2[iStave] = new TGraph();
        SetStyle(time_L2[iStave], col[iStave % 10], 24 + iStave / 10);
      }
      if (iLayer == 3) {
        time_L3[iStave] = new TGraph();
        SetStyle(time_L3[iStave], col[iStave % 10], 24 + iStave / 10);
      }
      if (iLayer == 4) {
        time_L4[iStave] = new TGraph();
        SetStyle(time_L4[iStave], col[iStave % 10], 24 + iStave / 10);
      }
      if (iLayer == 5) {
        time_L5[iStave] = new TGraph();
        SetStyle(time_L5[iStave], col[iStave % 10], 24 + iStave / 10);
      }
      if (iLayer == 6) {
        time_L6[iStave] = new TGraph();
        SetStyle(time_L6[iStave], col[iStave % 10], 24 + iStave / 10);
      }

      for (int iRun = 0; iRun < NRun; iRun++) {
        double num = 0;
        int den = 0;
        if (iLayer == 0) {
          for (int ifee = binedge_L0[iStave]; ifee < binedge_L0[iStave + 1]; ifee++)
            if (htimeratio[iRun]->GetBinContent(ifee + 1) > -1) {  // -1 meas that the corresponding number of triggers was not found. 
              den++;
              num += htimeratio[iRun]->GetBinContent(ifee + 1);
            }
          time_L0[iStave]->SetPoint(iRun, iRun, den == 0 ? 1 : num / den + 1.e-7); // adding extra 1.e-7 not to have zero and to keep the point displayed in the log-y graph. Set ratio to 1 in case the full stave has no trigger entries
        }
        if (iLayer == 1) {
          for (int ifee = binedge_L1[iStave]; ifee < binedge_L1[iStave + 1]; ifee++)
            if (htimeratio[iRun]->GetBinContent(ifee + 1) > -1) {
              den++;
              num += htimeratio[iRun]->GetBinContent(ifee + 1);
            }
          time_L1[iStave]->SetPoint(iRun, iRun, den == 0 ? 1 : num / den + 1.e-7);
        }
        if (iLayer == 2) {
          for (int ifee = binedge_L2[iStave]; ifee < binedge_L2[iStave + 1]; ifee++)
            if (htimeratio[iRun]->GetBinContent(ifee + 1) > -1) {
              den++;
              num += htimeratio[iRun]->GetBinContent(ifee + 1);
            }
          time_L2[iStave]->SetPoint(iRun, iRun, den == 0 ? 1 : num / den + 1.e-7);
        }
        if (iLayer == 3) {
          for (int ifee = binedge_L3[iStave]; ifee < binedge_L3[iStave + 1]; ifee++)
            if (htimeratio[iRun]->GetBinContent(ifee + 1) > -1) {
              den++;
              num += htimeratio[iRun]->GetBinContent(ifee + 1);
            }
          time_L3[iStave]->SetPoint(iRun, iRun, den == 0 ? 1 : num / den + 1.e-7);
        }
        if (iLayer == 4) {
          for (int ifee = binedge_L4[iStave]; ifee < binedge_L4[iStave + 1]; ifee++)
            if (htimeratio[iRun]->GetBinContent(ifee + 1) > -1) {
              den++;
              num += htimeratio[iRun]->GetBinContent(ifee + 1);
            }
          time_L4[iStave]->SetPoint(iRun, iRun, den == 0 ? 1 : num / den + 1.e-7);
        }
        if (iLayer == 5) {
          for (int ifee = binedge_L5[iStave]; ifee < binedge_L5[iStave + 1]; ifee++)
            if (htimeratio[iRun]->GetBinContent(ifee + 1) > -1) {
              den++;
              num += htimeratio[iRun]->GetBinContent(ifee + 1);
            }
          time_L5[iStave]->SetPoint(iRun, iRun, den == 0 ? 1 : num / den + 1.e-7);
        }
        if (iLayer == 6) {
          for (int ifee = binedge_L6[iStave]; ifee < binedge_L6[iStave + 1]; ifee++)
            if (htimeratio[iRun]->GetBinContent(ifee + 1) > -1) {
              den++;
              num += htimeratio[iRun]->GetBinContent(ifee + 1);
            }
          time_L6[iStave]->SetPoint(iRun, iRun, den == 0 ? 1 : num / den + 1.e-7);
        }


        cout << "Time average irun" << iRun << " ilayer" << iLayer << " istave"<< iStave << "  value:" << num / den << endl;
      }
    }
  }



  TCanvas *ctrend2[NStatus][NLayer];
  TCanvas *ctrend3[NLayer]; //time averaged plots

  TLegend *legLayer[NLayer];
  for (int iLayer = 0; iLayer < NLayer; iLayer++) {
    int tempstave = nstave[iLayer];
    legLayer[iLayer] = new TLegend(0.904, 0.197, 0.997, 0.898);
    //	legLayer[iLayer]->SetHeader(Form("Layer%d Stave",iLayer));
    legLayer[iLayer]->SetTextSize(0.04);
    if (iLayer > 3)
      legLayer[iLayer]->SetNColumns(2);
    for (int iStave = 0; iStave < tempstave; iStave++) {
      if (iLayer == 0)
        legLayer[iLayer]->AddEntry(grt_L0[0][iStave], Form(" %d", iStave), "p");
      if (iLayer == 1)
        legLayer[iLayer]->AddEntry(grt_L1[0][iStave], Form(" %d", iStave), "p");
      if (iLayer == 2)
        legLayer[iLayer]->AddEntry(grt_L2[0][iStave], Form(" %d", iStave), "p");
      if (iLayer == 3)
        legLayer[iLayer]->AddEntry(grt_L3[0][iStave], Form(" %d", iStave), "p");
      if (iLayer == 4)
        legLayer[iLayer]->AddEntry(grt_L4[0][iStave], Form(" %d", iStave), "p");
      if (iLayer == 5)
        legLayer[iLayer]->AddEntry(grt_L5[0][iStave], Form(" %d", iStave), "p");
      if (iLayer == 6)
        legLayer[iLayer]->AddEntry(grt_L6[0][iStave], Form(" %d", iStave), "p");
    }
  }

  TH1F *hblank[NStatus][NLayer];
  for (int iStatus = 0; iStatus < NStatus; iStatus++) {
    for (int iLayer = 0; iLayer < NLayer; iLayer++) {
      hblank[iStatus][iLayer] =
          new TH1F(Form("hblank_%d_%d", iStatus, iLayer),
                   Form("Layer %d - lanes into %s; Run; #Lanes into %s", iLayer,
                        StatusKind[iStatus].Data(), StatusKind[iStatus].Data()),
                   NRun, -0.5, (double)NRun - 0.5);
      for (int ir = 0; ir < (int)runnumbersErr.size(); ir++)
        hblank[iStatus][iLayer]->GetXaxis()->SetBinLabel(ir + 1, Form("%06d",stoi(runnumbersErr[runnumbersErr.size() - 1 - ir]))); // runnumbersErr is a descending order

      hblank[iStatus][iLayer]->GetYaxis()->SetRangeUser(0, 30); // for number of nonzero lanes in each FEEID
      hblank[iStatus][iLayer]->GetXaxis()->SetTitleOffset(2.8);
      hblank[iStatus][iLayer]->GetXaxis()->LabelsOption("v");
      hblank[iStatus][iLayer]->GetXaxis()->SetLabelOffset(0.01);
      ctrend2[iStatus][iLayer] = new TCanvas();
      ctrend2[iStatus][iLayer]->cd();
      ctrend2[iStatus][iLayer]->SetTickx();
      ctrend2[iStatus][iLayer]->SetTicky();
      // ctrend2[iStatus][iLayer]->SetLogy();
      ctrend2[iStatus][iLayer]->SetMargin(0.0988, 0.1, 0.194, 0.0993);
      hblank[iStatus][iLayer]->Draw();
      // gPad->SetLogy();//for total number of errors
      gPad->SetGridx();

      if (iLayer == 0) {
        for (int iStave = 0; iStave < nstave[0]; iStave++) {
          grt_L0[iStatus][iStave]->SetName(
              Form("L%d_St%d_F%s", iLayer, iStave, StatusKind[iStatus].Data()));
          grt_L0[iStatus][iStave]->Draw("P");
        }
      }
      if (iLayer == 1) {
        for (int iStave = 0; iStave < nstave[1]; iStave++) {
          grt_L1[iStatus][iStave]->SetName(
              Form("L%d_St%d_F%s", iLayer, iStave, StatusKind[iStatus].Data()));
          grt_L1[iStatus][iStave]->Draw("P");
        }
      }
      if (iLayer == 2) {
        for (int iStave = 0; iStave < nstave[2]; iStave++) {
          grt_L2[iStatus][iStave]->SetName(
              Form("L%d_St%d_F%s", iLayer, iStave, StatusKind[iStatus].Data()));
          grt_L2[iStatus][iStave]->Draw("P");
        }
      }
      if (iLayer == 3) {
        for (int iStave = 0; iStave < nstave[3]; iStave++) {
          grt_L3[iStatus][iStave]->SetName(
              Form("L%d_St%d_F%s", iLayer, iStave, StatusKind[iStatus].Data()));
          grt_L3[iStatus][iStave]->Draw("P");
        }
      }
      if (iLayer == 4) {
        for (int iStave = 0; iStave < nstave[4]; iStave++) {
          grt_L4[iStatus][iStave]->SetName(
              Form("L%d_St%d_F%s", iLayer, iStave, StatusKind[iStatus].Data()));
          grt_L4[iStatus][iStave]->Draw("P");
        }
      }
      if (iLayer == 5) {
        for (int iStave = 0; iStave < nstave[5]; iStave++) {
          grt_L5[iStatus][iStave]->SetName(
              Form("L%d_St%d_F%s", iLayer, iStave, StatusKind[iStatus].Data()));
          grt_L5[iStatus][iStave]->Draw("P");
        }
      }
      if (iLayer == 6) {
        for (int iStave = 0; iStave < nstave[6]; iStave++) {
          grt_L6[iStatus][iStave]->SetName(
              Form("L%d_St%d_F%s", iLayer, iStave, StatusKind[iStatus].Data()));
          grt_L6[iStatus][iStave]->Draw("P");
        }
      }
      legLayer[iLayer]->Draw();
      if (ccdb_upload) {
        string Runperiod =
            Form("%s", filepath.substr(filepath.find("from"), 27).c_str());
        ctrend2[iStatus][iLayer]->SetName(
            Form("Layer%d_Status_%s", iLayer, StatusKind[iStatus].Data()));
        auto mo = std::make_shared<o2::quality_control::core::MonitorObject>(
            ctrend2[iStatus][iLayer],
            TaskName + Form("/Status_%s", StatusKind[iStatus].Data()),
            TaskClass, DetectorName, 1, Runperiod);
        mo->setIsOwner(false);
        ccdb->storeMO(mo);
      }
      ctrend2[iStatus][iLayer]->SaveAs(
          Form("../Plots/LaneStatusFlag_%s.pdf",
               filepath
                   .substr(filepath.find("from"),
                           filepath.find(".root") - filepath.find("from"))
                   .c_str()));
      cout << "Just added to pdf file: Status:" << StatusKind[iStatus]
           << "  Layer:" << iLayer << endl;
      if (iStatus == NStatus - 1 && iLayer == NLayer - 1)
        ctrend2[iStatus][iLayer]->SaveAs(
            Form("../Plots/LaneStatusFlag_%s.pdf",
                 filepath
                     .substr(filepath.find("from"),
                             filepath.find(".root") - filepath.find("from"))
                     .c_str()));
    } // iLayer
  }   // iStatus




     // adding to the file trigger plot for suspicious runs
    for (int ir=0; ir<NRun; ir++)
      for (int il=0;il<7;il++){
	if (nOrbitDispersion[ir][il] > 0.10){
	  TCanvas *ctemp = new TCanvas;
	  int irev = NRun - ir -1;
	  htrg[irev]->SetTitle(Form("Trigger flags for sanity check of averaged time in run %s",runnumbersErr[irev].c_str()));
	  htrg[irev]->Draw("colz");
	  ctemp->SaveAs(
          Form("../Plots/LaneStatusFlag_%s.pdf",
               filepath
                   .substr(filepath.find("from"),
                           filepath.find(".root") - filepath.find("from"))
                   .c_str()));
	  delete ctemp;
	  break;
	}
      }


  TH1F *hblank2[NLayer];
  TLatex latex;
  TLatex latex2;
  latex.SetTextSize(0.024);
  latex2.SetTextSize(0.021);
  TLine *Avg[NLayer][(int)runnumbersErr.size()];
  for (int iLayer = 0; iLayer < NLayer; iLayer++) {
    hblank2[iLayer] = new TH1F(Form("hblank2_%d", iLayer),
                               Form("Layer %d - Percentage of orbits into "
                                    "ERROR or FAULT; Run; Time ratio",
                                    iLayer),
                               NRun, -0.5, (double)NRun - 0.5);
    for (int ir = 0; ir < (int)runnumbersErr.size(); ir++)
      hblank2[iLayer]->GetXaxis()->SetBinLabel(ir + 1,Form("%06d", stoi(runnumbersErr[runnumbersErr.size() - 1 - ir]))); // runnumbersErr is a descending order
     
    hblank2[iLayer]->GetXaxis()->LabelsOption("v");
    hblank2[iLayer]->GetYaxis()->SetRangeUser(1.e-7, 6);
    hblank2[iLayer]->GetXaxis()->SetTitleOffset(2.8);
    hblank2[iLayer]->GetXaxis()->SetLabelOffset(0.01);
    ctrend3[iLayer] = new TCanvas();
    ctrend3[iLayer]->cd();
    ctrend3[iLayer]->SetTickx();
    ctrend3[iLayer]->SetTicky();
    ctrend3[iLayer]->SetMargin(0.0988, 0.1, 0.194, 0.0993);
    hblank2[iLayer]->Draw();
    latex2.DrawLatex(-0.4,6.5,"Normalization uncertainty (%):");
    int yoffset=0;
    for (int ir = 0; ir < (int)runnumbersErr.size(); ir++){
      latex.DrawLatex(ir-0.95, 2.5-yoffset, Form("%4.0f",100.*nOrbitDispersion[ir][iLayer]));
      yoffset = 1 - yoffset;
    }
    gPad->SetLogy(); // for total number of errors
    gPad->SetGridx();
    gPad->SetGridy();

    // drawing averge over Layer
    
    for (int ir = 0; ir < (int)runnumbersErr.size(); ir++){
      double average=0;
      if (iLayer ==0) for (int ist = 0; ist<nstave[iLayer]; ist++) average += time_L0[ist]->GetPointY(ir) / nstave[iLayer];
      if (iLayer ==1) for (int ist = 0; ist<nstave[iLayer]; ist++) average += time_L1[ist]->GetPointY(ir) / nstave[iLayer];
      if (iLayer ==2) for (int ist = 0; ist<nstave[iLayer]; ist++) average += time_L2[ist]->GetPointY(ir) / nstave[iLayer];
      if (iLayer ==3) for (int ist = 0; ist<nstave[iLayer]; ist++) average += time_L3[ist]->GetPointY(ir) / nstave[iLayer];
      if (iLayer ==4) for (int ist = 0; ist<nstave[iLayer]; ist++) average += time_L4[ist]->GetPointY(ir) / nstave[iLayer];
      if (iLayer ==5) for (int ist = 0; ist<nstave[iLayer]; ist++) average += time_L5[ist]->GetPointY(ir) / nstave[iLayer];
      if (iLayer ==6) for (int ist = 0; ist<nstave[iLayer]; ist++) average += time_L6[ist]->GetPointY(ir) / nstave[iLayer];
      Avg[iLayer][ir] = new TLine(ir-0.5,average,ir+0.5,average);
      Avg[iLayer][ir]->SetLineColor(2);
      Avg[iLayer][ir]->SetLineWidth(2);
      Avg[iLayer][ir]->Draw("same");
    }
    

    if (iLayer == 0) {
      for (int iStave = 0; iStave < nstave[0]; iStave++) {
        time_L0[iStave]->SetName(Form("L%d_St%d", iLayer, iStave));
        time_L0[iStave]->Draw("P");
      }
    }
    if (iLayer == 1) {
      for (int iStave = 0; iStave < nstave[1]; iStave++) {
        time_L1[iStave]->SetName(Form("L%d_St%d", iLayer, iStave));
        time_L1[iStave]->Draw("P");
      }
    }
    if (iLayer == 2) {
      for (int iStave = 0; iStave < nstave[2]; iStave++) {
        time_L2[iStave]->SetName(Form("L%d_St%d", iLayer, iStave));
        time_L2[iStave]->Draw("P");
      }
    }
    if (iLayer == 3) {
      for (int iStave = 0; iStave < nstave[3]; iStave++) {
        time_L3[iStave]->SetName(Form("L%d_St%d", iLayer, iStave));
        time_L3[iStave]->Draw("P");
      }
    }
    if (iLayer == 4) {
      for (int iStave = 0; iStave < nstave[4]; iStave++) {
        time_L4[iStave]->SetName(Form("L%d_St%d", iLayer, iStave));
        time_L4[iStave]->Draw("P");
      }
    }
    if (iLayer == 5) {
      for (int iStave = 0; iStave < nstave[5]; iStave++) {
        time_L5[iStave]->SetName(Form("L%d_St%d", iLayer, iStave));
        time_L5[iStave]->Draw("P");
      }
    }
    if (iLayer == 6) {
      for (int iStave = 0; iStave < nstave[6]; iStave++) {
        time_L6[iStave]->SetName(Form("L%d_St%d", iLayer, iStave));
        time_L6[iStave]->Draw("P");
      }
    }
    legLayer[iLayer]->Draw();
    // if(ccdb_upload){ // TODO: CHECK THIS
    // string Runperiod =
    // Form("%s",filepath.substr(filepath.find("from"),27).c_str());
    // ctrend3[iLayer]->SetName(Form("Layer%d",iLayer));
    // auto mo2=
    // std::make_shared<o2::quality_control::core::MonitorObject>(ctrend2[iStatus][iLayer],
    // TaskName+Form("/time"), TaskClass, DetectorName,1,Runperiod);
    // mo->setIsOwner(false);
    // ccdb->storeMO(mo);
    // }


     
  
    ctrend3[iLayer]->SaveAs(
        Form("../Plots/LaneStatusFlag_%s.pdf",
             filepath
                 .substr(filepath.find("from"),
                         filepath.find(".root") - filepath.find("from"))
                 .c_str()));
    if (iLayer == NLayer - 1)
      ctrend3[iLayer]->SaveAs(
          Form("../Plots/LaneStatusFlag_%s.pdf]",
               filepath
                   .substr(filepath.find("from"),
                           filepath.find(".root") - filepath.find("from"))
                   .c_str()));

  
	 

	  
  } // iLayer

  ccdb->disconnect();
}

void SetCommonAxis(TH2D *hSummary) {
  // hSummary1[ilay]->GetXaxis()->SetNdivisions(530);
  // hSummary1[ilay]->GetYaxis()->SetNdivisions(516);
  hSummary->GetXaxis()->SetLabelSize(0.045);
  hSummary->GetYaxis()->SetLabelSize(0.045);
  hSummary->GetYaxis()->SetTitleOffset(1.3);
  hSummary->GetZaxis()->SetLabelSize(0.045);
  hSummary->GetXaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleOffset(0.7);
  hSummary->GetZaxis()->SetTitleSize(0.05);
  hSummary->GetZaxis()->SetTitleOffset(0.9);
}
