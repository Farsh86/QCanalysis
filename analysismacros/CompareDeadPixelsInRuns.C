#include "inc/itsAnalysis.hh"
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/CcdbDatabase.h"
#include "inc/ccdb.h"


void SetStyle(TGraphErrors *ge, Color_t col){
  ge->SetMarkerStyle(0);
  ge->SetMarkerColor(col);
  ge->SetFillColor(col);
  ge->SetLineColor(col);
}

std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2){
  std::array<long int,5> noisypix = {0, 0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int iybin=1; iybin<=512; iybin++){
    href->GetAxis(1)->SetRange(iybin,iybin); //select a row
    h2->GetAxis(1)->SetRange(iybin,iybin);//select a row
    TH1D *hrefproj = (TH1D*)href->Projection(0);
    TH1D *h2proj = (TH1D*)h2->Projection(0);
    for(int ixbin=1; ixbin<=9216; ixbin++){
      if(hrefproj->GetBinContent(ixbin)==1 && h2proj->GetBinContent(ixbin)==1){//dead in both runs
        noisypix[2]++;
      }
      else if(hrefproj->GetBinContent(ixbin)==1 && h2proj->GetBinContent(ixbin)==0){//dead only in ref run
        noisypix[0]++;
      }
      else if(hrefproj->GetBinContent(ixbin)==0 && h2proj->GetBinContent(ixbin)==1){//dead only in second run
        noisypix[1]++;
      }
      else continue;
    }
    delete hrefproj;
    delete h2proj;
  }
  return noisypix;
}

Int_t col[] = {810, 807, 797, 827, 417, 841, 868, 867, 860, 602, 921, 874};

// Main function
void CompareDeadPixelsInRuns(){

bool ccdb_upload;
string CCDB_up;

cout<<"Would you like to upload the output to ccdb? [y/n] ";
  cin>>CCDB_up;
  cout<<endl;
 if(CCDB_up =="y"||CCDB_up =="Y") ccdb_upload= true;
  else ccdb_upload= false;

if(ccdb_upload)SetTaskName(__func__);

std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");

auto ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());

  ccdb->connect(ccdbport.c_str(), "", "", "");
 

// itsAnalysis myAnalysis("Dead pixel Hits");
  itsAnalysis myAnalysis("Hits on Layer");

  auto laynums      = myAnalysis.Layers();      //vec of layers
  auto runNumbers   = myAnalysis.Runs();        //vec of run numbers
  auto hmaps        = myAnalysis.loadedHists(); // all histograms for layers and runs needed
  auto nRuns        = myAnalysis.nRuns();

  cout<<"Runs that will be used for comparisons: "<<endl;
  for (auto x: runNumbers){
    cout<<x<<endl;
  }
  string refrun;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>refrun;

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  long int first[7][100], second[7][100], both[7][100];
  bool filled[7][100];
  for (string layer : laynums){ // loop over layers
    int ilay=stoi(layer);
    for(int i=0; i<100; i++){
      first[ilay][i]=0; second[ilay][i]=0; both[ilay][i]=0;
      filled[ilay][i] = false;
    }
  }

  long int refHist[7][100];
  vector<array<long int,5>> noisypix;
  for (string layer : laynums){ // loop over layers
    auto hist = myAnalysis.loadLayerSparse(stoi(layer));
    for(long int i = 0; i<hist.size();i++ ){ //loop hist for layer
      string hname = hist[i]->GetName();
      string stvnum = hname.substr(hname.find("Stv")+3,2);
      if(stvnum.find("_")!=string::npos){
        stvnum = hname.substr(hname.find("Stv")+3,1);
      }
      if(myAnalysis.getRunNumber(hist[i])==refrun){
        refHist[stoi(layer)][stoi(stvnum)] = i;
      }
    }
    for(int i = 0; i<hist.size();i++ ){ //loop hist for layer
      if(hist[i]->GetEntries()>1e4){
        cout<<hist[i]->GetName()<<" skipped because has more than 10000 entries (probably a bad run)."<<endl;
        continue;
      }
      string hname = hist[i]->GetName();
      string stvnum = hname.substr(hname.find("Stv")+3,2);
      if(stvnum.find("_")!=string::npos){
        stvnum = hname.substr(hname.find("Stv")+3,1);
      }
      auto runn = myAnalysis.getRunNumber(hist[i]);
      auto itf = find(runNumbers.begin(), runNumbers.end(), runn);
      int irun = distance(runNumbers.begin(), itf);
      if(runn==refrun) continue;
      cout<<"Layer: "<<layer<<" stvnum: "<<stvnum<<" Run Number: "<<runNumbers[irun]<<endl;
      noisypix.push_back(CompareTwoRuns(hist[refHist[stoi(layer)][stoi(stvnum)]], hist[i]));
      noisypix[noisypix.size()-1][3] = -9999; // Not needed for now
      noisypix[noisypix.size()-1][4] = stol(layer);
      first[stoi(layer)][irun]+=noisypix[noisypix.size()-1][0];
      second[stoi(layer)][irun]+=noisypix[noisypix.size()-1][1];
      both[stoi(layer)][irun]+=noisypix[noisypix.size()-1][2];
      filled[stoi(layer)][irun] = true;
    }
  }

  //Make plot for each layer and for each stave in the root file
  TGraphErrors *ge_nref[7];
  TGraphErrors *ge_n2[7];
  TGraphErrors *ge_ncom1[7];
  TGraphErrors *ge_ncom2[7];

  double xshift = 3.;
  double max[7];
  double min[7];

  for (string layer : laynums){ // loop over layers
    int ilay=stoi(layer);

    ge_nref[ilay] = new TGraphErrors();
    ge_n2[ilay] = new TGraphErrors();
    ge_ncom1[ilay] = new TGraphErrors();
    ge_ncom2[ilay] = new TGraphErrors();

    max[ilay] = -1.;
    min[ilay] = 1e35;

    int ipoint = 0;
    for(int ir=0; ir<100; ir++){//first the older data and last the most recent
      if(!filled[ilay][ir]) continue; //skip if not filled
      if(!ipoint) xshift=1.;
      else xshift = 3.;
      ge_nref[ilay]->SetPoint(ipoint, ipoint*xshift, (double)both[ilay][ir]/2.+(double)first[ilay][ir]/2.);
      ge_nref[ilay]->SetPointError(ipoint, 0.5, (double)first[ilay][ir]/2.);
      ge_ncom1[ilay]->SetPoint(ipoint, ipoint*xshift, 0.);
      ge_ncom1[ilay]->SetPointError(ipoint, 0.5, (double)both[ilay][ir]/2.);
      if((double)both[ilay][ir]/2.+(double)first[ilay][ir] > max[ilay]) max[ilay] = (double)both[ilay][ir]/2.+(double)first[ilay][ir];

      //second couple of bar on the right
      ge_n2[ilay]->SetPoint(ipoint, ipoint*xshift+1, -(double)both[ilay][ir]/2.-(double)second[ilay][ir]/2.);
      ge_n2[ilay]->SetPointError(ipoint, 0.5, (double)second[ilay][ir]/2.);
      ge_ncom2[ilay]->SetPoint(ipoint, ipoint*xshift+1, 0.);
      ge_ncom2[ilay]->SetPointError(ipoint, 0.5, (double)both[ilay][ir]/2.);
      if(-(double)both[ilay][ir]/2.-(double)second[ilay][ir] < min[ilay]) min[ilay] = -(double)both[ilay][ir]/2.-(double)second[ilay][ir];
      ipoint++;
    }//end first loop on runs
    //Style
    SetStyle(ge_nref[ilay], kBlue);
    SetStyle(ge_ncom1[ilay], kBlack);
    SetStyle(ge_ncom2[ilay], kBlack);
    SetStyle(ge_n2[ilay], kRed+2);
  }//end loop on layers

  //Draw plot for each layer
  for (string layer : laynums){ // loop over layers
    int ilay=stoi(layer);

    TCanvas *canvas = new TCanvas(Form("mycanvas_%d",ilay), Form("mycanvas_%d",ilay), 1300, 800);
    canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
    canvas->cd();
    double x2,y2;
    ge_ncom2[ilay]->GetPoint(ge_ncom2[ilay]->GetN()-1, x2,y2);
    TH1F *hfake = new TH1F("hfake","hfake", (int)3*runNumbers.size()+1, -3, 3*runNumbers.size()-2);
    hfake->SetStats(0);

    ///draw labels on x axis
    int counter = 0;
    for(Int_t k=4;k<3*runNumbers.size();k+=3){
      if(runNumbers[counter]==refrun){
        k-=3;
        counter++;
        continue;
      }
      hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runNumbers[runNumbers.size()-1-counter].c_str()));
      counter++;
    }
    
    hfake->Draw();
    hfake->SetTitle(Form("Layer-%d - run%s compared to all",ilay, refrun.c_str()));
    ge_nref[ilay]->Draw("P E2 same");
    ge_ncom1[ilay]->Draw("E2 same");
    ge_ncom2[ilay]->Draw("E2 same");
    ge_n2[ilay]->Draw("E2 same");
    hfake->GetYaxis()->SetRangeUser(min[ilay]+0.1*min[ilay], max[ilay]+0.1*max[ilay]);
    hfake->GetYaxis()->SetTickLength(0.005);
    hfake->GetYaxis()->SetMaxDigits(4);
    TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
    lineref->SetLineColor(kGray-1);
    lineref->SetLineStyle(2);
    lineref->Draw("same");

    //Legend
    TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->AddEntry(ge_nref[0], "#splitline{#dead pix}{ref. run only}", "f");
    leg->AddEntry(ge_n2[0], "#splitline{#dead pix}{2nd run only}", "f");
    leg->AddEntry(ge_ncom1[0], "#splitline{#dead pix}{both}");
    leg->Draw("same");
if(ccdb_upload){
     string Runperiod = Form("from_run%s_to_run%s",runNumbers[0].c_str(),runNumbers[nRuns-1].c_str());
//   string Runperiod = Form("%s",filepath.substr(filepath.find("from"),27).c_str()); //This should be used for actual data 
    canvas->SetName(Form("Layer%s_Dead_pixels_Correlation",layer.c_str()));
    auto mo1= std::make_shared<o2::quality_control::core::MonitorObject>(canvas, TaskName+Form("/Layer%s",layer.c_str()), TaskClass, DetectorName,std::stoi(runNumbers.back()),Runperiod);
    mo1->addMetadata("ReferenceRunNumber",refrun.c_str());
    mo1->setIsOwner(false);
    ccdb->storeMO(mo1);
               }
    canvas->SaveAs(Form("../Plots/Layer%s_DeadPixComparison_run%s_compared_to_run%s_run%s.pdf", layer.c_str(),refrun.c_str(),runNumbers[0].c_str(),runNumbers[nRuns-1].c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_DeadPixComparison_run%s_compared_to_run%s_run%s.root", layer.c_str(),refrun.c_str(),runNumbers[0].c_str(),runNumbers[nRuns-1].c_str()));

    delete canvas;
    delete hfake;
    delete lineref;
  }//end loop on layers

//Disconnencting from the database
  ccdb->disconnect();

} // end of analyseLayerThresholds()
