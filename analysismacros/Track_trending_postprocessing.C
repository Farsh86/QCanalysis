#include <string>
#include <iostream>
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
#include <TGraphErrors.h>
#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/Reductor.h"
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/RootClassFactory.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/CcdbDatabase.h"
#include "inc/ccdb.h"

void Track_trending_postprocessing(){

string fpath;
cout<<"Available file(s) for the analysis of track task:"<<endl;
gSystem->Exec("ls ../Data/*TrackTask* -Art | tail -n 500");
cout<<"\nCopy the file name to analyse:"<<endl;
cin>>fpath;
cout<<endl;


TFile *file = TFile::Open(fpath.c_str());

   TList * keyslist=(TList *) file->GetListOfKeys();
   TIter next(keyslist);
   TKey *key;
   vector<string> runs, runs1, runs2, runs3, runs4, runsZ, runs5, runs_eta, runs_phi;
   vector<double>phi_ave, phicounts, phicounts1, phicounts2, relative_phi1, relative_phi2;
   vector<double>eta_ave, etacounts, relative_eta1, relative_eta2;
   vector<double>cycle,cycle1,cycle2, cycle3, cycle4, cycle5, cycle_eta, cycle_phi;
   vector<double>binminCenters, binminCenters1, binminCenters2, binminCenters3, binminCenters4, binminCenters5;
   vector<double>binmaxCenters, binmaxCenters1, binmaxCenters2, binmaxCenters3, binmaxCenters4, binmaxCenters5;
   vector<double>means, means1, means2, means3, means4, means5;
   vector<double>aves, aves_err, aves_nozero, aves_nozero_err, aves1, aves1_err, aves1_nozero, aves1_nozero_err, aves2, aves2_err, aves3, aves3_err;
   vector<double>rmss, rmss_err, rmss_nozero, rmss_nozero_err, rmss1, rmss1_err, rmss1_nozero, rmss1_nozero_err, rmss2, rmss2_err, rmss3, rmss3_err;
   vector<double>xpoints, xpoints_err, ypoints, ypoints_err, rmsx_vertex, rmsy_vertex, rmsx_vertex_err, rmsy_vertex_err;
   vector<TH2*>histos, histos1;
   vector<double>bx1,by1,xmaxs1,xmins1,ymaxs1,ymins1;
   double a,b,d,e,f,z,j,k=0;
   bool ccdb_upload;

   cout<<"List of runs:"<<endl;
   for(auto&& keyAsObj : *file->GetListOfKeys()){
       auto keysub = (TKey*) keyAsObj;
       string name= keysub->GetName();
       if(name.find("NClustersPerTrackEta")!=string::npos){
        string run = name.substr(name.find("run")+3, 6);
        cout<<run<<endl;
       }
    }

   string skipans, skipruns, CCDB_up;
   cout<<"Would you like to skip some run(s)? [y/n]"<<endl;
   cin>>skipans;
   if(skipans=="y" || skipans=="Y"){
     cout<<endl;
     cout<<"Specify run number(s) to skip separated by comma (no white spaces!):"<<endl;
     cin>>skipruns;
     cout<<endl;
   }
   else
     skipruns=" ";

 cout<<"Would you like to upload the output to ccdb? [y/n] ";
  cin>>CCDB_up;
  cout<<endl;
 if(CCDB_up =="y"||CCDB_up =="Y") ccdb_upload= true;
  else ccdb_upload= false;

if(ccdb_upload)SetTaskName(__func__);

string vtxRun3;
bool isVtxRun3 = false;
cout<<"Is this the vertex trend for full Run 3? [y/n] ";
cin>>vtxRun3;
cout<<endl;
if(vtxRun3 == "y" || vtxRun3 == "Y") isVtxRun3 = true;

//Setting up the connection to the ccdb database
std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");

auto ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());

  ccdb->connect(ccdbport.c_str(), "", "", "");

   while ((key=(TKey*)next())){

      string keyname = (string)key->GetName();
      string runnum = keyname.substr(keyname.find("run")+3, 6);
      if(skipruns.find(runnum)!=string::npos) continue;

//////////////////////////////////////////////////////////////////////////////////////////////

       if (keyname.find("EtaDistribution")!=string::npos){
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs_eta.push_back(runnum);
          cycle_eta.push_back(j);
          j++;

          TH1D *h= (TH1D *) key->ReadObj();

          //Integral in all range
          double integral = h->Integral();

          //Extract the mean value in x axis
          double x_mean = h->GetMean(1);
          eta_ave.push_back(x_mean);


          //Find the bin number for a value in x axis
          TAxis* xaxis = h->GetXaxis();
          int binx = xaxis->FindBin(x_mean);

          //Find the value in y axis
          double bin = h->GetBinContent(binx);
          //std::cout<<"etaCounts_y_value= "<< bin << std::endl;
          etacounts.push_back(bin);

          TH1D *hc1 = (TH1D *)h->Clone();
          hc1->GetXaxis()->SetRangeUser(-1.5, 0);

          //Integral in range (-pi/2)-(0)
          double integral1 = hc1->Integral();

          TH1D *hc2 = (TH1D *)h->Clone();
          hc2->GetXaxis()->SetRangeUser(0, 1.5);

          //Integral in range 0-pi/2
          double integral2 = hc2->Integral();

          //Relative values
          double rel_eta1;
          if (integral != 0){
            rel_eta1 = (integral1/integral);
          }
          else {rel_eta1 = 0;}

          relative_eta1.push_back(rel_eta1);

          double rel_eta2;
          if (integral != 0){
            rel_eta2 = (integral2/integral);
          }
          else {rel_eta2 = 0;}

          relative_eta2.push_back(rel_eta2);

     }

//////////////////////////////////////////////////////////////////////////////////////////////

       if (keyname.find("PhiDistribution")!=string::npos){
          //key->Print();
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs_phi.push_back(runnum);
          cycle_phi.push_back(k);
          k++;

          TH1D *h= (TH1D *) key->ReadObj();

          //Integral in range 0-2pi
          double integral = h->Integral();

          //Extract the mean value in x axis
          double x_mean = h->GetMean(1);
          phi_ave.push_back(x_mean);
          //Find the bin number for a value in x axis
          TAxis* xaxis = h->GetXaxis();
          int binx = xaxis->FindBin(x_mean);
          //Find the value in y axis
          double bin = h->GetBinContent(binx);
          //cout<<"phi_y_value= "<<bin<<endl;
          phicounts.push_back(bin);

          TH1D *hc1 = (TH1D *)h->Clone();
          hc1->GetXaxis()->SetRangeUser(0, 3);

          //Integral in range 0-pi
          double integral1 = hc1->Integral();

          //Extract the mean value in x axis
          double x_mean1 = hc1->GetMean(1);
          //Find the bin number for a value in x axis
          TAxis* xaxis1 = hc1->GetXaxis();
          int binx1 = xaxis1->FindBin(x_mean1);
          //Find the value in y axis
          double bin1 = hc1->GetBinContent(binx1);
          //cout<<"y value clone 1= "<<bin1<<endl;
          phicounts1.push_back(bin1);

          TH1D *hc2 = (TH1D *)h->Clone();
          hc2->GetXaxis()->SetRangeUser(3, 6);

          //Integral in range pi-2pi
          double integral2 = hc2->Integral();

          //Extract the mean value in x axis
          double x_mean2 = hc2->GetMean(1);
          //Find the bin number for a value in x axis
          TAxis* xaxis2 = hc2->GetXaxis();
          int binx2 = xaxis2->FindBin(x_mean2);
          //Find the value in y axis
          double bin2 = hc2->GetBinContent(binx2);
          //cout<<"y value clone 2= "<<bin2<<endl;
          phicounts2.push_back(bin2);

          //Relative values
          double rel_phi1;
          if (integral != 0){
            rel_phi1 = (integral1/integral);
          }
          else {rel_phi1 = 0;}

          relative_phi1.push_back(rel_phi1);

          double rel_phi2;
          if (integral != 0){
            rel_phi2 = (integral2/integral);
          }
          else {rel_phi2 = 0;}

          relative_phi2.push_back(rel_phi2);
      }

//////////////////////////////////////////////////////////////////////////////////
    if (keyname.find("VertexZ")!=string::npos){

         cout<<keyname<<endl;
         cout<<"Run number: "<<runnum<<endl;

         runs5.push_back(runnum);
         cycle5.push_back(z);
         z++;

         TH1D *hz1= (TH1D *) key->ReadObj();

         //Mean value in x axis
         double ave3 = hz1->GetMean(1);
         aves3.push_back(ave3);
         aves3_err.push_back(hz1->GetMeanError());

         //Extract RMS of the distribution
         double rms3 = hz1->GetRMS();
         rmss3.push_back(rms3);
         rmss3_err.push_back(hz1->GetRMSError());
      }

//////////////////////////////////////////////////////////////////////////////////
    if (keyname.find("VertexRvsZ")!=string::npos){
        cout<<keyname<<endl;
        cout<<"Run number: "<<runnum<<endl;

        runsZ.push_back(runnum);

        TH2D *hz= (TH2D *) key->ReadObj();
        histos1.push_back(hz);

        int binsx1 = hz->GetNbinsX();
        bx1.push_back(binsx1);
        int binsy1 = hz->GetNbinsY();
        by1.push_back(binsy1);
        double xmax1=hz->GetXaxis()->GetXmax();
        xmaxs1.push_back(xmax1);
        double xmin1=hz->GetXaxis()->GetXmin();
        xmins1.push_back(xmin1);
        double ymax1 = hz->GetYaxis()->GetBinLowEdge(hz->GetNbinsY()) + hz->GetYaxis()->GetBinWidth(hz->GetNbinsY());
        ymaxs1.push_back(ymax1);
        double ymin1 = hz->GetYaxis()->GetBinLowEdge(1);
        ymins1.push_back(ymin1);
     }
////////////////////////////////////////////////////////////////////////////////////
      if (keyname.find("VertexCoordinates")!=string::npos){
         cout<<keyname<<endl;
         cout<<"Run number: "<<runnum<<endl;

         runs4.push_back(runnum);
         cycle4.push_back(f);
         f++;

         TH2D *hf= (TH2D *) key->ReadObj();
         histos.push_back(hf);
         //x
         xpoints.push_back(hf->GetMean(1));
         xpoints_err.push_back(hf->GetMeanError(1));
         rmsx_vertex.push_back(hf->GetRMS(1));
         rmsx_vertex_err.push_back(hf->GetRMSError(1));
         //y
         ypoints.push_back(hf->GetMean(2));
         ypoints_err.push_back(hf->GetMeanError(2));
         rmsy_vertex.push_back(hf->GetRMS(2));
         rmsy_vertex_err.push_back(hf->GetRMSError(2));

      }

//////////////////////////////////////////////////////////////////////////////////
       if (keyname.find("NVertexContributors")!=string::npos){
         cout<<keyname<<endl;
         cout<<"Run number: "<<runnum<<endl;

         runs3.push_back(runnum);
         cycle3.push_back(e);
         e++;

         TH1D *he= (TH1D *) key->ReadObj();

         //Mean value in x axis
         aves2.push_back(he->GetMean());
         aves2_err.push_back(he->GetMeanError());

         //Extract RMS of the distribution
         rmss2.push_back(he->GetRMS());
         rmss2_err.push_back(he->GetRMSError());
      }

//////////////////////////////////////////////////////////////////////////////////
       if (keyname.find("Ntracks")!=string::npos){
         cout<<keyname<<endl;
         cout<<"Run number: "<<runnum<<endl;

         runs2.push_back(runnum);
         cycle2.push_back(d);
         d++;

         TH1D *ht= (TH1D *) key->ReadObj();

         //Mean value in x axis
         aves1.push_back(ht->GetMean(1));
         aves1_err.push_back(ht->GetMeanError(1));
         //Extract RMS of the distribution
         rmss1.push_back(ht->GetRMS(1));
         rmss1_err.push_back(ht->GetRMSError(1));
         //Mean value removing first bin with x = 0
         ht->SetBinContent(1,0);
         aves1_nozero.push_back(ht->GetMean(1));
         aves1_nozero_err.push_back(ht->GetMeanError(1));
         rmss1_nozero.push_back(ht->GetRMS());
         rmss1_nozero_err.push_back(ht->GetRMSError(1));

      }
//////////////////////////////////////////////////////////////////////////////////
      if (keyname.find("AssociatedClusterFraction")!=string::npos){
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs1.push_back(runnum);
          cycle1.push_back(b);
          b++;

          TH1D *hh= (TH1D *) key->ReadObj();

          //Mean value in x axis
          aves.push_back(hh->GetMean(1));
          aves_err.push_back(hh->GetMeanError(1));
          rmss.push_back(hh->GetRMS());
          rmss_err.push_back(hh->GetRMSError());

          //remove entries at x = 0
          hh->SetBinContent(1,0);

          aves_nozero.push_back(hh->GetMean());
          aves_nozero_err.push_back(hh->GetMeanError());
          rmss_nozero.push_back(hh->GetRMS());
          rmss_nozero_err.push_back(hh->GetRMSError());
      }

//////////////////////////////////////////////////////////////////////////////////
      if (keyname.find("NClustersPerTrackEta")!=string::npos){
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs.push_back(runnum);
          cycle.push_back(a);
          a++;

          TH2D *h= (TH2D *) key->ReadObj();

//////////////Projection for -1.5 < eta < 1.5////////////////////////////
          TH1D *hproj = h->ProjectionY();

          int nbinsx=hproj->GetNbinsX();
          //Average
          double x_mean = hproj->GetMean(1);
          means.push_back(x_mean);
          //Find the bin with minimum content != 0
          double min = 1e20;
          int minbin = 1;;
          for (int b=1; b<=hproj->GetNbinsX(); b++){
            double content = hproj->GetBinContent(b);
            if (content > 0 && content < min) {
              min = content;
              minbin = b;
            }
          }

          double binmincenter=hproj->GetBinCenter(minbin);
          binminCenters.push_back(x_mean<1e-15 ? 0 : binmincenter);

          //Find the bin with maximum content
          int max=hproj->GetMaximumBin();
          double binmaxcenter=hproj->GetBinCenter(max);
          binmaxCenters.push_back(x_mean<1e-15 ? 0 : binmaxcenter);

////////////////Projection for -1.2 < eta < 1.2////////////////////////
          TH2D *hclon1 = (TH2D *)h->Clone();
          hclon1->GetXaxis()->SetRangeUser(-1.2, 1.2);
          TH1D *hproj1 = hclon1->ProjectionY();

          int nbinsx1=hproj1->GetNbinsX();

          //Average
          double x_mean1 = hproj1->GetMean(1);
          means1.push_back(x_mean1);

          //Find the bin with minimum content != 0
          int minbin1 = 1;
          min = 1e20;
          for (int b=1; b<=hproj1->GetNbinsX(); b++){
            double content = hproj1->GetBinContent(b);
            if (content > 0 && content < min) {
              min = content;
              minbin1 = b;
            }
          }
          double binmincenter1=hproj1->GetBinCenter(minbin1);
          binminCenters1.push_back(x_mean1<1e-15 ? 0 : binmincenter1);

          //Find the bin with maximum content
          int max1=hproj1->GetMaximumBin();
          double binmaxcenter1=hproj1->GetBinCenter(max1);
          binmaxCenters1.push_back(x_mean1<1e-15 ? 0 : binmaxcenter1);

///////////////Projection for -0.8 < eta < 0.8////////////////////////
          TH2D *hclon2 = (TH2D *)h->Clone();
          hclon2->GetXaxis()->SetRangeUser(-0.8, 0.8);
          TH1D *hproj2 = hclon2->ProjectionY();

          int nbinsx2=hproj2->GetNbinsX();

          //Average
          double x_mean2 = hproj2->GetMean(1);
          means2.push_back(x_mean2);

          //Find the bin with minimum content != 0
          int minbin2 = 1;
          min = 1e20;
          for (int b=1; b<=hproj2->GetNbinsX(); b++){
            double content = hproj2->GetBinContent(b);
            if (content > 0 && content < min) {
              min = content;
              minbin2 = b;
            }
          }
          double binmincenter2=hproj2->GetBinCenter(minbin2);
          binminCenters2.push_back(x_mean2<1e-15 ? 0 : binmincenter2);

          //Find the bin with maximum content
          int max2=hproj2->GetMaximumBin();
          double binmaxcenter2=hproj2->GetBinCenter(max2);
          binmaxCenters2.push_back(x_mean2<1e-15 ? 0 : binmaxcenter2);

///////////////Projection for -1.5 < eta < -1.4///////////////////////////////
          TH2D *hclon3 = (TH2D *)h->Clone();
          hclon3->GetXaxis()->SetRangeUser(-1.5, -1.4);
          TH1D *hproj3 = hclon3->ProjectionY();

          int nbinsx3=hproj3->GetNbinsX();

          //Average
          double x_mean3 = hproj3->GetMean(1);
          means3.push_back(x_mean3);

          //Find the bin with minimum content != 0
          int minbin3 = 1;
          min = 1e20;
          for (int b=1; b<=hproj3->GetNbinsX(); b++){
            double content = hproj3->GetBinContent(b);
            if (content > 0 && content < min){
              min = content;
              minbin3 = b;
            }
          }
          double binmincenter3=hproj3->GetBinCenter(minbin3);
          binminCenters3.push_back(x_mean3<1e-15 ? 0 : binmincenter3);

          //Find the bin with maximum content
          int max3=hproj3->GetMaximumBin();
          double binmaxcenter3=hproj3->GetBinCenter(max3);
          binmaxCenters3.push_back(x_mean3<1e-15 ? 0 : binmaxcenter3);

////////////////Projection for 1.4 < eta < 1.5///////////////////////////////
          TH2D *hclon4 = (TH2D *)h->Clone();
          hclon4->GetXaxis()->SetRangeUser(1.4, 1.5);
          TH1D *hproj4 = hclon4->ProjectionY();

          int nbinsx4=hproj4->GetNbinsX();

          //Average
          double x_mean4 = hproj4->GetMean(1);
          means4.push_back(x_mean4);

          //Find the bin with minimum content != 0
          int minbin4 = 1;
          min = 1e20;
          for (int b=1; b<=hproj4->GetNbinsX(); b++){
            double content = hproj4->GetBinContent(b);
            if (content > 0 && content < min){
              min = content;
              minbin4 = b;
            }
          }
          double binmincenter4=hproj4->GetBinCenter(minbin4);
          binminCenters4.push_back(x_mean4<1e-15 ? 0 : binmincenter4);

          //Find the bin with maximum content
          int max4=hproj4->GetMaximumBin();
          double binmaxcenter4=hproj4->GetBinCenter(max4);
          binmaxCenters4.push_back(x_mean4<1e-15 ? 0 : binmaxcenter4);

////////////Projection for -0.2 < eta < 0.2/////////////////////////////////////
          TH2D *hclon5 = (TH2D *)h->Clone();
          hclon5->GetXaxis()->SetRangeUser(-0.2, 0.2);
          TH1D *hproj5 = hclon5->ProjectionY();

          int nbinsx5=hproj5->GetNbinsX();

          //Average
         double x_mean5 = hproj5->GetMean(1);
         means5.push_back(x_mean5);

          //Find the bin with minimum content != 0
          int minbin5 = 1;
          min = 1e20;
          for (int b=1; b<=hproj5->GetNbinsX(); b++){
            double content = hproj5->GetBinContent(b);
            if (content > 0 && content < min) {
              min = content;
              minbin5 = b;
            }
          }
          double binmincenter5=hproj5->GetBinCenter(minbin5);
          binminCenters5.push_back(x_mean5<1e-15 ? 0 : binmincenter5);

          //Find the bin with maximum content
          int max5=hproj5->GetMaximumBin();
          double binmaxcenter5=hproj5->GetBinCenter(max5);
          binmaxCenters5.push_back(x_mean5<1e-15 ? 0 : binmaxcenter5);
     }

    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(runs_eta.begin(),runs_eta.end());
  std::reverse(eta_ave.begin(),eta_ave.end());

  auto c2 = new TCanvas("Average of Eta distributions");
  c2->SetGrid(0,1);
  auto gr2 = new TGraph(cycle_eta.size(),&cycle_eta[0],&eta_ave[0]);
  double b2=0;

  gr2->GetXaxis()->SetNdivisions(cycle_eta.size());
  TH1F hfake("hfake","hfake", gr2->GetN(),-0.5,(double)gr2->GetN()-0.5);
  hfake.SetStats(0);
  for (int i=1;i<=(int)cycle_eta.size();i++) hfake.GetXaxis()->SetBinLabel(i,runs_eta[i-1].data());

  hfake.GetXaxis()->SetLabelSize(0.027);
  hfake.SetNameTitle("Average of Eta distributions", "Average of Eta distributions");
  hfake.GetXaxis()->SetTitle("");
  hfake.GetYaxis()->SetTitle("Mean angle (rad)");
  hfake.GetYaxis()->SetRangeUser(-0.9,0.9);
  gr2->SetMarkerStyle(20);
  hfake.Draw();
  gr2->Draw("PL same");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(relative_eta1.begin(),relative_eta1.end());
  std::reverse(relative_eta2.begin(),relative_eta2.end());

  auto c21 = new TCanvas("Relative number of tracks in Eta distributions");
  c21->SetGrid(0,1);

  auto gr21 = new TGraph(cycle_eta.size(),&cycle_eta[0],&relative_eta1[0]);
  auto gr22 = new TGraph(cycle_eta.size(),&cycle_eta[0],&relative_eta2[0]);
  double b21=0;

  gr21->GetXaxis()->SetNdivisions(cycle_eta.size());
  TH1F hfake2("hfake2","hfake2", gr21->GetN(),-0.5,(double)gr21->GetN()-0.5);
  hfake2.SetStats(0);
  for (int i=1;i<=(int)cycle_eta.size();i++) hfake2.GetXaxis()->SetBinLabel(i,runs_eta[i-1].data());


  hfake2.GetYaxis()->SetRangeUser(0,1.);
  hfake2.GetXaxis()->SetLabelSize(0.027);
  hfake2.SetNameTitle("Relative number of tracks in Eta distributions", "Relative number of tracks in Eta distributions");
  hfake2.GetXaxis()->SetTitle("");
  hfake2.GetYaxis()->SetTitle("Number of counts in subrange / Total counts");
  gr21->SetMarkerStyle(20);
  gr21->SetMarkerColor(kRed);
  gr21->SetLineColor(kRed);
  hfake2.Draw();
  gr21->Draw("PL same");
  gr22->SetMarkerStyle(20);
  gr22->SetMarkerColor(kBlue);
  gr22->SetLineColor(kBlue);
  gr22->Draw("PL SAME");

  auto legend=new TLegend(0.773381,0.85429,0.905276,0.909763);
  legend->SetNColumns(2);
  legend->SetFillColor(0);
  legend->SetHeader("#eta range","");
  legend->AddEntry(gr21,"-1.5-0","lp");
  legend->AddEntry(gr22,"0-1.5","lp");
  legend->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(runs_phi.begin(),runs_phi.end());
  std::reverse(phicounts.begin(),phicounts.end());
  std::reverse(phi_ave.begin(),phi_ave.end());

  auto c3 = new TCanvas("Average of Phi distributions (0 < Phi < 2#pi)");
  c3->SetGrid(0,1);

  auto gr3 = new TGraph(cycle_phi.size(),&cycle_phi[0],&phi_ave[0]);
  double b3=0;

  TH1F hfake3("hfake3","hfake3", gr3->GetN(),-0.5,(double)gr3->GetN()-0.5);
  hfake3.SetStats(0);
  gr3->GetXaxis()->SetNdivisions(cycle_phi.size());

  for (int i=1;i<=(int)cycle_phi.size();i++) hfake3.GetXaxis()->SetBinLabel(i,runs_phi[i-1].data());

  hfake3.GetXaxis()->SetLabelSize(0.027);
  hfake3.SetNameTitle("Average of Phi distributions (0 < Phi < 2#pi)", "Average of Phi distributions (0 < Phi < 2#pi)");
  hfake3.GetXaxis()->SetTitle("");
  hfake3.GetYaxis()->SetTitle("Mean angle (rad)");
  hfake3.GetYaxis()->SetRangeUser(0, 6);
  gr3->SetMarkerStyle(20);
  hfake3.Draw();
  gr3->Draw("PL same");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(phicounts1.begin(),phicounts1.end());
  std::reverse(phicounts2.begin(),phicounts2.end());
  std::reverse(relative_phi1.begin(),relative_phi1.end());
  std::reverse(relative_phi2.begin(),relative_phi2.end());

  auto c31 = new TCanvas("Relative number of tracks in Phi distributions");
  c31->SetGrid(0,1);

  auto gr31 = new TGraph(cycle_phi.size(),&cycle_phi[0],&relative_phi1[0]);
  auto gr32 = new TGraph(cycle_phi.size(),&cycle_phi[0],&relative_phi2[0]);
  double b31=0;

  gr31->GetXaxis()->SetNdivisions(cycle_phi.size());
  TH1F hfake4("hfake4","hfake4", gr31->GetN(),-0.5,(double)gr31->GetN()-0.5);
  hfake4.SetStats(0);
  for (int i=1;i<=(int)cycle_phi.size();i++) hfake4.GetXaxis()->SetBinLabel(i,runs_phi[i-1].data());

  hfake4.GetYaxis()->SetRangeUser(0,1);
  hfake4.GetXaxis()->SetLabelSize(0.027);
  hfake4.SetNameTitle("Relative number of tracks in Phi distributions", "Relative number of tracks in Phi distributions");
  hfake4.GetXaxis()->SetTitle("");
  hfake4.GetYaxis()->SetTitle("Number of counts in subrange / Total counts");
  gr31->SetMarkerStyle(20);
  gr31->SetMarkerColor(kRed);
  gr31->SetLineColor(kRed);
  hfake4.Draw();
  gr31->Draw("PL same");
  gr32->SetMarkerStyle(20);
  gr32->SetMarkerColor(kBlue);
  gr32->SetLineColor(kBlue);
  gr32->Draw("PL SAME");

  auto legend1=new TLegend(0.773381,0.85429,0.905276,0.909763);
  legend1->SetNColumns(2);
  legend1->SetFillColor(0);
  legend1->SetHeader("#phi range","");
  legend1->AddEntry(gr31,"0-#pi","lp");
  legend1->AddEntry(gr32,"#pi-2#pi","lp");
  legend1->Draw();

////////////////////////////////VertexZ Avg//////////////////////////////////////////////////

std::reverse(runs5.begin(),runs5.end());
std::reverse(aves3.begin(),aves3.end());
std::reverse(aves3_err.begin(),aves3_err.end());
std::reverse(rmss3.begin(),rmss3.end());
std::reverse(rmss3_err.begin(),rmss3_err.end());
auto cz1 = new TCanvas();
cz1->SetGrid(0,1);
auto gerr3 = new TGraphErrors(cycle5.size(),&cycle5[0],&aves3[0],0,&aves3_err[0]);
gerr3->GetXaxis()->SetNdivisions(cycle5.size());
TH1F hfake5("hfake5","hfake5", gerr3->GetN(),-0.5,(double)gerr3->GetN()-0.5);
hfake5.SetStats(0);
for (int i=1;i<=(int)cycle5.size();i++) hfake5.GetXaxis()->SetBinLabel(i,runs5[i-1].data());
hfake5.GetXaxis()->SetLabelSize(0.027);
gerr3->SetMarkerStyle(20);
hfake5.GetYaxis()->SetRangeUser(-5,5);
hfake5.SetTitle("Mean Z coordinate of vertices");
hfake5.GetXaxis()->SetTitle("");
hfake5.GetYaxis()->SetTitle("<Z coordinates> (cm)");
hfake5.Draw();
gerr3->Draw("PL same");

////////////////////////////////VertexZ RMS//////////////////////////////////////////////////
auto cz11 = new TCanvas();
cz11->SetGrid(0,1);
auto gerr31 = new TGraphErrors(cycle5.size(),&cycle5[0],&rmss3[0],0,&rmss3_err[0]);
gerr31->GetXaxis()->SetNdivisions(cycle5.size());
TH1F hfake51("hfake51","hfake51", gerr31->GetN(),-0.5,(double)gerr31->GetN()-0.5);
hfake51.SetStats(0);
for (int i=1;i<=(int)cycle5.size();i++) hfake51.GetXaxis()->SetBinLabel(i,runs5[i-1].data());
hfake51.GetXaxis()->SetLabelSize(0.027);
gerr31->SetMarkerStyle(20);
hfake51.GetYaxis()->SetRangeUser(0,10);
hfake51.SetTitle("RMS of Z coordinate of vertices");
hfake51.GetXaxis()->SetTitle("");
hfake51.GetYaxis()->SetTitle("RMS (cm)");
hfake51.Draw();
gerr31->Draw("PL same");

///////////////////////////// VertexRvsZ //////////////////////////////////////////////////////////

std::reverse(runsZ.begin(),runsZ.end());
auto c_summary1 = new TCanvas();
TH2D *hSummary1 =  (TH2D*)histos1[0]->Clone("hSummary1");
for(int iplot=1; iplot<(int)histos1.size(); iplot++){
  hSummary1->Add(histos1[iplot]);
}
c_summary1->cd();
c_summary1->SetLogz();
hSummary1->SetTitle("Distance in transverse plane vs Z. Summary plot for all runs");
hSummary1->Draw("colz");


///////////////////////////VertexCoordinates AVG///////////////////////////////////////////////////

std::reverse(runs4.begin(),runs4.end());
std::reverse(xpoints.begin(),xpoints.end());
std::reverse(ypoints.begin(),ypoints.end());
std::reverse(xpoints_err.begin(),xpoints_err.end());
std::reverse(ypoints_err.begin(),ypoints_err.end());
std::reverse(rmsx_vertex.begin(), rmsx_vertex.end());
std::reverse(rmsx_vertex_err.begin(), rmsx_vertex_err.end());
std::reverse(rmsy_vertex.begin(), rmsy_vertex.end());
std::reverse(rmsy_vertex_err.begin(), rmsy_vertex_err.end());
auto c_summary = new TCanvas();
TH2D *hSummary =  (TH2D*)histos[0]->Clone("hSummary");
for(int iplot=1; iplot<(int)histos.size(); iplot++){
  hSummary->Add(histos[iplot]);
}
c_summary->cd();
c_summary->SetLogz();
hSummary->SetTitle("Coordinates of track vertex. Summary plot for all runs");
hSummary->Draw("colz");

auto cxy = new TCanvas();
cxy->SetGrid(0,1);
auto grx = new TGraphErrors(cycle4.size(),&cycle4[0],&xpoints[0], NULL, &xpoints_err[0]);
grx->GetXaxis()->SetNdivisions(cycle4.size());
TH1F hfake6("hfake6","hfake6", grx->GetN(),-0.5,(double)grx->GetN()-0.5);
hfake6.SetStats(0);
for (int i=1;i<=(int)cycle4.size();i++) hfake6.GetXaxis()->SetBinLabel(i,runs4[i-1].data());
hfake6.GetXaxis()->SetLabelSize(0.027);
grx->SetMarkerStyle(20);
hfake6.SetTitle("Vertex X-Y average positions");
hfake6.GetXaxis()->SetTitle("");
hfake6.GetYaxis()->SetTitle("Average (cm)");
hfake6.GetYaxis()->SetRangeUser(-0.3,0.3);
hfake6.Draw();
grx->Draw("PL same");

auto gry = new TGraphErrors(cycle4.size(),&cycle4[0],&ypoints[0], NULL, &ypoints_err[0]);
gry->SetMarkerStyle(20);
gry->SetMarkerColor(kRed);
gry->SetLineColor(kRed);
gry->Draw("PL SAME");

auto legend_xy=new TLegend(0.837272,0.746976,0.899033,0.900202);
legend_xy->SetFillColor(0);
legend_xy->SetHeader("Vertex","");
legend_xy->AddEntry(grx,"X","lp");
legend_xy->AddEntry(gry,"Y","lp");
legend_xy->Draw();


///////////////////////VertexCoordinates RMS////////////////////////////////////////
auto cxyrms = new TCanvas();
cxyrms->SetGrid(0,1);
auto grxrms = new TGraphErrors(cycle4.size(),&cycle4[0],&rmsx_vertex[0], NULL, &rmsx_vertex_err[0]);
grxrms->GetXaxis()->SetNdivisions(cycle4.size());
TH1F hfake6rms("hfake6rms","hfake6rms", grxrms->GetN(),-0.5,(double)grxrms->GetN()-0.5);
hfake6rms.SetStats(0);
for (int i=1;i<=(int)cycle4.size();i++) hfake6rms.GetXaxis()->SetBinLabel(i,runs4[i-1].data());
hfake6rms.GetXaxis()->SetLabelSize(0.027);
grxrms->SetMarkerStyle(20);
hfake6rms.SetTitle("Vertex X-Y RMS");
hfake6rms.GetXaxis()->SetTitle("");
hfake6rms.GetYaxis()->SetTitle("RMS (cm)");
hfake6rms.GetYaxis()->SetRangeUser(0,0.3);
hfake6rms.Draw();
grxrms->Draw("PL same");

auto gryrms = new TGraphErrors(cycle4.size(),&cycle4[0],&rmsy_vertex[0], NULL, &rmsy_vertex_err[0]);
gryrms->SetMarkerStyle(20);
gryrms->SetMarkerColor(kRed);
gryrms->SetLineColor(kRed);
gryrms->Draw("PL SAME");

legend_xy->Draw();

///////////////////////NVertexContributors AVG//////////////////////////////////////////

std::reverse(runs3.begin(),runs3.end());
std::reverse(aves2.begin(),aves2.end());
std::reverse(rmss2.begin(),rmss2.end());
std::reverse(aves2_err.begin(),aves2_err.end());
std::reverse(rmss2_err.begin(),rmss2_err.end());
auto ce = new TCanvas();
ce->SetGrid(0,1);
auto gerr2 = new TGraphErrors(cycle3.size(),&cycle3[0],&aves2[0],0,&aves2_err[0]);
gerr2->GetXaxis()->SetNdivisions(cycle3.size());
TH1F hfake7("hfake7","hfake7", gerr2->GetN(),-0.5,(double)gerr2->GetN()-0.5);
hfake7.SetStats(0);
for (int i=1;i<=(int)cycle3.size();i++) hfake7.GetXaxis()->SetBinLabel(i,runs3[i-1].data());
hfake7.GetXaxis()->SetLabelSize(0.027);
gerr2->SetMarkerStyle(20);
hfake7.SetTitle("Mean NVertexContributors");
hfake7.GetXaxis()->SetTitle("");
hfake7.GetYaxis()->SetTitle("<# of contributors for vertex>");
hfake7.GetYaxis()->SetRangeUser(0,40);
hfake7.Draw();
gerr2->Draw("PL same");

/////////////////// NVertexContributors RMS /////////////////////////////////////////////
auto cerms = new TCanvas();
cerms->SetGrid(0,1);
auto gerr2rms = new TGraphErrors(cycle3.size(),&cycle3[0],&rmss2[0],0,&rmss2_err[0]);
gerr2rms->GetXaxis()->SetNdivisions(cycle3.size());
TH1F hfake7rms("hfake7rms","hfake7rms", gerr2rms->GetN(),-0.5,(double)gerr2rms->GetN()-0.5);
hfake7rms.SetStats(0);
for (int i=1;i<=(int)cycle3.size();i++) hfake7rms.GetXaxis()->SetBinLabel(i,runs3[i-1].data());
hfake7rms.GetXaxis()->SetLabelSize(0.027);
gerr2rms->SetMarkerStyle(20);
hfake7rms.SetTitle("RMS of NVertexContributors");
hfake7rms.GetXaxis()->SetTitle("");
hfake7rms.GetYaxis()->SetTitle("Contributors RMS");
hfake7rms.GetYaxis()->SetRangeUser(0,40);
hfake7rms.Draw();
gerr2rms->Draw("PL same");


///////////////////////////Ntracks averages/////////////////////////////////////////////

std::reverse(runs2.begin(),runs2.end());
std::reverse(aves1.begin(),aves1.end());
std::reverse(aves1_err.begin(),aves1_err.end());
std::reverse(aves1_nozero.begin(),aves1_nozero.end());
std::reverse(aves1_nozero_err.begin(),aves1_nozero_err.end());
std::reverse(rmss1.begin(),rmss1.end());
std::reverse(rmss1_err.begin(),rmss1_err.end());
std::reverse(rmss1_nozero.begin(),rmss1_nozero.end());
std::reverse(rmss1_nozero_err.begin(),rmss1_nozero_err.end());

auto ct = new TCanvas();
ct->SetGrid(0,1);
auto gerr1 = new TGraphErrors(cycle2.size(),&cycle2[0],&aves1[0],0,&aves1_err[0]);
auto gerr1_nozero = new TGraphErrors(cycle2.size(),&cycle2[0],&aves1_nozero[0],0,&aves1_nozero_err[0]);
gerr1->GetXaxis()->SetNdivisions(cycle2.size());
gerr1_nozero->GetXaxis()->SetNdivisions(cycle2.size());
TH1F hfake8("hfake8","hfake8", gerr1->GetN(),-0.5,(double)gerr1->GetN()-0.5);
hfake8.SetStats(0);
for (int i=1;i<=(int)cycle2.size();i++) hfake8.GetXaxis()->SetBinLabel(i,runs2[i-1].data());
hfake8.GetXaxis()->SetLabelSize(0.027);
gerr1->SetMarkerStyle(20);
gerr1->SetMarkerColor(kRed);
gerr1->SetLineColor(kRed);
gerr1_nozero->SetMarkerStyle(20);
gerr1_nozero->SetMarkerColor(kBlue);
gerr1_nozero->SetLineColor(kBlue);
hfake8.SetTitle("Mean number of tracks event by event (blue: removed bin at n_tracks = 0)");
hfake8.GetXaxis()->SetTitle("");
hfake8.GetYaxis()->SetTitle("<# tracks>");
hfake8.GetYaxis()->SetRangeUser(0, 100);
hfake8.Draw();
gerr1->Draw("PL same");
gerr1_nozero->Draw("PL same");

//////////////////////// Ntracks RMS /////////////////////////////////////////////////
auto ctrms = new TCanvas();
ctrms->SetGrid(0,1);
auto gerr1rms = new TGraphErrors(cycle2.size(),&cycle2[0],&rmss1[0],0,&rmss1_err[0]);
auto gerr1rms_nozero = new TGraphErrors(cycle2.size(),&cycle2[0],&rmss1_nozero[0],0,&rmss1_nozero_err[0]);
gerr1rms->GetXaxis()->SetNdivisions(cycle2.size());
gerr1rms_nozero->GetXaxis()->SetNdivisions(cycle2.size());
TH1F hfake8rms("hfake8rms","hfake8rms", gerr1rms->GetN(),-0.5,(double)gerr1rms->GetN()-0.5);
hfake8rms.SetStats(0);
for (int i=1;i<=(int)cycle2.size();i++) hfake8rms.GetXaxis()->SetBinLabel(i,runs2[i-1].data());
hfake8rms.GetXaxis()->SetLabelSize(0.027);
gerr1rms->SetMarkerStyle(20);
gerr1rms->SetMarkerColor(kRed);
gerr1rms->SetLineColor(kRed);
gerr1rms_nozero->SetMarkerStyle(20);
gerr1rms_nozero->SetMarkerColor(kBlue);
gerr1rms_nozero->SetLineColor(kBlue);
hfake8rms.SetTitle("Rms of number of tracks event by event (blue: removed bin at n_tracks = 0)");
hfake8rms.GetXaxis()->SetTitle("");
hfake8rms.GetYaxis()->SetTitle("# tracks RMS");
hfake8rms.GetYaxis()->SetRangeUser(0, 100);
hfake8rms.Draw();
gerr1rms->Draw("PL same");
gerr1rms_nozero->Draw("PL same");


////////////////////////AssociatedClusterFraction Avg/////////////////////////////////////

std::reverse(runs1.begin(),runs1.end());
std::reverse(aves.begin(),aves.end());
std::reverse(aves_err.begin(),aves_err.end());
std::reverse(rmss.begin(),rmss.end());
std::reverse(rmss_err.begin(),rmss_err.end());
std::reverse(aves_nozero.begin(),aves_nozero.end());
std::reverse(aves_nozero_err.begin(),aves_nozero_err.end());
std::reverse(rmss_nozero.begin(),rmss_nozero.end());
std::reverse(rmss_nozero_err.begin(),rmss_nozero_err.end());
auto c0 = new TCanvas();
c0->SetGrid(0,1);
auto gerr = new TGraphErrors(cycle1.size(),&cycle1[0],&aves[0],0,&aves_err[0]);
auto gerr_nozero = new TGraphErrors(cycle1.size(),&cycle1[0],&aves_nozero[0],0,&aves_nozero_err[0]);
gerr->GetXaxis()->SetNdivisions(cycle1.size());
gerr_nozero->GetXaxis()->SetNdivisions(cycle1.size());
TH1F hfake9("hfake9","hfake9", gerr->GetN(),-0.5,(double)gerr->GetN()-0.5);
hfake9.SetStats(0);
for (int i=1;i<=(int)cycle1.size();i++) hfake9.GetXaxis()->SetBinLabel(i,runs1[i-1].data());
hfake9.GetXaxis()->SetLabelSize(0.027);
gerr->SetMarkerStyle(20);
gerr->SetMarkerColor(kRed);
gerr->SetLineColor(kRed);
gerr_nozero->SetMarkerStyle(20);
gerr_nozero->SetMarkerColor(kBlue);
gerr_nozero->SetLineColor(kBlue);
hfake9.SetTitle("Mean of fraction of clusters into tracks event by event (blue: w/o entries at fraction = 0)");
hfake9.GetXaxis()->SetTitle("");
hfake9.GetYaxis()->SetTitle("<Clusters in tracks/All clusters>");
hfake9.GetYaxis()->SetRangeUser(0,1);
hfake9.Draw();
gerr->Draw("PL same");
gerr_nozero->Draw("PL same");

/////////////////////////AssociatedClusterFraction RMS////////////////////////////////
auto c0rms = new TCanvas();
c0rms->SetGrid(0,1);
auto gerrrms = new TGraphErrors(cycle1.size(),&cycle1[0],&rmss[0],0,&rmss_err[0]);
auto gerrrms_nozero = new TGraphErrors(cycle1.size(),&cycle1[0],&rmss_nozero[0],0,&rmss_nozero_err[0]);
gerrrms->GetXaxis()->SetNdivisions(cycle1.size());
gerrrms_nozero->GetXaxis()->SetNdivisions(cycle1.size());
TH1F hfake9rms("hfake9rms","hfake9rms", gerrrms->GetN(),-0.5,(double)gerrrms->GetN()-0.5);
hfake9rms.SetStats(0);
for (int i=1;i<=(int)cycle1.size();i++) hfake9rms.GetXaxis()->SetBinLabel(i,runs1[i-1].data());
hfake9rms.GetXaxis()->SetLabelSize(0.027);
gerrrms->SetMarkerStyle(20);
gerrrms->SetMarkerColor(kRed);
gerrrms->SetLineColor(kRed);
gerrrms_nozero->SetMarkerStyle(20);
gerrrms_nozero->SetMarkerColor(kBlue);
gerrrms_nozero->SetLineColor(kBlue);
hfake9rms.SetTitle("Rms of fraction of clusters into tracks event by event (blue: w/o entries at fraction = 0)");
hfake9rms.GetXaxis()->SetTitle("");
hfake9rms.GetYaxis()->SetTitle("RMS");
hfake9rms.GetYaxis()->SetRangeUser(0,2);
hfake9rms.Draw();
gerrrms->Draw("PL same");
gerrrms_nozero->Draw("PL same");

/////////////////////////NClustersPerTrackEta/////////////////////////////////////////
std::reverse(runs.begin(),runs.end());
std::reverse(means.begin(),means.end());
std::reverse(means1.begin(),means1.end());
std::reverse(means2.begin(),means2.end());
std::reverse(means3.begin(),means3.end());
std::reverse(means4.begin(),means4.end());
std::reverse(means5.begin(),means5.end());

auto c_ = new TCanvas();
c_->SetGrid(0,1);
auto gr_ = new TGraph(cycle.size(),&cycle[0],&means[0]);
gr_->GetXaxis()->SetNdivisions(cycle.size());
TH1F hfake10("hfake10","hfake10", gr_->GetN(),-0.5,(double)gr_->GetN()-0.5);
hfake10.SetStats(0);
for (int i=1;i<=(int)cycle.size();i++) hfake10.GetXaxis()->SetBinLabel(i,runs[i-1].data());
hfake10.GetXaxis()->SetLabelSize(0.027);
gr_->SetMarkerStyle(20);
hfake10.SetTitle("Mean number of clusters per track");
hfake10.GetXaxis()->SetTitle("");
hfake10.GetYaxis()->SetTitle("<# of clusters per track>");
hfake10.GetYaxis()->SetRangeUser(0,10);
hfake10.Draw();
gr_->Draw("PL same");

auto gr_1 = new TGraph(cycle.size(),&cycle[0],&means1[0]);
gr_1->SetMarkerStyle(20);
gr_1->SetMarkerColor(kRed);
gr_1->SetLineColor(kRed);
gr_1->Draw("PL SAME");

auto gr_2 = new TGraph(cycle.size(),&cycle[0],&means2[0]);
gr_2->SetMarkerStyle(20);
gr_2->SetMarkerColor(kBlue);
gr_2->SetLineColor(kBlue);
gr_2->Draw("PL SAME");

auto gr_3 = new TGraph(cycle.size(),&cycle[0],&means3[0]);
gr_3->SetMarkerStyle(20);
gr_3->SetMarkerColor(kGreen);
gr_3->SetLineColor(kGreen);
gr_3->Draw("PL SAME");

auto gr_4 = new TGraph(cycle.size(),&cycle[0],&means4[0]);
gr_4->SetMarkerStyle(20);
gr_4->SetMarkerColor(kOrange);
gr_4->SetLineColor(kOrange);
gr_4->Draw("PL SAME");

auto gr_5 = new TGraph(cycle.size(),&cycle[0],&means5[0]);
gr_5->SetMarkerStyle(20);
gr_5->SetMarkerColor(kMagenta);
gr_5->SetLineColor(kMagenta);
gr_5->Draw("PL SAME");

auto legend_=new TLegend(0.828142,0.698589,0.928034,0.933468);
legend_->SetFillColor(0);
legend_->SetHeader("#eta range","");
legend_->AddEntry(gr_,"|#eta|<1.5","lp");
legend_->AddEntry(gr_1,"|#eta|<1.2","lp");
legend_->AddEntry(gr_2,"|#eta|<0.8","lp");
legend_->AddEntry(gr_5,"|#eta|<0.2","lp");
legend_->AddEntry(gr_3,"-1.5<#eta<-1.4","lp");
legend_->AddEntry(gr_4,"1.4<#eta<1.5","lp");
legend_->Draw();

std::reverse(binminCenters.begin(),binminCenters.end());
std::reverse(binminCenters1.begin(),binminCenters1.end());
std::reverse(binminCenters2.begin(),binminCenters2.end());
std::reverse(binminCenters3.begin(),binminCenters3.end());
std::reverse(binminCenters4.begin(),binminCenters4.end());
std::reverse(binminCenters5.begin(),binminCenters5.end());
auto c = new TCanvas();
c->SetGrid(0,1);
auto gr = new TGraph(cycle.size(),&cycle[0],&binminCenters[0]);
gr->GetXaxis()->SetNdivisions(cycle.size());
TH1F hfake11("hfake11","hfake11", gr->GetN(),-0.5,(double)gr->GetN()-0.5);
hfake11.SetStats(0);
for (int i=1;i<=(int)cycle.size();i++) hfake11.GetXaxis()->SetBinLabel(i,runs[i-1].data());
hfake11.GetXaxis()->SetLabelSize(0.027);
gr->SetMarkerStyle(20);
hfake11.SetTitle("Less probable #clusters per track");
hfake11.GetXaxis()->SetTitle("");
hfake11.GetYaxis()->SetTitle("Less probable #clusters per track");
hfake11.GetYaxis()->SetRangeUser(0,10);
hfake11.Draw();
gr->Draw("PL same");

auto gr1 = new TGraph(cycle.size(),&cycle[0],&binminCenters1[0]);
gr1->SetMarkerStyle(20);
gr1->SetMarkerColor(kRed);
gr1->SetLineColor(kRed);
gr1->Draw("PL SAME");

auto gr20 = new TGraph(cycle.size(),&cycle[0],&binminCenters2[0]);
gr20->SetMarkerStyle(20);
gr20->SetMarkerColor(kBlue);
gr20->SetLineColor(kBlue);
gr20->Draw("PL SAME");

auto gr30 = new TGraph(cycle.size(),&cycle[0],&binminCenters3[0]);
gr30->SetMarkerStyle(20);
gr30->SetMarkerColor(kGreen);
gr30->SetLineColor(kGreen);
gr30->Draw("PL SAME");

auto gr4 = new TGraph(cycle.size(),&cycle[0],&binminCenters4[0]);
gr4->SetMarkerStyle(20);
gr4->SetMarkerColor(kOrange);
gr4->SetLineColor(kOrange);
gr4->Draw("PL SAME");

auto gr5 = new TGraph(cycle.size(),&cycle[0],&binminCenters5[0]);
gr5->SetMarkerStyle(20);
gr5->SetMarkerColor(kMagenta);
gr5->SetLineColor(kMagenta);
gr5->Draw("PL SAME");

auto legend2=new TLegend(0.828142,0.698589,0.928034,0.933468);
legend2->SetFillColor(0);
legend2->SetHeader("#eta range","");
legend2->AddEntry(gr,"|#eta|<1.5","lp");
legend2->AddEntry(gr1,"|#eta|<1.2","lp");
legend2->AddEntry(gr20,"|#eta|<0.8","lp");
legend2->AddEntry(gr5,"|#eta|<0.2","lp");
legend2->AddEntry(gr30,"-1.5<#eta<-1.4","lp");
legend2->AddEntry(gr4,"1.4<#eta<1.5","lp");
legend2->Draw();

std::reverse(binmaxCenters.begin(),binmaxCenters.end());
std::reverse(binmaxCenters1.begin(),binmaxCenters1.end());
std::reverse(binmaxCenters2.begin(),binmaxCenters2.end());
std::reverse(binmaxCenters3.begin(),binmaxCenters3.end());
std::reverse(binmaxCenters4.begin(),binmaxCenters4.end());
std::reverse(binmaxCenters5.begin(),binmaxCenters5.end());
auto c1 = new TCanvas();
c1->SetGrid(0,1);
auto gr0 = new TGraph(cycle.size(),&cycle[0],&binmaxCenters[0]);
gr0->GetXaxis()->SetNdivisions(cycle.size());
TH1F hfake12("hfake12","hfake12", gr0->GetN(),-0.5,(double)gr0->GetN()-0.5);
hfake12.SetStats(0);
for (int i=1;i<=(int)cycle.size();i++) hfake12.GetXaxis()->SetBinLabel(i,runs[i-1].data());
hfake12.GetXaxis()->SetLabelSize(0.027);
gr0->SetMarkerStyle(20);
hfake12.SetTitle("Most probable #clusters per track");
hfake12.GetXaxis()->SetTitle("");
hfake12.GetYaxis()->SetTitle("Most probable #clusters per track");
hfake12.GetYaxis()->SetRangeUser(0,10);
hfake12.Draw();
gr0->Draw("PL same");

auto gr01 = new TGraph(cycle.size(),&cycle[0],&binmaxCenters1[0]);
gr01->SetMarkerStyle(20);
gr01->SetMarkerColor(kRed);
gr01->SetLineColor(kRed);
gr01->Draw("PL SAME");

auto gr02 = new TGraph(cycle.size(),&cycle[0],&binmaxCenters2[0]);
gr02->SetMarkerStyle(20);
gr02->SetMarkerColor(kBlue);
gr02->SetLineColor(kBlue);
gr02->Draw("PL SAME");

auto gr03 = new TGraph(cycle.size(),&cycle[0],&binmaxCenters3[0]);
gr03->SetMarkerStyle(20);
gr03->SetMarkerColor(kGreen);
gr03->SetLineColor(kGreen);
gr03->Draw("PL SAME");

auto gr04 = new TGraph(cycle.size(),&cycle[0],&binmaxCenters4[0]);
gr04->SetMarkerStyle(20);
gr04->SetMarkerColor(kOrange);
gr04->SetLineColor(kOrange);
gr04->Draw("PL SAME");

auto gr05 = new TGraph(cycle.size(),&cycle[0],&binmaxCenters5[0]);
gr05->SetMarkerStyle(20);
gr05->SetMarkerColor(kMagenta);
gr05->SetLineColor(kMagenta);
gr05->Draw("PL SAME");

auto legend0=new TLegend(0.828142,0.698589,0.928034,0.933468);
legend0->SetFillColor(0);
legend0->SetHeader("#eta range","");
legend0->AddEntry(gr0,"|#eta|<1.5","lp");
legend0->AddEntry(gr01,"|#eta|<1.2","lp");
legend0->AddEntry(gr02,"|#eta|<0.8","lp");
legend0->AddEntry(gr05,"|#eta|<0.2","lp");
legend0->AddEntry(gr03,"-1.5<#eta<-1.4","lp");
legend0->AddEntry(gr04,"1.4<#eta<1.5","lp");
legend0->Draw();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

c2->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf[", runs.front().c_str(), runs.back().c_str()));
c2->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c21->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c3->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c31->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
cz1->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
cz11->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c_summary1->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c_summary->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
cxy->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
cxyrms->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
ce->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
cerms->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
ct->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
ctrms->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c0->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c0rms->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c_->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c1->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf", runs.front().c_str(), runs.back().c_str()));
c1->SaveAs(Form("../Plots/track_qc_from_run%s_to_run%s.pdf]", runs.front().c_str(), runs.back().c_str()));

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if(ccdb_upload){
string Runperiod = Form("run%s_to_run%s%s",runs.front().c_str(), runs.back().c_str(), isVtxRun3 ? "_vtxRun3" :"");

	c1->SetName("MostProbable_Clusters_per_track_vs_eta");
auto mo1 = std::make_shared<o2::quality_control::core::MonitorObject>(c1, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo1->setIsOwner(false);
        ccdb->storeMO(mo1);

	c2->SetName("Average_Eta_distribution");
auto mo2 = std::make_shared<o2::quality_control::core::MonitorObject>(c2, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo2->setIsOwner(false);
        ccdb->storeMO(mo2);

	c21->SetName("Rel_number_tracks_Eta_distribution");
auto mo21 = std::make_shared<o2::quality_control::core::MonitorObject>(c21, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo21->setIsOwner(false);
        ccdb->storeMO(mo21);

	c3->SetName("Average_Phi_distribution");
auto mo3 = std::make_shared<o2::quality_control::core::MonitorObject>(c3, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo3->setIsOwner(false);
        ccdb->storeMO(mo3);

	c31->SetName("Rel_number_tracks_Phi_distribution");
auto mo31 = std::make_shared<o2::quality_control::core::MonitorObject>(c31, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo31->setIsOwner(false);
        ccdb->storeMO(mo31);

	cz1->SetName("Mean_Zvertex_coordinate");
auto moz1 = std::make_shared<o2::quality_control::core::MonitorObject>(cz1, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        moz1->setIsOwner(false);
        ccdb->storeMO(moz1);

        cz11->SetName("Rms_Zvertex_coordinate");
auto moz11 = std::make_shared<o2::quality_control::core::MonitorObject>(cz11, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        moz11->setIsOwner(false);
        ccdb->storeMO(moz11);

	c_summary1->SetName("Summary_Distance_Primary_Vertex");
auto mosummary1 = std::make_shared<o2::quality_control::core::MonitorObject>(c_summary1, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mosummary1->setIsOwner(false);
        ccdb->storeMO(mosummary1);

        c_summary->SetName("Summary_Track_Vertex_Coordinates");
auto mosummary = std::make_shared<o2::quality_control::core::MonitorObject>(c_summary, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mosummary->setIsOwner(false);
        ccdb->storeMO(mosummary);

        cxy->SetName("Track_Vertex_Coordinates");
auto moxy = std::make_shared<o2::quality_control::core::MonitorObject>(cxy, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        moxy->setIsOwner(false);
        ccdb->storeMO(moxy);

cxyrms->SetName("Track_Vertex_Coordinates_Rms");
auto moxyrms = std::make_shared<o2::quality_control::core::MonitorObject>(cxyrms, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        moxyrms->setIsOwner(false);
        ccdb->storeMO(moxyrms);

        ce->SetName("Mean_NVertexContributors");
auto moe = std::make_shared<o2::quality_control::core::MonitorObject>(ce, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        moe->setIsOwner(false);
        ccdb->storeMO(moe);

        cerms->SetName("Rms_NVertexContributors");
auto moerms = std::make_shared<o2::quality_control::core::MonitorObject>(cerms, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        moerms->setIsOwner(false);
        ccdb->storeMO(moerms);

        ct->SetName("Mean_NTracks_EbyE");
auto mot = std::make_shared<o2::quality_control::core::MonitorObject>(ct, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mot->setIsOwner(false);
        ccdb->storeMO(mot);

        ctrms->SetName("Rms_NTracks_EbyE");
auto motrms = std::make_shared<o2::quality_control::core::MonitorObject>(ctrms, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        motrms->setIsOwner(false);
        ccdb->storeMO(motrms);

        c0->SetName("Mean_Fraction_Clusters_Tracks");
auto mo0 = std::make_shared<o2::quality_control::core::MonitorObject>(c0, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo0->setIsOwner(false);
        ccdb->storeMO(mo0);

        c0rms->SetName("Rms_Fraction_Clusters_Tracks");
auto mo0rms = std::make_shared<o2::quality_control::core::MonitorObject>(c0rms, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo0rms->setIsOwner(false);
        ccdb->storeMO(mo0rms);


        c_->SetName("Mean_NumberCluster_per_Track");
auto mo_ = std::make_shared<o2::quality_control::core::MonitorObject>(c_, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        mo_->setIsOwner(false);
        ccdb->storeMO(mo_);

        c->SetName("LessProbable_Clusters_per_track_vs_eta");
auto moc = std::make_shared<o2::quality_control::core::MonitorObject>(c, TaskName, TaskClass, DetectorName,std::stoi(runs.front()),Runperiod);
        moc->setIsOwner(false);
        ccdb->storeMO(moc);


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TFile *f1 = new TFile(Form("../Plots/TrackTask_Result_trends_from_run%s_to_run%s.root", runs.front().c_str(), runs.back().c_str()), "RECREATE");
f1->cd();
c2->Write();
c21->Write();
c3->Write();
c31->Write();
cz1->Write();
cz11->Write();
c_summary1->Write();
c_summary->Write();
cxy->Write();
cxyrms->Write();
ce->Write();
cerms->Write();
ct->Write();
ctrms->Write();
c0->Write();
c0rms->Write();
c_->Write();
c->Write();
c1->Write();
f1->Close();

delete gr2; delete gr21; delete gr22; delete legend; delete c2; delete c21;
delete gr3; delete gr31; delete gr32; delete legend1; delete c3; delete c31;
delete gerr3; delete cz1; delete cz11;
delete hSummary1; delete c_summary1;
delete hSummary; delete c_summary;
delete grx; delete gry; delete legend_xy; delete cxy; delete cxyrms;
delete gerr2; delete ce; delete cerms;
delete gerr1; delete ct; delete ctrms;
delete gerr; delete c0; delete c0rms;
delete gr_; delete gr_1; delete gr_2; delete gr_3; delete gr_4; delete gr_5; delete legend_; delete c_;
delete gr; delete gr1; delete gr20; delete gr30; delete gr4; delete gr5; delete legend2; delete c;
delete gr0; delete gr01; delete gr02; delete gr03; delete gr04; delete gr05; delete legend0; delete c1;

ccdb->disconnect();
}
