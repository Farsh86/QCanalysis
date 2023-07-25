//utilities and constants
#include "inc/constants.h"
#include "inc/utilities.h"

#include<algorithm>
#include <cstdio>

using namespace std;

void DoAnalysis(string filepath, const int nChips, bool isIB);

//
// MAIN
//
void ReportTHR_shifters(){
  //string fpath;
  int nchips=9;
  //cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  string fpath = (string)gSystem->GetFromPipe("ls ../Data/*THRMAPS_DEADPIXMAPS* -Art | tail -n 1"); //take most recent data file

  bool isIB;
  if(fpath.find("IB")!=string::npos){
    isIB = kTRUE;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    if(stoi(layernum)>=0 && stoi(layernum)<=2) nchips = 9;
    else if (stoi(layernum)==3 && stoi(layernum)==4) nchips = 54*2;
    else nchips = 98*2;
    if(nchips==9) isIB=kTRUE;
    else isIB=kFALSE;
  }

  //Call: do full analysis
  DoAnalysis(fpath, nchips, isIB);
}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, bool isIB){

  gStyle->SetOptStat(0000);

  string localdatetime = GetCurrentDateTime(1);

  std::freopen(Form("../logs/logTHR_%s.log",localdatetime.c_str()), "w", stdout);

  std::vector<TH2*> hmapsTHR;
  std::vector<TH2*> hmapsDEAD;
  std::vector<THnSparse*> hmapsDEADPIX;

  std::vector<string> timestamps, runnumbers, laynums, stavenums, runlabel;
  int nTimes=0, nRuns=1, cstv=0;

  //Read the file and the list of plots with entries that are necessary for all the analyses
  cout<<endl;
  cout<<"... Opening file: "<<filepath<<endl;
  cout<<endl;
  TFile *infile=new TFile(filepath.c_str());
  TList *list = infile->GetListOfKeys();
  TKey *key;
  TObject *obj;
  TIter next(list);
  TH2 *h2;
  THnSparse *hsparse;
  bool isfirst = true;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D or sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    string objtitle = (string)obj->GetTitle();

    if(objname.find("Stv")==string::npos) h2 = (TH2*)obj;
    else hsparse = (THnSparse*)obj;
    //if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    string timestamp = objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = "nostave";

    if(objname.find("Stv")!=string::npos){ //DEAD pix MAPs
      stvnum = objname.substr(objname.find("Stv")+3,2);
      if(stvnum.find("_")!=string::npos){
        stvnum = objname.substr(objname.find("Stv")+3,1);
      }
      stavenums.push_back(stvnum);
      hmapsDEADPIX.push_back(hsparse);
      if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-2] && isfirst){
        isfirst=false;
        continue;
      }
    }
    else if(objtitle.find("Threshold")!=string::npos){// THR maps
      hmapsTHR.push_back(h2);
      timestamps.push_back(timestamp);
      runnumbers.push_back(runnum);
      laynums.push_back(laynum);
      nTimes++;
      if(nTimes>1 && laynum==laynums[laynums.size()-2]){
        nRuns++;
      }
      else nRuns=1;

      if(nTimes<=1) runlabel.push_back(runnum);
      if(nTimes>1 && laynum==laynums[laynums.size()-2] && isfirst){
        runlabel.push_back(runnum);
      }
      if(nTimes>1 && laynum!=laynums[laynums.size()-2]){
        isfirst = false;
      }

    }
    else if(objtitle.find("DeadPixel")){
      hmapsDEAD.push_back(h2);
    }
    //cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
  }

  cout<<"#Runs: "<<nRuns<<endl;
  const int nLayers = (int)hmapsTHR.size()==nRuns ? 1 : stoi(laynums[laynums.size()-1])+1;

  //**************************************************************************************
  //************************* Avg threshold vs run number ********************************
  //**************************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* Average threshold vs run number ****************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;
  TGraph *trend[nLayers][100];
  int ilayer=nLayers-1;
  for(int ihist=(int)hmapsTHR.size()-1; ihist>=0; ihist--){
    for(int ibiny=1; ibiny<=hmapsTHR[ihist]->GetNbinsY(); ibiny++){
      trend[ilayer][ibiny-1] = new TGraph();
    }
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
  }
  ilayer=nLayers-1;
  string histname = hmapsTHR[0]->GetName();
  int irun=0;
  for(int ihist=(int)hmapsTHR.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    for(int ibiny=1; ibiny<=hmapsTHR[ihist]->GetNbinsY(); ibiny++){//loop on y bins (staves)
      TH1D *hproj = hmapsTHR[ihist]->ProjectionX("proj",ibiny,ibiny); //single stave
      trend[ilayer][ibiny-1]->SetName(Form("gr_L%s_stave%d",laynums[ihist].c_str(),ibiny-1));
      int deadchips = 0;
      for(int ibinx=1; ibinx<=hmapsTHR[ihist]->GetNbinsX(); ibinx++){//evaluate the number of disabled chips
        if(hmapsTHR[ihist]->GetBinContent(ibinx,ibiny)<1e-15)
          deadchips++;
      }
      if(deadchips>0)
        cout<<"Layer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" --> Chips active:"<<nChips-deadchips<<endl;

      if(deadchips!=nChips)
        trend[ilayer][ibiny-1]->SetPoint(irun, irun, hproj->Integral()/(nChips-deadchips));
      else
        trend[ilayer][ibiny-1]->SetPoint(irun, irun, 0.);

      if((ibiny-1)<hmapsTHR[ihist]->GetNbinsY()/2)
        SetStyle(trend[ilayer][ibiny-1], col[ibiny-1], 24);
      else
        SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-hmapsTHR[ihist]->GetNbinsY()/2], 26);
    }
    irun++;
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        irun=0;
        ilayer--;
      }
  }

  //Draw
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas canvas;
    canvas.cd();
    canvas.SetTickx();
    canvas.SetTicky();
    canvas.SetMargin(0.0988,0.1,0.194,0.0993);
    TLegend leg(0.904, 0.197,0.997,0.898);
    for(int istave=0; istave<hmapsTHR[ilay*nRuns]->GetNbinsY(); istave++)
      leg.AddEntry(trend[ilay][istave], Form("Stv%d",istave), "p");

    int npoints = trend[ilay][0]->GetN();
    TH1F *hfake = new TH1F("hfake", "; Run; Avg. Threshold (DAC)", npoints, -0.5, (double)npoints-0.5);
    int start = 0;
    for(int inum=0; inum<(int)laynums.size(); inum++){
      if(laynums[inum].find(to_string(ilay))!=string::npos){
        start = inum;
        break;
      }
    }
    for(int ir=0; ir<npoints; ir++){
        hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[start+npoints-1-ir])));
    }
    hfake->GetYaxis()->SetRangeUser(8.5, 14.5);
    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->SetTitle(Form("Layer-%s, %s",laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    hfake->Draw();

    for(int istave=0; istave<hmapsTHR[ilay*nRuns]->GetNbinsY(); istave++)
      trend[ilay][istave]->Draw("P same");
    leg.Draw("same");
    if(!ilay) canvas.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf[", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }

  //**************************************************************************************
  //************************* # Dead pixels vs run number ********************************
  //**************************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* # Dead Pixels vs run number ********************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  TGraph *trendDead[nLayers][100];
  ilayer=nLayers-1;
  for(int ihist=(int)hmapsDEAD.size()-1; ihist>=0; ihist--){
    for(int ibiny=1; ibiny<=hmapsDEAD[ihist]->GetNbinsY(); ibiny++){
      trendDead[ilayer][ibiny-1] = new TGraph();
    }
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
  }
  ilayer=nLayers-1;
  histname = hmapsDEAD[0]->GetName();
  irun=0;
  double maxtotdead=-1.;
  int staveswithdead[nLayers];
  for(int ilay=0; ilay<nLayers; ilay++)
    staveswithdead[ilay] = 0;
  for(int ihist=(int)hmapsDEAD.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    for(int ibiny=1; ibiny<=hmapsDEAD[ihist]->GetNbinsY(); ibiny++){//loop on y bins (staves)
      TH1D *hproj = hmapsDEAD[ihist]->ProjectionX("proj",ibiny,ibiny); //single stave
      trendDead[ilayer][ibiny-1]->SetName(Form("grDead_L%s_stave%d",laynums[ihist].c_str(),ibiny-1));
      int chipswithdeadpix = 0;
      for(int ibinx=1; ibinx<=hmapsDEAD[ihist]->GetNbinsX(); ibinx++){//evaluate the chips with dead pixels
        if(hmapsDEAD[ihist]->GetBinContent(ibinx,ibiny)>0)
          chipswithdeadpix++;
      }
      if(chipswithdeadpix>0){
        cout<<"Layer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" --> # Chips with dead pix: "<<chipswithdeadpix<<endl;
        staveswithdead[ilayer]++;
      }

      trendDead[ilayer][ibiny-1]->SetPoint(irun, irun, hproj->Integral());//total number of dead pix for this stave
      if(hproj->Integral()>maxtotdead)
        maxtotdead=hproj->Integral();

      if((ibiny-1)<hmapsDEAD[ihist]->GetNbinsY()/2)
        SetStyle(trendDead[ilayer][ibiny-1], col[ibiny-1], 24);
      else
        SetStyle(trendDead[ilayer][ibiny-1], col[ibiny-1-hmapsDEAD[ihist]->GetNbinsY()/2], 26);
    }
    irun++;
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        irun=0;
        ilayer--;
      }
  }

  //Draw
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas canvas;
    canvas.cd();
    canvas.SetLogy();
    canvas.SetTickx();
    canvas.SetTicky();
    canvas.SetMargin(0.0988,0.1,0.194,0.0993);
    TLegend leg(0.904, 0.197,0.997,0.898);
    for(int istave=0; istave<hmapsDEAD[ilay*nRuns]->GetNbinsY(); istave++)
      leg.AddEntry(trendDead[ilay][istave], Form("Stv%d",istave), "p");

    int start = 0;
    for(int inum=0; inum<(int)laynums.size(); inum++){
      if(laynums[inum].find(to_string(ilay))!=string::npos){
        start = inum;
        break;
      }
    }
    int npoints = trendDead[ilay][0]->GetN();
    TH1F *hfakeDead = new TH1F("hfakeDead", "; Run; # Dead Pixels", npoints, -0.5, (double)npoints-0.5);
    for(int ir=0; ir<npoints; ir++)
        hfakeDead->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[start+npoints-1-ir])));
    hfakeDead->GetYaxis()->SetRangeUser(8e-1, maxtotdead+0.5*maxtotdead);
    hfakeDead->GetXaxis()->SetTitleOffset(2.8);
    hfakeDead->SetTitle(Form("Layer-%s (%d Staves with dead pix), %s",laynums[ilay*nRuns].c_str(), (int)ceil(staveswithdead[ilay]/nRuns), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    hfakeDead->Draw();
    for(int istave=0; istave<hmapsDEAD[ilay*nRuns]->GetNbinsY(); istave++)
      trendDead[ilay][istave]->Draw("P same");
    leg.Draw("same");

    canvas.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }

  //***********************************************************************************
  //***************************** Threshold correlation *******************************
  //***********************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* Threshold correlation analysis *****************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  //count occurrences of runs listed in the runlabel vector
  vector<int> occurrences;
  for(int ilab=0; ilab<(int)runlabel.size(); ilab++){
    int count = 0;
    for(int iall=0; iall<(int)runnumbers.size(); iall++){
      if(runlabel[ilab]==runnumbers[iall])
        count++;
    }
    occurrences.push_back(count);
  }

  //Reference runs
  //Choose two random reference runs
  const int pos1 = (int)nRuns*0.25;//pos of ref run 1
  const int pos2 = (int)nRuns*0.75;//pos of ref run 2;
  long int refrun[2];
  vector<array<int,2>> posrefrun;
  array<int,2> tmparray;
  for(int ilab=0; ilab<(int)runlabel.size(); ilab++){
    if(occurrences[ilab]==nLayers){
      refrun[0] = stol(runlabel[ilab]);
      break;
    }
  }
  for(int ilab=(int)runlabel.size()-1; ilab>=0; ilab--){
    if(occurrences[ilab]==nLayers){
      refrun[1] = stol(runlabel[ilab]);
      break;
    }
  }
  int icnt = 0;
  for(int ihist=0; ihist<(int)hmapsTHR.size(); ihist++){
    string hname = hmapsTHR[ihist]->GetName();
    string runn =  hname.substr(hname.find("run")+3, 6);
    if(stol(runn)==refrun[icnt]){ //position of refence run
      tmparray[icnt] = ihist;
      icnt++;
      if(icnt==2) {
        posrefrun.push_back(tmparray);
        icnt=0;
      }
    }
  }
  cout<<endl;
  cout<<"Reference Run 1: "<<refrun[0]<<endl;
  cout<<"Reference Run 2: "<<refrun[1]<<endl;
  cout<<endl;

  //analysis

  TH2D *hCorrTHR[2][nLayers];

  for(int iref=0; iref<2; iref++)
    for(int ilay=0; ilay<nLayers; ilay++)
      hCorrTHR[iref][ilay] = new TH2D(Form("hCorrTHR_L%s_refrun_%ld",laynums[ilay*nRuns].c_str(),refrun[iref]), Form("Layer-%s - THR corr. %s - Ref. run: %ld; Chip Threshold (run%ld); Chip Threshold (runs)",laynums[ilay*nRuns].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str(),refrun[iref],refrun[iref]), 150, 7.5, 13.5, 150, 7.5, 13.5);

  for(int iref=0; iref<2; iref++){
    ilayer=nLayers-1;
    for(int ihist=(int)hmapsTHR.size()-1; ihist>=0; ihist--){
      if(stol(runnumbers[ihist])==refrun[iref]){
        if(ihist>0)
          if(laynums[ihist-1]!=laynums[ihist]){// in case the ref run is the first in the list of all layers
            ilayer--;
          }
        continue; //skip ref run
      }
      for(int ibinx=1; ibinx<=hmapsTHR[ihist]->GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<=hmapsTHR[ihist]->GetNbinsY(); ibiny++){
          double data_refrun = hmapsTHR[posrefrun[ilayer][iref]]->GetBinContent(ibinx,ibiny);
          double data_run = hmapsTHR[ihist]->GetBinContent(ibinx,ibiny);
          hCorrTHR[iref][ilayer]->Fill(data_refrun, data_run);
          //cout<<fhr_refrun<<"   "<<fhr_run<<endl;
        }
      }
      if(ihist>0)
        if(laynums[ihist-1]!=laynums[ihist]){
          ilayer--;
        }
    }
  }//end loop on reference runs

  gStyle->SetPalette(1);
  //Draw
  TLine *line = new TLine(7.5,7.5,13.5,13.5);
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++){
      TCanvas canvas;
      canvas.cd();
      canvas.SetTickx();
      canvas.SetTicky();

      hCorrTHR[iref][ilay]->Draw("COLZ");
      line->Draw("same");
      hCorrTHR[iref][ilay]->GetXaxis()->SetTitleOffset(1.2);
      canvas.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    }
  }

  //***********************************************************************************
  //***************************** Dead pixels correlation *****************************
  //***********************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* Dead pixels correlation analysis ***************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  //find max for plots (same for all staves)
  double maxdeadpixcorr = -1.;
  for(int ihist=0; ihist<(int)hmapsDEAD.size(); ihist++)
    for(int ix=1; ix<=hmapsDEAD[ihist]->GetNbinsX();ix++)
      for(int iy=1; iy<=hmapsDEAD[ihist]->GetNbinsY();iy++)
        if(hmapsDEAD[ihist]->GetBinContent(ix,iy)>maxdeadpixcorr){
          maxdeadpixcorr=hmapsDEAD[ihist]->GetBinContent(ix,iy);
        }
  if(maxdeadpixcorr>5000) maxdeadpixcorr=5000;
  double maxlimit = maxdeadpixcorr+0.5*maxdeadpixcorr;
  int nbins = (int)(maxdeadpixcorr+0.5*maxdeadpixcorr-0.9)*1.0;

  TH2I *hCorrDEAD[2][nLayers];

  for(int iref=0; iref<2; iref++)
    for(int ilay=0; ilay<nLayers; ilay++)
      hCorrDEAD[iref][ilay] = new TH2I(Form("hCorr_L%s_refrun_%ld",laynums[ilay*nRuns].c_str(),refrun[iref]), Form("Layer-%s - DeadPix corr. %s - Ref. run: %ld; # Dead Pixel per Chip (run%ld); # Dead Pixel per Chip (runs)",laynums[ilay*nRuns].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str(),refrun[iref],refrun[iref]), nbins, 0.9, maxlimit, nbins, 0.9, maxlimit);

  for(int iref=0; iref<2; iref++){
    for(int ihist=(int)hmapsDEAD.size()-1; ihist>=0; ihist--){
      string hname = hmapsDEAD[ihist]->GetName();
      int ilayer = stoi(hname.substr(hname.find("L")+1,1));
      if(stol(runnumbers[ihist])==refrun[iref]){
        continue; //skip ref run
      }
      for(int ibinx=1; ibinx<=hmapsDEAD[ihist]->GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<=hmapsDEAD[ihist]->GetNbinsY(); ibiny++){
          int data_refrun = (int)hmapsDEAD[posrefrun[ilayer][iref]]->GetBinContent(ibinx,ibiny);
          int data_run = (int)hmapsDEAD[ihist]->GetBinContent(ibinx,ibiny);
          if(data_run && data_refrun)
            hCorrDEAD[iref][ilayer]->Fill(data_refrun, data_run);
          //cout<<fhr_refrun<<"   "<<fhr_run<<endl;
        }
      }
    }
  }

  //Draw
  TLine *linedead = new TLine(1,1,maxlimit,maxlimit);
  linedead->SetLineStyle(2);
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++){
      TCanvas canvas;
      canvas.cd();
      canvas.SetLogx();
      canvas.SetLogy();
      canvas.SetTickx();
      canvas.SetTicky();

      hCorrDEAD[iref][ilay]->Draw("COLZ");
      linedead->Draw("same");
      hCorrDEAD[iref][ilay]->GetXaxis()->SetTitleOffset(1.2);
      canvas.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    }
  }

  //***********************************************************************************
  //***************************** Dead pixel comparison in runs ***********************
  //***********************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* Dead pixel comparison in runs ******************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  cout<<"Reference Run 1: "<<refrun[0]<<endl;
  cout<<"Reference Run 2: "<<refrun[1]<<endl;
  cout<<endl;

  //reset posrefrun run to find new positions
  posrefrun.clear();
  icnt = 0;
  int posrefrun2[nLayers>1 ? 48:nStavesInLay[stoi(laynums[0])]][2];
  for(int ir=0; ir<2; ir++){
    int stop = (nLayers>1) ? 48:nStavesInLay[stoi(laynums[0])];
    for(int is=0; is< stop; is++)
      posrefrun2[is][ir] = -1;
  }
  for(int ihist=0; ihist<(int)hmapsDEADPIX.size(); ihist++){
    string hname = hmapsDEADPIX[ihist]->GetName();
    string runn =  hname.substr(hname.find("run")+3, 6);
    string stvnum = hname.substr(hname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = hname.substr(hname.find("Stv")+3,1);
    }
    int snum = stoi(stvnum);
    int lnum = stoi(hname.substr(hname.find("L")+1,1));
    int staveindex = snum;
    if(nLayers>1){
      staveindex = lnum==0 ? snum: lnum==1? snum+nStavesInLay[0]:snum+nStavesInLay[0]+nStavesInLay[1];
    }
    if(stol(runn)==refrun[0]){ //position of refence run
      posrefrun2[staveindex][0] = ihist;
    }
    if(stol(runn)==refrun[1]){ //position of refence run
      posrefrun2[staveindex][1] = ihist;
    }
  }

  //Compare all the runs (non-empty ones) with the reference runs chosen by the user
  long int first[2][nLayers][100], second[2][nLayers][100], both[2][nLayers][100]; // 100 is just to put a large number of runs that will be never reached
  bool filled[2][nLayers][100];
  vector<array<long int,5>> noisypix;
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++)
      for(int i=0; i<100; i++){
        first[iref][ilay][i]=0; second[iref][ilay][i]=0; both[iref][ilay][i]=0;
        filled[iref][ilay][i] = false;
      }
  }

  for(int iref=0; iref<2; iref++){
    for(int ihist=(int)hmapsDEADPIX.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run

      string hname = hmapsDEADPIX[ihist]->GetName();
      string runn =  hname.substr(hname.find("run")+3, 6);
      string lnum =  hname.substr(hname.find("L")+1,1);
      string stvnum = hname.substr(hname.find("Stv")+3,2);
      if(stvnum.find("_")!=string::npos){
        stvnum = hname.substr(hname.find("Stv")+3,1);
      }
      int snum = stoi(stvnum);
      auto itf = find(runlabel.begin(), runlabel.end(), runn);
      int irun = distance(runlabel.begin(), itf);
      int ilayer = nLayers>1 ? stoi(lnum) : 0;
      int istave = snum;
      if(nLayers>1)
        istave = (lnum=="0") ? snum : (lnum=="1")? snum+nStavesInLay[0]:snum+nStavesInLay[0]+nStavesInLay[1];
      cout<<hmapsDEADPIX[ihist]->GetName()<<endl;
      if(hmapsDEADPIX[ihist]->GetEntries()>1e4){
        cout<<"L"<<lnum<<"_"<<snum<<" - run"<<runn<<" skipped because has more than 10000 entries (probably a bad run)."<<endl;
        continue;
      }
      if(posrefrun2[istave][iref]==-1){
        continue; //skip comparison if no run is found
      }
      if(hmapsDEADPIX[posrefrun2[istave][iref]]->GetEntries()>1e4){
        cout<<"L"<<lnum<<"_"<<snum<<" - run"<<runn<<" ref run skipped because has more than 10000 entries (probably a bad run)."<<endl;
        continue;
      }

      if(runn.find(std::to_string(refrun[iref]))!=string::npos){
        continue;
      }
      noisypix.push_back(CompareTwoRuns(hmapsDEADPIX[posrefrun2[istave][iref]], hmapsDEADPIX[ihist]));
      noisypix[noisypix.size()-1][3] = stol(stavenums[ihist]);
      noisypix[noisypix.size()-1][4] = stol(lnum);
      first[iref][ilayer][irun]+=noisypix[noisypix.size()-1][0];
      second[iref][ilayer][irun]+=noisypix[noisypix.size()-1][1];
      both[iref][ilayer][irun]+=noisypix[noisypix.size()-1][2];
      filled[iref][ilayer][irun] = true;
    }//end loop on histograms
  }//end loop on reference runs

  //Make plot for each layer and for each stave in the root file
  TGraphErrors *ge_nref[2][nLayers];
  TGraphErrors *ge_n2[2][nLayers];
  TGraphErrors *ge_ncom1[2][nLayers];
  TGraphErrors *ge_ncom2[2][nLayers];

  double xshift = 3.;
  double maxbar[2][nLayers];
  double minbar[2][nLayers];

  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++){
      ge_nref[iref][ilay] = new TGraphErrors();
      ge_n2[iref][ilay] = new TGraphErrors();
      ge_ncom1[iref][ilay] = new TGraphErrors();
      ge_ncom2[iref][ilay] = new TGraphErrors();
      maxbar[iref][ilay] = -1.;
      minbar[iref][ilay] = 1e35;

      int ipoint = 0;
      for(int ir=0; ir<100; ir++){//first the oldest data and last the most recent
        if(!filled[iref][ilay][ir]) continue; //skip if not filled
        //first couple of bar on the left
        //int ipoint = (int)noisypix.size()-icomp-1;
        if(!ipoint) xshift=1.;
        else xshift = 3.;
        ge_nref[iref][ilay]->SetPoint(ipoint, ipoint*xshift, (double)both[iref][ilay][ir]/2.+(double)first[iref][ilay][ir]/2.);
        ge_nref[iref][ilay]->SetPointError(ipoint, 0.5, (double)first[iref][ilay][ir]/2.);
        ge_ncom1[iref][ilay]->SetPoint(ipoint, ipoint*xshift, 0.);
        ge_ncom1[iref][ilay]->SetPointError(ipoint, 0.5, (double)both[iref][ilay][ir]/2.);
        if((double)both[iref][ilay][ir]/2.+(double)first[iref][ilay][ir] > maxbar[iref][ilay]) maxbar[iref][ilay] = (double)both[iref][ilay][ir]/2.+(double)first[iref][ilay][ir];

        //second couple of bar on the right
        ge_n2[iref][ilay]->SetPoint(ipoint, ipoint*xshift+1, -(double)both[iref][ilay][ir]/2.-(double)second[iref][ilay][ir]/2.);
        ge_n2[iref][ilay]->SetPointError(ipoint, 0.5, (double)second[iref][ilay][ir]/2.);
        ge_ncom2[iref][ilay]->SetPoint(ipoint, ipoint*xshift+1, 0.);
        ge_ncom2[iref][ilay]->SetPointError(ipoint, 0.5, (double)both[iref][ilay][ir]/2.);
        if(-(double)both[iref][ilay][ir]/2.-(double)second[iref][ilay][ir] < minbar[iref][ilay]) minbar[iref][ilay] = -(double)both[iref][ilay][ir]/2.-(double)second[iref][ilay][ir];

        ipoint++;
      }//end first loop on runs

      //Style
      SetStyle(ge_nref[iref][ilay], kBlue);
      SetStyle(ge_ncom1[iref][ilay], kBlack);
      SetStyle(ge_ncom2[iref][ilay], kBlack);
      SetStyle(ge_n2[iref][ilay], kRed+2);
    }//end loop on layers
  }//end loop on reference runs

  //Legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref[0][0], "#splitline{#dead pix}{ref. run only}", "f");
  leg->AddEntry(ge_n2[0][0], "#splitline{#dead pix}{2nd run only}", "f");
  leg->AddEntry(ge_ncom1[0][0], "#splitline{#dead pix}{both}");

  //Draw plot for each layer
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++){
      TCanvas canvas(Form("mycanvas_%d_%d",ilay,iref), Form("mycanvas_%d_%d",ilay,iref), 1300, 800);
      canvas.SetMargin(0.08, 0.1271, 0.1759, 0.0996);
      canvas.cd();

      //fake histo (just for the axes)
      double x2,y2;
      int npoints = ge_ncom2[iref][ilay]->GetN();
      ge_ncom2[iref][ilay]->GetPoint(ge_ncom2[iref][ilay]->GetN()-1, x2,y2);
      TH1F *hfake = new TH1F("hfakedeadpix","hfakedeadpix", (int)x2+6, -3, x2+3);
      //draw labels on x axis
      int counter = runlabel.size()-1;
      for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
        if(stol(runlabel[counter])==refrun[iref]){
          k-=3;
          counter--;
          continue;
        }
        hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runlabel[counter].c_str()));
        counter--;
      }

      hfake->Draw();
      //canvas->SetLogy();
      hfake->SetTitle(Form("Layer-%s - %s%06ld compared to all",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(), filepath.find("run")==string::npos? "":"run",refrun[iref]));
      ge_nref[iref][ilay]->Draw("P E2 same");
      ge_ncom1[iref][ilay]->Draw("E2 same");
      ge_ncom2[iref][ilay]->Draw("E2 same");
      ge_n2[iref][ilay]->Draw("E2 same");
      hfake->GetYaxis()->SetRangeUser(minbar[iref][ilay]+0.1*minbar[iref][ilay], maxbar[iref][ilay]+0.1*maxbar[iref][ilay]);
      //hfake->GetYaxis()->SetLabelColor(kWhite);
      hfake->GetYaxis()->SetTickLength(0.005);
      hfake->GetYaxis()->SetMaxDigits(4);
      TLine lineref(-0.5, 0, x2+0.5, 0);
      lineref.SetLineColor(kGray-1);
      lineref.SetLineStyle(2);
      lineref.Draw("same");

      //draw legend
      leg->Draw("same");

      canvas.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      delete hfake;
    }//end loop on layers
  }//end loop on reference runs

  //***********************************************************************************
  //*****************************Dead pixel maps **************************************
  //***********************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************************ Dead pixel maps *********************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;
  //Draw hot pixel maps for each layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay));
    cnv.SetTopMargin(0.4);
    cnv.Divide(1,nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])],0,0);
    for(int istv=0; istv<nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])]; istv++){
      TH2I *hHotMap = new TH2I(Form("hHotMap_L%s_Stv%d",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(), istv), "; ; ", 9216,0.5,9216.5,512,0.5,512.5);
      int cnt = 0;
      string badrunalert = "";
      for(int ihist=0; ihist<(int)hmapsDEADPIX.size(); ihist++){
        string hname = hmapsDEADPIX[ihist]->GetName();
        string lnum =  hname.substr(hname.find("L")+1,1);
        bool istherebadrun = false;
        if(nLayers>1){
          if(stoi(lnum)==ilay && stoi(stavenums[ihist]) != istv) continue;
          if(stoi(lnum)!=ilay) continue;
        }
        else if(stoi(stavenums[ihist]) != istv) {
          continue;
        }
        if(hmapsDEADPIX[ihist]->GetEntries()>5000){
          cout<<"run: "<<hname.substr(hname.find("run")+3, 6)<<" skipped for "<<"L"<<lnum<<"_"<<istv<<" because it has more than 5k entries"<<endl;
          istherebadrun = true;
          badrunalert = "There is a bad run! Check logs";
        }
        if(!cnt && !istherebadrun) {hHotMap = (TH2I*)hmapsDEADPIX[ihist]->Projection(1,0);}
        else if(cnt>0 && !istherebadrun){
          TH2I *htemp = (TH2I*)hmapsDEADPIX[ihist]->Projection(1,0);
          hHotMap->Add(htemp);
          delete htemp;
        }

        if(ihist>0 && stavenums[ihist+1]!=stavenums[ihist]) break;
        if(!istherebadrun) cnt++;
      }
      cnv.cd(istv+1);
      cnv.GetPad(istv+1)->SetTickx();
      cnv.GetPad(istv+1)->SetTicky();
      cnv.GetPad(istv+1)->SetRightMargin(0.01);
      if(!istv) cnv.GetPad(istv+1)->SetTopMargin(0.1);

      hHotMap->SetTitle(" ");
      hHotMap->SetMarkerStyle(20);
      hHotMap->SetMarkerSize(0.6);
      hHotMap->SetMarkerColor(kRed);
      hHotMap->SetLineColor(kRed);

      hHotMap->GetXaxis()->SetRangeUser(0.,9216.);
      hHotMap->GetYaxis()->SetRangeUser(0.,512.);

      hHotMap->GetXaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetLabelSize(0.13);
      hHotMap->GetXaxis()->SetLabelSize(0.13);
      if(istv>0){
        hHotMap->GetXaxis()->SetLabelOffset(999);
        hHotMap->GetXaxis()->SetTickLength(0.05);
        hHotMap->GetXaxis()->SetNdivisions(530);
      }
      else{
        hHotMap->GetXaxis()->SetLabelOffset(0.003);
        hHotMap->GetXaxis()->SetNdivisions(530);
        hHotMap->GetXaxis()->SetTickLength(0.05);
      }

      hHotMap->DrawCopy("P X+");

      TLatex lat0;
      lat0.SetNDC();
      lat0.SetTextSize(0.2);
      lat0.DrawLatex(0.5,0.5,badrunalert.c_str());

      TLatex lat;
      lat.SetTextAngle(90);
      lat.SetNDC();
      lat.SetTextSize(0.15);
      lat.DrawLatex(0.04,0.3,Form("Stv%d",istv));

      delete hHotMap;
    }//end loop on staves
    cnv.cd();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.DrawLatex(0.01,0.98,Form("L%s",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str()));

    cnv.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    if(ilay==nLayers-1) cnv.SaveAs(Form("../Plots/ShiftReport24h_THR_%s_%s.pdf]", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

  }

  //_____________________________________________
  //scp to copy the shift report (only on flp6)
  /*string user = (string)gSystem->GetFromPipe("whoami");
  if(user=="its"){
    cout<<endl;
    cout<<"... Copying the Report on eos"<<endl;
    cout<<endl;
    cout<<"Insert the password of user itsshift (ask shift leader if you do not know!):"<<endl;
    gSystem->Exec(Form("scp ../Plots/ShiftReport24h_THR_%s_%s.pdf itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/THR", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }*/


}
