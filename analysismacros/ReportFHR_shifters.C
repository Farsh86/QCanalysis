//utilities and constants
#include "inc/constants.h"
#include "inc/utilities.h"

#include <cstdio>

using namespace std;

void DoAnalysis(string filepath, const int nChips, bool isIB);

//
// MAIN
//
void ReportFHR_shifters(){
  //string fpath;
  int nchips=9;
  //cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  string fpath = (string)gSystem->GetFromPipe("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 1"); //take most recent data file

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

  std::freopen(Form("../logs/logFHR_%s.log",localdatetime.c_str()), "w", stdout);

  std::vector<TH2*> hmapsFHR;
  std::vector<THnSparse*> hmapsHIT;
  std::vector<TH2*> hmapsERR;
  std::vector<TH2*> hmapsTRG;

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
  TH2 *h2 = NULL;
  THnSparse *hsparse = NULL;
  bool isfirst = true;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
         && (!obj->InheritsFrom("THnSparse"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D or sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();

    if(objname.find("Stv")==string::npos) h2 = (TH2*)obj;
    else hsparse = (THnSparse*)obj;

    //if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    string timestamp = objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = "nostave";
    if(objname.find("Stv")!=string::npos){ //Hitmaps
      stvnum = objname.substr(objname.find("Stv")+3,2);
      if(stvnum.find("_")!=string::npos){
        stvnum = objname.substr(objname.find("Stv")+3,1);
      }
      stavenums.push_back(stvnum);
      hmapsHIT.push_back(hsparse);
      if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-2] && isfirst){
        isfirst=false;
        continue;
      }
      runlabel.push_back(runnum);
    }
    else if(objname.find("err")!=string::npos){
      hmapsERR.push_back(h2);
    }
    else if(objname.find("trg")!=string::npos){
      hmapsTRG.push_back(h2);
    }
    else{// FHR maps
      hmapsFHR.push_back(h2);
      timestamps.push_back(timestamp);
      runnumbers.push_back(runnum);
      laynums.push_back(laynum);
      nTimes++;
      if(nTimes>1 && laynum==laynums[laynums.size()-2])
        nRuns++;
      else nRuns=1;
    }
    //cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
  }

  const int nLayers = (int)hmapsFHR.size()==nRuns ? 1 : stoi(laynums[laynums.size()-1])+1;

  //**************************************************************************************
  //************************* Fake-hit rate vs run number ********************************
  //**************************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* Fake-hit vs run number *************************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;
  TGraph *trend[nLayers][100];
  int ilayer=nLayers-1;
  for(int ihist=(int)hmapsFHR.size()-1; ihist>=0; ihist--){
    for(int ibiny=1; ibiny<=hmapsFHR[ihist]->GetNbinsY(); ibiny++){
      trend[ilayer][ibiny-1] = new TGraph();
    }
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
  }
  ilayer=nLayers-1;
  //TH1F *hproj = new TH1F();
  TH1F *hproj;
  string histname = hmapsFHR[0]->GetName();
  int irun=0;
  for(int ihist=(int)hmapsFHR.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    for(int ibiny=1; ibiny<=hmapsFHR[ihist]->GetNbinsY(); ibiny++){//loop on y bins (staves)
      TH1D *hproj = (TH1D*)hmapsFHR[ihist]->ProjectionX("proj",ibiny,ibiny); //single stave
      trend[ilayer][ibiny-1]->SetName(Form("gr_L%s_stave%d",laynums[ihist].c_str(),ibiny-1));
      int deadchips = 0;
      for(int ibinx=1; ibinx<=hmapsFHR[ihist]->GetNbinsX(); ibinx++){//evaluate the number of disabled chips
        if(hmapsFHR[ihist]->GetBinContent(ibinx,ibiny)<1e-20)
          deadchips++;
      }
      if(deadchips>0)
        cout<<"Layer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" --> Chips active:"<<nChips-deadchips<<endl;

      if(deadchips!=nChips)
        trend[ilayer][ibiny-1]->SetPoint(irun, irun, hproj->Integral()/(nChips-deadchips));
      else
        trend[ilayer][ibiny-1]->SetPoint(irun, irun, 0.);

      if((ibiny-1)<hmapsFHR[ihist]->GetNbinsY()/2)
        SetStyle(trend[ilayer][ibiny-1], col[ibiny-1], 24);
      else
        SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-hmapsFHR[ihist]->GetNbinsY()/2], 26);
    }
    irun++;
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        irun=0;
        ilayer--;
      }
  }

  int npoints = trend[0][0]->GetN();
  TH1F *hfake = new TH1F("hfake", "; Run; Fake-hit Rate (/event/pixel)", npoints, -0.5, (double)npoints-0.5);

  for(int ir=0; ir<(int)runnumbers.size()/nLayers; ir++)
      hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[(int)runnumbers.size()/nLayers-1-ir])));


  //Draw
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogy();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.0988,0.1,0.194,0.0993);
    TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
    for(int istave=0; istave<hmapsFHR[ilay*nRuns]->GetNbinsY(); istave++)
      leg->AddEntry(trend[ilay][istave], Form("Stv%d",istave), "p");
    hfake->GetYaxis()->SetRangeUser(1e-14, 1e-3);
    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->SetTitle(Form("Layer-%s, %s",laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
    hfake->Draw();
    for(int istave=0; istave<hmapsFHR[ilay*nRuns]->GetNbinsY(); istave++)
      trend[ilay][istave]->Draw("P same");
    leg->Draw("same");

    if(!ilay) canvas->SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf[", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));

    delete canvas;
    delete leg;
  }

  //**************************************************************************************
  //************************* Fake-hit rate vs masked hot pixels *************************
  //**************************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* Fake-hit rate vs masked hot pixels *************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;
  //Calculate number of triggers for each run
  cout<<endl;
  cout<<"... Extract the number of triggers (events) from for each run"<<endl;
  vector<double> ntrig;
  for(int ir=0; ir<nRuns; ir++){
    double fhr_run = hmapsFHR[ir]->GetBinContent(1,1);
    int stavefound = 0;
    int chipfound = 0;
    if(fhr_run<1e-20){
      for(int ibinx=1; ibinx<=hmapsFHR[ir]->GetNbinsX(); ibinx++){
        for(int ibiny=hmapsFHR[ir]->GetNbinsY(); ibiny>=1; ibiny--){
          fhr_run = hmapsFHR[ir]->GetBinContent(ibinx,ibiny);
          if(fhr_run>1e-20) {
            stavefound = ibiny-1;
            chipfound = ibinx-1;
            break;
          }
        }
        if(fhr_run>1e-20) break;
      }
    }
    if(fhr_run<1e-20){
      fhr_run = -1.;
      ntrig.push_back(-1.);
      cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers (ignored in the calculation of the average fhr)"<<endl;
      continue;
      //cout<<"INVALID FHR... setting it to -1"<<endl;
    }
    hmapsHIT[stavefound*nRuns+ir] -> GetAxis(0) ->SetRange(1+1024*chipfound, 1024+1024*chipfound);
    TH2F *hprojsparse = (TH2F*)hmapsHIT[stavefound*nRuns+ir]->Projection(1,0);
    hmapsHIT[stavefound*nRuns+ir] -> GetAxis(0) ->SetRange(1, 9216);//reset the range

    double hits_chip = hprojsparse->Integral(1,1024,1,512);
    if(hits_chip/(512.*1024.*fhr_run) < 1e-15){//to avoid bad runs
      ntrig.push_back(-1.);
      cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers (ignored in the calculation of the average fhr)"<<endl;
    }
    else{
      ntrig.push_back(hits_chip/(512.*1024.*fhr_run));
      cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers"<<endl;
    }
    delete hprojsparse;
  }

  //Start masking hottest pixels for each stave in each run, Fill also the histo with the hot pixel maps for each stave
  cout<<endl;
  cout<<"... Analysing FHR with (hot) pixel masking (Making also hot pixel map)"<<endl;
  vector<array<float,nMasked+1>> fhrall;
  TH2F *hHotMap[nLayers][20];
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int istave=0; istave<20; istave++)
      hHotMap[ilay][istave] = new TH2F(Form("hHotMap_L%s_Stv%d",laynums[ilay*nRuns].c_str(), istave), "; ; ", 2304,-0.5,9215.5, 128,-0.5,511.5);//4x4 pixel cells
  irun = nRuns-1;
  for(int ihist=(int)hmapsHIT.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    int nchipsactive = GetNchipsActive(hmapsHIT[ihist],nChips);
    string hname = hmapsHIT[ihist]->GetName();
    string layn = hname.substr(hname.find("L")+1,1);
    string runnum =  hname.substr(hname.find("run")+3, 6);
    if(nchipsactive<nChips)
      cout<<"Layer "<<layn<<" Stave "<<stavenums[ihist]<<" Run: "<<runnum<<" --> Chips active:"<<nchipsactive<<endl;
    fhrall.push_back(GetFHRwithMasking(hmapsHIT[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : stoi(layn)][stoi(stavenums[ihist])]));
    irun--;
    if(ihist>0){
      if(stavenums[ihist-1]!=stavenums[ihist]){
        irun=nRuns-1;
      }
    }
  }

  //special binning
  double binstart = 0.4;
  double binsmasked[nMasked+2];
  binsmasked[0] = 0.1;
  binsmasked[1] = 0.5;
  for(int i=2; i<=nMasked+1; i++){
    binsmasked[i] = binsmasked[i-1]+1.;
  }

  TH2F *hFhrStv[nLayers][100];
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int is=0; is<nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])]; is++)
      hFhrStv[ilay][is] = new TH2F(Form("h2FhrStv_%s_%d", laynums[ilay*nRuns].c_str(),is), Form("Layer-%s - Stave-%d; # Hot Pixel Masked;Run", laynums[ilay*nRuns].c_str(),is),nMasked+1, binsmasked, nRuns, 0.5, nRuns+0.5);
  //Fill histogram
  ilayer = nLayers-1;
  cout<<"nLayers: "<<nLayers<<endl;
  int istave = nStavesInLay[nLayers>1 ? ilayer:stoi(laynums[0])]-1;
  irun=0;
  for(int i=0; i<(int)fhrall.size(); i++){

    for(int ifhr=0; ifhr<(int)fhrall[i].size(); ifhr++){
      hFhrStv[ilayer][istave]->SetBinContent(ifhr+1, irun+1, fhrall[i][ifhr]);
    }
    irun++;
    if(i<(int)fhrall.size()-1){
      if(stavenums[fhrall.size()-i-2]!=stavenums[fhrall.size()-i-1]){
        istave--;
        irun=0;
      }
      if(istave<0){
        ilayer--;
        istave = nStavesInLay[nLayers>1 ? ilayer:stoi(laynums[0])]-1;
      }
    }
  }

  //Make FHR (averaged on all runs) vs #masked pix for all staves in a layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay));
    cnv.cd();
    cnv.SetLogy();
    cnv.SetLogx();
    cnv.SetTickx();
    cnv.SetTicky();
    cnv.SetMargin(0.0988,0.1,0.1,0.0993);
    TH1F *hframe = cnv.DrawFrame(0.1,7e-15,3*(nMasked),1e-3,Form("Layer %s - Average FHR %s; # Hot Pixel Clusters masked ; FHR (/event/pixel)",laynums[ilay*nRuns].c_str(),filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
    //legend
    TLegend leg(0.904, 0.127,0.997,0.898);
    for(int is=0; is<nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])];is++){
      TH1F *proj = (TH1F*)hFhrStv[ilay][is]->ProjectionX(Form("proj_%d%d",ilay,is));
      int runswohits = GetNrunsWOhits(hFhrStv[ilay][is]);
      //cout<<runswohits<<endl;
      proj->Scale(1./(nRuns-runswohits)); //Divide by the number of runs minus the ones without hits
      SetStyle(proj, col[is<nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])]/2 ? is : is-nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])]/2],is<nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])]/2 ? 24:26);
      proj->Draw("PL same");
      leg.AddEntry(proj, Form("Stv%d",is),"pl");
    }
    leg.Draw("same");
    cnv.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  }

  //Draw hot pixel maps for each layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay));
    cnv.SetTopMargin(0.4);
    cnv.Divide(1,nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])],0,0);
    for(int istave=0; istave<nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])]; istave++){
      hHotMap[ilay][istave]->SetMarkerStyle(20);
      hHotMap[ilay][istave]->SetMarkerSize(0.6);
      hHotMap[ilay][istave]->SetMarkerColor(kRed);
      hHotMap[ilay][istave]->SetLineColor(kRed);

      cnv.cd(istave+1);
      cnv.GetPad(istave+1)->SetTickx();
      cnv.GetPad(istave+1)->SetTicky();
      cnv.GetPad(istave+1)->SetRightMargin(0.01);
      if(!istave) cnv.GetPad(istave+1)->SetTopMargin(0.1);

      hHotMap[ilay][istave]->Draw("P X+");
      hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.005);
      hHotMap[ilay][istave]->GetYaxis()->SetTickLength(0.005);
      hHotMap[ilay][istave]->GetYaxis()->SetLabelSize(0.13);
      hHotMap[ilay][istave]->GetXaxis()->SetLabelSize(0.13);
      if(istave>0){
        hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(999);
        hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
        hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
      }
      else{
        hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(0.003);
        hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
        hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
      }

      TLatex lat;
      lat.SetTextAngle(90);
      lat.SetNDC();
      lat.SetTextSize(0.15);
      lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
    }
    cnv.cd();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.DrawLatex(0.01,0.98,Form("L%s",laynums[ilay*nRuns].c_str()));

    cnv.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  }

  //**************************************************************************************
  //************************* Noisy pixel comparison  ************************************
  //**************************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"************************* Noisy pixel comparison *************************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  //Choose two random reference runs
  const int pos1 = (int)nRuns*0.25;//pos of ref run 1
  const int pos2 = (int)nRuns*0.75;//pos of ref run 2;
  long int refrun[2];
  vector<array<int,2>> posrefrun;
  array<int,2> tmparray;
  refrun[0] = stol(runlabel[pos1]);
  refrun[1] = stol(runlabel[pos2]);
  int icnt = 0;
  for(int ihist=0; ihist<(int)hmapsHIT.size(); ihist++){
    string hname = hmapsHIT[ihist]->GetName();
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
  //Compare all the runs (non-empty ones) with the reference runs chosen above
  istave = 0;
  for(int ilay=0; ilay<nLayers; ilay++)
    istave+=nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])];
  istave--;
  int stavetotal = istave;
  ilayer = nLayers-1;
  long int first[2][nLayers][nRuns], second[2][nLayers][nRuns], both[2][nLayers][nRuns]; //2 because we have 2 ref runs
  irun=0;
  vector<array<long int,5>> noisypix;
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++)
      for(int i=0; i<nRuns; i++){
        first[iref][ilay][i]=0; second[iref][ilay][i]=0; both[iref][ilay][i]=0;
      }
  }
  for(int iref=0; iref<2; iref++){
    irun=0;
    istave = 0;
    for(int ilay=0; ilay<nLayers; ilay++)
      istave+=nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])];
    istave--;
    //loop on histos
    for(int ihist=(int)hmapsHIT.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
      string hname = hmapsHIT[ihist]->GetName();
      string runn =  hname.substr(hname.find("run")+3, 6);
      string lnum =  hname.substr(hname.find("L")+1,1);
      ilayer = stoi(lnum);
      if(runn.find(std::to_string(refrun[iref]))!=string::npos){
        if(ihist>0){// in case ref run is the first into the list of runs
          if(stavenums[ihist-1]!=stavenums[ihist]){
            istave--;
            irun=0;
          }
        }
        continue;
      }
      //if(!hmaps[ihist]->GetEntries()) continue; // do not compare ref run with empty run (= empty maps)

      noisypix.push_back(CompareNoisyTwoRuns(hmapsHIT[posrefrun[istave][iref]], hmapsHIT[ihist]));
      noisypix[noisypix.size()-1][3] = stol(stavenums[ihist]);
      noisypix[noisypix.size()-1][4] = stol(lnum);
      //cout<<noisypix[noisypix.size()-1][0]<<"  "<<noisypix[noisypix.size()-1][1]<<"  "<<noisypix[noisypix.size()-1][2]<<"  "<<noisypix[noisypix.size()-1][3]<<"  "<<noisypix[noisypix.size()-1][4]<<endl;
      first[iref][ilayer][irun]+=noisypix[noisypix.size()-1][0];
      second[iref][ilayer][irun]+=noisypix[noisypix.size()-1][1];
      both[iref][ilayer][irun]+=noisypix[noisypix.size()-1][2];
      //cout<<"ilayer: "<<ilayer<<"  irun: "<<irun<<"  #noisyincommon: "<<noisypix[noisypix.size()-1][2]<<endl;
      irun++;

      if(ihist>0){
        if(stavenums[ihist-1]!=stavenums[ihist]){
          istave--;
          irun=0;
        }
      }

    }//end loop on histograms
  }//end loop on ref runs


  //Make plot for each layer and for each stave in the root file
  TGraphErrors *ge_nref[2][nLayers];
  TGraphErrors *ge_n2[2][nLayers];
  TGraphErrors *ge_ncom1[2][nLayers];
  TGraphErrors *ge_ncom2[2][nLayers];
  double xshift = 3.;
  double max[2][nLayers];
  double min[2][nLayers];
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++){
      ge_nref[iref][ilay] = new TGraphErrors();
      ge_n2[iref][ilay] = new TGraphErrors();
      ge_ncom1[iref][ilay] = new TGraphErrors();
      ge_ncom2[iref][ilay] = new TGraphErrors();
      max[iref][ilay] = -1.;
      min[iref][ilay] = 1e35;

      for(int ir=0; ir<nRuns-1; ir++){//first the older data and last the most recent
        //first couple of bar on the left
        //int ipoint = (int)noisypix.size()-icomp-1;
        if(!ir) xshift=1.;
        else xshift = 3.;
        ge_nref[iref][ilay]->SetPoint(ir, ir*xshift, (double)both[iref][ilay][ir]/2.+(double)first[iref][ilay][ir]/2.);
        ge_nref[iref][ilay]->SetPointError(ir, 0.5, (double)first[iref][ilay][ir]/2.);
        ge_ncom1[iref][ilay]->SetPoint(ir, ir*xshift, 0.);
        ge_ncom1[iref][ilay]->SetPointError(ir, 0.5, (double)both[iref][ilay][ir]/2.);
        if((double)both[iref][ilay][ir]/2.+(double)first[iref][ilay][ir] > max[iref][ilay]) max[iref][ilay] = (double)both[iref][ilay][ir]/2.+(double)first[iref][ilay][ir];

        //second couple of bar on the right
        ge_n2[iref][ilay]->SetPoint(ir, ir*xshift+1, -(double)both[iref][ilay][ir]/2.-(double)second[iref][ilay][ir]/2.);
        ge_n2[iref][ilay]->SetPointError(ir, 0.5, (double)second[iref][ilay][ir]/2.);
        ge_ncom2[iref][ilay]->SetPoint(ir, ir*xshift+1, 0.);
        ge_ncom2[iref][ilay]->SetPointError(ir, 0.5, (double)both[iref][ilay][ir]/2.);
        if(-(double)both[iref][ilay][ir]/2.-(double)second[iref][ilay][ir] < min[iref][ilay]) min[iref][ilay] = -(double)both[iref][ilay][ir]/2.-(double)second[iref][ilay][ir];
      }//end first loop on runs

      //Style
      SetStyle(ge_nref[iref][ilay], kBlue);
      SetStyle(ge_ncom1[iref][ilay], kBlack);
      SetStyle(ge_ncom2[iref][ilay], kBlack);
      SetStyle(ge_n2[iref][ilay], kRed+2);
    }//end loop on layers
  }//end loop on iref

  //Legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref[0][0], "#splitline{#noisy pix}{ref. run only}", "f");
  leg->AddEntry(ge_n2[0][0], "#splitline{#noisy pix}{2nd run only}", "f");
  leg->AddEntry(ge_ncom1[0][0], "#splitline{#noisy pix}{both}");

  //Draw plot for each layer
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++){
      TCanvas canvas(Form("mycanvas_%d",ilay), Form("mycanvas_%d",ilay), 1300, 800);
      canvas.SetMargin(0.08, 0.1271, 0.1759, 0.0996);
      canvas.cd();

      //fake histo (just for the axes)
      double x2,y2;
      ge_ncom2[iref][ilay]->GetPoint(ge_ncom2[iref][ilay]->GetN()-1, x2,y2);
      TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
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
      hfake->SetTitle(Form("Layer-%s - %s%06ld compared to all (#hits/pixel>2)",laynums[ilay*nRuns].c_str(), filepath.find("run")==string::npos? "":"run",refrun[iref]));
      ge_nref[iref][ilay]->Draw("P E2 same");
      ge_ncom1[iref][ilay]->Draw("E2 same");
      ge_ncom2[iref][ilay]->Draw("E2 same");
      ge_n2[iref][ilay]->Draw("E2 same");
      hfake->GetYaxis()->SetRangeUser(min[iref][ilay]+0.1*min[iref][ilay], max[iref][ilay]+0.1*max[iref][ilay]);
      //hfake->GetYaxis()->SetLabelColor(kWhite);
      hfake->GetYaxis()->SetTickLength(0.005);
      hfake->GetYaxis()->SetMaxDigits(4);
      TLine lineref(-0.5, 0, x2+0.5, 0);
      lineref.SetLineColor(kGray-1);
      lineref.SetLineStyle(2);
      lineref.Draw("same");

      //draw legend
      leg->Draw("same");

      //Save
      canvas.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
      delete hfake;
    }//end loop on layers
  }//end loop on iref

  //**************************************************************************************
  //************************* Fake-hit rate correlation **********************************
  //**************************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"******************** Chip-by-chip fake-hit rate correlation study ********************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  cout<<"Reference Run 1: "<<refrun[0]<<endl;
  cout<<"Reference Run 2: "<<refrun[1]<<endl;
  cout<<endl;
  icnt = 0;
  vector<array<int,2>> posrefrunfhr;
  for(int ihist=0; ihist<(int)hmapsFHR.size(); ihist++){
    string hname = hmapsFHR[ihist]->GetName();
    string runn =  hname.substr(hname.find("run")+3, 6);
    if(stol(runn)==refrun[icnt]){ //position of refence run
      tmparray[icnt] = ihist;
      icnt++;
      if(icnt==2) {
        posrefrunfhr.push_back(tmparray);
        icnt=0;
      }
    }
  }

  TH2D *hCorr[2][nLayers];// 2 because we have 2 reference runs here
  //bins
  double exponent = -14.;
  double base = 10;
  double bins[100];
  bins[0] = 1e-14;
  for(int i=1; i<100; i++){
    bins[i] = i%9==0 ?TMath::Power(base, exponent) : bins[i-1] + TMath::Power(base, exponent);
    if((i+1)%9==0) {
      exponent++;
    }
  }
  for(int iref=0; iref<2; iref++)
    for(int ilay=0; ilay<nLayers; ilay++)
      hCorr[iref][ilay] = new TH2D(Form("hCorr_L%s_refrun_%ld",laynums[ilay*nRuns].c_str(),refrun[iref]), Form("Layer-%s - FHR corr. %s - Ref. run: %ld; FHR (run%ld); FHR (runs)",laynums[ilay*nRuns].c_str(),filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str(),refrun[iref],refrun[iref]), 99, bins, 99, bins);

  for(int iref=0; iref<2; iref++){
    ilayer=nLayers-1;
    for(int ihist=(int)hmapsFHR.size()-1; ihist>=0; ihist--){
      if(stol(runnumbers[ihist])==refrun[iref]){
        if(ihist>0)
          if(laynums[ihist-1]!=laynums[ihist]){// in case the ref run is the first in the list of all layers
            ilayer--;
          }
        continue; //skip ref run
      }
      for(int ibinx=1; ibinx<=hmapsFHR[ihist]->GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<=hmapsFHR[ihist]->GetNbinsY(); ibiny++){
          double fhr_refrun = hmapsFHR[posrefrunfhr[ilayer][iref]]->GetBinContent(ibinx,ibiny);
          double fhr_run = hmapsFHR[ihist]->GetBinContent(ibinx,ibiny);
          hCorr[iref][ilayer]->Fill(fhr_refrun, fhr_run);
          //cout<<fhr_refrun<<"   "<<fhr_run<<endl;
        }
      }
      if(ihist>0)
        if(laynums[ihist-1]!=laynums[ihist]){
          ilayer--;
        }
    }
  }//end loop on iref


  gStyle->SetPalette(1);
  //Draw
  TLine line2(1e-14,1e-14,1e-3,1e-3);
  for(int iref=0; iref<2; iref++){
    for(int ilay=0; ilay<nLayers; ilay++){
      TCanvas canvas;
      canvas.cd();
      canvas.SetLogy();
      canvas.SetLogx();
      canvas.SetTickx();
      canvas.SetTicky();

      hCorr[iref][ilay]->Draw("COLZ");
      line2.Draw("same");
      hCorr[iref][ilay]->GetXaxis()->SetTitleOffset(1.2);
      canvas.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
    }
  }

  //***********************************************************************
  //***************** Error plots *****************************************
  //***********************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"*********************************** Error plots **************************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  //sum all the histos in a single histogram (for summary plot) for each layer
  TH2D *hSummaryErr = (TH2D*)hmapsERR[0]->Clone("hSummaryErr");
  hSummaryErr->Reset();
  for(int iplot=0; iplot<(int)hmapsERR.size(); iplot++){
    if(hmapsERR[iplot]->GetEntries()) hSummaryErr->Add(hmapsERR[iplot]);
  }

  //Make plots with Error IDs vs Run for each layer
  TGraph *trendErr[hSummaryErr->GetNbinsY()];

  int ir = 0;
  double maxErr = -1.;
  for(int iplot=(int)hmapsERR.size()-1; iplot>=0; iplot--){
    TH1D *hproj = (TH1D*)hmapsERR[iplot]->ProjectionY(Form("hmapsERR_%d",iplot));
    for(int ibin=1; ibin<=hproj->GetNbinsX(); ibin++){
      if(ir==0){
        trendErr[ibin-1] = new TGraph();
        trendErr[ibin-1]->SetName(Form("gr_errID%d",ibin));
        SetStyle(trendErr[ibin-1], col[ibin<=10?ibin-1:ibin<=20?ibin-11:ibin-21], ibin<=10?24:ibin<=20?25:26);
      }
      trendErr[ibin-1]->SetPoint(ir,ir, hproj->GetBinContent(ibin));
      if(hproj->GetBinContent(ibin)>maxErr)
        maxErr=hproj->GetBinContent(ibin);
    }
    delete hproj;
    ir++;
  }

  //Draw summary plot for errors
  TCanvas canvasErr;
  canvasErr.cd();
  canvasErr.SetTickx();
  canvasErr.SetTicky();
  canvasErr.SetLogz();
  canvasErr.SetMargin(0.0988,0.2,0.194,0.0993);
  canvasErr.SetRightMargin(0.15);
  hSummaryErr->SetTitle(Form("Errors IB, %s",filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  hSummaryErr->Draw("COLZ");
  hSummaryErr->SetMinimum(1);
  hSummaryErr->GetXaxis()->SetLabelSize(0.045);
  hSummaryErr->GetYaxis()->SetLabelSize(0.045);
  hSummaryErr->GetZaxis()->SetLabelSize(0.045);
  hSummaryErr->GetXaxis()->SetTitleSize(0.05);
  hSummaryErr->GetYaxis()->SetTitleSize(0.05);
  hSummaryErr->GetYaxis()->SetTitleOffset(0.7);
  hSummaryErr->GetZaxis()->SetTitleSize(0.05);
  hSummaryErr->GetZaxis()->SetTitleOffset(0.9);
  canvasErr.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  //Draw trends for errors
  int npointsErr = trendErr[0]->GetN();
  TH1F *hfakeErr = new TH1F("hfakeErr", "; Run; # Errors", npointsErr, -0.5, (double)npointsErr-0.5);
  for(int ir=0; ir<(int)runnumbers.size()/nLayers; ir++)
      hfakeErr->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[(int)runnumbers.size()/nLayers-1-ir])));

  TLegend *legErr = new TLegend(0.904, 0.197,0.997,0.898);
  legErr->SetHeader("Error IDs");
  legErr->SetNColumns(2);
  for(int iid=1; iid<=hSummaryErr->GetNbinsY();iid++)
    legErr->AddEntry(trendErr[iid-1], Form("%d",iid), "p");

  TCanvas canvasErr2;
  canvasErr2.cd();
  canvasErr2.SetTickx();
  canvasErr2.SetTicky();
  canvasErr2.SetLogy();
  canvasErr2.SetMargin(0.0988,0.1,0.194,0.0993);

  hfakeErr->GetXaxis()->SetTitleOffset(2.8);
  hfakeErr->SetTitle(Form("IB, Error trends %s", filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  hfakeErr->GetYaxis()->SetRangeUser(1, 10*maxErr);
  hfakeErr->GetXaxis()->SetTitleOffset(2.8);
  hfakeErr->Draw();
  for(int iid=1; iid<=hSummaryErr->GetNbinsY();iid++){
    trendErr[iid-1]->Draw("P same");
  }
  legErr->Draw("same");

  canvasErr2.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  cout<<"here 4"<<endl;
  //*******************************************************************************************
  //********************************* Trigger and Flags plots *********************************
  //*******************************************************************************************
  cout<<endl; cout<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<"********************************** Trigger&flags plots *******************************"<<endl;
  cout<<"**************************************************************************************"<<endl;
  cout<<endl;

  //sum all the histos in a single histogram (for summary plot) for each layer --> Trigger summary
  TH2D *hSummaryTrg = (TH2D*)hmapsTRG[0]->Clone("hSummaryTrg");
  for(int iplot=1; iplot<(int)hmapsTRG.size(); iplot++){
    hSummaryTrg->Add(hmapsTRG[iplot]);
  }

  //Make plots with Error IDs vs Run for each layer
  TGraph *trendTrg[hSummaryTrg->GetNbinsY()];

  ir = 0;
  double maxTrg = -1.;
  for(int iplot=(int)hmapsTRG.size()-1; iplot>=0; iplot--){
    TH1D *hproj = (TH1D*)hmapsTRG[iplot]->ProjectionY(Form("herrTrg_%d",iplot));
    for(int ibin=1; ibin<=hSummaryTrg->GetNbinsY(); ibin++){
      if(ir==0){
        trendTrg[ibin-1] = new TGraph();
        trendTrg[ibin-1]->SetName(Form("gr_trgID%d",ibin));
        SetStyle(trendTrg[ibin-1], col[ibin<=10?ibin-1:ibin-11], ibin<=10?24:25);
      }
      trendTrg[ibin-1]->SetPoint(ir,ir, hproj->GetBinContent(ibin));
      if(hproj->GetBinContent(ibin)>maxTrg)
        maxTrg=hproj->GetBinContent(ibin);
    }
    delete hproj;
    ir++;
  }

  //Draw summary plot
  TCanvas canvasTrg;
  canvasTrg.cd();
  canvasTrg.SetTickx();
  canvasTrg.SetTicky();
  canvasTrg.SetLogz();
  canvasTrg.SetMargin(0.18,0.2,0.194,0.0993);
  canvasTrg.SetRightMargin(0.15);
  hSummaryTrg->SetTitle(Form("Trigger & Flags IB, %s", filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  hSummaryTrg->Draw("colz");
  //hSummary[ilay]->GetXaxis()->SetNdivisions(530);
  //hSummary[ilay]->GetYaxis()->SetNdivisions(516);
  hSummaryTrg->GetXaxis()->SetLabelSize(0.045);
  hSummaryTrg->GetYaxis()->SetLabelSize(0.045);
  hSummaryTrg->GetZaxis()->SetLabelSize(0.045);
  hSummaryTrg->GetXaxis()->SetTitleSize(0.05);
  hSummaryTrg->GetYaxis()->SetTitleSize(0.05);
  hSummaryTrg->GetYaxis()->SetTitleOffset(0.9);
  hSummaryTrg->GetZaxis()->SetTitleSize(0.05);
  hSummaryTrg->GetZaxis()->SetTitleOffset(0.9);

  canvasTrg.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));

  //Draw trends
  int npointstrg = trendTrg[0]->GetN();
  TH1F *hfakeTrg = new TH1F("hfakeTrg", "; Run; # Errors", npointstrg, -0.5, (double)npointstrg-0.5);
  for(int ir=0; ir<(int)runnumbers.size()/nLayers; ir++)
      hfakeTrg->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[(int)runnumbers.size()/nLayers-1-ir])));

  TLegend *legtrg = new TLegend(0.904, 0.197,0.997,0.898);
  legtrg->SetHeader("IDs");
  for(int iid=1; iid<=hSummaryTrg->GetNbinsY();iid++)
    legtrg->AddEntry(trendTrg[iid-1], Form("%s",hSummaryTrg->GetYaxis()->GetBinLabel(iid)), "p");

    TCanvas canvasTrg2;
    canvasTrg2.cd();
    canvasTrg2.SetTickx();
    canvasTrg2.SetTicky();
    canvasTrg2.SetLogy();
    canvasTrg2.SetMargin(0.0988,0.1,0.194,0.0993);

    hfakeTrg->GetXaxis()->SetTitleOffset(2.8);
    hfakeTrg->SetTitle(Form("IB, Trigger & Flag trends %s", filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
    hfakeTrg->GetYaxis()->SetRangeUser(1, 10*maxTrg);
    hfakeTrg->GetXaxis()->SetTitleOffset(2.8);
    hfakeTrg->Draw();
    for(int iid=1; iid<=hSummaryTrg->GetNbinsY();iid++){
      trendTrg[iid-1]->Draw("P same");
    }
    legtrg->Draw("same");
    canvasTrg2.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
    canvasTrg2.SaveAs(Form("../Plots/ShiftReport24h_FHR_%s_%s.pdf]", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));

  //scp to copy the shift report (only on flp6)
  /*string user = (string)gSystem->GetFromPipe("whoami");
  if(user=="its"){
    cout<<endl;
    cout<<"... Copying the Report on eos"<<endl;
    cout<<endl;
    cout<<"Insert the password of user itsshift (ask shift leader if you do not know!):"<<endl;
    gSystem->Exec(Form("scp ../Plots/ShiftReport24h_FHR_%s_%s.pdf itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/FHR", localdatetime.c_str(), filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  }*/


}
