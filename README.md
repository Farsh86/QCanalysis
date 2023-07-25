# QC analysis instructions
In the following some instructions about how to use the scripts are presented. 

## Purpose of this repository
The purpose of the code developed in this repository is to analyse data coming from the commissioning runs on the available ITS layers. The data are downloaded from the CCDB database or read directly from the computer (FLP) on which data are saved after each run. 

## Prepare the environment
In order to be able to run the scripts you will need to build the QualityControl (origin/master) and O2 (origin/dev) software. Be sure to **build the mentioned versions**. The installation is based on *aliBuild*. For installing aliBuild, see the pre-requisites [here](https://alice-doc.github.io/alice-analysis-tutorial/building/custom.html#prerequisites) and follow [these](https://alice-doc.github.io/alice-analysis-tutorial/building/custom.html#get-or-upgrade-alibuild) instructions for the installation. Then, you need to build the packages. To do so you can follow the standard guide [here](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#%F0%9F%9B%A0-build-the-packages). The main things to consider are:
1. [Preparation](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#prepare-your-source-code) of the souce code. Consider the command for O2. For QualityControl you will need to run the command: 
```bash
aliBuild init QualityControl@master --defaults o2
```
2. Check the [pre-requisites](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#check-your-prerequisites-skip-if-using-alidock) with aliDoctor following the command for O2. 
3. [Build](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#build-and-rebuild) the O2 software following the specific command line for O2. For QualityControl software use the following command:
```bash
aliBuild build QualityControl --defaults o2
```
The script *startanalysis.sh* will load the QualityControl environment with:
```bash
alienv load QualityControl/latest
```
Before thisi, make sure that you have set the WORKDIR path of aliBuild with the path of the `sw` folder:
```
export ALIBUILD_WORK_DIR=$HOME/<path-to-sw>/sw
```
and:
```
eval "`alienv shell-helper`"
```

**Be sure that this command works on your machine**. Then, you need to have permissions to execute the script *startanalysis.sh*. To do so, just do on your terminal (in the folder where the script is):
```bash
chmod +x startanalysis.sh
``` 

## Start the analysis
To start the download (optional) and analysis of the data you simply need to run the script *startanalysis.sh*:
```bash
./startanalysis.sh
```
Then you need to follow the instructions that appear on the terminal window. The workflow is described in the following. **Remember that you can exit the program with Ctrl+C at any moment**. 

### Software update
Everytime you run the script *startanalysis.sh*, a pull of this repository is done automatically. In case you performed local modifications (as for the Makefile), they will be kept (git stash).

### Choose what to do (main menu)
When starting the script, it will automatically load the QualityControl modules (after the git repository update). After this, a first menu will appear on the terminal allowing the choice of the following options:

1. Download and analyse data from CCDB
2. Analyse data on flp

You just need to enter the option number (1 or 2). Any other character that will be inserted will be recognized as invalid and the option will have to be retyped. **IMPORTANT! Option 2 doesn't work at present**. 

#### Option 1 - Download and analyse data from CCDB
In case option 1 is chosen, two options are presented
1. Download data from CCDB and perform the analysis
2. Analyse data only (I have already a data sample)
##### Option 1 - Download data from CCDB and perform the analysis
If option 1 (**Download data from CCDB and perform the analysis**) is chosen, the script compiles automatically the code for database access (getObject.cxx) provided that the MakeFile has been properly edited. Then, you need to choose what to download:

1. Fake-hit scan data

only this option is available at present. 
Now you need to enter the layer number. At present only the Inner Barrel data are saved into the database so, the possibile layer numbers are 0, 1 or 2. If "-1" is typed, the data of all the three IB layers are downloaded at the same time. 
Later, you need to choose whether to download Inner Barrel Bottom (IBB) or Inner Barrel Top (IBT) data. Simply type a B (or b) or a T (or t), respectively. The next thing to choose is related to the error plots: if you type "y", the error plots in CCDB path qc/ITS/ITSRawTask/General/ErrorFile (ITSRawTask can be also ITSRawTaskIBB1 or 2 in case IBB is chosen) are also dowloaded and saved in the same file of the other data. 

Then menu with the following options will appear:

1. Enter run numbers
2. Enter timestamps

You need to choose one of the options by typing the corresponding number (any other character will be recognized as invalid and the correct option will have to be retyped).
If option 1 (**Enter run numbers**) is chosen, you will need to type an interval of runs from which to download data (i.e. if the interval is 100 - 200, all the data related to runs having a number between 100 and 200 will be downloaded). **Note that the run number is available in the database from July 3rd, 2019 only!**. For older data you will need to choose option 2 (**Enter timestamps**).  
Instead, if option 2 (**Enter timestamps**) is chosen, you will need to type a time interval from which to download data: the time interval will include a date (DD MM YY) and a time in H24 format (HH MM SS). 
In both cases, the final result will be a *.root* file saved in the directory *Data/* containing the histrograms of each run or timestamp. This is the file to use for the analysis (next step). **IMPORTANT: It's strongly suggested to type run numbers instead of timestamps**. 

After the downloading of the data, a menu with the available analyses is shown (**!IMPORTANT! Other analyses will be added soon**):

1. Fake-hit rate run by run
2. Compare number of noisy pixels between runs
3. Fake-hit rate study with hot pixels masking
4. Fake-hit rate correlation analysis (reference run to be chosen)
5. Error analysis for all runs

If option 1 (**Fake-hit rate run by run**) is chosen, the software will show you the fake-hit rate as a function of the run number for all the chips in a Stave. The ROOT macro is run automatically and a list of (good) files for the analysis is shown. Copy and paste the file name in the input line that will appear. If needed, there is also the possibilty to skip runs within the run interval selected (specify runs separated by comma and without white spaces!). The final plot(s) is(are) automatically saved in *pdf* and *root* format in the repository *Plots/*. Also *gifs* are created (run-by-run animation). 
If option 2 (**Compare noisy pixels between runs**) is chosen, the software will compare the number of noisy pixels run by run. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. Then, a reference run number (**choose a run among the ones in the input file!**) has to be typed. All the runs in the input file will be compared with the chosen reference run. If needed, there is also the possibilty to skip runs within the run interval selected (specify runs separated by comma and without white spaces!). The plot(s) is(are) automatically saved in *pdf* and *root* format in the repository *Plots/*. In particular a summary plot for the layer(s) is created and then, a pdf containing plots of single staves is also added to the output in case of deeper investigations. 
If option 3 (**Fake-hit rate study with hot pixels masking**) is chosen, the fake-hit rate is studied as a function of the number of hot pixel masked. An average of all runs is shown. By default the 100 hottest pixels are masked. Note that in this case a single pixel corresponds to a cluster of 4x4 pixels. Furthermore, the pixels masked in each run are shown in a map for each layer and each stave. If needed, there is also the possibilty to skip runs within the run interval selected (specify runs separated by comma and without white spaces!). The plot(s) is(are) automatically saved in *pdf* and *root* format in the repository *Plots/*. 
If option 4 (**Fake-hit rate correlation analysis (reference run to be chosen**) is chosen, a correlation of the fake-hit rate is studied: a reference run has to be chosen and the chip fake-hit rates of this run are correlated to the rates in all the other runs within the selected interval. If needed, there is also the possibilty to skip runs within the run interval selected (specify runs separated by comma and without white spaces!). The plot(s) is(are) automatically saved in *pdf* and *root* format in the repository *Plots/*. 
For all the options, in case of multiple layers, the plots are created layer by layer automatically. 
If option 5 (**Error analysis for all runs**) is chosen, the software makes a single plot summing also the errors found in the set of runs chosen for the analysis. The colored scale represents the frequency of the different errors. A legend shows the meening of each error ID. The x-axis shows the file ID that is linked to the 24 staves of IBT or IBB. If needed, there is also the possibilty to skip runs within the run interval selected (specify runs separated by comma and without white spaces!). The plot(s) is(are) automatically saved in *pdf* and *root* format in the repository *Plots/*. 

##### Option 2 - Analyse data only (I have already a data sample)
For this option you will need to have already a file containing a set of histograms to analyse (downloaded from CCDB!). If this is the case, choosing this option, you will access directly the analysis menu. Choose a single analysis as explained in the previous paragraph. 

#### Option 2 - Analyse data on flp 
Not available at present. 
<!---This option allows you to analyse data directly on FLP skipping the download from the database. This option can be chosen both if you are working on your computer and on FLP. First, you need to choose on which FLP the data you want to analyse are. **Type only the number of the FLP**. Depending on your choice, if you are not working already on the selected FLP, the script will ask you to type your CERN username. This is used to make a connection via ssh to the selected flp (through lxplus). You just need to insert first your CERN account password and then the one of the "its" account on the selected flp. 
Then, the script asks to type the layer number. After this, an update of the repository is automatically done keeping any local (on flp) modifications. 
Later, the QualityControl environment is loaded on the flp and the script asks to type the path in which data are saved (tipycally: */data/shifts/* or */data/L0_shifts/*). When the path is selected, the script asks you to choose an analysis type:
1. Average stave thresholds run by run
2. Compare dead pixels between runs
3. Option 1 and 2 together
4. Fake-hit rate run by run
5. Compare noisy pixels between runs
6. Option 4 and 5 together
In general, once you have selected one of the options above the script asks you whether you want to prepare and analyse a data sample or you want to analyse an already existing sample. You just need to enter the option you want to perform. 
When you choose to prepare and analyse a data sample, for all the options 1-6 above, the script will ask you to type a run interval. If you choose "starting run = 100" and "final run = 200", all the runs between 100 and 200 will be included in the sample and hence in the analysis. Then, for options 1, 3, 4, 6 you need to type also the number of the first stave in the layer (for example "6" if you are analysing the layer 0 currently under commissioning -> stave numbered from 6 to 11). Now, the details of each option (1-6) will be described. 
##### Option 1 - Average stave thresholds run by run
If this option is chosen, the software will show you the average threshold (stave by stave on the layer) as a function of the run number. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 2 - Compare dead pixels between runs
If this option is chosen, the software will compare the number of dead pixels (i.e. having null threshold) run by run. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. Then, a reference run number (**choose a run among the ones in the input file!**) has to be typed. All the runs in the input file will be compared with the reference run. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 3 - Option 1 and 2 together 
This option will simply do the option 1 and option 2 at the same time since the data sample for the analysis can be the same. 
##### Option 4 - Fake-hit rate run by run
If this option is chosen, the software will show you the fake-hit rate as a function of the run number for all the staves composing the layer. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. The final plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 5 - Compare noisy pixels between runs
If this option is chosen, the software will compare the number of noisy pixels run by run on all the layer (no stave-by-stave analysis at present). The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. Then, a reference run number (**choose a run among the ones in the input file!**) has to be typed. All the runs in the input file will be compared with the reference run. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 6 - Option 4 and 5 together
This option will simply do the option 4 and option 5 at the same time since the data sample for the analysis can be the same.---> 
## Contacts
For any issue, please contact <ivan.ravasenga@cern.ch>.
