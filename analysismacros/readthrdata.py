import numpy
import re
import math
import gzip
import sys
from array import array
from ROOT import TTree, TFile

first_arg = sys.argv[1]
second_arg = sys.argv[2]

def readdata(inputfile):
    fl = gzip.GzipFile(inputfile,"r")
    data = numpy.load(fl)
    return data

def insertrun(word):
    rn = input('Type the {} run number of a threshold scan: '.format(word))
    return rn

def askruns():
    run1 = 0
    run2 = 0
    #input run1 and run2
    rn1 = insertrun('starting')
    try:
        run1 = int(rn1)
    except ValueError:
        print("Invalid number")
        rn1 = insertrun('starting')
        run1 = int(rn1)

    rn2 = insertrun('final')
    try:
        run2 = int(rn2)
    except ValueError:
        print("Invalid number")
        rn2 = insertrun('starting')
        run2 = int(rn2)

    if(run1>run2):
        print("Error: Starting run number must be lower than final run number, retype a correct interval")
        askruns(run1, run2)
    return [run1, run2]


########
#MAIN
########
def main(barrel=first_arg, layer=second_arg):
    runs = askruns()
    staven = input('Type the number of the first stave in the (half-)layer: ')

    runnum = array('i', [0])
    currstave = array('i', [0])
    currchip = array('i', [0])
    avgthr = array('f', [0.])
    col = array('i', [0])
    row = array('i', [0])

    f = open("datatoanalyse.txt", "r") ## file with all paths

    #open file with all the paths
    ftree = TFile.Open("../Data/{}_Layer{}_thresholds_tree_from_run{}_to_run{}.root".format(barrel, layer, runs[0], runs[1]), "recreate") ## file containing the tree
    roottree1 = TTree("thrscan_thr", "thrscan_thr");
    roottree1.Branch("runnum", runnum, "runnum/I")
    roottree1.Branch("stavenum", currstave, "stavenum/I")
    roottree1.Branch("chipnum", currchip, "chipnum/I")
    roottree1.Branch("avgchipthr", avgthr, "avgchipthr/F")
    roottree2 = TTree("thrscan_deadpix", "thrscan_deadpix");
    roottree2.Branch("runnum", runnum, "runnum/I")
    roottree2.Branch("stavenum", currstave, "stavenum/I")
    roottree2.Branch("chipnum", currchip, "chipnum/I")
    roottree2.Branch("col", col, "col/I")
    roottree2.Branch("row", row, "row/I")
    sum = [0] * 200
    counter = [0] * 200
    for xline in f: #loop on file lines (paths)
        run = re.search('run(.+?)/thr', xline)
        if(run):
            runnum[0] = int(run.group(1)) # run number
            if(runs[0] <= runnum[0] <= runs[1]):
                #read data
                datathr = readdata(str(xline.rstrip()))
                # calculate the average thr chip by chip
                for i in range(len(datathr)):## loop on rows
                    currstave[0] = int(staven) + math.floor(i / 512.);

                    for j in range(len(datathr[i])): ##loop on columns
                        currchip[0] = math.floor(j / 1024)
                        if(datathr[i][j]!=0):
                            sum[currchip[0]]+=datathr[i][j]
                            counter[currchip[0]]+=1
                        if(datathr[i][j]==0):#dead pixels
                            row[0] = i;
                            col[0] = j;
                            roottree2.Fill()
                    if(int(staven) + math.floor((i+1) / 512.) > currstave[0]):## when changing stave, write into the tree and reset the counters
                        for ichip in range(len(counter)):
                            if(counter[ichip]!=0):
                                avgthr[0] = sum[ichip]/counter[ichip]
                                currchip[0] = ichip
                                roottree1.Fill()
                        sum = [0] * 200
                        counter = [0] * 200
    roottree1.Write()
    roottree2.Write()
    ftree.Close()
    f.close()


main(first_arg, second_arg)
