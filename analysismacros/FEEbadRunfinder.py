# FEEbadRunfinder.py macro
# Author: G. Mandaglio - 8 August 2022

# The macro looks for the bad runs by using the criteria described above
# The macro analyze the txt files (error, warning and fault) produced in the general analysis program (lane_dump_....txt) 
#[Analyses on FEE]
#	 13. Dump in txt file of lanes into ERROR, FAULT, WARNING
# Run the macro with the following command line (put the macro in the same folder of the txt files)
# python FEEbadRunfinder.py FirstRunNumber  LastRunNumber

#Inner Barrel: L0, L1, L2    #9 lanes x staves
#L0 -> 12 staves -> 108 lanes
#L1 -> 16 staves -> 144 lanes
#L2 -> 20 staves -> 180 lanes

#Outer Barrel - Middle layers: L3, L4 #16 lanes x staves
#L3 -> 24 staves -> 384 lanes
#L4 -> 30 staves -> 480 lanes

#Per l'Outer Barrel - Outer Layers: L5, L6  #28 lanes per stave
#L5 -> 42 staves -> 1176 lanes
#L6 -> 48 staves -> 1344 lanes

#total 3816 lanes

#  Bad finder criteria
#  there is at least 1 layer with >25% staves with >80% of the lanes (per stave) in any status (warning/fault/error)
#  there are >10% of the lanes in any status overall (full detector).



import re       #regular expression in string finder package

def StripLaneFile(file_name):
    File = open(file_name, 'r')
    Lines = File.readlines()

    lista = []

    for line in Lines:
        s = re.findall(r'\d+', line)    #list of number in the string
        s.pop(1)                        #FEEID number removed
        s_int = [int(x) for x in s]     #string to integer
        lista.append(s_int)             #list of array = matrix
    if not lista:
        lista = [[0, 0, 0, 0],[0, 0, 0, 0]]  #it removes the list mismatch issue in case of an empty file

    return lista



import numpy as np
import sys
#Bad tresholds
too_many_lanes = 3816 * 0.1

too_many_staves = np.array([12,16,20,24,30,42,48]) # too_many_staves[NLayer]
too_many_staves = np.multiply(too_many_staves,0.25) #25% staves per layer
lanes_limit_x_staves = np.array([9,9,9,16,16,28,28])#lanes_limit_x_staves[NLayer]
lanes_limit_x_staves = np.multiply(lanes_limit_x_staves,0.8) #80% staves per layer

errorfile ='lane_dump_error_run'+sys.argv[1]+'_to_run'+sys.argv[2]+'.txt'
faultfile ='lane_dump_fault_run'+sys.argv[1]+'_to_run'+sys.argv[2]+'.txt'
warningfile ='lane_dump_warning_run'+sys.argv[1]+'_to_run'+sys.argv[2]+'.txt'
error = StripLaneFile(errorfile)
fault = StripLaneFile(faultfile)
warning = StripLaneFile(warningfile)
#we merge all problem sources
data = np.concatenate((error, warning), axis =0)
data = np.concatenate((data, fault), axis =0)
#we remove the multiple problem for the same lane
data = np.unique(data, axis=0)

spy_lane = 1
spy_stave = 1
lane_repetition = 1
stave_repetition = 0
for i in range(len(data[:,0])-1):
    if data[i,0] == data[i+1,0]: #same run
        if data[i,1] == data[i+1,1]: #same layer
            if data[i,2] == data[i+1,2]:#same stave
                lane_repetition = lane_repetition +1
            else:
                lane_repetition =1
                spy_lane = 1
            if lane_repetition > lanes_limit_x_staves[data[i,1]] and spy_lane ==1:
                stave_repetition = stave_repetition +1
                spy_lane = 0
            if  stave_repetition > too_many_staves[data[i,1]] and spy_stave == 1:
                print("run ",data[i,0]," bad, too many staves in the layer ", data[i,1])
                spy_stave = 0
        else:
            stave_repetition =0
            lane_repetition = 1
            spy_lane = 1
    else:
        stave_repetition =0
        lane_repetition = 1
        spy_lane = 1
        spy_stave = 1


############check numbers of lans in efw in the whole detector##########
#we catch the run list
runs = data[:,0]
runs_lane_efw = data[:,0]
runs = np.unique(runs,axis=0)
for a in runs:
#    print ('run ',a,' numero di lane in efw',  np.count_nonzero(runs_lane_efw == a))
    if np.count_nonzero(runs_lane_efw == a) > too_many_lanes :
        print ('run ',a,' number of lanes in e_or_f_or_w',  np.count_nonzero(runs_lane_efw == a))
######################################################
