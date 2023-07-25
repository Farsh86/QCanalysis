import yaml
import argparse
import ROOT
import array as a
import numpy as np

def main():

    #HS-L B0
    coltochipid0l = {0:16, 1:17, 2:18, 3:19, 4:20, 5:21, 6:22,
                    7:32, 8:33, 9:34, 10:35, 11:36, 12:37, 13:38,
                    14:48, 15:49, 16:50, 17:51, 18:52, 19:53, 20:54,
                    21:64, 22:65, 23:66, 24:67, 25:68, 26:69, 27:70,
                    28:80, 29:81, 30:82, 31:83, 32:84, 33:85, 34:86,
                    35:96, 36:97, 37:98, 38:99, 39:100, 40:101, 41:102,
                    42:112, 43:113, 44:114, 45:115, 46:116, 47:117, 48:118}
    #HS-L A8
    coltochipid8l = {0:30, 1:29, 2:28, 3:27, 4:26, 5:25, 6:24,
                    7:46, 8:45, 9:44, 10:43, 11:42, 12:41, 13:40,
                    14:62, 15:61, 16:60, 17:59, 18:58, 19:57, 20:56,
                    21:78, 22:77, 23:76, 24:75, 25:74, 26:73, 27:72,
                    28:94, 29:93, 30:92, 31:91, 32:90, 33:89, 34:88,
                    35:110, 36:109, 37:108, 38:107, 39:106, 40:105, 41:104,
                    42:126, 43:125, 44:124, 45:123, 46:122, 47:121, 48:120}
    #HS-U B0
    coltochipid0u = {0:144, 1:145, 2:146, 3:147, 4:148, 5:149, 6:150,
                    7:160, 8:161, 9:162, 10:163, 11:164, 12:165, 13:166,
                    14:176, 15:177, 16:178, 17:179, 18:180, 19:181, 20:182,
                    21:192, 22:193, 23:194, 24:195, 25:196, 26:197, 27:198,
                    28:208, 29:209, 30:210, 31:211, 32:212, 33:213, 34:214,
                    35:224, 36:225, 37:226, 38:227, 39:228, 40:229, 41:230,
                    42:240, 43:241, 44:242, 45:243, 46:244, 47:245, 48:246}

    #HS-U A8
    coltochipid8u = {0:158, 1:157, 2:156, 3:155, 4:154, 5:153, 6:152,
                    7:174, 8:173, 9:172, 10:171, 11:170, 12:169, 13:168,
                    14:190, 15:189, 16:188, 17:187, 18:186, 19:185, 20:184,
                    21:206, 22:205, 23:204, 24:203, 25:202, 26:201, 27:200,
                    28:222, 29:221, 30:220, 31:219, 32:218, 33:217, 34:216,
                    35:238, 36:237, 37:236, 38:235, 39:234, 40:233, 41:232,
                    42:254, 43:253, 44:252, 45:251, 46:250, 47:249, 48:248}

    rowselector = {0: coltochipid0l, 1: coltochipid8l, 2: coltochipid0u, 3: coltochipid8u}
    '''
    noisyrowstave = {"L4_16":[0, 2, 177], "L6_04":[1,313,201], "L6_06":[2,250,252], "L6_20":[3,111,211],
    "L4_03":[4,360,160], "L5_26":[5,319,25],"L5_09":[6,473,196],"L5_13":[7,198,148],"L5_32":[8,251,18], "L5_38":[9,227,19],
    "L5_04":[10,189,169], "L5_07":[11,169,166]}# "stave": {id, row, chipid}
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, help="Input file to be analysed")
    parser.add_argument("-m", "--merge", default=False, choices=('True','False'), help="Merge with existing noise mask")
    args = parser.parse_args()
    print(f"Analysing file: {args.file}")
    print(f"Merging set to: {args.merge}")
    runnum = args.file[50:56]
    #open file
    infl = ROOT.TFile.Open(args.file, "READ")
    ntriggers = GetTriggers(infl)
    print(f"Number of triggers: {ntriggers}")

    file1 = open(f"../yaml/noise_masks/number_of_noisy_pix_{runnum}.txt","a")
    #Loop over all THnSparse and prepare yaml file with noisy pixels
    NOISECUT_IB = 1e-2
    NOISECUT_OB = 2e-6

    for key in infl.GetListOfKeys():
        obj=key.ReadObj()
        if obj.InheritsFrom("THnSparse"):
            if obj.GetEntries()-1 < 0: #skip empty THnSparse
                continue
            ## Loop of pixels which fired
            layer = obj.GetName()[9:10]
            dict = {}
            npix = 0
            for ibin in range(obj.GetNbins()):
                coord = np.array([0,0], dtype=np.int32)
                pixelhits = obj.GetBinContent(ibin, np.asarray(coord))
                fhr = pixelhits/ntriggers #fhr of the pixel
                NOISECUT = NOISECUT_IB if int(layer) < 3 else NOISECUT_OB
                if fhr < NOISECUT: ##noise cut
                    continue
                npix = npix+1
                chipid = 0
                if int(layer)<3: #ib
                    chipid = int((coord[0]-1)/1024)
                    if chipid not in dict:
                        dict.update({chipid:[[int(coord[0]-1-chipid*1024),int(coord[1]-1),fhr]]})
                    else:
                        dict[chipid].append([int(coord[0]-1-chipid*1024),int(coord[1]-1),fhr])
                else: #ob
                    rowidx = int(int(coord[1]-1) / 512)
                    colidx = int(int(coord[0]-1) / 1024)
                    chipid = rowselector[rowidx][colidx]
                    colfinal = 0
                    rowfinal = 0
                    if rowidx==0 or rowidx==2:
                        rowfinal = 511*(rowidx+1) + rowidx - int(coord[1]-1)
                        colfinal = int(coord[0]-1) - 1024*colidx
                    if rowidx==1 or rowidx==3:
                        rowfinal = int(coord[1]-1) - 512*(rowidx)
                        colfinal = 1023*(colidx+1) + colidx - int(coord[0]-1)
                    if chipid not in dict:
                        dict.update({chipid:[[colfinal,rowfinal,fhr]]})
                    else:
                        dict[chipid].append([colfinal,rowfinal,fhr])

            stavenum = obj.GetName()[14:15] if obj.GetName()[15:16] == "_" else obj.GetName()[14:16]

            #Add the entire noisy rows for some chip in some staves
            staveid = f"L{layer}_{int(stavenum):02d}"
            '''
            if staveid in noisyrowstave: ## if the stave is the one with noisy rows -> MASK FULL ROW
                for icol in range(0,1024):
                    check = 0
                    if noisyrowstave[staveid][2] in dict:
                        for el in dict[noisyrowstave[staveid][2]]:
                            if el[0]==icol and el[1] == noisyrowstave[staveid][1]:
                                check = 1
                                break
                        if check==0:
                            dict[noisyrowstave[staveid][2]].append([icol, noisyrowstave[staveid][1], 1.0]) ##dummy FHR = 1.0
                    else:
                        dict.update({noisyrowstave[staveid][2]:[[icol,noisyrowstave[staveid][1],1.0]]})
            '''

            #Add stuck pix in IB and OB (NOTE: VERY IMPORTANT IN IB BECAUSE WE GET MORE PIXELS FUNCTIONING)
            #OB
            if layer=="4" and stavenum=="3":
                if 198 not in dict:
                    dict.update({198:[[988,2,1.0]]})
                else:
                    dict[198].append([988,2,1.0]) ## this is a stuck pixel found during data taking (fhr is random)
            if layer=="6" and stavenum=="22":
                if 85 not in dict:
                    dict.update({85:[[226,248,1.0]]})
                else:
                    dict[85].append([226,248,1.0]) ## this is a stuck pixel found during data taking (fhr is random)
            if layer=="3" and stavenum=="7":
                if 65 not in dict:
                    dict.update({65:[[631,69,1.0]]})
                else:
                    dict[65].append([631,69,1.0]) ## this is a stuck pixel found during data taking (fhr is random)
            #IB
            if layer=="1" and stavenum=="15":
                if 4 not in dict:
                    dict.update({4:[[732,14,1.0]]})
                else:
                    dict[4].append([732,14,1.0]) ## this is a stuck pixel found during surf comm (fhr is random)
            if layer=="2" and stavenum=="9":
                if 3 not in dict:
                    dict.update({3:[[637,3,1.0]]})
                else:
                    dict[3].append([637,3,1.0]) ## this is a stuck pixel found during surf comm (fhr is random)
            if layer=="2" and stavenum=="18":
                if 6 not in dict:
                    dict.update({6:[[526,290,1.0]]})
                else:
                    dict[6].append([526,290,1.0]) ## this is a stuck pixel found during surf comm (fhr is random)
            if layer=="1" and stavenum=="3":
                if 4 not in dict:
                    dict.update({4:[[422,433,1.0],[423,433,1.0]]})
                else:
                    dict[4].append([422,433,1.0]) ## this is a stuck pixel found during surf comm (fhr is random)
                    dict[4].append([423,433,1.0])
            if layer=="2" and stavenum=="14":
                if 0 not in dict:
                    dict.update({0:[[362,173,1.0]]})
                else:
                    dict[0].append([362,173,1.0]) ## this is a stuck pixel found during surf comm (fhr is random)

            print(f"L{layer}_{int(stavenum):02d}: {npix} hot pixels above cut")

            #Merge with existing mask (if specified in the option)
            merge = dict.copy()
            if args.merge == "True":
                with open(f"../yaml/noise_masks/L{layer}_{int(stavenum):02d}.yml", 'r') as f:
                    dataold=yaml.load(f, Loader=yaml.FullLoader) or {}
                merge.update(dataold)
                for i in dict:
                    if dict[i] != merge[i]:
                        for j in dict[i]:
                            check = 0
                            for k in merge[i]:
                                if j[0]==k[0] and j[1]==k[1]:
                                    check = 1
                                    break
                            if j not in merge[i] and check==0:
                                merge[i].append(j)

            ##save yaml
            file1.write(f"L{layer}_{int(stavenum):02d} {npix}\n")
            with open(f"../yaml/noise_masks/L{layer}_{int(stavenum):02d}.yml", 'w') as f:
                yaml.dump(merge, f)
    file1.close()

## Function to get number of triggers
def GetTriggers(infl):
    fhr_chip_ib = 0.
    nhits_chip_ib = 0.
    for key in infl.GetListOfKeys():
        obj=key.ReadObj()
        name = obj.GetName()
        if obj.InheritsFrom("TH2") and ("L0" in name): ##just take L0
            fhr_chip_ib = obj.GetBinContent(1,2)
            if fhr_chip_ib<1e-15:
                print("Please change chip, this bin is empty")
        if obj.InheritsFrom("THnSparse") and ("L0_Stv1" in name): ##just take L0_00
            h2 = obj.Projection(1,0)
            nhits_chip_ib = h2.Integral(1,1024,1,512)
            del h2
            break
    return nhits_chip_ib / (512.*1024.*fhr_chip_ib)

if __name__ == '__main__':
    main()
