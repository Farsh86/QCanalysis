import yaml

def main():

    maxstave = [11, 15, 19, 23, 29, 41, 47]
    #output file
    fileout = open(f"../Data/noisescan_P2_run518567.txt","w")
    fileout.write("Layer Stave Hs HIC Chip Row Col\n")
    for ilay in range(0,3): ## IB only for now
        for istave in range(0,20):
            if istave > maxstave[ilay]:
                continue;
            with open(f"../yaml/noise_masks/L{ilay}_{istave:02d}.yml", 'r') as f:
                yamldata=yaml.load(f, Loader=yaml.FullLoader) or {}

            ##dump in the txt file
            for (ichip, ivalue) in enumerate(yamldata.items()):
                for ipix in ivalue[1]:
                    if ilay<3: ##IB
                        fileout.write(f"{ilay} {istave} L 0 {ichip} {ipix[1]} {ipix[0]}\n")


if __name__ == '__main__':
    main()
