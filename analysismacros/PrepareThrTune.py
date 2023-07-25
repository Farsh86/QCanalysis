import yaml
import json
import argparse
import ROOT
import array as a
import numpy as np
'''
Prepare threshold tuning yamls from json file in ib-commissioning-tools (on-surface commissioning data)
'''

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, help="Input file to be analysed")
    args = parser.parse_args()
    print(f"Analysing file: {args.file}")
    bb = "0v" if ("0v" in args.file) else "3v"

    #open json file
    f = open(args.file, 'r')
    jsondata = json.load(f)
    dict = {}
    dict2 = {}
    #loop on staves, dump ITHR values
    for stave in jsondata:
        for ic, chip in enumerate(jsondata[stave]):
            dict.update({ic:jsondata[stave][chip]["ITHR"]})
            dict2.update({ic:jsondata[stave][chip]["VCASN"]})
        with open(f"../yaml/ithr/{bb}/{stave}.yml", 'w') as f:
            yaml.dump(dict, f)
        with open(f"../yaml/vcasn/{bb}/{stave}.yml", 'w') as f:
            yaml.dump(dict2, f)

if __name__ == '__main__':
    main()
