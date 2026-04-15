#!/bin/env python3

import pandas as pd
import argparse

def args_parse_function(args_in):

    all_args = args_in

    all_args.add_argument("-a", "--path_anarci", default = [], required = True, help = "path to the input ANANRCI heavy chain csvfile")
    all_args.add_argument("-o", "--out_file", default = [], required = True, help = "name of the outfile")
    all_args.add_argument("-c", "--chain", choices = ["heavy","light"], required = True, help = "the expected chain of the sequence")
    args_return = all_args.parse_args()
    return args_return

def main():

    # get the arguments using argparse
    args_return = args_parse_function(argparse.ArgumentParser())

    # load the RS file
    anarci_pd = pd.read_csv(args_return.path_anarci, sep = ",", header = 0, index_col = 0)

    # subset the file
    anarci_filtered = anarci_pd.loc[:,["23","104","118"]]
    

    # find the all three ccw
    ccw = []
    for i,r in anarci_filtered.iterrows():
        
        if args_return.chain == "heavy":
            ccw.append(True if (r[0]=="C" and r[1]=="C" and r[2]=="W") else False)
        else:
            ccw.append(True if (r[0]=="C" and r[1]=="C" and r[2]=="F") else False)
    
    anarci_filtered["C_C_WF"] = ccw

    anarci_filtered.rename(columns = {"23":"IMGT_C1","104":"IMGT_C2","118":"IMGT_WF"}, inplace = True)

    print("The number of input sequences with the expected CCW/F amino acids is: " + str(sum(anarci_filtered["C_C_WF"]==True)))
    print("The number of input sequences without the expected CCW/F amino acids is: " +  str(sum(anarci_filtered["C_C_WF"]==False)))

    anarci_filtered.to_csv(args_return.out_file)

if __name__ == "__main__":
    main()

