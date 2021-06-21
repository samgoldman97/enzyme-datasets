#!/usr/bin/env python
# coding: utf-8
""" Extract kinase dataset from Bradley et al. 

Sample usage: 

Extraction steps:
- Take raw MSA file and convert it into actual sequences
"""


import os
import pandas as pd
import numpy as np
import json
import time
import argparse
from collections import defaultdict

from enzpred.utils import parse_utils

LEN_CUTOFF = 700

def main(outdir, data_folder, enzyme_file):
    """ Main method"""


    # Shift columns
    COL_SHIFT = -8
    # Only keep +5 and -5 regions
    UB = 5
    LB = -5


    name_to_seq = dict()
    for name, seq in parse_utils.fasta_iter(enzyme_file): 
        #new_name =name.split("|")[2].split(" ")[0]
        org = name.split("|")[2].split(" ")[0].split("_")[1]
        # Kinases are named by group
        gn = name.split("GN=")[1].split()[0]
        new_name = f"{gn}_{org}".upper()
        new_seq = seq.replace("-", "")


        if len(new_seq) < LEN_CUTOFF:
            name_to_seq[new_name] = new_seq

    kinase_dicts = {}
    for root, dirs, files in os.walk(data_folder, topdown=False):
        for file_ in files:
            out_dict = defaultdict(lambda : dict())
            if file_.endswith('.txt'): 
                with open(os.path.join(root, file_), "r") as fp: 
                    lines = [i.strip() for i in fp.readlines()]
                    cols = lines[0].split(" ")
                    for line in lines[1:]:
                        entries = line.split()
                        aa = entries[0]
                        # Verify zip has equal length
                        assert(len(entries) - 1 == len(cols))
                        for col, entry in zip(cols, entries[1:]): 
                            col = int(col)
                            col += COL_SHIFT
                            out_dict[col][aa] = float(entry )
                kinase_dicts[file_] = out_dict

    # Now fill all entries and PWM's
    aas = "ACDEFGHIKLMNPQRSTVWY"
    df = []
    for kinase_file, kinase_pwm  in kinase_dicts.items(): 
        kinase_key = kinase_file.split(".txt")[0]
        kinase_key = "_".join(np.array(kinase_key.split("_"))[[0,2]]).upper()
        if kinase_key in name_to_seq:    
            seq = name_to_seq[kinase_key]
            for pos in range(1,16): 
                pos = pos + COL_SHIFT
                if pos <= UB and pos >= LB: 
                    for aa in aas:
                        value = kinase_pwm[pos][aa] 
                        substrate = f"{aa}{pos}"
                        df.append({"SEQ" : seq,
                                   "SUBSTRATES" : substrate,
                                   "Probability" : value})


    df = pd.DataFrame(df)

    print(f"Number of kinases added: {len(pd.unique(df['SEQ']))}")

    ## Make heatmap 
    out = os.path.join(outdir, f'kinases.csv')
    df.to_csv(out)

    import matplotlib.pyplot as plt
    import seaborn as sns

    pivoted = df.pivot_table(index="SEQ",
                             columns="SUBSTRATES",
                             values="Probability")
    plt.figure()
    sns.clustermap(pivoted.values)
    plt.xlabel("(Position,amino)")
    plt.ylabel("Sequences")
    plt.title(f"Kinases sequence_logo_probability")
    plt.savefig(os.path.join(outdir, f"kinases_heatmap.png"), 
                bbox_inches="tight")

if __name__=="__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir",
                        type=str,
                        help="outdir",
                        default="/Users/samgoldman/Desktop/kinases")
    parser.add_argument("--data-folder",
                        type=str,
                        help="Folder containing the PWM files for each kinase",
                        default='data/raw/kinases/mmc5')
    parser.add_argument("--enzyme-file",
                        type=str,
                        help="Name of the file containing enzyme alignment",
                        default='data/raw/kinases/kinases_aligned.txt')
    args = parser.parse_args()

    outdir = args.outdir
    data_folder = args.data_folder
    enzyme_file = args.enzyme_file 

    os.makedirs(outdir, exist_ok=True)
    main(outdir, data_folder, 
         enzyme_file)



