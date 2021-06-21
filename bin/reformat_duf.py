#!/usr/bin/env python
# coding: utf-8
""" Extract DUF dataset from Bastard et al. (2014)

Sample usage: 
    python enzpred/scripts/reformat_duf.py --outdir ~/Desktop/duf --data-file data/raw/duf/BinaryActivities.DB --enzyme-file data/raw/duf/duf_msa.txt --manually-annotated data/raw/duf/duf_smiles.xlsx

Extraction steps:
- Use R code supplied by Karine Bastard to export a binary classification of
  her data
- Manually map from compond (whether it was tested in forward or reverse
  direction) into the _substrate_ smiles string for simplification. This is for
  convenience 
- Take the _max_ of all substrate, enzyme experimental pairs to get some ground
  truth value


TODO: 
- Change the input data file if Karine replies. 

"""


import os
import pandas as pd
import numpy as np
import json
import time
import argparse
from collections import defaultdict

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolops import RemoveStereochemistry

from enzpred.utils import parse_utils




def main(outdir, data_file, manually_annotated, 
         enzyme_file):
    """ Main method"""

    # Get all substrate mappings
    new_df = pd.read_excel(manually_annotated)
    mapping = dict(zip(new_df['Annotated item'],
                       new_df['substrate_smiles']))


    # Only select where we map the substrate out
    entries = [i.strip().split("$") 
               for i in open(data_file, "r").readlines() 
               if i.strip().split("$")[2] in mapping]
    enzyme, replicate, substrate, value = [np.array(i) for i in zip(*entries)]
    print("Unique substrates: ", set(substrate))

    # Handle sequence mapping
    seq_mapping = {name.replace("BKACE", "KCE") : seq.replace("-", "") 
                   for name,seq in  parse_utils.fasta_iter(enzyme_file)}
    not_found = set()
    for e in set(enzyme):
        if e not in seq_mapping:
            not_found.add(e)

    print(f"Could not find {len(not_found)} enzyme ids in the msa: {len(not_found)}")
    print(f"Found {len(set(enzyme)) - len(not_found)} enzymes")

    out_dicts_bin = []
    for enz_id, replicate_num, substrate_name, val in zip(enzyme, replicate,
                                                          substrate, value):

        if enz_id not in seq_mapping: continue

        aa_seq = seq_mapping[enz_id]
        smiles = mapping[substrate_name]

        out_dict_bin = {
            "SEQ" : aa_seq,
            "SUBSTRATES": smiles,
            "Activity": int(val), 
            #"Replicate": replicate_num 
        }

        out_dicts_bin.append(out_dict_bin)

    bin_df = pd.DataFrame(out_dicts_bin)
    # Take the maximum value for equivalent seq,smiles pairs
    bin_df = bin_df.groupby(["SEQ", "SUBSTRATES"]).max().reset_index()
    out_bin = os.path.join(outdir, f'duf_binary.csv')
    bin_df.to_csv(out_bin)

    import matplotlib.pyplot as plt
    import seaborn as sns

    pivoted = bin_df.pivot_table(index="SEQ",
                                 columns="SUBSTRATES",
                                 values="Activity")
    plt.figure()
    sns.clustermap(pivoted.values)
    plt.xlabel("Substrates")
    plt.ylabel("Sequences")
    plt.title(f"DUF Enzymes")
    plt.savefig(os.path.join(outdir, f"duf_heatmap.png"), 
                bbox_inches="tight")

if __name__=="__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir",
                        type=str,
                        help="outdir",
                        default="/Users/samgoldman/Desktop/duf")
    parser.add_argument("--data-file",
                        type=str,
                        help="Data file as input. This is the binary outfile from Karine's scripts",
                        default='data/raw/duf/BinaryActivities.DB')
    parser.add_argument("--manually-annotated",
                        type=str,
                        help="Manually annotated excel sheet mapping names to smiles",
                        default='data/raw/duf/duf_smiles.xlsx')
    parser.add_argument("--enzyme-file",
                        type=str,
                        help="Name of the file containing enzyme annotations",
                        default='data/raw/duf/duf_msa.txt')
    args = parser.parse_args()

    outdir = args.outdir
    data_file = args.data_file 
    manually_annotated = args.manually_annotated

    os.makedirs(outdir, exist_ok=True)
    main(outdir, data_file, manually_annotated,
         args.enzyme_file)



