"""Reformat the rubisco dataset from Davidi et al. 2020
    python enzpred/scripts/reformat_rubisco.py --outdir ~/Desktop/rubiscos --data-file data/raw/rubiscos/rubisco_data.xlsx
"""

import os
import pandas as pd
import json
import time
import argparse
from collections import defaultdict

def main(outdir, data_file): 
    """ Main method"""

    df  = pd.read_excel(data_file, "spectroscopic assay")
    proteins = df["Protein seq"].values
    rate_means = df["Rate mean [s-1] "].values
    selection = ~pd.isna(rate_means)

    proteins = proteins[selection]
    rate_means = rate_means[selection]

    # Make matrix of binarized conversion 
    entries = []
    for protein, rate in zip(proteins, rate_means): 
        smiles_string = "C(=O)=O"
        new_entry = {"SEQ" : protein.replace("*", ""), 
                     "SUBSTRATES": smiles_string, 
                     "Rate" : rate}
        entries.append(new_entry)

    reg_df = pd.DataFrame(entries)

    # Take the maximum value for equivalent seq,smiles pairs
    # bin_df = bin_df.groupby(["SEQ", "SUBSTRATES"]).max().reset_index()
    out_reg = os.path.join(outdir, f'rubisco.csv')
    reg_df.to_csv(out_reg)
    print(f"Length of dataset {len(reg_df)}")

    import matplotlib.pyplot as plt
    import seaborn as sns
    rate_vals= reg_df["Rate"].values    
    sns.distplot(rate_vals, bins=20)
    plt.xlabel("Rate ($s^{-1}$)")
    plt.ylabel("Frequency")
    plt.title("Rubiscos")
    plt.savefig(os.path.join(outdir, "Rubiscos.png"), bbox_inches="tight")

if __name__=="__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir",
                        type=str,
                        help="outdir",
                        default="/Users/samgoldman/Desktop/rubisco")
    parser.add_argument("--data-file",
                        type=str,
                        help="Data file as input")
    args = parser.parse_args()


    outdir = args.outdir
    data_file = args.data_file 

    os.makedirs(outdir, exist_ok=True)

    main(outdir, data_file)







