""" Reformat amino transferases

    python enzpred/scripts/reformat_aminotransferases.py --outdir ~/Desktop/aminotransferase --data-file data/raw/aminotransferase/aminotransferase_data.xlsx
"""

import os
import pandas as pd
import json
import time
import argparse
from collections import defaultdict
import numpy as np

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize

from enzpred.utils import download_utils


# Unlogged threshold 
THRESH = 0.01
# Motivated by 1e-2 as the lower bound of detection
MIN_FACTOR = 1e-3

def main(outdir, data_file): 
    """ Main method"""

    df_conversion = pd.read_excel(data_file,
                                  header=1)

    prot_ids = [i.strip() for i in df_conversion["SeqID"].values]

    uniprot_file = os.path.join(outdir, "aminotransferase_sequence_file.txt")
    if os.path.exists(uniprot_file):
        os.remove(uniprot_file)
    download_utils.download_uniprot_uniparc_list(uniprot_file, prot_ids)

    df = pd.read_csv(uniprot_file, sep="\t")
    id_to_seq = dict(
        zip(df['Entry'].values.tolist(), df['Sequence'].values.tolist())
    )


    # Make matrix of binarized conversion 
    entries = []
    for num, row in df_conversion.iterrows(): 
        seq_id = row["SeqID"].strip()
        for substrate, val in row.items():
            if substrate not in {"Unnamed: 0", "SeqID"}:
                new_entry = {"SEQ" : id_to_seq[seq_id], 
                             "SUBSTRATES": substrate,  
                             "LogSpActivity" : np.log(val+MIN_FACTOR)}
                entries.append(new_entry)

    df = pd.DataFrame(entries)
    print(f"Unique substrates: {len(pd.unique(df['SUBSTRATES']))}")
    print(f"Unique sequences: {len(pd.unique(df['SEQ']))}")
    print(f"Length of dataset {len(df)}")

    # Take the maximum value for equivalent seq,smiles pairs
    # Export!
    # First is regression
    for dataset_type in ['', '_binary']:  
        out_name = f"aminotransferase{dataset_type}"
        out_file = os.path.join(outdir, f'{out_name}.csv')

        if dataset_type == "_binary": 
            df = df.copy()
            df["LogSpActivity"] = [0 if i < np.log(THRESH) else 1 
                                for i in df["LogSpActivity"]]

        df.to_csv(out_file)

        import matplotlib.pyplot as plt
        import seaborn as sns
        pivoted = df.pivot_table(index="SEQ",
                                     columns="SUBSTRATES",
                                     values="LogSpActivity")
        plt.figure()
        sns.heatmap(pivoted.values)
        plt.xlabel("Substrates")
        plt.ylabel("Sequences")
        plt.title(f"Aminotransferases{dataset_type}, Thresh = {THRESH}")
        plt.savefig(os.path.join(outdir, f"{out_name}.png"), 
                    bbox_inches="tight")

if __name__=="__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir",
                        type=str,
                        help="outdir",)
    parser.add_argument("--data-file",
                        type=str,
                        help="Data file as input (excel file extracted from si)")
    args = parser.parse_args()


    outdir = args.outdir
    data_file = args.data_file 

    os.makedirs(outdir, exist_ok=True)

    main(outdir, data_file)





