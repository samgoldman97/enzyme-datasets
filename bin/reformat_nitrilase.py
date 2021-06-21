""" Reformat nitrilase dataa from black et al.

    python enzpred/scripts/reformat_nitrilase.py --outdir ~/Desktop/nitrilase --data-file data/raw/nitrilase/nitrilase_data.xlsx
"""

import os
import pandas as pd
import json
import time
import argparse
from collections import defaultdict

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize

from enzpred.utils import download_utils

def main(outdir, data_file): 
    """ Main method"""

    df_conversion = pd.read_excel(data_file)
    protein_keys = set(df_conversion.keys())
    prot_ids = protein_keys.difference({"Number", "Substrate", "Smiles"})
    uniprot_file = os.path.join(outdir, "nitrilase_sequence_file.txt")
    if os.path.exists(uniprot_file):
        os.remove(uniprot_file)
    download_utils.download_uniprot_uniparc_list(uniprot_file, prot_ids)
    df = pd.read_csv(uniprot_file, sep="\t")
    id_to_seq = dict(
        zip(df['Entry'].values.tolist(), df['Sequence'].values.tolist()))

    # Make matrix of binarized conversion 
    entries = []
    for num, row in df_conversion.iterrows(): 
        smiles_string = row["Smiles"]
        for seq_id in id_to_seq.keys(): 
            ammonia_conversion = row[seq_id]
            conversion = 1 if ammonia_conversion > 0 else 0
            new_entry = {"SEQ" : id_to_seq[seq_id], 
                         "SUBSTRATES": smiles_string, 
                         "Conversion" : conversion}
            entries.append(new_entry)

    bin_df = pd.DataFrame(entries)

    # Take the maximum value for equivalent seq,smiles pairs
    # bin_df = bin_df.groupby(["SEQ", "SUBSTRATES"]).max().reset_index()
    out_bin = os.path.join(outdir, f'nitrilase_binary.csv')
    bin_df.to_csv(out_bin)
    print(f"Length of dataset {len(bin_df)}")

    import matplotlib.pyplot as plt
    import seaborn as sns
    pivoted = bin_df.pivot_table(index="SEQ",
                                 columns="SUBSTRATES",
                                 values="Conversion")
    sns.heatmap(pivoted.values)
    plt.xlabel("Substrates")
    plt.ylabel("Sequences")
    plt.title("Nitrilases")
    plt.savefig(os.path.join(outdir, "Nitrilases.png"), bbox_inches="tight")

if __name__=="__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir",
                        type=str,
                        help="outdir",
                        default="/Users/samgoldman/Desktop/nitrilase")
    parser.add_argument("--data-file",
                        type=str,
                        help="Data file as input.")
    args = parser.parse_args()



    outdir = args.outdir
    data_file = args.data_file 

    os.makedirs(outdir, exist_ok=True)

    main(outdir, data_file)





