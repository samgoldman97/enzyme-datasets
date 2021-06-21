""" Reformat PafA high throughput screen from Markin et al., 2020 

    python enzpred/scripts/reformat_pafa.py --outdir ~/Desktop/pafa --data-file data/raw/pafa/pafa_data.csv
"""

import os
import re 
import pandas as pd
import json
import time
import argparse
from collections import defaultdict

#import rdkit
#from rdkit import Chem
#from rdkit.Chem import Draw
#from rdkit.Chem.MolStandardize import rdMolStandardize

import numpy as np


smiles_dict = {
    "cMUP" : "OC(CC(C1=C(O2)C=C(OP(O)(O)=O)C=C1)=CC2=O)=O", 
    "MecMUP" : "OC(CC(C1=C(O2)C=C(OP(OC)(O)=O)C=C1)=CC2=O)=O"
}


# ID Q9KJX5
uniprot_seq = ("MLTPKKWLLGVLVVSGMLGAQKTNAVPRPKLVVGLVVDQMRWDYLYRYYSKYGEGGFKRM"
               "LNTGYSLNNVHIDYVPTVTAIGHTSIFTGSVPSIHGIAGNDWYDKELGKSVYCTSDETVQ"
               "PVGTTSNSVGQHSPRNLWSTTVTDQLGLATNFTSKVVGVSLKDRASILPAGHNPTGAFWF"
               "DDTTGKFITSTYYTKELPKWVNDFNNKNVPAQLVANGWNTLLPINQYTESSEDNVEWEGL"
               "LGSKKTPTFPYTDLAKDYEAKKGLIRTTPFGNTLTLQMADAAIDGNQMGVDDITDFLTVN"
               "LASTDYVGHNFGPNSIEVEDTYLRLDRDLADFFNNLDKKVGKGNYLVFLSADHGAAHSVG"
               "FMQAHKMPTGFFVEDMKKEMNAKLKQKFGADNIIAAAMNYQVYFDRKVLADSKLELDDVR"
               "DYVMTELKKEPSVLYVLSTDEIWESSIPEPIKSRVINGYNWKRSGDIQIISKDGYLSAYS"
               "KKGTTHSVWNSYDSHIPLLFMGWGIKQGESNQPYHMTDIAPTVSSLLKIQFPSGAVGKPI"
               "TEVIGR")

# Regex for parsing variants
MUT_RE = "([A-Z])([0-9]*)([A-Z])"

# Names of all columns to extract
extract_values = {"cMUP": "kcatOverKM_cMUP_M-1s-1", 
                  "MecMUP" : "kcatOverKM_MecMUP_M-1s-1"}

def main(outdir, data_file): 
    """ Main method"""
    df_conversion = pd.read_csv(data_file)

    entries = []
    for row_index, row in df_conversion.iterrows(): 
        variant = row['variant']
        if variant == "WT": 
            new_seq = "".join(list(uniprot_seq))
        else:
            new_seq = list(uniprot_seq)
            regex_matches = re.findall(MUT_RE, variant)
            # There's 1 multimutant
            for index, (start, pos, change) in enumerate(regex_matches):
                pos = int(pos)
                if uniprot_seq[pos - 1] != start: 
                    raise ValueError()
                new_seq[pos - 1] = change

            new_seq = "".join(new_seq)
        for extract_name, extract_key in extract_values.items(): 
            new_entry = {"SEQ" : new_seq, 
                         "SUBSTRATES": smiles_dict[extract_name],
                         "LogKcatKm" : np.log(row[extract_key])}
            entries.append(new_entry)

    df = pd.DataFrame(entries)
    out_dir = os.path.join(outdir, f'pafa.csv')

    df.to_csv(out_dir)
    print(f"Length of dataset {len(df)}")
    print(f"Num unique prots: {len(pd.unique(df['SEQ']))}")

    import matplotlib.pyplot as plt
    import seaborn as sns
    for sub_df_key, sub_df in df.groupby("SUBSTRATES"): 
        param_dist = sub_df["LogKcatKm"]
        plt.figure()
        sns.distplot(param_dist.values)
        plt.xlabel("Sequences log(KcatKm)")
        plt.ylabel("Frequency")
        plt.title(f"PafA: {sub_df_key}")
        plt.savefig(os.path.join(outdir, f"{sub_df_key}.png"), bbox_inches="tight")

if __name__=="__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir",
                        type=str,
                        help="outdir",
                        default="/Users/samgoldman/Desktop/pafa")
    parser.add_argument("--data-file",
                        type=str,
                        help="Data file as input (sd02).xlsx",)
    args = parser.parse_args()


    outdir = args.outdir
    data_file = args.data_file 

    os.makedirs(outdir, exist_ok=True)

    main(outdir, data_file)







