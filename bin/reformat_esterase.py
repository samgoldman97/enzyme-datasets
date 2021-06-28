""" Reformat esterase data from martinez-martinez et al. 

All values >0 are considered hits as in the original paper. 

Remove one enzyme over 1000 residues long (second highest is 694)

"""
import os
from itertools import groupby
import re
import pandas as pd

RAW_FOLDER = "data/raw/esterase"
OUT_FOLDER = "data/processed/"

path_wrap = lambda x: os.path.join(RAW_FOLDER, x)

DATA_FILE = path_wrap("cb7b00996_si_002.xls")
OUT_NAME = path_wrap("esterase_binary.csv")


def create_new_table():
    """Create new table"""

    ### Grab sequences
    seq_df = pd.read_excel(DATA_FILE, sheet_name="Supplementary Table S1")
    ids = seq_df["ID enzyme"].values

    ## Remove trailing character
    seqs = seq_df["AMINO ACID SEQUENCE"].values
    seqs = [i.replace("*", "").upper() for i in seqs]
    seqs = [i.replace(" ", "") for i in seqs]
    id_to_seq = dict(zip(ids, seqs))
    seq_to_id = dict(zip(seqs, ids))

    ### Grab smiles and values
    val_df = pd.read_excel(DATA_FILE, sheet_name="Supplementary Table S3")

    seq_cols = [j for j in val_df.columns if j in id_to_seq]

    out_dicts = []
    for row_num, entry in val_df.iterrows():
        smiles = entry["Smiles code"]
        for seq_id in seq_cols:
            seq = id_to_seq[seq_id]
            conversion = entry[seq_id]
            new_entry = {
                "SEQ": seq,
                "SUBSTRATES": smiles.strip(),
                "activity": conversion,
            }
            out_dicts.append(new_entry)

    df_out = pd.DataFrame(out_dicts)

    ### Confirm that zero values are non hits according to them
    # import numpy as np
    # res = df_out.groupby("SEQ").apply(lambda x: np.sum(x["activity"] > 0))
    # out = {seq_to_id[i]: j for i, j in res.items()}

    # Binarize data
    pos_activity = df_out["activity"] > 0
    df_out.loc[pos_activity, "activity"] = 1
    df_out.loc[~pos_activity, "activity"] = 0

    # Filter by sequence length
    df_out = df_out[[len(i) < 1000 for i in df_out["SEQ"]]].reset_index(drop=True)

    ### Grab data
    df_out.to_csv(OUT_NAME)


def main():
    os.makedirs(OUT_FOLDER, exist_ok=True)
    create_new_table()


if __name__ == "__main__":
    main()
