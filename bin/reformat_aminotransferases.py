""" Reformat aminotransferases extracted from """

import os
import pandas as pd
import argparse
import numpy as np

import download_utils

# Unlogged threshold
THRESH = 0.01

# Motivated by 1e-2 as the lower bound of detection
MIN_FACTOR = 1e-3


def main(outdir, rawdir, data_file):
    """Main method"""

    df_conversion = pd.read_excel(data_file, header=1, engine="openpyxl")

    bool_select = [isinstance(i, str) for i in df_conversion["SeqID"].values]
    df_conversion = df_conversion[bool_select]

    prot_ids = [i.strip() for i in df_conversion["SeqID"].values]

    uniprot_file = os.path.join(rawdir, "aminotransferase_sequence_file.txt")
    if os.path.exists(uniprot_file):
        os.remove(uniprot_file)
    download_utils.download_uniprot_uniparc_list(uniprot_file, prot_ids)

    df = pd.read_csv(uniprot_file, sep="\t")
    id_to_seq = dict(zip(df["Entry"].values.tolist(), df["Sequence"].values.tolist()))

    # Make matrix of binarized conversion
    entries = []
    for num, row in df_conversion.iterrows():
        seq_id = row["SeqID"].strip()
        for substrate, val in row.items():
            if substrate not in {"Unnamed: 0", "SeqID"}:
                new_entry = {
                    "SEQ": id_to_seq[seq_id],
                    "SUBSTRATES": substrate,
                    "LogSpActivity": np.log(val + MIN_FACTOR),
                }
                entries.append(new_entry)

    df = pd.DataFrame(entries)
    print(f"Unique substrates: {len(pd.unique(df['SUBSTRATES']))}")
    print(f"Unique sequences: {len(pd.unique(df['SEQ']))}")
    print(f"Length of dataset {len(df)}")

    # Take the maximum value for equivalent seq,smiles pairs
    # Export both regression and binary class.
    for dataset_type in ["", "_binary"]:
        out_name = f"aminotransferase{dataset_type}"
        out_file = os.path.join(outdir, f"{out_name}.csv")

        if dataset_type == "_binary":
            df = df.copy()
            df["LogSpActivity"] = [
                0 if i < np.log(THRESH) else 1 for i in df["LogSpActivity"]
            ]

        df.to_csv(out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="outdir", default="data/processed")
    parser.add_argument(
        "--rawdir", type=str, help="raw directory", default="data/raw/aminotransferase"
    )
    parser.add_argument(
        "--data-file",
        type=str,
        help="Data file as input (excel file extracted from si)",
        default="data/raw/aminotransferase/aminotransferase_data.xlsx",
    )
    args = parser.parse_args()

    outdir = args.outdir
    rawdir = args.rawdir
    data_file = args.data_file

    os.makedirs(outdir, exist_ok=True)

    main(outdir, rawdir, data_file)
