""" Reformat nitrilase dataa from Black et al. 2014"""

import os
import pandas as pd
import argparse

import download_utils


def main(outdir, rawdir, data_file):
    """Main method"""

    df_conversion = pd.read_excel(data_file, engine="openpyxl")

    not_na = [isinstance(i, str) for i in df_conversion["Substrate"]]
    df_conversion = df_conversion[not_na]

    protein_keys = set(df_conversion.keys())
    prot_ids = protein_keys.difference({"Number", "Substrate", "Smiles"})
    uniprot_file = os.path.join(rawdir, "nitrilase_sequence_file.txt")
    if os.path.exists(uniprot_file):
        os.remove(uniprot_file)
    download_utils.download_uniprot_uniparc_list(uniprot_file, prot_ids)
    df = pd.read_csv(uniprot_file, sep="\t")
    id_to_seq = dict(zip(df["Entry"].values.tolist(), df["Sequence"].values.tolist()))

    # Make matrix of binarized conversion
    entries = []
    for num, row in df_conversion.iterrows():
        smiles_string = row["Smiles"]
        for seq_id in id_to_seq.keys():
            ammonia_conversion = row[seq_id]
            conversion = 1 if ammonia_conversion > 0 else 0
            new_entry = {
                "SEQ": id_to_seq[seq_id],
                "SUBSTRATES": smiles_string,
                "Conversion": conversion,
            }
            entries.append(new_entry)

    bin_df = pd.DataFrame(entries)

    # Take the maximum value for equivalent seq,smiles pairs
    # bin_df = bin_df.groupby(["SEQ", "SUBSTRATES"]).max().reset_index()
    out_bin = os.path.join(outdir, f"nitrilase_binary.csv")
    bin_df.to_csv(out_bin)
    print(f"Length of dataset {len(bin_df)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="outdir", default="data/processed")
    parser.add_argument(
        "--rawdir", type=str, help="rawdir", default="data/raw/nitrilase"
    )

    parser.add_argument(
        "--data-file",
        type=str,
        default="data/raw/nitrilase/nitrilase_data.xlsx",
        help="Data file as input.",
    )
    args = parser.parse_args()

    outdir = args.outdir
    rawdir = args.rawdir
    data_file = args.data_file

    os.makedirs(outdir, exist_ok=True)

    main(outdir, rawdir, data_file)
