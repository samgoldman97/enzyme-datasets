""" Extract DUF dataset from Bastard et al. (2014) """
import os
import pandas as pd
import numpy as np
import argparse

import parse_utils


def main(outdir, data_file, manually_annotated, enzyme_file):
    """Main method"""

    # Get all substrate mappings
    new_df = pd.read_excel(manually_annotated, engine="openpyxl")

    # Remove NaN
    bool_indices = [isinstance(i, str) for i in new_df["Annotated item"]]

    mapping = dict(
        zip(
            new_df["Annotated item"][bool_indices],
            new_df["substrate_smiles"][bool_indices],
        )
    )

    # Only select where we map the substrate out
    entries = [
        i.strip().split("$")
        for i in open(data_file, "r").readlines()
        if i.strip().split("$")[2] in mapping
    ]
    enzyme, replicate, substrate, value = [np.array(i) for i in zip(*entries)]
    print("Unique substrates: ", set(substrate))

    # Handle sequence mapping
    seq_mapping = {
        name.replace("BKACE", "KCE"): seq.replace("-", "")
        for name, seq in parse_utils.fasta_iter(enzyme_file)
    }
    not_found = set()
    for e in set(enzyme):
        if e not in seq_mapping:
            not_found.add(e)

    print(f"Could not find {len(not_found)} enzyme ids in the msa: {len(not_found)}")
    print(f"Found {len(set(enzyme)) - len(not_found)} enzymes")

    out_dicts_bin = []
    for enz_id, replicate_num, substrate_name, val in zip(
        enzyme, replicate, substrate, value
    ):

        if enz_id not in seq_mapping:
            continue

        aa_seq = seq_mapping[enz_id]
        smiles = mapping[substrate_name]

        out_dict_bin = {
            "SEQ": aa_seq,
            "SUBSTRATES": smiles,
            "Activity": int(val),
            # "Replicate": replicate_num
        }

        out_dicts_bin.append(out_dict_bin)

    bin_df = pd.DataFrame(out_dicts_bin)

    # Take the maximum value for equivalent seq,smiles pairs (replicates)
    bin_df = bin_df.groupby(["SEQ", "SUBSTRATES"]).max().reset_index()
    out_bin = os.path.join(outdir, f"duf_binary.csv")
    bin_df.to_csv(out_bin)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="outdir", default="data/processed/")
    parser.add_argument(
        "--data-file",
        type=str,
        help="Data file as input. This is the binary outfile from Karine's scripts",
        default="data/raw/duf/BinaryActivities.DB",
    )
    parser.add_argument(
        "--manually-annotated",
        type=str,
        help="Manually annotated excel sheet mapping names to smiles",
        default="data/raw/duf/duf_smiles.xlsx",
    )
    parser.add_argument(
        "--enzyme-file",
        type=str,
        help="Name of the file containing enzyme annotations",
        default="data/raw/duf/duf_msa.txt",
    )
    args = parser.parse_args()

    outdir = args.outdir
    data_file = args.data_file
    manually_annotated = args.manually_annotated

    os.makedirs(outdir, exist_ok=True)
    main(outdir, data_file, manually_annotated, args.enzyme_file)
