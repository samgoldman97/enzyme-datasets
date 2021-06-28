"""Reformat the glycosyltransferase dataset from Yang et al. 2017
"""

import os
import argparse

import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdmolops import RemoveStereochemistry

import parse_utils


def main(outdir, data_file_acceptors, data_file_donors, seq_file):
    """Main method"""

    # Convert to upper
    seq_dict = {
        i: j.replace("*", "").strip().upper()
        for i, j in parse_utils.fasta_iter(seq_file)
    }

    ### Start with data_file acceptrs
    for chiral in [True, False]:
        for input_file, base_name in zip(
            [data_file_acceptors, data_file_donors], ["gt_acceptors", "gt_donors"]
        ):

            base_name = f"{base_name}_chiral" if chiral else f"{base_name}_achiral"

            df = pd.read_excel(
                input_file, header=0, engine="openpyxl", index_col=0
            ).rename({"Unnamed: 1": "SUBSTRATE"}, axis=1)
            df = df.rename(lambda x: x.upper(), axis=1)

            entries = []
            for row_num, row in df.iterrows():
                sub = row["SMILES"]

                # Break if we don't have a substrate in this row
                if not isinstance(sub, str):
                    print(f"Moving onto next row from substrate: {sub}")
                    continue

                # If not chiral, remove stereochemistry
                if not chiral:
                    m = Chem.MolFromSmiles(sub)
                    # Make an unchiral version
                    RemoveStereochemistry(m)
                    sub = Chem.MolToSmiles(m)

                for key, value in row.items():
                    # Make sure it's an actual protein entry
                    if key != "SUBSTRATE" and key != "SMILES" and int(value) != 0:
                        protein_seq = seq_dict[key]
                        entry = {
                            "SEQ": protein_seq,
                            "SUBSTRATES": sub,
                            "Activity": int(value),
                        }
                        entries.append(entry)

            img_dir = os.path.join(outdir, "plots")
            os.makedirs(img_dir, exist_ok=True)

            # Send these to outfile
            # For categorical:
            # 3 is good, 2 is eh, 1 is bad
            for output_type in ["categorical", "binary"]:
                # Write these out to a categorical version
                base_df = pd.DataFrame(entries)
                # Take the maximum of equivalently named psespecies
                base_df = base_df.groupby(["SEQ", "SUBSTRATES"]).max().reset_index()

                # If we have a binary matrix
                if output_type == "binary":
                    # Write these out to a binary version
                    base_df["Activity"] = [
                        0 if i == 1 else 1 for i in base_df["Activity"]
                    ]

                name = f"{base_name}_{output_type}"
                out = os.path.join(outdir, f"{name}.csv")
                base_df.to_csv(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="outdir", default="data/processed/")
    parser.add_argument(
        "--sequence-file",
        type=str,
        help="Sequence name file as a fasta doc",
        default="data/raw/gts/sequence_file.txt",
    )
    parser.add_argument(
        "--data-file-acceptors",
        type=str,
        help="Data file containing the acceptors info",
        default="data/raw/gts/gt_acceptor_data.xlsx",
    )
    parser.add_argument(
        "--data-file-donors",
        type=str,
        help="Data file containing donors info",
        default="data/raw/gts/gt_donor_data.xlsx",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    main(
        args.outdir, args.data_file_acceptors, args.data_file_donors, args.sequence_file
    )
