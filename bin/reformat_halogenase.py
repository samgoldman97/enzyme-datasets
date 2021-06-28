""" Extract phosphatase dataset from Fisher et al. (2019)
"""
import os
import argparse
import itertools
from collections import defaultdict
import pandas as pd
import numpy as np

# Set in the paper
BINARY_THRESH = 0.08
SEQ_CUTOFF = 1000


def main(outdir, conversion_file, smiles_mapping_file, ssn_file, soluble_proteins_file):
    """Execute program logic"""

    ### Extract info from smiles_mapping_file
    df_subs = pd.read_excel(smiles_mapping_file, engine="openpyxl")

    unhalogenated_smiles = []
    name_to_smiles = {}
    all_smiles = set()
    for sub_name, smiles, halogenated in zip(
        df_subs["SubstrateName"], df_subs["Smiles"], df_subs["Halogenated"]
    ):
        name_to_smiles[sub_name] = smiles
        unhalogenated_smiles.append(smiles)
        all_smiles.add(smiles)

    ### Extract sequences
    df_seqs = pd.read_excel(ssn_file, engine="openpyxl")
    id_to_seq = dict(zip(df_seqs["fdh_id"].values, df_seqs["sequence"].values))

    # Pop this item
    # Get soluble protein sequence list
    df_seqs_sol = pd.read_excel(soluble_proteins_file, engine="openpyxl")
    all_seqs = set(list([id_to_seq[i] for i in df_seqs_sol["soluble_ids"].values]))
    all_seq_ids = set(list([i for i in df_seqs_sol["soluble_ids"].values]))

    ### Open conversion file and writeout
    df_conversion = pd.read_excel(conversion_file, engine="openpyxl")
    all_halides = set(df_conversion["halide"].values)
    all_seq_ids_in_conversion = set(i for i in df_conversion["fdh_id"].values)

    data_points_to_add = set(itertools.product(all_smiles, all_seqs, all_halides))

    true_data_points = []
    for idx, fdh_id, substrate, halide, conversion in df_conversion.itertuples():
        seq = id_to_seq[fdh_id]
        smiles = name_to_smiles[substrate]
        true_data_points.append(
            {
                "SEQ": seq,
                "SUBSTRATES": smiles,
                "Halide": halide,
                "Conversion": conversion,
            }
        )
        data_points_to_add.discard((smiles, seq, halide))

    # Now add zeros for all the others
    for smiles, seq, halide in data_points_to_add:
        new_entry = {
            "SEQ": seq,
            "SUBSTRATES": smiles,
            "Halide": halide,
            "Conversion": 0,
        }
        true_data_points.append(new_entry)

    #### Add unlisted examples
    # Let's make a dataframe
    df = pd.DataFrame(list(true_data_points))

    # Filter out by length
    filter_out_set = set([j for j in pd.unique(df["SEQ"]) if len(j) > SEQ_CUTOFF])

    # Num filtered out
    print("Num seqs filtered out by length: ", len(filter_out_set))
    bools = [i not in filter_out_set for i in df["SEQ"]]

    # Reset df
    df = df[bools].reset_index(drop=True)

    # Filter out sequences if they have _never_ been predicted positively > 5%.
    # This will control for poor solubility
    # (Could be a solubility problem)
    filter_set = set()
    for enzyme, subdf in df.groupby("SEQ"):
        max_conv = np.max(subdf["Conversion"])
        if max_conv <= 0.05:
            filter_set.add(enzyme)

    print("Num seqs filtered out by heatmap criteria: ", len(filter_set))
    bools = [i not in filter_set for i in df["SEQ"]]

    # Reset df
    df = df[bools].reset_index(drop=True)

    for halide, df_sub in df.groupby("Halide"):
        output_df = pd.DataFrame(
            [
                {"SEQ": seq, "SUBSTRATES": sub, f"Conversion_{halide}": conversion}
                for _, seq, sub, halide_, conversion in df_sub.itertuples()
            ]
        )

        # Average these
        output_df = output_df.groupby(["SEQ", "SUBSTRATES"]).mean().reset_index()
        output_df_binary = output_df.copy()
        output_df_binary[f"Conversion_{halide}"] = [
            1 if i > BINARY_THRESH else 0 for i in output_df[f"Conversion_{halide}"]
        ]

        # Export to data frame
        out_regr = os.path.join(outdir, f"halogenase_{halide}.csv")
        out_bin = os.path.join(outdir, f"halogenase_{halide}_binary.csv")

        # temp_vals = output_df_binary[f"Conversion_{halide}"].values
        # temporary_val = output_df_binary[(temp_vals != 0) & (temp_vals != 1]

        output_df_binary.to_csv(out_bin)
        output_df.to_csv(out_regr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="outdir", default="data/processed/")
    parser.add_argument(
        "--conversion-file",
        type=str,
        help="Name of conversion file",
        default="data/raw/halogenation/oc9b00835_si_004.xlsx",
    )
    parser.add_argument(
        "--smiles-mapping-file",
        type=str,
        help="Manually annotated common name to smiles string, tab separated. This was collected and exported from the provided chemdraw files in the halogenase paper",
        default="data/raw/halogenation/cdx_to_smiles.xlsx",
    )
    parser.add_argument(
        "--ssn-file",
        type=str,
        help="Name of the ssn file that conatins the actual sequences tested",
        default="data/raw/halogenation/oc9b00835_si_005.xlsx",
    )
    parser.add_argument(
        "--soluble-proteins-file",
        type=str,
        help="""This is a personally curated list of the FULL 87
                        proteins that are expressed as soluble. I mined this
                        from Fig S19.""",
        default="data/raw/halogenation/soluble_protein_ids.xlsx",
    )
    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    main(**args.__dict__)
