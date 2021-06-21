""" Reformat oleA (thiolase) dataset """
import os
import argparse
import re
from collections import defaultdict

import pandas as pd

import parse_utils
import file_utils

ORG_RE = r"^.+\.1"  # \_(.*$)"
XC_RE = r"4KU5_"

MANUAL_SUBSTITUTIONS = {
    "Xanthomonas campestris": "XC",
    "Kocuria varians strain G6": "Kocuria varians",
    "Halobacteriovorax marinus strain BE01": "Halobacteriovorax marinus",
    "Xanthomonas translucens pv. graminis": "Xanthomonas translucens",
    "Enhygromyxa salina strain SWB007": "Enhygromyxa salina",
    "Pseudoxanthomonas sp. NML171200": "Pseudoxanthomonas",
}

TAXON_SUBS = {"XC": "Xanthomonas campestris"}

COMPOUND_SUBS = {"oxidazole": "oxadiazole"}

THRESH = 1e-8

# Account for His tag
SEQ_SUBS = {
    "MGSSHHHHHHSSGLVPRGSHMLFQNVSIAGLAHIDAPHTLTSKEINERLQPTYDRLGIKTDVLGDVAGIHARRLWDQDVQASDAATQAARKALIDANIGIEKIGLLINTSVSRDYLEPSTASIVSGNLGVSDHCMTFDVANASLAFINGMDIAARMLERGEIDYALVVDGETANLVYEKTLERMTSPDVTEEEFRNELAALTLGCGAAAMVMARSELVPDAPRYKGGVTRSATEWNKLCRGNLDRMVTDTRLLLIEGIKLAQKTFVAAKQVLGWAVEELDQFVIHQVSRPHTAAFVKSFGIDPAKVMTIFGEHGNIGPASVPIVLSKLKELGRLKKGDRIALLGIGSGLNCSMAEVVW": "MLFQNVSIAGLAHIDAPHTLTSKEINERLQPTYDRLGIKTDVLGDVAGIHARRLWDQDVQASDAATQAARKALIDANIGIEKIGLLINTSVSRDYLEPSTASIVSGNLGVSDHCMTFDVANACLAFINGMDIAARMLERGEIDYALVVDGETANLVYEKTLERMTSPDVTEEEFRNELAALTLGCGAAAMVMARSELVPDAPRYKGGVTRSATEWNKLCRGNLDRMVTDTRLLLIEGIKLAQKTFVAAKQVLGWAVEELDQFVIHQVSRPHTAAFVKSFGIDPAKVMTIFGEHGNIGPASVPIVLSKLKELGRLKKGDRIALLGIGSGLNCSMAEVVW",
}


def get_args():
    """Get args"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compound-file",
        help="Compound file name",
        default="data/raw/OleA/final_15_cmpnds_tested.csv",
    )
    parser.add_argument(
        "--activity-file",
        help="Activity file",
        default="data/raw/OleA/enzyme_substrate_activity_data.csv",
    )
    parser.add_argument(
        "--protein-file",
        help="Protein file",
        default="data/raw/OleA/73_OleA_JGI_unaligned.fasta.txt",
    )
    parser.add_argument(
        "--compound-feature-file",
        help="Feature file",
        default="data/raw/OleA/all_molecular_properties.csv",
    )
    parser.add_argument(
        "--prot-feature-file",
        help="Feature file",
        default="data/raw/OleA/73_overall_calculated_protein_properties.csv",
    )
    parser.add_argument(
        "--taxon-file",
        help="Taxon containing file -- use for extracting class",
        default="data/raw/OleA/OleA_taxonomic_classification.csv",
    )
    parser.add_argument("--out-dir", help="Out prefix", default="data/processed/")
    parser.add_argument(
        "--raw-dir", help="Raw directory prefix", default="data/raw/OleA/"
    )

    return parser.parse_args()


def main(args):
    """Main method"""
    compounds = pd.read_csv(args.compound_file)
    compound_to_smiles = dict(
        zip(compounds["cmpnd_abbrev"].values, compounds["SMILES"].values)
    )
    for compound_sub, compound_targ in COMPOUND_SUBS.items():
        compound_to_smiles[compound_targ] = compound_to_smiles[compound_sub]
        # Remove this because it's mislabeled
        compound_to_smiles.pop(compound_sub)

    activity_df = pd.read_csv(args.activity_file)

    prot_dict = {}
    for header, seq in parse_utils.fasta_iter(args.protein_file):
        if "\ufeff" in header:
            header = header.replace("\ufeff", "")

        header = re.sub(ORG_RE, "", header)
        header = re.sub(XC_RE, "", header)
        header = re.sub("\_", " ", header)
        header = header.strip()

        # Manual substitutions:
        for sub, val in MANUAL_SUBSTITUTIONS.items():
            header = re.sub(sub, val, header)
        prot_dict[header] = seq
        if seq in SEQ_SUBS:
            prot_dict[header] = SEQ_SUBS[seq]

    # Create taxon mapping
    tax_df = pd.read_csv(args.taxon_file)
    # Map taxon to sequence
    tax_mapping = defaultdict(lambda: set())
    for genus, class_ in tax_df[["genus", "class"]].values:
        for prot_key in prot_dict.keys():
            # Get substitution with default being the queyr
            if genus in TAXON_SUBS.get(prot_key, prot_key):
                sequence = prot_dict[prot_key]
                tax_mapping[class_].add(sequence)

    sequences = set([j for i in tax_mapping.values() for j in i])
    unlabeled_taxon_seq = set(prot_dict.values()).difference(sequences)
    unlabeled_taxons = [k for k, v in prot_dict.items() if v in unlabeled_taxon_seq]

    file_utils.pickle_obj(
        {i: list(j) for i, j in tax_mapping.items()},
        os.path.join(args.raw_dir, f"OleA_taxon_split.p"),
    )

    org_keys = set(pd.unique(activity_df["org"]))
    prot_keys = set(prot_dict.keys())
    # print(org_keys.difference(prot_keys))
    # print(prot_keys.difference(org_keys))

    smiles_keys = set(compound_to_smiles.keys())
    cmpd_keys = set(pd.unique(activity_df["cmpnd"]))
    no_activity_compounds = smiles_keys.difference(cmpd_keys)
    no_activity_smiles = [compound_to_smiles[i] for i in no_activity_compounds]

    # print(cmpd_keys.difference(smiles_keys))

    # Follow chemprop nomenclature
    # seq, rxn, [predict_columns]
    rows = []
    for index, row in activity_df.iterrows():
        new_row = {
            "SEQ": prot_dict[row["org"]].strip(),
            "SUBSTRATES": compound_to_smiles[row["cmpnd"]].strip(),
            "log_slope": row["log_slope"],
        }
        rows.append(new_row)
    df = pd.DataFrame(rows)
    # Average
    df = df.groupby(["SEQ", "SUBSTRATES"]).mean().reset_index()
    df.to_csv(os.path.join(args.out_dir, f"olea.csv"))

    df["log_slope"] = [1 if i > THRESH else 0 for i in df["log_slope"]]
    df["log_slope"] = [1 if i > THRESH else 0 for i in df["log_slope"]]

    # Create classification table
    pivoted = df.pivot_table(
        values="log_slope", index="SEQ", columns="SUBSTRATES", fill_value=0
    )
    pivoted[pivoted > THRESH] = 1
    for j in no_activity_smiles:
        pivoted[j] = 0

    unpivoted = pd.melt(pivoted.reset_index(), id_vars="SEQ", value_name="Activity")
    unpivoted.to_csv(os.path.join(args.out_dir, f"olea_binary.csv"))

    ## Making feature file
    if args.compound_feature_file is not None:
        df = pd.read_csv(args.compound_feature_file)
        feat_dict = {}
        for index, row in df.iterrows():
            smiles = row.iloc[2]
            feats = row.iloc[3:].values
            feat_dict[smiles] = feats
        file_utils.pickle_obj(
            feat_dict, os.path.join(args.raw_dir, f"olea_compound_feats.p")
        )

    ## Making feature file
    if args.prot_feature_file is not None:
        df = pd.read_csv(args.prot_feature_file)
        feat_dict = {}
        for index, row in df.iterrows():
            seq = row.iloc[1]
            feats = row.iloc[4:].values
            feat_dict[seq] = feats
        file_utils.pickle_obj(
            feat_dict, os.path.join(args.raw_dir, f"olea_prot_feats.p")
        )


if __name__ == "__main__":
    args = get_args()
    file_utils.make_dir(args.out_dir)
    main(args)
