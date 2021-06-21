""" Reformat davis kinase data as prepared in hie et al. cell systems paper 

Pipeline to trim domains:

1. Download pfam seed PF000069_seed.txt
2. hmmbuild ~/Desktop/output_hmm_00069 PF00069_seed.txt
3. hmmsearch --domtblout domains_found_00069.txt ~/Desktop/output_hmm_00069 ~/Desktop/uniprot_sequences_processed.fasta
4. Manually edit these results such that if they sepcify domain 1 or domain 2,
we only keep one
5. If it's clear that htere's some insertion inside the envelopes, we should
include this insertion as well


"""
import os
import pandas as pd
from itertools import groupby
import re
import numpy as np

from enzpred.utils import parse_utils


raw_folder = "data/raw/davis/"

path_wrap = lambda x : os.path.join(raw_folder, x)
data_file = path_wrap("nbt.1990-S4.csv")
smiles_file = path_wrap("chem_smiles.csv")
domains_out = path_wrap("domains_found_00069_edited_manually.txt")
seq_file_preprocessed = path_wrap("uniprot_sequences.fasta")
seq_file_processed = path_wrap("uniprot_sequences_processed.fasta")
out_name = "data/processed/dense/davis.csv"
out_name_filter = "data/processed/dense/davis_filtered.csv"

MAX_VAL = 10000. 
SUB_RE = r"([A-Z]+)([0-9]+)([A-Z]+)"
DEL_RE = r"((?:[A-Z]+[0-9]+-*)+)del"
#DEL_RE = r"([A-Z]+)([0-9]+)del"
#DEL_RE = r"(([A-Z]+)([0-9]+)-*)*del"

def convert_kinase_to_seq(genes, proteins, seq_dict):
    """ Convert the proteins passed into amino acid sequences"""

    output_dict = {}
    for gene, prot  in zip(genes, proteins): 
        seq_returned = seq_dict.get(gene, None)
        if seq_returned is None: 
            seq_returned = seq_dict.get(prot, None)
        if seq_returned is None:
            raise ValueError(gene, prot)

        initial_seq = seq_returned
        split_prot = prot.split("(", 1)
        if len(split_prot) > 1:
            # Look for substitutions
            sub_str = split_prot[1]
            groups = re.findall(SUB_RE, sub_str)
            if groups is not None and len(groups) > 0:
                for group in groups: 
                    start, pos, sub = group
                    pos = int(pos)
                    pos -= 1
                    assert(seq_returned[pos] == start)
                    temp_seq = list(seq_returned)
                    temp_seq[pos] = sub
                    seq_returned = "".join(temp_seq)

            # Look for substitutions
            groups = re.search(DEL_RE, sub_str)
            if groups is not None:
                group = groups.groups()[0]
                new_indices = list(range(len(seq_returned)))
                seq_ar = np.array(list(seq_returned))
                for del_instance in group.split("-"):   
                    start = del_instance[0]
                    pos = int(del_instance[1:])
                    pos -= 1
                    assert(seq_ar[pos] == start)
                    new_indices.remove(pos)
                seq_returned = "".join(seq_ar[new_indices])
        # Remove handling of CDK 
        if prot.startswith("CDK4"): 
            pass
            #if prot.endswith('cyclinD1'): 
            #    seq_returned = "".join([seq_returned, seq_dict['CCND1']])
            #elif prot.endswith("cyclinD1"): 
            #    seq_returned = "".join([seq_returned, seq_dict['CCND3']])
            #else:
            #    raise ValueError()
        output_dict[prot] = seq_returned

    return output_dict

def create_processed_seqs():

    data_df = pd.read_csv(data_file)
    seq_dict = {i : j for i,j in parse_utils.fasta_iter(seq_file_preprocessed)}
    extract_header = lambda  x : re.search("GN=([0-9,A-Z,a-z,\.]*) ", x ).groups()[0]

    seq_dict = {extract_header(k) : v for k,v  in seq_dict.items()}

    genes = data_df["Entrez Gene Symbol"].values
    proteins = data_df["Kinase"].values

    output_dict = convert_kinase_to_seq(genes, proteins, seq_dict)
    out_fasta = "\n".join([f">{prot}\n{seq}" for prot,seq in output_dict.items()])
    with open(seq_file_processed, "w") as fp:
        fp.write(out_fasta)


def trim_domains(seq_dict, domains_out):
    """ Trim all the domains down given a domain mapping table"""

    new_dict = {}
    with open(domains_out, "r") as fp:
        for line in fp: 
            line_orig = line
            if line.startswith("#"): 
                continue
            line = line.strip().split()
            env_start = int(line[-4])
            env_end = int(line[-3])
            prot_name = line[0]

            old_len = len(seq_dict[prot_name])

            new_dict[prot_name] = seq_dict[prot_name][env_start : env_end+1]
            new_len = len(new_dict[prot_name])
            ## Pop
            if new_len < 200:
                new_dict.pop(prot_name)
    return new_dict

def create_new_table():
    """ Create new table"""


    seq_dict = {i : j for i,j in parse_utils.fasta_iter(seq_file_processed)}
    seq_dict = trim_domains(seq_dict, domains_out)

    #len_dist = [len(i) for i in seq_dict.values()]
    #print(f"Avg length: {np.mean(len_dist)}")
    #print(f"Num entries: {len(len_dist)}")

    #import matplotlib.pyplot as plt
    #plt.hist(len_dist, bins = 50)
    #plt.savefig("temp.png")

    for filter_data in [False, True]: 


        has_sub = lambda x: len(re.findall(SUB_RE, x)) > 0
        has_del = lambda x: len(re.findall(DEL_RE, x)) > 0

        if filter_data: 
            seq_dict = {key : val for key, val in seq_dict.items() 
                        if ((not has_sub(key)) and (not has_del(key)))}

        #from Levenshtein import distance
        #import matplotlib.pyplot as plt
        #all_seqs = list(seq_dict.values())
        #dist_mat = [distance(i,j) for index, i in enumerate(all_seqs[:-1]) for j in all_seqs[index + 1:]]
        #plt.hist(dist_mat, bins=50)
        #plt.savefig(f"/Users/samgoldman/Desktop/temp_{filter_data}.png")
        #pdb.set_trace()

        data_df = pd.read_csv(data_file)
        smiles_df = pd.read_csv(smiles_file)


        comp_to_smi = dict(zip(smiles_df["Compound Name"].values, 
                               smiles_df["SMILES"].values))
        proteins = data_df["Kinase"].values

        substrate_columns = set(data_df.columns)
        # Keep "Kinase" in columns
        substrate_columns = substrate_columns.difference(["Entrez Gene Symbol",
                                                          "Accession Number"])

        entries = []
        melted = data_df.loc[:, substrate_columns].melt(id_vars=["Kinase"], var_name="Compound", value_name="Val")
        for kinase, compound, val in melted[["Kinase", "Compound", "Val"]].values:

            seq = seq_dict.get(kinase, None)
            if seq is None: 
                continue

            smi = comp_to_smi[compound]
            if pd.isna(val):
                val = MAX_VAL 
            val = float(val)

            val = -np.log10(val)
            entry = {"SEQ": seq, "SUBSTRATES" : smi, "-Log10Kd" : val}
            entries.append(entry)

        df_out = pd.DataFrame(entries)
        print(f"Num entries: {len(df_out)}")
        # Get best measured kd by collapsing
        df_out = df_out.groupby(["SEQ", "SUBSTRATES"]).min().reset_index()

        print(f"Mean: {np.mean(df_out['-Log10Kd'])}")
        print(f"Std: {np.std(df_out['-Log10Kd'])}")
        print(f"Num entries: {len(df_out)}")
        if filter_data: 
            df_out.to_csv(out_name_filter)
        else:
            df_out.to_csv(out_name)

def main(): 
    # Create processed uniprot seqs
    # Once this is created, use hmmer profile to manually extract desired
    # sequences and post process
    #create_processed_seqs()

    # Now extract values from table
    create_new_table()

if __name__=="__main__": 
    main()
