""" parse_utils.py

Utils for parsing file formats

"""

import re
import numpy as np 
from collections import defaultdict 
from scipy import sparse

from itertools import groupby
import typing
from typing import Iterable, Tuple
from tqdm import tqdm

def fasta_iter(fasta_name: str) -> Iterable[Tuple[str, str]]:
    """fasta_iter.

     modified from Brent Pedersen
     Correct Way To Parse A Fasta File In Python
     given a fasta file. yield tuples of header, sequence
     
    Args:
        fasta_name (str): fasta_name

    Returns: 
        Iterable[tuple] of (headerStr, seq)
    """
    # open fasta file
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)


def write_fasta(seq_dict: dict, fasta_name: str) -> None:
    """write_fasta.

     modified from Brent Pedersen
     Correct Way To Parse A Fasta File In Python
     given a fasta file. yield tuples of header, sequence
     
    Args:
        seq_dict (dict): name to sequence map
        fasta_name (str): fasta_name

    Returns: 
    """
    # open fasta file
    fh = open(fasta_name, "w")
    for k, v in seq_dict.items():
        fh.write(f">{k.strip()}\n{v.strip().upper()}\n")
    fh.close()

def get_blast_tsv_dict(tsv_name): 
    """ Extract the blast TSV dict.

    Note: Make this symmetric by average

    """
    # open fasta file
    seq_to_seq_dist= defaultdict(lambda : {})

    print("Starting to read in tsv")
    # Use min to "symmetrize" blast distance 
    with open(tsv_name, "r") as fp: 
        for line_num, line in enumerate(fp): 
            cols = line.split("\t")

            # Percent identity
            e_val = float(cols[-2])
            query = cols[0]
            subject = cols[1]

            # Keep lowest distance!
            if subject in seq_to_seq_dist[query]: 
                old_val = seq_to_seq_dist[query][subject]
                e_val = min(e_val, old_val)

            if query in seq_to_seq_dist[subject]: 
                old_val = seq_to_seq_dist[subject][query]
                e_val = min(e_val, old_val)

            seq_to_seq_dist[query][subject] = e_val
            seq_to_seq_dist[subject][query] = e_val

    return seq_to_seq_dist

def construct_sparse_matrix(seq_to_seq_dist, MIN_N = 5):
    """ Build a sparse matrix """
    ### MIN NEIGHBORS
    total_seqs_read = len(list(seq_to_seq_dist.keys()))
    keep_filtering = True 
    loop_num = 0
    while keep_filtering: 
        print(f"Filter loop number: {loop_num}")
        loop_num += 1
        filter_seqs = [seq for seq in seq_to_seq_dist 
                       if len(seq_to_seq_dist[seq]) < MIN_N]

        if len(filter_seqs) == 0: 
            keep_filtering = False
        else:
            print(f"Num_seqs_removed: {len(filter_seqs)}")
            # First pop top level
            [seq_to_seq_dist.pop(key) for key in filter_seqs]

            # Now pop bottom level 
            print("Popping bottom level")
            [top_dict.pop(filter_key) 
             for key,top_dict in tqdm(seq_to_seq_dist.items())
             for filter_key in filter_seqs if filter_key in top_dict] 

    seq_ids = list(seq_to_seq_dist.keys())
    num_seqs = len(seq_ids)

    print(f"Num total seqs: {total_seqs_read}")
    print(f"Num keep seqs: {num_seqs}")
    seq_to_num = dict(zip(seq_ids, 
                          np.arange(num_seqs).tolist()))

    ### Enforce Zero dist to self
    for seq in seq_ids: 
        seq_to_seq_dist[seq][seq] = 0 

    ### Now construct CSR matrix
    ## Taken from scipy docs for incremental construction
    data = []
    row_ind = []
    col_ind = []
    min_n = num_seqs
    print("Constructing sparse matrix")
    for seq_1 in tqdm(seq_to_seq_dist): 
        num_n = 0 
        for seq_2 in seq_to_seq_dist[seq_1]: 
            num_n += 1
            e_val = seq_to_seq_dist[seq_1][seq_2]
            id_1 = seq_to_num[seq_1]
            id_2 = seq_to_num[seq_2]
            row_ind.append(id_1)
            col_ind.append(id_2)
            data.append(e_val)
        if num_n < min_n: min_n = num_n

    # data = np.array(data)
    # row_ind = np.array(row_ind)
    # col_ind = np.array(col_ind)
    ret_mat = sparse.csr_matrix((data, (row_ind, col_ind)), 
                                shape=(len(seq_ids), 
                                       len(seq_ids)))

    # Return ret_mat, mapping from seq to ind, and 
    # the min number of neighbors per seq
    return ret_mat, seq_to_num, min_n 

def parse_blast_tsv_res(tsv_name : str, MIN_N = 5):
    """parse_blast_tsv_res.
     
    Args:
        tsv_name (str): Name of output tsv file

    Returns: 
    """
    seq_to_seq_dist = get_blast_tsv_dict(tsv_name)
    return construct_sparse_matrix(seq_to_seq_dist, MIN_N)



def parse_ssa_reference(file_name : str) -> Tuple[str, np.ndarray]: 
    """ parse_ssa_reference.
   
    Parse the reference file. 
    
    Args: 
        file_name (str): Name of file to parse

    Return:
        Tuple[str, np.ndarray]: seq and positions

    """

    tuple_regex = "\(([0-9]*),\s*\'([A-Z,?])\'\)"

    with open(file_name, "r") as fp: 
        active_pos = []
        for line_num, line in enumerate(fp): 
            if line_num == 0: continue
            if line_num == 1: 
                seq = line.strip() 
            else: 
                line = line.strip()
                position, residue = re.search(tuple_regex, line).groups()

                # Make this zero indexed
                position = int(position) 
                position -= 1

                assert(residue == seq[position])
                active_pos.append(position)

    return seq, np.array(active_pos)


def map_aligned(aligned_seq : str) -> dict: 
    """ map_aligned.

    Map positions in the unlaigned seq to positions in the aligned seq.

    Args: 
        aligned_seq (str): Seq with gaps

    Return: 
        dict: mapping unaligned seq positions to their corresponding location in
            the aligned seq
    """

    unaligned_index = 0
    mapping = {}
    for index, i in enumerate(aligned_seq):
        if i != "-":
            mapping[unaligned_index] = index
            unaligned_index += 1

    return mapping

def map_unaligned(aligned_seq : str) -> dict: 
    """ map_unaligned.

    Map positions in the aligned seq to positions in the unaligned seq.

    Args: 
        aligned_seq (str): Seq with gaps

    Return: 
        dict: mapping aligned seq positions to their corresponding location in
            the aligned seq. Gaps are not included
    """

    unaligned_index = 0
    mapping = {}
    for index, i in enumerate(aligned_seq):
        if i != "-":
            mapping[index] = unaligned_index
            unaligned_index += 1

    return mapping


def extract_pool_residue_dict(seq_msa : str, ref_seq : str,
                              pool_residues : list, 
                              rand : bool = False, 
                              pool_num : int = 1) -> dict: 
    """extract_pool_residue_dict.

    Args:
        seq_msa (str): Seq MSA file
        ref_seq (str): Unaligned ref seq that must be in seq msa
        pool_residues (str): Residues in the unaligned ref seq to extract
        rand (bool): If true, get random msa positions
        pool_num (int): If rand, this sets how many random positions to choose

    Return: 
        Dict mapping every unaligned sequence in the MSA to the positions in
        its sequence we should pool over based on the reference pooled
        residues.
    """

    # Map seq to seq with msa
    msa_map = {seq.replace('-', '') : seq 
               for seq_name, seq  in fasta_iter(seq_msa)}

    # Set of all positions in the alignment we should pool over
    if rand:
        temp_seq = list(msa_map.values())[0]
        msa_pos_list = np.arange(len(temp_seq))
        alignment_pool_positions = np.random.choice(msa_pos_list, pool_num, 
                                                    replace=False)
    else:
        ref_aligned_mapping = map_aligned(msa_map[ref_seq])
        alignment_pool_positions = [ref_aligned_mapping[j] for j in pool_residues]

    # Now extract these for _all_ sequences in alignment 
    pool_residues_mapping = dict()
    for unaligned_seq, aligned_seq in msa_map.items(): 
        unaligned_mapping = map_unaligned(aligned_seq) 

        # Make sure it's not a gap in sequence of interest!
        pool_set = [unaligned_mapping[pos] for pos in alignment_pool_positions 
                    if pos in unaligned_mapping]
        pool_residues_mapping[unaligned_seq] = pool_set

    return pool_residues_mapping 


def get_col_freqs(msa_map : dict):
    """ Get column frqeuencies"""

    aligned_seqs = list(msa_map.values() )
    align_length = len(aligned_seqs[0])
    # Get frequencies in each column; use 100 rather than converting characters
    # to an alphabet. We can ignore unused positions
    frequencies = np.zeros((align_length, 100))
    for align_seq in msa_map.values():
        for pos_index, pos in enumerate(align_seq): 
            frequencies[pos_index, ord(pos)] += 1
    return frequencies

def extract_coverage_residue_dict(seq_msa : str,
                                  pool_num : int = 1) -> dict: 
    """extract_coverage_residue_dict.

    Args:
        seq_msa (str): Seq MSA file
        pool_num (int): If rand, this sets how many random positions to choose

    Return: 
        Dict mapping every unaligned sequence in the MSA to the positions in
        its sequence we should pool over based on coverage ranking
        residues.
    """

    # Map seq to seq with msa
    msa_map = {seq.replace('-', '') : seq 
               for seq_name, seq  in fasta_iter(seq_msa)}

    frequencies = get_col_freqs(msa_map)

    # Now tally them up and find the columns with fewest gaps
    gap_pos = ord('-')
    sequence_coverage = frequencies[:, gap_pos]
    align_position_ordering = np.argsort(sequence_coverage)

    # get the proper positions
    position_cap = sequence_coverage[align_position_ordering[pool_num - 1]]
    potential_positions = np.argwhere(sequence_coverage <= position_cap).flatten()
    alignment_pool_positions = np.random.choice(potential_positions, pool_num, replace = False)

    #alignment_pool_positions  = align_position_ordering[:pool_num]


    # Now extract these for _all_ sequences in alignment 
    pool_residues_mapping = dict()
    for unaligned_seq, aligned_seq in msa_map.items(): 
        unaligned_mapping = map_unaligned(aligned_seq) 

        # Make sure it's not a gap in sequence of interest!
        pool_set = [unaligned_mapping[pos] for pos in alignment_pool_positions 
                    if pos in unaligned_mapping]
        pool_residues_mapping[unaligned_seq] = pool_set

    return pool_residues_mapping 

def extract_conserve_residue_dict(seq_msa : str,
                                  pool_num : int = 1) -> dict: 
    """extract_conserve_residue_dict.

    Get the residues that have the most conservation

    Args:
        seq_msa (str): Seq MSA file
        pool_num (int): If rand, this sets how many random positions to choose

    Return: 
        Dict mapping every unaligned sequence in the MSA to the positions in
        its sequence we should pool over based on coverage ranking
        residues.
    """

    # Map seq to seq with msa
    msa_map = {seq.replace('-', '') : seq 
               for seq_name, seq  in fasta_iter(seq_msa)}

    frequencies = get_col_freqs(msa_map)

    # Now tally them up and find the columns with fewest gaps
    # Get all frequencies _not_ gaps 
    gap_pos = ord('-')
    bool_ar = np.ones(frequencies.shape[1]) 
    bool_ar[gap_pos] = 0
    frequencies = frequencies[:, bool_ar.astype(bool)]

    # Get the max conservation in any one column
    positional_max = frequencies.max(1)

    align_position_ordering = np.argsort(positional_max)[::-1]

    # get the proper positions
    #alignment_pool_positions  = align_position_ordering[:pool_num]

    # get the proper positions
    position_cap = positional_max[align_position_ordering[pool_num - 1]]
    potential_positions = np.argwhere(positional_max >= position_cap).flatten()
    alignment_pool_positions = np.random.choice(potential_positions, pool_num, replace = False)

    # Now extract these for _all_ sequences in alignment 
    pool_residues_mapping = dict()
    for unaligned_seq, aligned_seq in msa_map.items(): 
        unaligned_mapping = map_unaligned(aligned_seq) 

        # Make sure it's not a gap in sequence of interest!
        pool_set = [unaligned_mapping[pos] for pos in alignment_pool_positions 
                    if pos in unaligned_mapping]
        pool_residues_mapping[unaligned_seq] = pool_set

    return pool_residues_mapping 

