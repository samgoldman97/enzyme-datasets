"""Extract structures and alignments. 

In our analysis, we are interested in using a reference structure and also
alignments for each sequence. This script extracts these alignments and
structure references for the datasets considered.

Note: MuscleCommandLine must be installed and enabled 
"""

import os
import sys
from Bio import PDB
from Bio.Align.Applications import MuscleCommandline
import tempfile

import pandas as pd
import numpy as np

from structure_utils import *


IMPT_RESIDUES = [
    {
        "dataset_name": "esterase",
        "pdb_code": "5a6v",
        "references": [
            {"structure": 0, "chain": "A", "key": (" ", 105, " ")},
            {"structure": 0, "chain": "A", "key": (" ", 187, " ")},
            {"structure": 0, "chain": "A", "key": (" ", 224, " ")},
        ],
        "exclude_chains": ["B", "C", "D"],
        "substitutions": [],
    },
    {
        "dataset_name": "davis",
        "pdb_code": "2CN5",
        "references": [
            {"structure": 0, "chain": "A", "key": ("H_ADP", 1509, " ")},
        ],
        "seq": "HMSVYPKALRDEYIMSKTLGSGACGEVKLAFERKTCKKVAIKIISKRKFAIGSAREADPALNVETEIEILKKLNHPCIIKIKNFFDAEDYYIVLELMEGGELFDKVVGNKRLKEATCKLYFYQMLLAVQYLHENGIIHRDLKPENVLLSSQEEDCLIKITDFGHSKILGETSLMRTLCGTPTYLAPEVLVSVGTAGYNRAVDCWSLGVILFICLSGYPPFSEHRTQVSLKDQITSGKYNFIPEVWAEVSEKALDLVKKLLVVDPKARFTTEEALRHPWLQDEDMKRKFQDLLSEENESTALPQVLAQPSTSRKRPREGEAEGAE",
        "seq_offset": 207,
        "substitutions": [],
    },
    {
        "dataset_name": "aminotransferase",
        "pdb_code": "3QPG",
        "references": [
            {
                "structure": 0,  # C169, E53, K135
                "chain": "A",
                "key": ("H_3QP", 422, " "),
            },
        ],
        "exclude_chains": ["B"],  # Exclude substrate chain B
        "substitutions": [],
        "seq": "HHHHHHHHHHHHMFENITAAPADPILGLADLFRADERPGKINLGIGVYKDETGKTPVLTSVKKAEQYLLENETTKNYLGIDGIPEFGRCTQELLFGKGSALINDKRARTAQTPGGTGALRVAADFLAKNTSVKRVWVSNPSWPNHKSVFNSAGLEVREYAYYDAENHTLDFDALINSLNEAQAGDVVLFHGCCHNPTGIDPTLEQWQTLAQLSVEKGWLPLFDFAYQGFARGLEEDAEGLRAFAAMHKELIVASSYSKNFGLYNERVGACTLVAADSETVDRAFSQMKAAIRANYSNPPAHGASVVATILSNDALRAIWEQELTDMRQRIQRMRQLFVNTLQEKGANRDFSFIIKQNGMFSFSGLTKEQVLRLREEFGVYAVASGRVNVAGMTPDNMAPLCEAIVAVL",
    },
    {
        "dataset_name": "nitrilase",
        "pdb_code": "3WUY",
        "references": [
            {"structure": 0, "chain": "A", "key": (" ", 53, " ")},  # C169, E53, K135
            {"structure": 0, "chain": "A", "key": (" ", 135, " ")},
            {"structure": 0, "chain": "A", "key": (" ", 169, " ")},
        ],
        "exclude_chains": [],  # Exclude substrate chain B
        "substitutions": [],
        "seq": "MLGKIMLNYTKNIRAAAAQISPVLFSQQGTMEKVLDAIANAAKKGVELIVFPETFVPYYPYFSFVEPPVLMGKSHLKLYQEAVTVPGKVTQAIAQAAKTHGMVVVLGVNEREEGSLYNTQLIFDADGALVLKRRKITPTYHERMVWGQGDGAGLRTVDTTVGRLGALACWEHYNPLARYALMAQHEQIHCGQFPGSMVGQIFADQMEVTMRHHALESGCFVINATGWLTAEQKLQITTDEKMHQALSGGCYTAIISPEGKHLCEPIAEGEGLAIADLDFSLIAKRKRMMDSVGHYARPDLLQLTLNNQPWSALEANPVTPNAIPAVSDPELTETIEALPNNPIFSH",
    },
    {
        "dataset_name": "phosphatase",
        "pdb_code": "3l8e",
        "references": [
            {"structure": 0, "chain": "A", "key": ("H_ACY", 901, " ")},
        ],
        "seq": "MAKSVPAIFLDRDGTINVDHGYVHEIDNFEFIDGVIDAMRELKKMGFALVVVTNQSGIARGKFTEAQFETLTEWMDWSLADRDVDLDGIYYCPHHPQGSVEEFRQVCDCRKPHPGMLLSARDYLHIDMAASYMVGDKLEDMQAAVAANVGTKVLVRTGKPITPEAENAADWVLNSLADLPQAIKKQQ",
        "substitutions": [],
    },
    {
        "dataset_name": "halogenase",
        "pdb_code": "2AR8",
        "references": [{"structure": 0, "chain": "A", "key": ("H_CTE", 650, " ")}],
        "substitutions": [],
    },
    {
        "dataset_name": "olea",
        "pdb_code": "4KU5",
        "references": [{"structure": 0, "chain": "A", "key": (" ", 143, " ")}],
        "substitutions": [("S", 143, "C")],
    },
    {
        "dataset_name": "duf",
        "pdb_code": "2Y7F",
        "references": [{"structure": 0, "chain": "A", "key": ("H_KMH", 1276, " ")}],
        "substitutions": [],
        "seq": "HEPLILTAAITGAETTRADQPNLPITPEEQAKEAKACFEAGARVIHLHIREDDGRPSQRLDRFQEAISAIREVVPEIIIQISTGGAVGESFDKRLAPLALKPEMATLNAGTLNFGDDIFINHPADIIRLAEAFKQYNVVPEVEVYESGMVDAVARLIKKGIITQNPLHIQFVLGVPGGMSGKPKNLMYMMEHLKEEIPTATWAVAGIGRWHIPTSLIAMVTGGHIRCGFEDNIFYHKGVIAESNAQLVARLARIAKEIGRPLATPEQAREILALN",
    },
    {
        "dataset_name": "gt",
        "pdb_code": "3HBF",
        "references": [
            {"structure": 0, "chain": "A", "key": ("H_UDP", 900, " ")},
            {"structure": 0, "chain": "A", "key": ("H_MYC", 901, " ")},
        ],
        "substitutions": [],
    },
]


def get_temp_fasta(seqs, ref_seq):
    """Get a temporary fasta file to align with a list of sequences and a ref seq"""
    # Make fasta file to align
    seqs_fasta = tempfile.NamedTemporaryFile(mode="w+b")
    fasta_text = "\n".join([f">seq{index}\n{seq}" for index, seq in enumerate(seqs)])
    fasta_text = f">reference\n{ref_seq}\n" + fasta_text
    seqs_fasta.write(fasta_text.encode())
    seqs_fasta.seek(0)
    return seqs_fasta


def create_MSA(dataset_dir):
    """create_MSAs."""
    alignment_savedir = os.path.join(dataset_dir, "alignments")
    os.makedirs(alignment_savedir, exist_ok=True)

    for impt_residue in IMPT_RESIDUES:
        data_name = impt_residue["dataset_name"]

        # Get ref seq
        lines = open(
            os.path.join(
                dataset_dir, f"structure_references/{data_name}_reference_1.txt"
            ),
            "r",
        ).readlines()
        ref_seq = lines[1].strip()

        # Get all sequences
        data_files = [
            os.path.join(dataset_dir, i)
            for i in os.listdir(dataset_dir)
            if (f"{data_name}.csv" in i or f"{data_name}_" in i)
            and not i.startswith(".")
        ]

        # sort them to get a key
        seqs = sorted(list(pd.read_csv(data_files[0])["SEQ"].unique()))

        seqs_fasta = get_temp_fasta(seqs, ref_seq)

        cline = MuscleCommandline(
            cmd="muscle",
            input=seqs_fasta.name,
            out=os.path.join(alignment_savedir, f"{data_name}_alignment.fasta"),
        )
        cline()


def create_markdown_table():
    """Create a markdown table with results"""

    df = pd.DataFrame(IMPT_RESIDUES)
    df["Dataset"] = df["dataset_name"]
    df["PDB ID"] = df["pdb_code"]
    markdown = df[["Dataset", "PDB ID"]].to_markdown()
    open("markdown_temp.md", "w").write(markdown)


def extract_structure_references(savedir):
    """extract_structure_references."""
    angstrom_dists = np.arange(1, 50)
    os.makedirs(savedir, exist_ok=True)
    for impt_residue in IMPT_RESIDUES:
        structure, outputs = run_extraction(impt_residue, angstrom_dists, savedir)


if __name__ == "__main__":
    create_markdown_table()
    dataset_dir = "data/processed/"
    extract_structure_references(dataset_dir)
    create_MSA(dataset_dir)
