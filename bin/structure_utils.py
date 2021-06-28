"""structure_utils.py

Helper functions to extract active site shells
"""
import os
import tempfile
from Bio import PDB
from urllib import request

AA3 = [
    "ALA",
    "CYS",
    "ASP",
    "GLU",
    "PHE",
    "GLY",
    "HIS",
    "ILE",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "VAL",
    "TRP",
    "TYR",
    "MSE",  # Add selenomethionine
    "KCX",  # Add carboxymalated lysine
    "SEP",  # Add phosphoserine
    "TPO",  # Add phosphothreonine
]


AA1 = list("ACDEFGHIKLMNPQRSTVWYMKST")
THREE_TO_ONE = dict(zip(AA3, AA1))


def get_residue_pairs(
    structure,
    center_atoms: list,
    angstrom_dist: int,
    seq_end: int,
    exclude_chains: list = [],
):
    """Get a list containing tuples (pos, residue) for all residues within angstrom_dist of
    any of the atoms in center_atoms.

    Args:
        structure: Protein structure object
        center_atoms (list) : List of center atoms
        angstrom_dist (int): Angstrom radius from center atoms
        seq_end (int): Last residue position of the sequence in the chain
            Use this to filter out non amino acid residues
        exclude_chains (list): List of chains NOT to include in pairings. use
            this to avoid taking the kinase substrates

    Return:
        List[Tuple[int, str]]
    """

    atoms = PDB.Selection.unfold_entities(structure, "A")
    ns = PDB.NeighborSearch(atoms)
    close_atoms = []
    for center_atom in center_atoms:
        close_atoms.extend(ns.search(center_atom.coord, angstrom_dist, level="R"))

    # Filter for right seq len
    filtered_close = [
        i
        for i in close_atoms
        if i.id[1] < seq_end
        and i.id[1] > 0
        and i.get_full_id()[2] not in exclude_chains
    ]

    # Export residue pairs
    pos_residue_pairs = [
        (i.id[1], THREE_TO_ONE[i.get_resname()])
        for i in filtered_close
        if i.get_resname() != "HOH"
    ]
    pos_residue_pairs = list(set(pos_residue_pairs))

    return pos_residue_pairs


def get_structure(pdb_code: str, savedir: str):
    """Get PDB Structure.

    Args:
        pdb_code (str): PDB code
        savedir (str) : Location to save

    Return:
        Structure object
    """
    cif_parser = PDB.MMCIFParser()
    pdbl = PDB.PDBList()
    pdbl.retrieve_pdb_file(pdb_code, pdir=savedir)
    structure = cif_parser.get_structure(
        pdb_code, filename=os.path.join(savedir, f"{pdb_code}.cif")
    )
    return structure


# Code for helping to select this in pymol separately
# Takes in as input list of tuples
build_query_str = lambda x: "(" + " or ".join([f"resi {i[0]}" for i in x]) + ")"


def get_pdb_sequence(pdb_code: str):
    """Download pdb sequence.

    Download the pdb sequence in fasta form.
    Args:
        pdb_code (str) : PDB Code

    Return:
        str, str, int, int: fasta text, sequence, sequence start, sequence end

    """

    output = request.urlretrieve(
        f"https://www.rcsb.org/fasta/entry/{pdb_code}/download"
    )
    fasta_text = "\n".join([i.strip() for i in open(output[0], "r").readlines()])

    # Split up again into text and sequence
    # Only take first two
    fasta_text, seq = fasta_text.split("\n")[:2]
    seq_start, seq_end = 1, len(seq) + 1
    return fasta_text, seq, seq_start, seq_end


def get_shells(structure, angstrom_dists: list, impt_residue: dict, seq_end: int):
    """get_shells.

    Get the shells of residues in the structure beyond a certain radius from the impt_residue.

    Args:
        structure: Structure
        angstrom_dists (list): List of angstrom distances to test
        impt_residue (dict): Dictionary containing all the important residues to select
        seq_end (int): Sequence end position to use to filter out waters and things

    Return:
        List of num_residues in each shell
        List of List of pairings of (position, residue)
    """
    center_atom_set = set()
    for reference in impt_residue["references"]:
        substrate = structure[reference["structure"]][reference["chain"]][
            reference["key"]
        ]
        center_atoms = PDB.Selection.unfold_entities(substrate, "A")
        center_atom_set.update(set(center_atoms))

    num_residues = []
    outputs = []
    for i in angstrom_dists:
        output = get_residue_pairs(
            structure,
            center_atom_set,
            i,
            seq_end,
            exclude_chains=impt_residue.get("exclude_chains", []),
        )
        outputs.append(output)
        num_residues.append(len(output))
    return num_residues, outputs


def make_shell_residues_plot(
    angstrom_dist: list, num_residues: list, dataset_name: str, savedir: str
):
    """Save shells"""
    import matplotlib.pyplot as plt

    savedir = os.path.join(savedir, "output_plots")
    os.makedirs(savedir, exist_ok=True)
    save_name = os.path.join(savedir, f"{dataset_name}_residues_plot.png")
    plt.plot(angstrom_dist, num_residues)
    plt.xlabel("Angstrom Radius from References")
    plt.ylabel("Number of Residues")
    plt.title(f"{dataset_name}")
    plt.savefig(save_name, bbox_inches="tight")
    plt.close()


def verify_shell(outputs, seq, substitutions):
    """Verify that all the residues match the sequence downloaded

    Args:
        outputs: List of export outputs
        seq: Sequence to compare against
        substitutions: List of tuples containing (old residue, pos, new residue)
            Note: these are 1 indexed
        offset (int): Integer offset to subtract from all shell positions

    """

    # Force this to iterate at least once
    if len(substitutions) == 0:
        substitutions = [(None, None, None)]

    # New seq
    seq_list = list(seq)
    # Use 1 indexed positions
    for old_residue, pos, new_residue in substitutions:
        for output_index in range(len(outputs)):
            angstrom_shell = outputs[output_index]

            # i is position with 1 index
            # j is residue
            for shell_index in range(len(angstrom_shell)):
                i, j = angstrom_shell[shell_index]
                try:
                    if i == pos:

                        assert j == old_residue
                        assert j == seq[i - 1]

                        # Now make substitution
                        # At both sequence position and pairing
                        seq_list[i - 1] = new_residue
                        angstrom_shell[shell_index] = (pos, new_residue)
                    else:
                        assert j == seq[i - 1]
                except:
                    raise ValueError(f"Failed on residue pair: {i,j}")
    return outputs, "".join(seq_list)


def export_shells(
    output_shells: list,
    angstrom_dists: list,
    output_header: str,
    seq: str,
    dataset_name: str,
    savedir: str,
):
    """export_shells.

    Save the output shells to a file.
    Args:
        output_shells (list): List of lists containing tuples of
            residue num, residue symbol
        angstrom_dists (list): List of distances corresponding to output shell exports
        output_header (str) : Fasta header for the reference export
        seq (str) : Fasta seq for the reference export
        dataset_name (str) : Name of dataset
        savedir (str): Save directory
    """
    savedir = os.path.join(savedir, "structure_references")
    os.makedirs(savedir, exist_ok=True)
    for output_shell, angstrom_dist in zip(output_shells, angstrom_dists):
        out_file = os.path.join(
            savedir, f"{dataset_name}_reference_{int(angstrom_dist)}.txt"
        )
        with open(out_file, "w") as fp:
            fp.write(output_header + "\n")
            fp.write(seq + "\n")
            for i, j in output_shell:
                fp.write(f"({i},'{j}')\n")


def run_extraction(
    impt_residue: dict, angstrom_dists: list, savedir: str, plot: bool = False
):
    """run_extraction.

    Args:
        impt_residue (dict): Dictionary containing key info about
            the shells to extract. See main for an example
        angstrom_dists (list): Angstrom distances to compute shells for
        savedir (str): Directory to save results to
        plot (bool): If true, create a plot showing num residues at each ang
            dist.


    Return:
        Returns structure and outputs list
    """

    with tempfile.TemporaryDirectory() as tempdir:
        structure = get_structure(impt_residue["pdb_code"], tempdir)

    fasta_text, seq, seq_start, seq_end = get_pdb_sequence(impt_residue["pdb_code"])

    # Optionally add an override to sequence in case the pdb file is messed up
    if impt_residue.get("seq", False):
        seq = impt_residue.get("seq")

    # Get seq offset
    seq_offset = impt_residue.get("seq_offset", 0)

    # Add seq offset
    seq_end += seq_offset

    num_residues, outputs = get_shells(structure, angstrom_dists, impt_residue, seq_end)

    ### Shift all outputs by subtracting the seq offset
    outputs = [[(i - seq_offset, j) for i, j in output] for output in outputs]

    if plot:
        make_shell_residues_plot(
            angstrom_dists, num_residues, impt_residue["dataset_name"], savedir
        )
    outputs, seq = verify_shell(outputs, seq, impt_residue["substitutions"])
    export_shells(
        outputs, angstrom_dists, fasta_text, seq, impt_residue["dataset_name"], savedir
    )
    return structure, outputs


### Example
if __name__ == "__main__":
    impt_residue = {
        "dataset_name": "phosphatase",
        "pdb_code": "3l8e",
        "references": [
            {"structure": 0, "chain": "A", "key": ("H_ACY", 901, " ")},
        ],
        "seq": "MAKSVPAIFLDRDGTINVDHGYVHEIDNFEFIDGVIDAMRELKKMGFALVVVTNQSGIARGKFTEAQFETLTEWMDWSLADRDVDLDGIYYCPHHPQGSVEEFRQVCDCRKPHPGMLLSARDYLHIDMAASYMVGDKLEDMQAAVAANVGTKVLVRTGKPITPEAENAADWVLNSLADLPQAIKKQQ",
        "substitutions": [],
    }

    angstrom_dists = np.arange(1, 30)
    savedir = "pdb_exports"
    os.makedirs(savedir, exist_ok=True)

    structure, outputs = run_extraction(impt_residue, angstrom_dists, savedir)
