""" parse_utils.py

Utils for parsing file formats

"""
from itertools import groupby
from typing import Iterable, Tuple


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
