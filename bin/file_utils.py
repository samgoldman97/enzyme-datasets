"""Module containing some helper scripts"""
import pickle
from typing import Any
import os


def pickle_obj(obj: Any, outfile: str = "temp_dir/temp.p") -> None:
    """pickle_obj.

    Helper fn to pickle object

    Args:
        obj (Any): obj
        outfile (str): outfile

    Returns:
        None
    """
    with open(outfile, "wb") as fp:
        pickle.dump(obj, fp)


def pickle_load(infile: str = "temp_dir/temp.p") -> Any:
    """pickle_load.

    Args:
        infile (str): infile, the name of input object

    Returns:
        Any: the object loaded from pickled file

    """
    with open(infile, "rb") as fp:
        return pickle.load(fp)


def make_dir(filename: str) -> None:
    """make_dir.

    Makes the directory that should contain this file

    Args:
        filename (str): filename

    Returns:
        None
    """
    # Make outdir if it doesn't exist
    out_folder = os.path.dirname(filename)
    os.makedirs(out_folder, exist_ok=True)
