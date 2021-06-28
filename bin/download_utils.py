""" Download data from internet. 

Helper utils to pull down biology data

"""
import os
import numpy as np
import logging
from typing import Tuple, Optional, List
import tempfile
import copy

import re
import requests
import shutil
import urllib.request as request
from contextlib import closing
from lxml import etree as ET
import time
from collections import defaultdict

from tqdm import tqdm

MAX_BATCH_SIZE = 199


# Download from uniprot
def num_lines(f_name: str) -> int:
    """Count num lines in a file"""

    count = 0
    for line in open(f_name):
        count += 1
    return count


# Download from uniprot
def has_lines(f_name: str) -> bool:
    """Return true if there's more than 1 line in file"""

    count = 0
    ret = False
    for line in open(f_name):
        ret = True
        break

    return ret


def download_enzymes(
    f_name: str,
    step_size: int = 1000,
    total_lines_to_read: Optional[int] = None,
    read_timeout: int = 240,
    columns: List[str] = ["id", "sequence"],
    query: str = "ec:*",
):
    """download_enzymes.

    Args:
        f_name (str) : File name
        step_size (int) : How many lines to raed in each request
        total_lines_to_read (Optionial[int]) : Number of lines to read
        read_timeout (int) : Seceonds to wait before timing out
        columns (List[str]): List of columns to download
        query (str): String for query to search. To pull down all seqs in
            uniprot labeled with an EC num, use "ec:*" to pull down everything
            with a cross referencee to BRENDA, use "database:(type:brenda)"
    """

    # https://www.uniprot.org/uniprot/?query=ec%3A*&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2Clength%2Cec%2Cfeature(DOMAIN%20EXTENT)%2Cdatabase(PDB)%2Csequence%2Ccomment(CATALYTIC%20ACTIVITY)%2Cfeature(ACTIVE%20SITE)&sort=score
    base_url = "https://www.uniprot.org/uniprot/"

    # See if we've already read part of the file
    num_lines_read = 0
    write_code = "w"
    if os.path.exists(f_name):
        # Subtract 1 for header
        num_lines_read = max(num_lines(f_name) - 1, 0)
        write_code = "a"
    params = {
        "query": query,  # "database:(type:brenda) OR ec:*",
        "columns": ",".join(columns),
        "format": "tab",
        "offset": num_lines_read,
        "limit": step_size,
    }

    try:
        with open(f_name, write_code) as f:

            # Continue to place requests until we read all the lines
            num_requests = 0
            pbar = None

            while total_lines_to_read is None or num_lines_read < total_lines_to_read:
                r = requests.get(base_url, params=params, timeout=read_timeout)
                text = r.text

                # First pass in this download
                if num_requests == 0:
                    header_lines = int(r.headers["X-Total-Results"])
                    if total_lines_to_read is None:
                        total_lines_to_read = header_lines
                    else:
                        total_lines_to_read = min(header_lines, total_lines_to_read)

                    # Make download bar
                    pbar = tqdm(total=total_lines_to_read)
                    # Set at the starting position
                    pbar.update(num_lines_read)

                # If starting position is 0, add header
                if num_lines_read == 0:
                    f.write(text)
                else:
                    # Remove the header (first line)
                    f.write(text.split("\n", 1)[1])

                # Increment num lines read and the offset
                num_lines_read += step_size
                params["offset"] += step_size
                pbar.update(step_size)
                num_requests += 1
            if pbar is not None:
                pbar.close()

    except (requests.exceptions.ReadTimeout, requests.exceptions.ChunkedEncodingError):
        # We got a requests.exceptions.ChunkedEncodingError last time???
        if pbar is not None:
            pbar.close()

        logging.warning(f"Timeout Exception! Starting download_enzymes again.")

        # Recursive call
        download_enzymes(
            f_name, step_size, total_lines_to_read, read_timeout * 2, columns, query
        )


UNIPROT_COLUMNS = [
    "id",
    "organism",
    "ec",
    "sequence",
    "feature(ACTIVE SITE)",
    "organism-id",
    "length",
    "annotation score",
    "protein names",
]


def download_enzyme_list(
    f_name: str,
    enzyme_list: list,
    step_size: int = 1000,
    read_timeout: int = 240,
    columns: List[str] = ["id", "sequence"],
    database: str = "uniprot",
):
    """download_enzyme_list.

    Args:
        f_name (str) : File name
        enzyme_list (list) : List of uniprot ids to query
            NOTE: If we want to search speciifically ids or uniprot ids ini
            uniparc, this liist should contain that full string. E.g.
            ["id:uniprotaccession",..] or ["uniprot":accession,...,]
        step_size (int) : How many lines to raed in each request
        total_lines_to_read (Optionial[int]) : Number of lines to read
        read_timeout (int) : Seceonds to wait before timing out
        columns (List[str]): List of columns to download
        database (str): Name of database to be querying
    """

    if len(enzyme_list) == 0:
        logging.warning("No enzymes given, returning from download_enzyme_list")
        return
    base_url = f"https://www.uniprot.org/{database}/"

    # See if we've already read part of the file
    has_header = False
    write_code = "w"
    if os.path.exists(f_name):
        # Subtract 1 for header
        has_header = has_lines(f_name)
        write_code = "a"

    num_partitions = (len(enzyme_list) // step_size) + 1
    splits = np.array_split(enzyme_list, num_partitions)

    # Calculate which column corresponds to an alternate name
    alt_name_index = None
    for index, j in enumerate(columns):
        if j == "protein names":
            alt_name_index = index

    with open(f_name, write_code) as f:
        for split_index, split in tqdm(enumerate(splits)):
            params = {
                "query": " OR ".join(split),
                "columns": ",".join(columns),
                "format": "tab",
            }
            try:
                # Continue to place requests until we read all the lines

                r = requests.get(base_url, params=params, timeout=read_timeout)
                text = r.text
                if text == "":
                    continue

                # If starting position is 0, add header
                if not has_header:
                    f.write(text)
                else:
                    # Remove the header (first line)
                    f.write(text.split("\n", 1)[1])

                    # Increment num lines read and the offset
                    # num_lines_read += len(split)

            except (
                requests.exceptions.ReadTimeout,
                requests.exceptions.ChunkedEncodingError,
            ):
                # We got a requests.exceptions.ChunkedEncodingError last time???

                logging.warning(f"Timeout Exception! Starting again.")

                # Recursive call
                download_enzyme_list(
                    f_name,
                    list(np.concatenate(splits[split_index:])),
                    step_size,
                    read_timeout * 2,
                    columns,
                )


def replace_uniparc_ids(
    f_name: str,
    step_size: int = 1000,
    read_timeout: int = 240,
    columns: List[str] = ["id", "sequence"],
):
    """replace_uniparc_ids.

    Place a call into all uniprot to pick up all entries that have been merged

    Args:
        f_name (str) : File name that has columns already.
        step_size (int) : How many lines to raed in each request
        total_lines_to_read (Optionial[int]) : Number of lines to read
        read_timeout (int) : Seceonds to wait before timing out
        columns (List[str]): List of columns to download
    """

    required_mapping = {}
    deleted_entries = set()
    with open(f_name, "r") as fp:
        alt_name_index = None
        id_index = None
        seq_index = None
        for line, j in enumerate(open(f_name, "r")):
            # Parse header
            if line == 0:
                headers = j.strip().split("\t")
                if len(headers) != len(columns):
                    raise RuntimeError(
                        """Require the file being opened to
                                       have the same # of headers as the
                                       columns passed as argument"""
                    )
                for header_index, header in enumerate(headers):
                    if header.strip() == "Protein names":
                        alt_name_index = header_index
                    if header.strip() == "Entry":
                        id_index = header_index
                    if header.strip() == "Sequence":
                        seq_index = header_index
                if alt_name_index is None or id_index is None or seq_index is None:
                    raise RuntimeError(
                        f"File {f_name} has no col Protein names or Entry or sequence "
                    )
            # Parse body to get sequences that have been merged
            else:
                cols = j.strip("\n").split("\t")
                alt_name = cols[alt_name_index].strip()
                uniprot_id = cols[id_index].strip()
                merge_regex = "Merged into (.*?)\."
                # Only capture the first demerged enzyme word
                demerge_regex = "Demerged into (\w+).*$"
                if "Merged" in alt_name:
                    r = re.search(merge_regex, alt_name)
                    if r is None or len(r.groups()) != 1:
                        raise RuntimeError(
                            f"Found merged but failed to match single regex to {alt_name}"
                        )
                    else:
                        merged_id = r.groups()[0]
                        required_mapping[uniprot_id] = merged_id
                elif "Deleted." in alt_name:
                    deleted_entries.add(uniprot_id)
                elif "Demerged" in alt_name:
                    r = re.search(demerge_regex, alt_name)
                    if r is None or len(r.groups()) != 1:
                        raise RuntimeError(
                            f"Found merged but failed to match single regex to {alt_name}"
                        )
                    else:
                        merged_id = r.groups()[0]
                        # Drop these into the same required mapping
                        required_mapping[uniprot_id] = merged_id

    # Query Uniparc for all "Deleted" ids and
    # turn them into mappings from old id
    # to needed id
    uniparc_cols = ["id", "kb", "sequence"]

    # Get all uniparc ids
    # Load all remapped into memory directly
    deleted_ids = list(deleted_entries)
    temp_uniparc = tempfile.NamedTemporaryFile(delete=False)
    temp_uniparc_name = temp_uniparc.name
    # Download enzymes from this mapping into tf
    download_enzyme_list(
        temp_uniparc_name, deleted_ids, step_size, read_timeout, uniparc_cols, "uniparc"
    )

    # Map old id to sequence if we don't find a currently valid annotation for
    # this organism
    uniparc_mapping = {}
    num_uniparc_added = 0

    # Process temp uniparc file to extract new downloaded sequences
    for row_index, entry in enumerate(temp_uniparc):
        entry = entry.decode()
        if row_index == 0:
            pass
        else:
            split_entry = entry.strip("\n").split("\t")
            uniprot_ids = split_entry[1].strip().split(";")
            sequence = split_entry[2].strip()

            valid_id = None
            found_deleted = set()
            # Loop over every item stored in this column entry. Find a valid
            # mapping in the current uniprot database and find the deleted
            # entry
            for uniprot_id in uniprot_ids:
                if uniprot_id == "":
                    break
                uniprot_id = uniprot_id.strip()
                # Just extract top match that is a valid id
                if "obsolete" not in uniprot_id and valid_id is None:
                    valid_id = uniprot_id
                else:
                    # Split them at the dot to get the old ids
                    temp_id = uniprot_id.split(".")[0]
                    if temp_id in deleted_ids:
                        found_deleted.add(temp_id)

            for deleted_id in found_deleted:
                if valid_id is not None:
                    required_mapping[deleted_id] = valid_id
                    num_uniparc_added += 1
                elif sequence != "":
                    # If we can't get a sequence
                    uniparc_mapping[deleted_id] = sequence
                    num_uniparc_added += 1
                else:
                    raise RuntimeError("Unexpected loop branch")

    logging.info(
        f"""Looking through the uniprot downloads, we found:
                 {num_uniparc_added} sequences to query in uniparc"""
    )
    temp_uniparc.close()
    os.remove(temp_uniparc_name)

    # Get all MERGED id's
    # Load all remapped into memory directly
    new_values = list(required_mapping.values())
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf_name = tf.name
    # Download enzymes from this mapping into tempfiile, tf
    download_enzyme_list(tf_name, new_values, step_size, read_timeout, columns)

    id_to_new_row = {}
    for line_index, j in enumerate(open(tf_name, "r")):
        if line_index > 0:
            cols = j.strip("\n").split("\t")
            new_id = cols[id_index].strip()
            # Map new id to the old column list
            id_to_new_row[new_id] = cols
    tf.close()
    # Remove the tempfile
    os.remove(tf_name)

    # Now rewrite the old file into a new file
    # Copy the old file into a temporary file
    tf = tempfile.NamedTemporaryFile(delete=False)
    with open(f_name, "r+b") as f:
        shutil.copyfileobj(f, tf)

    # Rewrite old file
    # Go to top of tf
    tf.seek(0)
    # Now read from tempfile and rewrite into f_name replacing the info
    # the alt name rows with the "merged into" informatoin
    with open(f_name, "w") as fp:
        for line_index, line in enumerate(tf):
            # Convert to str
            line = line.decode()

            # Copy header
            if line_index == 0:
                fp.write(line)
            else:
                cols = line.strip().split("\t")
                uniprot_id = cols[id_index].strip()

                if uniprot_id in required_mapping:

                    # If alt id, insert the alt id into the right position
                    new_line = copy.copy(id_to_new_row[required_mapping[uniprot_id]])

                    # Replace the mapped to id with the correct uniprot id
                    new_line[id_index] = uniprot_id
                    new_line = "\t".join(new_line)

                    # Add next line token
                    new_line = f"{new_line}\n"
                    fp.write(new_line)
                elif uniprot_id in uniparc_mapping:
                    # If we couldn't find a currently dated uniprot entry,
                    # check uniparc for the sequence alone

                    # Replace the sequence entry to the correct one
                    cols[seq_index] = uniparc_mapping[uniprot_id]
                    new_line = "\t".join(cols)

                    # Add next line token
                    new_line = f"{new_line}\n"
                    fp.write(new_line)

                # If this doesn't have an alt id, just write it
                else:
                    fp.write(line)

    tf.close()
    os.remove(tf.name)


def download_uniprot_uniparc_list(
    f_name: str, enzyme_list: list, step_size: int = 200, read_timeout: int = 240
):
    """download_uniprot_uniparc_list.

    Consolidate the methods for pulling a list of sequences such that both
    merged and deleted sequences are still pulled w/ this method.

    Args:
        f_name (str) : File name
        enzyme_list (list) : List of uniprot ids
        step_size (int) : How many lines to raed in each request
        total_lines_to_read (Optionial[int]) : Number of lines to read
        read_timeout (int) : Seceonds to wait before timing out
    """
    download_enzyme_list(
        f_name,
        [f"id:{i}" for i in enzyme_list],
        database="uniprot",
        columns=UNIPROT_COLUMNS,
        step_size=step_size,
        read_timeout=read_timeout,
    )
    replace_uniparc_ids(
        f_name, step_size=step_size, read_timeout=read_timeout, columns=UNIPROT_COLUMNS
    )


######


def download_ftp(ftp_link: str, outfile: str):
    """download_ftp.

    Args:
        ftp_link (str): ftp_link
        outfile (str): outfile
    """

    with closing(request.urlopen(ftp_link)) as r:
        with open(outfile, "wb") as f:
            shutil.copyfileobj(r, f)


def query_pubchem(
    ids: list,
    query_type: str = "inchi",
    save_file: str = "pubchem_save.txt",
    rm_save: bool = True,
    encoding: str = "utf8",
) -> dict:
    """query_pubchem.

    Args:
        ids (list):
        query_type (str): inchi, chebi, or synonym
        save_file (str):
        rm_save (bool): If true, delete the file saved afterward
        encoding (str): Encoding to send request, defaults to utf-8

    Return:
        dict mapping ids to smiles lists
    """
    # Add options for query_type

    # 60 seconds
    WAIT_TIME = 60
    DOCHEADER = (
        '<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "NCBI_PCTools.dtd">'
    )
    URL = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
    REQUIRED_HEADINGS = [
        "PCT-Data_input",
        "PCT-InputData",
        "PCT-InputData_query",
        "PCT-Query",
        "PCT-Query_type",
        "PCT-QueryType",
        "PCT-QueryType_id-exchange",
        "PCT-QueryIDExchange",
    ]

    QUERY_SUBTREE_NAMES = ["PCT-QueryIDExchange_input", "PCT-QueryUids"]
    OUTPUT_HEADERS = [
        "PCT-QueryIDExchange_operation-type",
        "PCT-QueryIDExchange_output-type",
        "PCT-QueryIDExchange_output-method",
        "PCT-QueryIDExchange_compression",
    ]
    OUTPUT_VALUES = ["same", "smiles", "file-pair", "none"]

    # Start building the query tree
    root = ET.Element("PCT-Data")
    cur_pos = root

    for new_heading in REQUIRED_HEADINGS:
        cur_pos = ET.SubElement(cur_pos, new_heading)

    # Navigate down to where we add the inchis
    query_subtree = cur_pos
    for query_subtree_name in QUERY_SUBTREE_NAMES:
        query_subtree = ET.SubElement(query_subtree, query_subtree_name)

    # Now add the things SPECIFIC to inchi
    if query_type == "inchi":
        query_root, query_name = "PCT-QueryUids_inchis", "PCT-QueryUids_inchis_E"
        query_subtree = ET.SubElement(query_subtree, query_root)
        for id_ in ids:
            new_id = ET.SubElement(query_subtree, query_name)
            # give this the id text
            try:
                new_id.text = id_
            except ValueError:
                logging.warning(f"Couldn't query {id_} due to bad encoding")

    elif query_type == "synonym":
        query_root, query_name = "PCT-QueryUids_synonyms", "PCT-QueryUids_synonyms_E"
        query_subtree = ET.SubElement(query_subtree, query_root)
        for id_ in ids:
            new_id = ET.SubElement(query_subtree, query_name)
            # give this the id text
            try:
                new_id.text = id_
            except ValueError:
                logging.warning(f"Couldn't query {id_} due to bad encoding")

    elif query_type == "chebi":
        for i in ["PCT-QueryUids_source-ids", "PCT-RegistryIDs"]:
            query_subtree = ET.SubElement(query_subtree, i)
        source_id_name = ET.SubElement(query_subtree, "PCT-RegistryIDs_source-name")
        source_id_name.text = "ChEBI"

        query_subtree = ET.SubElement(query_subtree, "PCT-RegistryIDs_source-ids")
        for id_ in ids:
            new_id = ET.SubElement(query_subtree, "PCT-RegistryIDs_source-ids_E")
            # give this the id text
            try:
                new_id.text = id_
            except ValueError:
                logging.warning(f"Couldn't query {id_} due to bad encoding")
    else:
        raise NotImplementedError()

    # Go back up to to current position holder
    # Add the output specification
    for output_header, output_value in zip(OUTPUT_HEADERS, OUTPUT_VALUES):
        output_xml = ET.SubElement(cur_pos, output_header)
        output_xml.set("value", output_value)

    out_xml = ET.tostring(
        root, encoding=encoding, method="xml", xml_declaration=True, doctype=DOCHEADER
    ).decode()

    # Post the request!
    resp = requests.post(URL, data=out_xml.encode("utf-8"))

    # Handle response and build a request to check on status
    resp_tree = ET.fromstring(resp.text)
    waiting_id = resp_tree.xpath("//PCT-Waiting_reqid")
    waiting_id = waiting_id[0].text if waiting_id else None

    STATUS_CHECK_HEADERS = [
        "PCT-Data_input",
        "PCT-InputData",
        "PCT-InputData_request",
        "PCT-Request",
    ]

    root = ET.Element("PCT-Data")
    cur_pos = root
    for header in STATUS_CHECK_HEADERS:
        cur_pos = ET.SubElement(cur_pos, header)
    req_id = ET.SubElement(cur_pos, "PCT-Request_reqid")
    req_id.text = waiting_id
    req_type = ET.SubElement(cur_pos, "PCT-Request_type")
    req_type.set("value", "status")
    query_xml = ET.tostring(
        root, encoding=encoding, method="xml", xml_declaration=True, doctype=DOCHEADER
    ).decode()

    download_link = None
    waiting_time = 0

    # TODO: Add stop timeout condition?
    # Repeatedly query to see if the results are done, then sleep for WAITIME
    # in case they aren't
    while not download_link:
        resp = requests.post(URL, data=query_xml.encode("utf-8"))
        resp_tree = ET.fromstring(resp.text)
        download_link = resp_tree.xpath("//PCT-Download-URL_url")
        download_link = download_link[0].text if download_link else None
        time.sleep(WAIT_TIME)
        waiting_time += WAIT_TIME
        logging.warning(f"Waiting time: {waiting_time} seconds")

    # At conclusion, download the ftp file
    download_ftp(download_link, save_file)

    # Also parse this
    ret_dict = defaultdict(lambda: [])
    with open(save_file, "r") as fp:
        for linenum, line in enumerate(fp):
            line = line.strip()
            split_line = line.split("\t")
            if len(split_line) == 2:
                mol, smiles = split_line
                mol = mol.strip()
                smiles = smiles.strip()
                ret_dict[mol].append(smiles)
            else:
                logging.debug(f"No smiles mol found for {line}")

    # Remove temp file
    if os.path.exists(save_file) and rm_save:
        os.remove(save_file)

    return ret_dict


def debug_pubchem():
    """Helper fn for debugging pubchem"""
    test_inchis = [
        "InChI=1S/I2/c1-2",
        "InChI=1S/C30H46O3/",
        "InChI=1S/C6H11NO2/c8-6(9)5-3-1-2-4-7-5/h5,7H,1-4H2,(H,8,9)/t5-/m0/s1",
    ]
    test_chebis = ["chebi:61185"]
    query_syns = ["glucose", "NADH"]

    mapping = query_pubchem(
        test_inchis, query_type="inchi", save_file="pubchem_inchi.txt"
    )
    print(mapping)
    mapping = query_pubchem(
        test_chebis, query_type="chebi", save_file="pubchem_chebi.txt"
    )
    mapping = query_pubchem(
        query_syns, query_type="synonym", save_file="pubchem_syns.txt"
    )


if __name__ == "__main__":
    debug_pubchem()
