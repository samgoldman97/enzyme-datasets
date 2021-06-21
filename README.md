# enzyme-datasets
Enzyme datasets used to benchmark enzyme-substrate promiscuity models

Enzyme datasets figured currently

## Dataset details

TODO:   
1. [Done] `reformat_esterase.py`   
2. [Done] `reformat_gts.py`    
3. [Done] `reformat_halogenase.py`    
4. [Done] `reformat_phosphatase.py`    
5. [Done] `reformat_olea.py`    
6. [Done] `reformat_duf.py`    
7. [Done] `reformat_davis.py`    

Not benchmarked 

8. [Done] `reformat_ainotransfrase.py`    
9. [Done] `reformat_nitrilase.py`    

After compiling this, add list of dataset processed files and the amount of enzymes, substrates, and pairings for each

Auto generate MSA? 

Also auto generate structure references? 

Make distinction between nitrilase and aminotransferase datasets

## Install

Creating an env: 

`conda create -c conda-forge -n enz-datasets rdkit python=3.6`

other packages: xlrd, scipy, tqdm, openpyxl, cirpy, Biopython, requests


## Dataset descriptions

### Esterase

**Source:** Martinez-Martinez et al. *ACS Chem. Biol*. 2018.  

Raw data can be found in the SI and includes listings of smiles strings and enzyme sequences

### Glycosyltransferases (gts)

**Source:** Yang et al. *Nature Chem. Bio.* 2017. 

Raw data is extracted as an excel spreadsheet with sequences. Numbers are manually added to each spreadsheet to reflect the low, medium, high activity scoring color equivalent to the green, amber, red screen. Some combinations were not tested. Results from both the acceptor and donor screen are extracted and parsed into categorical and binary data. In the binary setting, medium activity is considered to be "active." 

### Halogenases 

**Source:** Lewis et al. *ACS Cent. Sci.* 2019. 

Conversion, chemdraw files, the ssn file with sequences tested, and solubility files are extracted from the SI of the halogenase paper. Data was processed into regression and binarized prediction tasks using a cutoff of 0.08 for binary thresholding. Sequences were also cutoff at 1,000 amino acids in length, removing a single sequence. Lastly, to remove sequences that have never been seen to perform halogenatoin reactions, all sequences that never achieve maximum conversion of 0.05 were filtered. 

Separate data files were created for bromination and chlorination reactions.

### Phosphatases

**Source:** Huang et al. *PNAS.* 2015. 

Raw phosphatase data is extracted from the SI and stored as two excel notebooks. Each sheet corresponds to a different protein screened. Each cell's corresponding compound is listed as a comment on the cell. These are extracted programmatically. A corresponding manually annotated smiles data file was created by redrawing chemical structures from the SI of the original paper, in addition to a series of other name to smiles mappings (pubchem and cirpy). The proper smiles is chosen from this dat file in order of manual mapping, pubchem mapping, and cirpy mapping. 

Smiles are standardized by uncharging followed by uncharging. An achiral version of each molecule is also created for a second, achiral version of the dataset. 

Sequences are resolved from uniprot ID's by programmatically querying Uniprot. In the case where uniprot ID's are not found due to merged metagenomic entries, uniparc is further queried to find the original sequences tested. 

All activity is binarized at the suggested 0.2 threshold and also exported with the original value for regression tasks.

### OleA (thiolase) 

**Source:** Robinson et al. *Synthetic Biology*. 2020. 

Data input files are extracted from the corresponding github repository provided by Robinson et al. 

Binarized and regression versions of the dataset are created as in the original analaysis. 

### DUF (BKACE) 

**Source:** Bastard et al. *Nature Chem. Bio*. 2014. 

Data files are taken from Bastard et al in the form of BinaryActivities.DB. Smiles strings are manually redrawn from the apper and converted into a mapping between chemical name and smiles string. Sequences are further extracted from an MSA. 

### Davis (kinase inhibitors) 

**Source:** Davis et al. *Nature Biotech.* 2011. 

Unlike other enzyme datasets, this dataset contains a kinase inhibitor profile of 72 inhibitors against 442 kinases. We process this dataset in two ways. First, using the gene symbols that have deletions and insertions, we modify each sequence to have the appropriately listed insertions and deletions. Certain entries are given for both domains or a single domain of the kinase tested. To make sure we only model the kinase domain actually tested and stay within a single family, we use HMMER to trim each sequence to its corresponding domain PF000069. Further details can be found in the parse file `bin/reformat_davis.py`. 

All kinase Kd pairs not listed are assumed to have been tested and ascribed a Kd value of 10,000. We create a second version of the dataset, "filtered," that contain no substitutions, insertions, or deletions.

### Nitrilase 

**Source:** Black et al. *RSC Chem. Commun..* 2014. 

Nitrilase data is extracted directly from Black et al. and uniprot ids are converted to corresponding sequences.

### Aminotransferase

**Source:** Li et al. *ACS Catal.* 2020. 

Aminotransferase data is extracted directly from Li et al. and uniprot ids are converted to corresponding sequences.
