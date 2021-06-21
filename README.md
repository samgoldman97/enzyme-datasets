# enzyme-datasets
Enzyme datasets used to benchmark enzyme-substrate promiscuity models

Enzyme datasets figured currently

## Dataset details

TODO:   
1. [Done] `reformat_esterase.py`   
2. [Done] `reformat_gts.py`    
3. [Done] `reformat_halogenase.py`    
4. [Done] `reformat_phosphatase.py`    
5. `reformat_olea.py`    
6. `reformat_duf.py`    
7. `reformat_davis.py`    
8. [OPTIONAL] `reformat_ainotransfrase.py`    
9. [OPTIONAL] `reformat_nitrilase.py`    
10. [OPTIONAL] `reformat_pafa.py`    
11. [OPTIONAL] `reformat_rubisco.py`    

After compiling this, add list of dataset processed files and the amount of enzymes, substrates, and pairings for each

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






