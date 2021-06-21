# enzyme-datasets
Enzyme datasets used to benchmark enzyme-substrate promiscuity models

## Dataset details

| Dataset Class    | Dataset Name                         | Dataset Type   |   # Enz. |   # Sub. |   # Pairs |
|:-----------------|:-------------------------------------|:---------------|---------:|---------:|----------:|
| aminotransferase | aminotransferase.csv                 | Regression     |       25 |       18 |       450 |
| aminotransferase | aminotransferase_binary.csv          | Binary         |       25 |       18 |       450 |
| olea             | olea.csv                             | Regression     |       73 |       14 |       550 |
| gt               | gt_donors_achiral_binary.csv         | Binary         |       55 |       10 |       514 |
| gt               | gt_donors_achiral_categorical.csv    | Categorical    |       55 |       10 |       514 |
| gt               | gt_donors_chiral_categorical.csv     | Categorical    |       55 |       13 |       667 |
| gt               | gt_donors_chiral_binary.csv          | Binary         |       55 |       13 |       667 |
| nitrilase        | nitrilase_binary.csv                 | Binary         |       18 |       38 |       684 |
| olea             | olea_binary.csv                      | Binary         |       73 |       15 |      1095 |
| halogenase       | halogenase_NaBr_binary.csv           | Binary         |       42 |       62 |      2604 |
| halogenase       | halogenase_NaBr.csv                  | Regression     |       42 |       62 |      2604 |
| halogenase       | halogenase_NaCl.csv                  | Regression     |       42 |       62 |      2604 |
| halogenase       | halogenase_NaCl_binary.csv           | Binary         |       42 |       62 |      2604 |
| duf              | duf_binary.csv                       | Binary         |      161 |       17 |      2737 |
| gt               | gt_acceptors_achiral_categorical.csv | Categorical    |       54 |       90 |      4298 |
| gt               | gt_acceptors_achiral_binary.csv      | Binary         |       54 |       90 |      4298 |
| gt               | gt_acceptors_chiral_binary.csv       | Binary         |       54 |       91 |      4347 |
| gt               | gt_acceptors_chiral_categorical.csv  | Categorical    |       54 |       91 |      4347 |
| esterase         | esterase_binary.csv                  | Binary         |      146 |       96 |     14016 |
| davis            | davis_filtered.csv                   | Regression     |      318 |       72 |     22896 |
| davis            | davis.csv                            | Regression     |      405 |       72 |     29160 |
| phosphatase      | phosphatase_achiral.csv              | Regression     |      218 |      108 |     23544 |
| phosphatase      | phosphatase_achiral_binary.csv       | Binary         |      218 |      108 |     23544 |
| phosphatase      | phosphatase_chiral.csv               | Regression     |      218 |      165 |     35970 |
| phosphatase      | phosphatase_chiral_binary.csv        | Binary         |      218 |      165 |     35970 |




## Install

Creating an env: 

`conda create -c conda-forge -n enz-datasets rdkit python=3.6`

other packages: xlrd, scipy, tqdm, openpyxl, cirpy, Biopython, requests, tabulate


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
