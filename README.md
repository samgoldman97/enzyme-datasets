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

### Citations

If you use any of these datasets, please cite the following, respective datasets: 

1.  Bastard, K. et al. Revealing the hidden functional diversity of an enzyme family. Nature Chemical Biology 10, 42–49 (2014).    

2.  Black, G. W. et al. A high-throughput screening method for determining the substrate scope of nitrilases. Chemical Communications 51, 2660–2662 (2015).   

3.  Davis, M. I. et al. Comprehensive analysis of kinase inhibitor selectivity. Nature biotechnology 29, 1046–1051 (2011).    

4.  Hie, B., Bryson, B. D. & Berger, B. Leveraging uncertainty in machine learning accelerates biological discovery and design. Cell Systems 11, 461-477. e9 (2020).    

5.  Huang, H. et al. Panoramic view of a superfamily of phosphatases through substrate profiling. PNAS 112, E1974–E1983 (2015).    

6.  Li, T. et al. Exploration of transaminase diversity for the oxidative conversion of natural amino acids into 2-ketoacids and high-value chemicals. ACS Catalysis 10, 7950–7957 (2020).    

7.  Martínez-Martínez, M. et al. Determinants and Prediction of Esterase Substrate Promiscuity Patterns. ACS Chem. Biol. 13, 225–234 (2018).    

8.  Robinson, S. L., Smith, M. D., Richman, J. E., Aukema, K. G. & Wackett, L. P. Machine learning-based prediction of activity and substrate specificity for OleA enzymes in the thiolase superfamily. Synth Biol 5, (2020).    

9.  Yang, M. et al. Functional and informatics analysis enables glycosyltransferase activity prediction. Nature Chemical Biology 14, 1109–1117 (2018).    



## Install

Creating an env: 

`conda create -c conda-forge -n enz-datasets rdkit python=3.6`

other packages: xlrd, scipy, tqdm, openpyxl, cirpy, Biopython, requests, tabulate


## Dataset descriptions

### Esterase

**Source:** Martinez-Martinez et al. *ACS Chem. Biol*. 2018.  
**Parser file:** `bin/reformat_esterase.py`  


Raw data is extracted from the paper supplement.

### Glycosyltransferases (gts)

**Source:** Yang et al. *Nature Chem. Bio.* 2017.   
**Parser file:** `bin/reformat_gts.py`  


Raw data is extracted as an excel spreadsheet with sequences. Numbers are manually added to each spreadsheet to reflect the low, medium, high activity scoring color equivalent to the green, amber, red screen. Some combinations were not tested and labeled as 0. Results from both the acceptor and donor screen are extracted and parsed into categorical and binary data. In the binary setting, medium activity is considered to be "active." 

### Halogenases 

**Source:** Lewis et al. *ACS Cent. Sci.* 2019.    
**Parser file:** `bin/reformat_halogenase.py`  


Conversion, chemdraw files, the sequence similarity network file with sequences, and solubility files are extracted from the paper supplement. Data was processed into regression and binarized prediction tasks using a cutoff of 0.08 for binary thresholding. Sequences were also cutoff at 1,000 amino acids in length, removing a single sequence. Lastly, to remove sequences that may not be halogenases, all sequences that never achieve a minimum of 0.05 conversion were filtered. 

Separate data files were created for bromination and chlorination reactions, measured separatley by Lewis et al.

### Phosphatases

**Source:** Huang et al. *PNAS.* 2015.    
**Parser file:** `bin/reformat_phosphatase.py`  


Raw phosphatase data is extracted from paper supplement and stored as two excel notebooks. Each sheet corresponds to a different protein screened. Each cell's corresponding compound is listed as a comment on the cell. These are extracted programmatically. A corresponding, manually annotated smiles data file was created by redrawing chemical structures from the original paper supplement, in addition to a series of other name to smiles mappings (Pubchem and Cirpy). The proper smiles is chosen from this data file in order of priority (1. manual mapping 2.pubchem mapping, and 3. Cirpy mapping). 

Smiles are standardized by uncharging. An achiral version of each molecule is also created for a second, achiral version of the dataset. 

Sequences are resolved from uniprot ID's by programmatically querying Uniprot. In the case where Uniprot ID's are not found due to merged metagenomic entries, Uniparc is further queried to find the original sequences tested. 

All activity is binarized at the suggested 0.2 threshold and also exported with the original value for regression tasks.

### OleA (thiolase) 

**Source:** Robinson et al. *Synthetic Biology*. 2020.    
**Parser file:** `bin/reformat_olea.py`  


Data input files are extracted from the corresponding github repository provided by Robinson et al. 

Binarized and regression versions of the dataset are created as in the original analysis. 

### DUF (BKACE) 

**Source:** Bastard et al. *Nature Chem. Bio*. 2014.     
**Parser file:** `bin/reformat_duf.py`  


Data files are taken from Bastard et al in the form of BinaryActivities.DB. Smiles strings are manually redrawn from the paper and converted into a mapping between chemical name and smiles string. Sequences are further extracted from an MSA. 

### Davis (kinase inhibitors) 

**Source:** Davis et al. *Nature Biotech.* 2011.     
**Parser file:** `bin/reformat_davis.py`  


Unlike other enzyme datasets, this dataset contains a kinase inhibitor profile of 72 inhibitors against 442 kinases. We process this dataset in two ways. First, using the gene symbols that have deletions and insertions, we modify each sequence to have the appropriately listed insertions and deletions. Certain entries are given for both domains or a single domain of the kinase tested. To make sure we only model the kinase domain actually tested and stay within a single family, we use HMMER to trim each sequence to its corresponding domain PF000069. Further details can be found in the parse file `bin/reformat_davis.py`. 

All kinase Kd pairs not listed are assumed to have been tested and ascribed a Kd value of 10,000. We create a second version of the dataset, "filtered," that contains only sequences without substitutions, insertions, or deletions.

### Nitrilase 

**Source:** Black et al. *RSC Chem. Commun..* 2014.    
**Parser file:** `bin/reformat_nitrilase.py`  


Nitrilase data is extracted directly from Black et al. and uniprot ids are converted to corresponding sequences.

### Aminotransferase

**Source:** Li et al. *ACS Catal.* 2020.     
**Parser file:** `bin/reformat_aminotransferases.py`  


Aminotransferase data is extracted directly from Li et al. and uniprot ids are converted to corresponding sequences.


## Structure and MSA Extraction


|    | Dataset          | PDB ID   |
|---:|:-----------------|:---------|
|  0 | esterase         | 5a6v     |
|  1 | davis            | 2CN5     |
|  2 | aminotransferase | 3QPG     |
|  3 | nitrilase        | 3WUY     |
|  4 | phosphatase      | 3l8e     |
|  5 | halogenase       | 2AR8     |
|  6 | olea             | 4KU5     |
|  7 | duf              | 2Y7F     |
|  8 | gt               | 3HBF     |

For certain pipelines, we may be interested in having access to one protein crystal structure representative of the dataset tested or an alignment of allt he proteins in the dataset. We extract these crystal structure references as well as alignments using the file: `create_ref_and_aligns.sh.`  The table above includes all the PDB ID's used as reference structures.  

Note that MuscleCommandline must be installed to do this. 