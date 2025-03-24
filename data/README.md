# Data sources

## PanRes v1.0.0
First metadata table for PanRes v.1.0.0 retrieved from: https://zenodo.org/records/8055116

## ResFinder
Phenotypes for ResFinder genes retrieved from: https://bitbucket.org/genomicepidemiology/resfinder_db

## ResFinderFG
The acronyms for the ResFinderFG genes are translated both with a manual lookup of the pubmed IDs, as well as using the file retrieved like this:
```
wget https://bitbucket.org/genomicepidemiology/resfinderfg_db/raw/cec72dd864faa11ae1301354cecd6cb71880c593/additional_info.txt
mv additional_info.txt resfinderfg_anno.txt
```

## CARD
Annotations for CARD genes retrieved from: https://card.mcmaster.ca/latest/data

## MegaRes
Information about the genes is extracted from the FASTA headers.

## AMRFinderPlus
The gene annotations for AMRFinderPlus retrieved from: https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt

## ARG-ANNOT
Information about the genes is extracted from the FASTA headers.

## BacMet
Annoations retrieved from http://bacmet.biomedicine.gu.se/download_temporary.html