# PLOS ONE data repository for nitrofurantoin and IPF analysis

This repository contains data files, analysis outputs, and scripts supporting the manuscript:

**Integrated network toxicology and single-cell transcriptomic analyses identify MMP13 as a putative link between nitrofurantoin and IPF-related aberrant epithelial remodeling**

## Repository contents

- `data/06ML/`: bulk transcriptomic machine-learning data package copied from the local `06ML` analysis directory.
- `DATA_MANIFEST.csv`: file-level manifest with relative paths, file sizes, and SHA-256 checksums.
- `CITATION.cff`: preferred citation metadata for this repository.
- `.zenodo.json`: metadata used by Zenodo when archiving a GitHub release.

The uploaded `06ML` package contains 74 files, including normalized expression matrices, model input files, trained R model objects, prediction outputs, ROC/confusion-matrix summaries, analysis scripts, and PDF figures. The total file size is approximately 118.79 MB.

## Public source data

The bulk transcriptomic datasets used in the study were obtained from the Gene Expression Omnibus (GEO):

- GSE10667
- GSE35145
- GSE24206
- GSE53845

The single-cell RNA-seq dataset used in the study was obtained from GEO:

- GSE128033

Additional public database inputs were obtained from:

- PubChem, for the nitrofurantoin chemical structure
- ChEMBL, SwissTargetPrediction, and SuperPred, for drug-related target retrieval
- UniProt, for gene/protein identifier standardization
- STRING, for protein-protein interaction information
- Protein Data Bank, for the MMP13 protein structure, accession ID 5B5P

## Data availability statement

All data underlying the findings of this study are publicly available without restriction. The bulk transcriptomic datasets used in this study were obtained from GEO under accession numbers GSE10667, GSE35145, GSE24206, and GSE53845, and the single-cell RNA-seq dataset was obtained from GEO under accession number GSE128033. Public database inputs were obtained from PubChem, ChEMBL, SwissTargetPrediction, SuperPred, UniProt, STRING, and the Protein Data Bank under accession ID 5B5P. The processed data, model inputs and outputs, analysis scripts, and supporting files from the bulk transcriptomic machine-learning analyses are available in this repository.

## Reuse

This repository is intended to support peer review, reproducibility, and reuse of the manuscript data package. Original third-party datasets and database records should also be cited according to the requirements of their source repositories.

