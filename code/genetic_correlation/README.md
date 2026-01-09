# Genetic correlation calculation using LDSC

This repository contains the `nextflow` code used to run **genetic correlation** analyses using `ldsc` between studies or between metabolites from different studies. Relies heavily on [this tutorial](https://cloufield.github.io/GWASTutorial/08_LDSC/#cross-trait-ld-score-regression)

---
## Requirements
- input file `studies.tsv` contains the name of the study and the location of the study folder.
- input file `metabolites.csv` contains a single row of comma-separated metabolite names.
- GWAS analyses are files with columns `SNP`, `A1`, `A2`, `FRQ`, `INFO`, `N`, `BETA`, `SE`, `P`. 
- the GWAS results files have the names of the metabolites.
---

## Nextflow script

`main.nf` is started by running submit_nextflow.sh. User should specify `-entry between_metabolites` or `-entry between_studies` as needed.
