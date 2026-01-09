# Fine mapping and LD matrix calculation

This repository contains the code used to run **LD matrix calculation** using [`LDstore`](http://www.christianbenner.com/#) and **fine mapping** analyses using [`susieR`](https://cran.r-project.org/web/packages/susieR/index.html).

---

## Scripts
- `01_loci_select.R` is used to select regions in which to calculate LD matrices.
- `02_ldstore_prep.R` is used to prepare input files for LDstore.
- `03_submit_to_dnanexus.txt` is the command used to submit LDstore jobs to DNAnexus.
- `04_finemap_locus.R` is the R script used to fine map the regions. It relies on [code](https://github.com/AlasooLab/reGSusie/blob/main/bin/susieR.R) written by Ida Rahu.

