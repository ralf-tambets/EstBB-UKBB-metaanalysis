# Genome-wide association study for circulating metabolic traits in 619,372 individuals

This repository contains code from the meta-analysis of 249 metabolic traits across the Estonian Biobank and UK Biobank ([preprint](https://doi.org/10.1101/2024.10.15.24315557))

## Browse results
* [meta_EUR PheWeb browser](https://nmrmeta.gi.ut.ee/)
* [meta_EUR colocalisation browser](https://elixir.ut.ee/eqtl/nmr-coloc)

## Download data
* [Full summary statistics](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/data/sumstats_paths.tsv)
* [Lead variants](https://dx.doi.org/10.5281/zenodo.13937265)
* [Colocalisation results](https://zenodo.org/records/17945143)
* [Fine-mapping results](https://dx.doi.org/10.5281/zenodo.18132538)


## Code overview
### Genetic colocalisation
* [gpu-coloc software](https://github.com/mjesse-github/gpu-coloc)
* [Prepare input data for gpu-coloc](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/tree/main/code/coloc/gpu-coloc)
* [Merge meta_EUR colocalisation results across data sources](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/coloc/meta_EUR_merge_gpu-coloc.R)

### Genetic correlation
* [Calculate genetic correlation between studies or datasets](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/genetic_correlation/)

### Fine mapping
* [Calculate LD matrices and use them for fine mapping](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/fine_mapping/)
* [Analyse fine-mapped credible sets](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/explore_finemapped_cs.r)

### Mendelian randomization
* [Run *cis*-MR analyses on loci of interest](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/MR/)
