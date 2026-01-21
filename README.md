# Meta-analysis on data from the Estonian Biobank and the UK Biobank 

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
* [Merge meta_EUR colocalisation results across data sources](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/coloc/meta_EUR_merge_gpu-coloc.R)

### Genetic correlation
* [Calculate genetic correlation between studies or datasets](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/genetic_correlation/)

### Fine mapping
* [Calculate LD matrices and use them for fine mapping](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/fine_mapping/)

### Mendelian randomization
* [Run *cis*-MR analyses on loci of interest](https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/blob/main/code/MR/)
