# Colocalisation with gpu-coloc

This repository contains the code used to run **colocalisation** using `gpu-coloc` between the following datasets:  
- EstBB + UKBB meta-analysis  
- FinnGen (fine-mapped and ABFs, reused from [doi.org/10.1101/2025.08.25.672103](https://doi.org/10.1101/2025.08.25.672103))  
- Million Veterans Programme (MVP)  
- PANUKBB  
- FinnGen + UKBB + MVP meta-analysis  

---

## Requirements

- **PanUKBB manifest**  
  [Google Sheet](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit?gid=1450719288#gid=1450719288)

- **Liftover chain file** `GRCh37_to_GRCh38.chain.gz`  
  [UCSC Download](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/)

- **Reference genome** `hg38.fa`  
  [UCSC Download](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

---

## Python scripts

Each Python script has a paired `.sh` file for submitting `sbatch` jobs.

- `download_panukbb.py` — download all PanUKBB summary stats  
- `panukbb_liftover.py` — liftover PanUKBB from GRCh37 → GRCh38 and check ref/alt alleles  
- `panukbb_signals.py` — pre-format lifted PanUKBB for gpu-coloc  
- `FINNGEN+UKBB+MVP.py` — pre-format FinnGen + UKBB + MVP for gpu-coloc  
- `meta_signals.py` — pre-format EstBB + UKBB meta-analysis (by population)  
- `meta_meta_eur.py` — pre-format EstBB + UKBB meta-analysis (meta_EUR)  
- `MVP_signals.py` - pre-format MVP for gpu-coloc (checks ref/alt alleles and computes -log10p, otherwise there will be p-val = 0)
- `annotate_metabolite_summaries.py` — to annotate colocalisation results with the -log10p, maf and nearest genes of the metabolite

---

## Running gpu-coloc

Other `.sh` scripts are used to run `gpu-coloc` directly via `sbatch`.
