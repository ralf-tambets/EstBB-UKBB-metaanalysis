library("data.table")
library("dplyr")

all_signals = readr::read_tsv("data/ge_vs_met_results.tsv")

#Microrarray signals
microarray = dplyr::filter(all_signals, signal1 %like% "ILMN")
other_signals = dplyr::filter(all_signals, !(signal1 %like% "ILMN"))

#Separate out dataset ids and molecular trait ids
sep_signals = tidyr::separate(other_signals, signal1, c("dataset_id", "molecular_trait_id", "cs_index"), "_")

#Split micrarray signals
split_signal = tidyr::separate(microarray, signal1, c("dataset_id", "rest"), sep = "_", extra = "merge") %>%
  tidyr::separate(rest, c("molecular_trait_id", "cs_index"), sep = "_L") %>%
  dplyr::mutate(cs_index = paste0("L", cs_index))

all_split_signals = dplyr::bind_rows(sep_signals, split_signal)

#Import dataset metadata
dataset_meta = readr::read_tsv("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/refs/heads/master/data_tables/dataset_metadata.tsv")

#Import gene meta
ge_meta = readr::read_tsv("https://zenodo.org/records/7808390/files/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz") %>%
  dplyr::transmute(molecular_trait_id = phenotype_id, gene_id, gene_name)

#Import microarray_meta
microarray_meta = readr::read_tsv("https://zenodo.org/records/7808390/files/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::transmute(molecular_trait_id = phenotype_id, gene_id, gene_name)

gene_meta = dplyr::bind_rows(ge_meta, microarray_meta) %>% dplyr::distinct()

joined_data = dplyr::left_join(all_split_signals, dataset_meta, by = "dataset_id") %>%
  dplyr::left_join(gene_meta, by = "molecular_trait_id")

write.table(joined_data, "data/sup_tables/ge_vs_met_coloc_results.tsv", sep = "\t", row.names = F, quote = F)

joined_data = readr::read_tsv("data/sup_tables/ge_vs_met_coloc_results.tsv")

#Filter LIPA colocalisations
lipa_coloc = dplyr::filter(joined_data, gene_name == "LIPA") %>% 
  dplyr::arrange(-PP.H4) %>% 
  dplyr::filter(signal2 == "GlycA_chr10:88243047-90243047")

#Filter SORT1 colocalisations
sort_coloc = dplyr::filter(joined_data, gene_name == "SORT1") %>% 
  dplyr::arrange(-PP.H4) %>% 
  dplyr::filter(signal2 == "LDL_C_chr1:108274570-110274570")

hmgcr_coloc = dplyr::filter(joined_data, gene_name == "HMGCR") %>% 
  dplyr::arrange(-PP.H4) 

