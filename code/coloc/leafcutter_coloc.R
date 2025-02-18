library(data.table)
library("dplyr")
library("stringr")
library("purrr")

all_signals = readr::read_tsv("data/lc_vs_met_results.tsv")

#Import leafcutter meta
file_paths = list.files("data/leafcutter_meta/", full.names = T)

dataset_labels = list.files("data/leafcutter_meta/") %>%
stringr::str_remove("leafcutter_") %>%
  stringr::str_remove("_Ensembl_105_phenotype_metadata.tsv.gz")

file_list = setNames(as.list(file_paths), dataset_labels)
lc_meta = purrr::map_df(file_list, readr::read_tsv, .id = "dataset_id") %>%
  dplyr::transmute(dataset_id, molecular_trait_id = phenotype_id, gene_id, gene_name)

split_signal = tidyr::separate(all_signals, signal1, c("dataset_id", "rest"), sep = "_", extra = "merge") %>%
  tidyr::separate(rest, c("molecular_trait_id", "cs_index"), sep = "_L") %>%
  dplyr::mutate(cs_index = paste0("L", cs_index))

#Import dataset metadata
dataset_meta = readr::read_tsv("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/refs/heads/master/data_tables/dataset_metadata.tsv")

#Join data together
joined_data = dplyr::left_join(split_signal, lc_meta, by = c("dataset_id", "molecular_trait_id")) %>%
  dplyr::left_join(dataset_meta, by = "dataset_id") 

write.table(joined_data, "data/sup_tables/leafcutter_vs_met_coloc_results.tsv", sep = "\t", row.names = F, quote = F)

joined_data = readr::read_tsv("data/sup_tables/leafcutter_vs_met_coloc_results.tsv")

hmgcr_coloc = dplyr::filter(joined_data, gene_name == "HMGCR") %>% 
  dplyr::filter(signal2 == "LDL_C_chr5:74360714-76360714") %>%
  dplyr::arrange(-PP.H4) 




