
#Import QTL coloc
lc_coloc = readr::read_tsv("data/sup_tables/leafcutter_vs_met_coloc_results.tsv")
ge_coloc = readr::read_tsv("data/sup_tables/ge_vs_met_coloc_results.tsv")

#Import GWAS
gwas_coloc = readr::read_tsv("data/met_vs_diasease_results.tsv") %>%
  dplyr::transmute(gwas_signal = signal1, gwas_lead = lead1, signal2, gwas_PP4 = PP.H4)

shared_eqtl = dplyr::left_join(gwas_coloc, ge_coloc, by = "signal2", relationship = "many-to-many")
shared_lc = dplyr::left_join(gwas_coloc, lc_coloc, by = "signal2", relationship = "many-to-many")
all_shared_coloc = dplyr::bind_rows(shared_eqtl, shared_lc) %>% dplyr::filter(!is.na(study_label))

#IL6R
dplyr::filter(all_shared_coloc, gene_name == "IL6R") %>% View()

#HMGCR
dplyr::filter(all_shared_coloc, gene_name == "HMGCR", signal2 == "LDL_C_chr5:74360714-76360714") %>% View()