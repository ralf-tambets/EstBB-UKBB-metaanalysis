#Import QTL coloc
lc_coloc = readr::read_tsv("data/sup_tables/leafcutter_vs_met_coloc_results.tsv")
ge_coloc = readr::read_tsv("data/sup_tables/ge_vs_met_coloc_results.tsv")

#Merge colocs
all_coloc = dplyr::bind_rows(ge_coloc, lc_coloc)

#Import GWAS
gwas_coloc = readr::read_tsv("data/met_vs_diasease_results.tsv") %>%
  dplyr::transmute(gwas_signal = signal1, gwas_lead = lead1, signal2, gwas_PP4 = PP.H4)

shared_eqtl = dplyr::left_join(gwas_coloc, ge_coloc, by = "signal2", relationship = "many-to-many")
shared_lc = dplyr::left_join(gwas_coloc, lc_coloc, by = "signal2", relationship = "many-to-many")
all_shared_coloc = dplyr::bind_rows(shared_eqtl, shared_lc) %>% dplyr::filter(!is.na(study_label))

#IL6R
dplyr::filter(all_shared_coloc, gene_name == "IL6R") %>% View()

#HMGCR
dplyr::filter(all_shared_coloc, gene_name == "HMGCR") %>% View()


#LIPA
lipa_coloc = dplyr::filter(all_shared_coloc, gene_name == "LIPA")

lipa_aragam = dplyr::filter(all_shared_coloc, gwas_signal == "GCST90132314_buildGRCh37_chr10:88247603-90247603") %>%
  dplyr::mutate(gwas_signal = "Aragam_2022_chr10:88247603-90247603") %>%
  dplyr::rename(metabolite_signal = signal2, met_gwas_PP4 = gwas_PP4, met_eqtl_PP4 = PP.H4)
write.table(lipa_aragam, "data/sup_tables/LIPA_locus_colocs.tsv", quote = F, row.names = F, sep = "\t")

#HMGCR
hmgcr_coloc = dplyr::filter(all_shared_coloc, gwas_signal == "GCST90132314_buildGRCh37_chr5:74360714-76360714")  %>%
  dplyr::mutate(gwas_signal = "Aragam_2022_chr5:74360714-76360714") %>%
  dplyr::rename(metabolite_signal = signal2, met_gwas_PP4 = gwas_PP4, met_eqtl_PP4 = PP.H4)
write.table(hmgcr_coloc, "data/sup_tables/HMGCR_locus_colocs.tsv", quote = F, row.names = F, sep = "\t")
