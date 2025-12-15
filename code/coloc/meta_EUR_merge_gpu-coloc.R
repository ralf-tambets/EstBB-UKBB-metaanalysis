library("data.table")
library("dplyr")
library("igraph")

#Import coloc results
biobank_meta = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/FinnGen+MVP+UKBB_coloc/meta_EUR_FinnGen+MVP+UKBB_coloc_results.tsv") %>%
  dplyr::distinct()
ukbppp = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/UKBPPP_coloc/UKBPPP_Combined_meta_EUR_results.tsv") %>%
  dplyr::distinct()
panukbb = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/PANUKBB_coloc/PANUKBB_EUR_meta_EUR_results.tsv") %>%
  dplyr::distinct()
interval = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/INTERVAL_coloc/INTERVAL_meta_EUR_results.tsv") %>%
  dplyr::distinct()
FinnGen = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/FinnGen_lbf_coloc/meta_EUR_FinnGen_lbf_results.tsv") %>%
  dplyr::distinct()
MVP = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/MVP_coloc/MVP_META_meta_EUR_results.tsv") %>%
  dplyr::distinct()
eqtl_catalogue_ge = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/eQTL_Catalogue_coloc/ge_microarray_aptamer_meta_EUR_results.tsv") %>%
  dplyr::distinct()
eqtl_catalogue_lc = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/eQTL_Catalogue_coloc/leafcutter_meta_EUR_results.tsv") %>%
  dplyr::distinct()
eqtl_catalogue = dplyr::bind_rows(eqtl_catalogue_ge, eqtl_catalogue_lc)
suzuki_aragam = readr::read_tsv("data/big_data/annotated_gpu-coloc_results/Suzuki_Aragam_coloc/meta_EUR_Suzuki_Aragam_coloc_results.tsv") %>%
  dplyr::distinct()
suzuki = dplyr::filter(suzuki_aragam, signal1 %like% "All_Metal")
aragam = dplyr::filter(suzuki_aragam, signal1 %like% "GCST90132314")

#Explore GP6 locus
biobank_gp6 = dplyr::filter(biobank_meta, gene_symbols %like% "GP6")
ukbpp_gp6 = dplyr::filter(ukbpp, HGNC.symbol == "GP6")
panukbb_gp6 = dplyr::filter(panukbb, gene_symbols %like% "GP6")
interval_gp6 = dplyr::filter(interval, gene_symbols %like% "GP6")
eqtl_catalogue_gp6 = dplyr::filter(eqtl_catalogue, gene_symbols %like% "GP6")
finngen_gp6 = dplyr::filter(FinnGen, gene_symbols %like% "GP6")
mvp_gp6 = dplyr::filter(MVP, gene_symbols %like% "GP6")

#Merge all pairs together in an unform format
ukbppp_selected = dplyr::transmute(ukbppp, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, UKBPPP_ProteinID, signal1_trait = HGNC.symbol, source = "UKB-PPP")
biobank_meta_selected = dplyr::transmute(biobank_meta, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = trait, source = "FinnGen+MVP+UKBB")
panukbb_selected = dplyr::transmute(panukbb, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = description, source = "Pan-UKBB")
mvp_selected = dplyr::transmute(MVP, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = MVP_trait, source = "MVP")
eqtl_catalogue_selected = dplyr::transmute(eqtl_catalogue, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = eqtl_catalogue_gene_name, source = "eQTLCatalogue", quant_method = eqtl_catalogue_quant_method, molecular_trait_id = eqtl_catalogue_molecular_trait_id)
interval_selected = dplyr::transmute(interval, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = gene_name, source = "INTEREVAL_eQTL", quant_method = "ge", molecular_trait_id = gene_id)
finngen_selected = dplyr::transmute(FinnGen, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = diagnosis, source = "FinnGen")
suzuki_selected = dplyr::transmute(suzuki, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = "T2D", source = "Suzuki_2024")
aragam_selected = dplyr::transmute(aragam, PP.H3, PP.H4, signal1, lead1, signal2, lead2, metabolite, maf2, log10p2, gene_symbols, signal1_trait = "CAD", source = "Aragam_2022")

merged = dplyr::bind_rows(ukbppp_selected, biobank_meta_selected, panukbb_selected, mvp_selected, eqtl_catalogue_selected, interval_selected, finngen_selected, suzuki_selected, aragam_selected)

#Identify colocalisation clusters
high_conf = dplyr::filter(merged, PP.H4 > 0.9)

# make the graph of connected components
edges = dplyr::select(high_conf, signal1, signal2) %>% as.matrix()
g <- igraph::graph_from_edgelist(edges, directed = F)
g_cc <- igraph::components(g)
members_df = tibble::tibble(signal2 = names(g_cc$membership), cluster = g_cc$membership)

clustered_signals = dplyr::left_join(high_conf, members_df)

#Count unique metabolites in each cluster
metabolite_counts = dplyr::select(clustered_signals, cluster, metabolite) %>% 
  dplyr::distinct() %>% dplyr::group_by(cluster) %>% dplyr::summarise(n_metabolites = n())

counted_signals = dplyr::left_join(clustered_signals, metabolite_counts)
write.table(counted_signals, "data/big_data/annotated_gpu-coloc_results/meta_EUR_big_coloc_151025.tsv", sep = "\t", row.names = F, quote = F)