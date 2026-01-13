library("dplyr")
library("arrow")
library("igraph")

#Import fine mapped credible sets
cs_df = arrow::read_parquet("data/big_data/finemapping/UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet") %>%
  dplyr::as_tibble()

#Import low freq lead variants
low_maf_leads = readr::read_tsv("data/big_data/finemapping/low_MAF_cluster_leaders.tsv") %>% dplyr::mutate(variant = snp_alternate)

#Explore Lactate credible set at the GP6 locus
dplyr::filter(cs_df, cs_id == "Lactate_chr19:54028085-56027610_L1") %>% View()

#GRK5
dplyr::filter(cs_df, variant == "10_119250744_A_G") %>% View()


#Identify variant present in multiple credible sets for the same trait
double_count = dplyr::group_by(cs_df, molecular_trait_id, variant) %>% dplyr::mutate(variant_count = n())

#Idntify credible sets that share at least one variant with another credible set for the same trait
shared_cs = dplyr::filter(double_count, variant_count > 1) %>% dplyr::ungroup() %>% dplyr::select(cs_id) %>% dplyr::distinct()

#Keep these credible sets
unique_cs = dplyr::anti_join(cs_df, shared_cs)
unqiue_cs_count = dplyr::select(unique_cs, cs_id) %>% dplyr::distinct() %>% nrow()
overlapping_cs_count = nrow(shared_cs)
total_count = unqiue_cs_count + overlapping_cs_count/2

#Identify unique fine mapped variants for each trait
#Distinct variant per credible set
dplyr::filter(cs_df, pip > 0.5) %>% dplyr::select(molecular_trait_id, variant) %>% dplyr::distinct()
#Distinct variants overall
dplyr::filter(cs_df, pip > 0.5) %>% dplyr::select(variant) %>% dplyr::distinct()


#Identify lead variants for credible sets calculate their size
low_maf_cs = dplyr::group_by(cs_df, cs_id) %>% 
  dplyr::arrange(-pip) %>% 
  dplyr::mutate(cs_size = n()) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(maf < 0.01) %>%
  dplyr::ungroup()

#Identify all variants belonging to low MAF cs
low_maf_cs_variants = dplyr::semi_join(cs_df, low_maf_cs, by = "cs_id") %>% dplyr::select(cs_id, variant)

#Identify non-overlapping cs
pairs = dplyr::left_join(low_maf_cs_variants, low_maf_cs_variants, by = "variant")

# make the graph of connected components
edges = dplyr::select(pairs,  cs_id.x,  cs_id.y) %>% dplyr::distinct() %>% as.matrix()
g <- igraph::graph_from_edgelist(edges, directed = F)
g_cc <- igraph::components(g)
members_df = tibble::tibble(signal2 = names(g_cc$membership), cluster = g_cc$membership)

#Count the number of clusters
length(g_cc$csize)

#How many low freq leads are contained in at least one credible set
cs_overlap = dplyr::semi_join(low_maf_cs, low_maf_leads, by = "variant")
length(unique(sort(cs_overlap$variant))) #n=70 variants are shared, 264 leads in that freq bracket

#Total variants
nrow(low_maf_leads) # n = 987 variants

#Below 0.1% frequency
rare_mad_leads = dplyr::filter(low_maf_leads, maf < 0.001) 
nrow(rare_mad_leads) #n = 723


#Explore how many low MAF credible sets contain missense variants
missense = readr::read_tsv("data/big_data/finemapping/meta_EUR_clpp_pip_filtered_missense_snps.tsv.gz")
missense_vars = dplyr::select(missense, chromosome, position, ref, alt, SYMBOL, Gene) %>% 
  dplyr::mutate(variant = paste(chromosome, position, ref, alt, sep = "_")) %>% 
  dplyr::distinct()

#Explore low MAF finemapped variants that are missense
dplyr::left_join(missense_vars, low_maf_cs) %>% 
  dplyr::select(variant, SYMBOL, molecular_trait_id, pip, maf) %>% 
  dplyr::filter(pip > 0.8) %>% View()

#Look at unique genes only
dplyr::left_join(missense_vars, low_maf_cs) %>% 
  dplyr::filter(pip > 0.8) %>%
  dplyr::select(variant, SYMBOL) %>%
  dplyr::distinct()

#Which of these low-freq variants are missense or splice?
s6 = readr::read_tsv("data/big_data/finemapping/Table_S6.tsv") %>%
  dplyr::mutate(CHR = ifelse(CHR == 23, "X",CHR)) %>%
  dplyr::mutate(variant = paste(CHR, POS, ALL0, ALL1, sep = "_")) %>%
  dplyr::left_join(dplyr::select(low_maf_leads, maf, LOG10P, variant))

#Leads with more stringent p-value threshold
dplyr::filter(low_maf_leads, LOG10P > 10, maf < 0.001)


