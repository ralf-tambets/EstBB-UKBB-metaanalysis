library("dplyr")
library("arrow")
library("igraph")

#Import fine mapped credible sets
complete_cs_df = arrow::read_parquet("data/big_data/finemapping/UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet") %>%
  dplyr::as_tibble()

#Import fine mapped credible sets (with overlapping cs removed)
unique_cs_df = arrow::read_parquet("data/big_data/finemapping/clpp_combined_with_MAF_no_overlaps.parquet") %>%
  dplyr::as_tibble()

#Count unique non-overlapping cs
length(unique(unique_cs_df$cs_id)) #116467
 
#Count the number of non_overlapping cs where PIP > 0.8
high_pip_unique_cs = dplyr::filter(unique_cs_df, pip > 0.8)
length(unique(high_pip_unique_cs$cs_id)) #31392

#Count the number of variants where PIP > 0.8 (including overlapping CS)
high_pip_all_cs = dplyr::filter(complete_cs_df, pip > 0.8)
high_pip_vars = dplyr::select(high_pip_all_cs, variant) %>% 
  dplyr::distinct()
length(unique(high_pip_cs$variant)) #2974

#Count how many of these variants are either missense or splice variants
splicing_scores = readr::read_tsv("data/big_data/finemapping/UKBB_EUR_fine_mapping_splicing_predictions.tsv.gz")
scores_selected = dplyr::mutate(splicing_scores, spliceai_max = pmax(abs(acceptor_score_spliceai_plus), 
                                                            abs(donor_score_spliceai_plus), 
                                                            abs(acceptor_score_spliceai_minus), 
                                                            abs(donor_score_spliceai_minus), na.rm = T)) %>%
  dplyr::mutate(alphagenome_max = pmax(abs(acceptor_score_alphagenome_plus), 
                                       abs(donor_score_alphagenome_plus), 
                                       abs(acceptor_score_alphagenome_minus), 
                                       abs(donor_score_alphagenome_minus), na.rm = T)) %>%
  dplyr::filter(spliceai_max > 0.1 | alphagenome_max > 0.1)

high_pip_splice = dplyr::left_join(high_pip_vars, scores_selected) %>% dplyr::filter(!is.na(region))
splice_variants = dplyr::select(high_pip_splice, gene_name, variant) %>% dplyr::distinct()

#Missense variants
#Explore missense variants
missense = readr::read_tsv("data/big_data/finemapping/UKBB_EUR_fine_mapping_missense_snps.tsv.gz")
missense_vars = dplyr::select(missense, chromosome, position, ref, alt, SYMBOL, Gene) %>% 
  dplyr::mutate(variant = paste(chromosome, position, ref, alt, sep = "_")) %>% 
  dplyr::semi_join(high_pip_vars) %>%
  dplyr::distinct()

#Merge missense and splice annotations into the fine mapping table
missense_df = dplyr::transmute(missense_vars, variant, missense_gene = SYMBOL, is_missense = TRUE)
splice_df = dplyr::group_by(splice_variants, variant) %>%
  dplyr::reframe(variant, splice_gene = paste(gene_name, collapse = ",")) %>%
  dplyr::distinct() %>%
  dplyr::mutate(is_splice = TRUE)
annotated_high_pip = dplyr::left_join(high_pip_all_cs, missense_df, by = "variant") %>%
  dplyr::left_join(splice_df, by = "variant")
write.table(annotated_high_pip,"data/big_data/finemapping/UKBB_EUR_annotated_high_PIP_variants.tsv",sep = "\t", row.names = F, quote = F)

#Count splice or missense variants
splice_or_missense = dplyr::filter(annotated_high_pip, is_missense | is_splice) %>% dplyr::select(variant, is_missense, is_splice, missense_gene, splice_gene, maf) %>% distinct()


# Explore low frequency variants
#Identify lead variants for low freq credible sets calculate their size
low_maf_cs = dplyr::group_by(unique_cs_df, cs_id) %>% 
  dplyr::arrange(-pip) %>% 
  dplyr::mutate(cs_size = n()) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(maf < 0.01) %>%
  dplyr::ungroup()
nrow(low_maf_cs) #10016

#Calculate percentage
10016/116467

#Identify overlapping cs
small_cs = dplyr::filter(low_maf_cs, cs_size < 100)

#Identify all variants belonging to low MAF cs
low_maf_cs_variants = dplyr::semi_join(unique_cs_df, small_cs, by = "cs_id") %>% dplyr::select(cs_id, variant)

#Identify non-overlapping cs
pairs = dplyr::left_join(low_maf_cs_variants, low_maf_cs_variants, by = "variant")

# make the graph of connected components
edges = dplyr::select(pairs,  cs_id.x,  cs_id.y) %>% dplyr::distinct() %>% as.matrix()
g <- igraph::graph_from_edgelist(edges, directed = F)
g_cc <- igraph::components(g)
members_df = tibble::tibble(signal2 = names(g_cc$membership), cluster = g_cc$membership)

#Count the number of clusters
length(g_cc$csize) #786

#Explore low-freq high PIP variants
dplyr::filter(annotated_high_pip, maf < 0.01) %>% 
  dplyr::select(variant) %>%
  dplyr::distinct() #583

#Count splice or missense variants with low maf
dplyr::filter(splice_or_missense, maf < 0.01) #135
#Fraction of low-MAF finemapped variants that are missense or splice-altering: 135/583
#Fraction of high-MAF finemapped variants that are missense or splice-altering: (415-135)/(3000-583)


#Import low freq lead variants
low_maf_leads = readr::read_tsv("data/big_data/finemapping/low_MAF_cluster_leaders.tsv") %>% 
  dplyr::mutate(variant = snp_alternate)

#How many low freq leads are contained in at least one credible set
cs_overlap = dplyr::semi_join(complete_cs_df, low_maf_leads, by = "variant")
length(unique(sort(cs_overlap$variant))) #n=96 variants are shared, 324 leads in that freq bracket

#Total variants
nrow(low_maf_leads) # n = 480 variants

#Below 0.1% frequency
rare_maf_leads = dplyr::filter(low_maf_leads, maf < 0.001) 
nrow(rare_maf_leads) #n = 156

#Explore how many low MAF credible sets contain missense variants
vep <- read_tsv("data/big_data/finemapping/low_MAF_cluster_leaders_VEP_response.txt") %>%
  rename(SNP = 1) %>%
  filter(grepl("missense_variant", Consequence) | grepl("splice_", Consequence))

missense_and_splice_variants <- low_maf_leads %>%
  select(CHR, SNP = snp_alternate, POS, ALL0, ALL1, MAF = maf, rsid, meta_EUR_associated_metabolites = metabolites_influenced_by_cluster, n_meta_EUR_associated_metabolites = n_metabolites_influenced_by_cluster) %>%
  inner_join(vep, by = "SNP")

unique_missense_splice = dplyr::select(missense_and_splice_variants, SNP, rsid, SYMBOL, MAF) %>% dplyr::distinct()

dplyr::filter(unique_missense_splice, MAF < .001) #19/156
dplyr::filter(unique_missense_splice, MAF > .001) #59/324
