scores = readr::read_tsv("data/big_data/finemapping/meta_EUR_clpp_pip_filtered_spliceai_alphagenome_scores.tsv")

scores_selected = dplyr::mutate(scores, spliceai_max = pmax(abs(acceptor_score_spliceai_plus), 
                                          abs(donor_score_spliceai_plus), 
                                          abs(acceptor_score_spliceai_minus), 
                                          abs(donor_score_spliceai_minus), na.rm = T)) %>%
  dplyr::mutate(alphagenome_max = pmax(abs(acceptor_score_alphagenome_plus), 
                                       abs(donor_score_alphagenome_plus), 
                                       abs(acceptor_score_alphagenome_minus), 
                                       abs(donor_score_alphagenome_minus), na.rm = T)) %>%
  dplyr::filter(spliceai_max > 0.1 | alphagenome_max > 0.1)
write.table(scores_selected, "data/big_data/finemapping/splicing_cs_predictions.tsv", sep = "\t", quote = F, row.names = F)

#Count the number of unique variants that disrupt splicing
splicing_count = dplyr::select(scores_selected, variant) %>% dplyr::distinct()

#Explore missense variants
missense = readr::read_tsv("data/big_data/meta_EUR_clpp_pip_filtered_missense_snps.tsv.gz")
missense_vars = dplyr::select(missense, chromosome, position, ref, alt, SYMBOL, Gene) %>% 
  dplyr::mutate(variant = paste(chromosome, position, ref, alt, sep = "_")) %>% 
  dplyr::distinct() %>% 
  dplyr::select(variant)

#Count misense and/or splice variants
shared = dplyr::bind_rows(missense_vars, splicing_count) %>% dplyr::group_by(variant) %>% summarise(count = n())
