library("dplyr")

#Import coloc results
coloc_table = readr::read_tsv("https://zenodo.org/records/17945143/files/meta_EUR_big_coloc_151025.tsv.gz")

#Identify clusters that contain GlycA
glyca_clusters = dplyr::filter(coloc_table, metabolite == "GlycA")$cluster
glyca_coloc_clusters = dplyr::filter(coloc_table, cluster %in% glyca_clusters)

#Number of unique clusters
length(unique(glyca_coloc_clusters$cluster))

#Number of clusters that also contain CRP
length(unique(dplyr::filter(glyca_coloc_clusters, signal1_trait == "C-reactive protein")$cluster))

#Count the number of colocalising loci per trait in FinnGen+MVP+UKBB meta-analysis
dplyr::filter(glyca_coloc_clusters, source %in% c("FinnGen+MVP+UKBB"), metabolite == "GlycA") %>% 
  dplyr::group_by(signal1_trait) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(-count) %>% View()