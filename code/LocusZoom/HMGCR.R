library("dplyr")
library("locuszoomr")
library("EnsDb.Hsapiens.v86")

#Import HMGCR spliceQTL
sqtl = readr::read_tsv("data/gwas_sumstats/Alasoo_2018_QTD000005.cc.tsv.gz")
hmgcr_sqtl = dplyr::filter(sqtl, molecular_trait_id == "5:75354697:75355365:clu_9027_-") #%>%
dplyr::select(-rsid) %>%
  dplyr::distinct()

hmgcr_df = dplyr::transmute(hmgcr_sqtl, chrom = chromosome, pos = position, rsid, other_allele = ref, effect_allele = alt, p = pvalue, beta, se, variant = rsid) %>%
  as.data.frame()
loc <- locus(data = hmgcr_df, gene = 'HMGCR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs12916")
summary(loc)

pdf("figures/HMGCR_Alasoo_2018_sQTL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Import CAD
aragam = arrow::read_parquet("data/gwas_sumstats/Aragam_2022_GCST90132314_harmonized.parquet")

aragam_df = dplyr::transmute(aragam, chrom = CHROM, pos = GENPOS, rsid = rsid_UKBB, other_allele = REF, effect_allele = ALT, pvalue = 10^-LOG10P, beta = BETA, se = SE, variant = rsid_UKBB) %>%
  as.data.frame()

loc <- locus(data = aragam_df, gene = 'HMGCR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs12916")
summary(loc)

pdf("figures/HMGCR_Aragam_CAD_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Import LDL
ds = arrow::open_dataset("data/gwas_sumstats/LDL_C_EstBB_UKBB_EUR_metaanalysis.parquet")
chr5_data = dplyr::filter(ds, CHROM == 5) %>% collect()

LDL_df = dplyr::transmute(chr5_data, chrom = CHROM, pos = GENPOS, 
                          other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
LDL_df = dplyr::mutate(LDL_df, rsid = variant)

loc <- locus(data = LDL_df, gene = 'HMGCR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "5_75360714_T_C")
summary(loc)

pdf("figures/HMGCR_LDL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()


#Import GlycA
ds = arrow::open_dataset("data/gwas_sumstats/GlycA_EstBB_UKBB_EUR_metaanalysis.parquet")
chr5_data = dplyr::filter(ds, CHROM == 5) %>% collect()

LDL_df = dplyr::transmute(chr5_data, chrom = CHROM, pos = GENPOS, 
                          other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
LDL_df = dplyr::mutate(LDL_df, rsid = variant)

loc <- locus(data = LDL_df, gene = 'HMGCR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "5_75360714_T_C")
summary(loc)

pdf("figures/HMGCR_GlycA_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()
