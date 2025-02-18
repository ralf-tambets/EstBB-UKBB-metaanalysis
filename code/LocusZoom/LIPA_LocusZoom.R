library("dplyr")
library("locuszoomr")
library("EnsDb.Hsapiens.v86")

#Import LIPA expression
eqtl = readr::read_tsv("data/gwas_sumstats/Schmidel_2018_QTD000504.cc.tsv.gz")
lipa_eqtl = dplyr::filter(eqtl, molecular_trait_id == "ENSG00000107798") #%>%
  dplyr::select(-rsid) %>%
  dplyr::distinct()

lipa_eqtl_df = dplyr::transmute(lipa_eqtl, chrom = chromosome, pos = position, rsid, other_allele = ref, effect_allele = alt, p = pvalue, beta, se, variant = rsid) %>%
  as.data.frame()

lipa_loc <- locus(data = lipa_eqtl_df, gene = 'LIPA', flank = 5e4, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs1412445")
lipa_loc$data = lipa_loc$data[!is.na(lipa_loc$data$variant),]
lipa_ld = link_LD(lipa_loc, token = "97af801007eb", method = "matrix")
summary(lipa_loc)

pdf("figures/LIPA_Schmiedel_2018_mono__LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(lipa_loc)
dev.off()

#Import CAD
aragam = arrow::read_parquet("data/gwas_sumstats/Aragam_2022_GCST90132314_harmonized.parquet")

aragam_df = dplyr::transmute(aragam, chrom = CHROM, pos = GENPOS, rsid = rsid_UKBB, other_allele = REF, effect_allele = ALT, pvalue = 10^-LOG10P, beta = BETA, se = SE, variant = rsid_UKBB) %>%
  as.data.frame()

loc <- locus(data = aragam_df, gene = 'LIPA', flank = 5e4, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs1412445")
summary(loc)

pdf("figures/LIPA_Aragam_CAD_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()


#Import GlycA
ds = arrow::open_dataset("data/gwas_sumstats/GlycA_EstBB_UKBB_EUR_metaanalysis.parquet")
chr10_data = dplyr::filter(ds, CHROM == 10) %>% collect()

glyca_df = dplyr::transmute(chr10_data, chrom = CHROM, pos = GENPOS, 
                            other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
glyca_df = dplyr::mutate(glyca_df, rsid = variant)

loc <- locus(data = glyca_df, gene = 'LIPA', flank = 5e4, ens_db = "EnsDb.Hsapiens.v86", index_snp = "10_89243047_C_T")
summary(loc)

pdf("figures/LIPA_GlycA_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Import LDL
ds = arrow::open_dataset("data/gwas_sumstats/LDL_C_EstBB_UKBB_EUR_metaanalysis.parquet")
chr10_data = dplyr::filter(ds, CHROM == 10) %>% collect()

LDL_df = dplyr::transmute(chr10_data, chrom = CHROM, pos = GENPOS, 
                            other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
LDL_df = dplyr::mutate(LDL_df, rsid = variant)

loc <- locus(data = LDL_df, gene = 'LIPA', flank = 5e4, ens_db = "EnsDb.Hsapiens.v86", index_snp = "10_89243047_C_T")
summary(loc)

pdf("figures/LIPA_LDL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()


