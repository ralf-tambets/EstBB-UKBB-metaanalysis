library("dplyr")
library("locuszoomr")
library("EnsDb.Hsapiens.v86")

#Import T2D
gwas = arrow::read_parquet("data/gwas_sumstats/Suzuki_2024_All_Metal_LDSC-CORR_Neff.v2_harmonized.parquet")

gwas_df = dplyr::transmute(gwas, chrom = CHROM, pos = GENPOS, rsid = rsid_UKBB, other_allele = REF, effect_allele = ALT, pvalue = 10^-LOG10P, beta = BETA, se = SE, variant = rsid_UKBB) %>%
  as.data.frame()

loc <- locus(data = gwas_df, gene = 'BCAT2', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86") #, index_snp = "rs1412445"
summary(loc)

pdf("figures/BCAT2_Suzuki_T2D_LocusZoom.pdf", width = 8, height = 4)
locus_plot(loc)
dev.off()

#Import Valine
ds = arrow::open_dataset("data/gwas_sumstats/Val_EstBB_UKBB_EUR_metaanalysis.parquet")
chr7_data = dplyr::filter(ds, CHROM == 19) %>% collect()

metabolite_df = dplyr::transmute(chr7_data, chrom = CHROM, pos = GENPOS, 
                                 other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
metabolite_df = dplyr::mutate(metabolite_df, rsid = variant)

loc <- locus(data = metabolite_df, gene = 'BCAT2', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86")
summary(loc)

pdf("figures/BCAT2_Val_LocusZoom.pdf", width = 8, height = 4)
locus_plot(loc)
dev.off()

#Import IDL_TG
ds = arrow::open_dataset("data/gwas_sumstats/IDL_TG_EstBB_UKBB_EUR_metaanalysis.parquet")
chr7_data = dplyr::filter(ds, CHROM == 19) %>% collect()

metabolite_df = dplyr::transmute(chr7_data, chrom = CHROM, pos = GENPOS, 
                                 other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
metabolite_df = dplyr::mutate(metabolite_df, rsid = variant)

loc <- locus(data = metabolite_df, gene = 'BCAT2', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86")
summary(loc)

pdf("figures/BCAT2_IDL_TG_LocusZoom.pdf", width = 8, height = 4)
locus_plot(loc)
dev.off()


#Make PPM1K plot
#Import T2D
gwas = arrow::read_parquet("data/gwas_sumstats/Suzuki_2024_All_Metal_LDSC-CORR_Neff.v2_harmonized.parquet")

gwas_df = dplyr::transmute(gwas, chrom = CHROM, pos = GENPOS, rsid = rsid_UKBB, other_allele = REF, effect_allele = ALT, pvalue = 10^-LOG10P, beta = BETA, se = SE, variant = rsid_UKBB) %>%
  as.data.frame()

loc <- locus(data = gwas_df, gene = 'PPM1K', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs1440581") #, index_snp = "rs1412445"
summary(loc)

pdf("figures/PPM1K_Suzuki_T2D_LocusZoom.pdf", width = 8, height = 4)
locus_plot(loc)
dev.off()

#Import Valine
ds = arrow::open_dataset("data/gwas_sumstats/Val_EstBB_UKBB_EUR_metaanalysis.parquet")
chr7_data = dplyr::filter(ds, CHROM == 4) %>% collect()

metabolite_df = dplyr::transmute(chr7_data, chrom = CHROM, pos = GENPOS, 
                                 other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
metabolite_df = dplyr::mutate(metabolite_df, rsid = variant)

loc <- locus(data = metabolite_df, gene = 'PPM1K', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "4_88305270_T_C")
summary(loc)

pdf("figures/PPM1K_Val_LocusZoom.pdf", width = 8, height = 4)
locus_plot(loc)
dev.off()

#Make DBT plot
#Import T2D
gwas = arrow::read_parquet("data/gwas_sumstats/Suzuki_2024_All_Metal_LDSC-CORR_Neff.v2_harmonized.parquet")

gwas_df = dplyr::transmute(gwas, chrom = CHROM, pos = GENPOS, rsid = rsid_UKBB, other_allele = REF, effect_allele = ALT, pvalue = 10^-LOG10P, beta = BETA, se = SE, variant = rsid_UKBB) %>%
  as.data.frame()

loc <- locus(data = gwas_df, gene = 'DBT', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "1_100218726_TCAA_T") #, index_snp = "rs1412445"
summary(loc)

pdf("figures/DBT_Suzuki_T2D_LocusZoom.pdf", width = 8, height = 4)
locus_plot(loc)
dev.off()

#Import Valine
ds = arrow::open_dataset("data/gwas_sumstats/Val_EstBB_UKBB_EUR_metaanalysis.parquet")
chr7_data = dplyr::filter(ds, CHROM == 1) %>% collect()

metabolite_df = dplyr::transmute(chr7_data, chrom = CHROM, pos = GENPOS, 
                                 other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
metabolite_df = dplyr::mutate(metabolite_df, rsid = variant)

loc <- locus(data = metabolite_df, gene = 'DBT', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86")
summary(loc)

pdf("figures/DBT_Val_LocusZoom.pdf", width = 8, height = 4)
locus_plot(loc)
dev.off()