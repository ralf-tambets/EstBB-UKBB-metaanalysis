library("dplyr")
library("locuszoomr")
library("EnsDb.Hsapiens.v86")

#Import ADCY5 expression
eqtl = readr::read_tsv("data/gwas_sumstats/PISA_QTD000574.cc.tsv.gz")
adcy5_eqtl = dplyr::filter(eqtl, molecular_trait_id == "ENSG00000173175") #%>%
  #dplyr::select(-rsid) %>%
  #dplyr::distinct()

eqtl_df = dplyr::transmute(adcy5_eqtl, chrom = chromosome, pos = position, rsid, other_allele = ref, effect_allele = alt, p = pvalue, beta, se, variant = rsid) %>%
  as.data.frame()

loc <- locus(data = eqtl_df, gene = 'ADCY5', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs11708067")
#loc$data = loc$data[!is.na(loc$data$variant),]

pdf("figures/ADCY5_PISA_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Import T2D
gwas = arrow::read_parquet("data/gwas_sumstats/Suzuki_2024_All_Metal_LDSC-CORR_Neff.v2_harmonized.parquet")

gwas_df = dplyr::transmute(gwas, chrom = CHROM, pos = GENPOS, rsid = rsid_UKBB, other_allele = REF, effect_allele = ALT, pvalue = 10^-LOG10P, beta = BETA, se = SE, variant = rsid_UKBB) %>%
  as.data.frame()

loc <- locus(data = gwas_df, gene = 'ADCY5', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs11708067") #, index_snp = "rs1412445"
summary(loc)

pdf("figures/ADCY5_Suzuki_T2D_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()


#Import Glucose
ds = arrow::open_dataset("data/gwas_sumstats/Glucose_EstBB_UKBB_EUR_metaanalysis.parquet")
chr7_data = dplyr::filter(ds, CHROM == 3) %>% collect()

metabolite_df = dplyr::transmute(chr7_data, chrom = CHROM, pos = GENPOS, 
                            other_allele = ALLELE0, effect_allele = ALLELE1, pvalue = 10^-LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
metabolite_df = dplyr::mutate(metabolite_df, rsid = variant)

loc <- locus(data = metabolite_df, gene = 'ADCY5', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "3_123346931_A_G")
summary(loc)

pdf("figures/ADCY5_Glucose_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()



#Perform MR
library("MendelianRandomization")
eQTL_beta = dplyr::filter(eqtl_df, variant == "rs11708067")
gwas_beta = dplyr::filter(gwas_df, variant == "rs11708067")
metabolite_meta = dplyr::filter(metabolite_df, variant == "3_123346931_A_G")

input = mr_input(bx = c(0.89443), bxse = c(0.0870403), by = c(-0.0812), , byse = c(0.0039))
mr_ivw(input)

pdf("figures/ADCY5_MR.pdf", width = 4, height = 4)
mr_plot(input, error = TRUE, orientate = FALSE, line = "ivw")
dev.off()

