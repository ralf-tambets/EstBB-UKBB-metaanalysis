library(arrow)
library(dplyr)
library(readr)

options(scipen=999)

z_folder <- "files_to_upload/z_files/"
master_folder_dnanexus <- "files_to_upload/master_files_dnanexus/"
incl_folder <- "files_to_upload/"

maf_limit = 0.001
incl_file = "EUR_samples.incl"

dir.create(z_folder, recursive = T)
dir.create(master_folder, recursive = T)

loci <- read_tsv("meta_EUR_loci_DNAnexus.tsv", col_names = c("chromosome", "start", "end", "width")) %>%
  filter(chromosome < 23)

samples <- read_delim("GEL_EUR_chr19.sample", delim = " ") %>%
  select(1:2)
samples_to_keep <- read_tsv("../../UKBB_EUR_all_metabolites_500k.tsv") %>%
  select(1:2) %>%
  filter(IID %in% samples$ID_1)

n_samples <- nrow(samples_to_keep)

write_tsv(x = samples_to_keep %>% select(IID),
          file = paste0(incl_folder, incl_file),
          col_names = F)

rm(samples)

for (i in 1:nrow(loci)){
  chr = loci$chromosome[i]
  locus_start = loci$start[i] %>% as.integer()
  locus_end = loci$end[i] %>% as.integer()
  
  region_name = paste0("chr", chr, "_", locus_start, "_", locus_end)
  
  message(region_name)
  
  #get EUR variants with maf > 0.1% in each locus
  EUR_high_maf_variants <- open_dataset("../Total_BCAA_EUR.parquet") %>%
    filter(CHROM == chr) %>%
    filter(between(GENPOS, locus_start, locus_end)) %>%
    collect() %>%
    mutate(maf = pmin(A1FREQ, 1-A1FREQ)) %>%
    filter(maf > maf_limit) %>%
    arrange(CHROM, GENPOS, ALLELE0, ALLELE1) %>%
    select(rsid = ID, CHROM, GENPOS, ALLELE0, ALLELE1)
  
  #reformat into z file
  z <- EUR_high_maf_variants %>%
    rename(chromosome = CHROM, position = GENPOS, allele1 = ALLELE0, allele2 = ALLELE1) %>%
    mutate(chromosome = paste0("chr", chromosome)) %>%
    mutate(rsid = ifelse(test = is.na(rsid) | !startsWith(rsid, "rs"), 
                         yes = paste0(chromosome, ":", position, "_", allele1, "_", allele2),
                         no = rsid)) %>%
    na.omit() %>%
    arrange(chromosome, position, allele1, allele2) %>%
    distinct()
  
  message(nrow(z))
  
  write_delim(x = z,
              file = paste0(z_folder, "EUR_", region_name, ".z"),
              delim = " ")
  
  # create corresponding master file for DNAnexus
  master_dnanexus <- data.frame(
    z = paste0("EUR_", region_name, ".z"),
    bgen = paste0("ukb21008_c", chr,"_b0_v1.bgen"),
    bgi = paste0("ukb21008_c", chr,"_b0_v1.bgen.bgi"),
    sample = paste0("ukb21008_c", chr,"_b0_v1.sample"),
    bcor = paste0("EUR_", region_name, ".bcor"),
    ld = paste0("EUR_", region_name, ".ld"),
    n_samples = n_samples,
    incl = incl_file
  )
  
  write_delim(x = master_dnanexus,
              file = paste0(master_folder_dnanexus, "EUR_", region_name, ".master"),
              delim = ";",
              na = "")
}
