#!/usr/bin/env Rscript
#based on work by Urmo Võsa
args = commandArgs(trailingOnly = TRUE)

print(paste("Chunk:", args[1]))
print(paste("Results folder", args[2]))
print(paste("Output path", args[3]))
print(paste("Chunk size", args[4]))
print(paste("Phenotype list", args[5]))

library(dplyr)
library(readr)
library(arrow)

# paths
gwas_result_path <- args[2]
output_path <- args[3]

# functions
IdentifyLeadSNPs <- function(data,
                             window = 1000000,
                             Pthresh = 5e-8,
                             snp_id_col = "snp",
                             snp_chr_col = "chr",
                             snp_pos_col = "pos",
                             eff_all_col = "ea",
                             other_all_col = "nea",
                             beta_col = "beta",
                             se_col = "se",
                             a1_freq_col = "A1FREQ",
                             log10p_col = NULL) {
  
  log10Pthresh = -log10(Pthresh)
  
  data <- data.frame(
    SNP = data[[snp_id_col]],
    chr = data[[snp_chr_col]],
    pos = data[[snp_pos_col]],
    A1FREQ = data[[a1_freq_col]],
    ea = data[[eff_all_col]],
    nea = data[[other_all_col]],
    beta = data[[beta_col]],
    se = data[[se_col]],
    LOG10P = data[[log10p_col]]
  )
  data$Z <- data$beta / data$se
  
  data_f <- data[data$LOG10P > as.numeric(log10Pthresh), ]
  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]
  
  while (max(data_f$LOG10P) >= log10Pthresh) {
    lead_snp <- data_f[abs(data_f$Z) == max(abs(data_f$Z)), ] #leia suurima Z skooriga variandid
    if (nrow(lead_snp) > 1) { #kui võrdseid on on mitu, vaata esimest
      lead_snp <- lead_snp[1, ]
    }
    res <- rbind(res, lead_snp) #lisa see res'i
    data_f <- data_f[!(data_f$chr == lead_snp$chr & data_f$pos > lead_snp$pos - window & data_f$pos < lead_snp$pos + window), ] #jäta alles need, mis on window'st väljas (kas teisel kromosoomil või piisavalt kaugel (window on tegelikult 2*window))
    # message(paste("Added:", lead_snp$snp_chr, lead_snp$snp_pos))
    if (nrow(data_f) == 0) {
      break
    }
  }
  return(res)
}

# Get phenotypes

all_metabolites <- colnames(readr::read_csv(args[5]))

# Define the chunk size
chunk_size <- as.numeric(args[4])

# Calculate the number of chunks
num_chunks <- length(all_metabolites) %/% chunk_size

# Create a grouping vector for splitting
group_vector <- rep(1:num_chunks, each = chunk_size)

# Split the vector into chunks
chunks <- split(all_metabolites[1:length(group_vector)], group_vector)

# Add the incomplete last chunk
remaining_elements <- all_metabolites[((num_chunks * chunk_size + 1):length(all_metabolites))]
if (length(remaining_elements) > 0) {
  chunks <- c(chunks, list(remaining_elements))
}

chunk_idx <- as.numeric(args[1])

metabolites <- unlist(chunks[[chunk_idx]])
colnames_written = FALSE
for (i in 1:length(metabolites)){
  metabolite <- metabolites[i]
  
  file_to_analyse <- list.files(path = paste0(gwas_result_path, metabolite), pattern = ".parquet$", full.names = TRUE)[1]
  
  inp <- arrow::open_dataset(file_to_analyse) %>%
    collect()

  message(paste(metabolite, "read in!"))
  
  # find the number of lead variants
  inp_f <- inp #%>% 
    #dplyr::filter(A1FREQ > 0.01 & A1FREQ < 0.99)

  
  if (inp_f %>% dplyr::filter(LOG10P >= -log10(5e-8)) %>% nrow() > 0){
    lead_var <- IdentifyLeadSNPs(inp_f, 
                                 window = 1000000,
                                 Pthresh = 5e-8,
                                 snp_id_col = "ID",
                                 snp_chr_col = "CHROM",
                                 snp_pos_col = "GENPOS",
                                 eff_all_col = "ALLELE0",
                                 other_all_col = "ALLELE1",
                                 a1_freq_col = "A1FREQ",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 log10p_col = "LOG10P")
    
    readr::write_tsv(
      x = data.frame(metabolite = metabolite, lead_var), 
      file = paste0(output_path, "/chunk", chunk_idx, "_lead_variants.txt"), 
      append = TRUE, 
      col_names = !colnames_written
      )
    colnames_written = TRUE
  } else {lead_var <- data.frame(NULL)}
  message(paste(metabolite, "lead variants found!"))
  
  
  # Calculate lambda1000
  
  lambda_001 <- median(inp_f$CHISQ) / qchisq(0.5, df = 1)
  

  # write out the results
  temp <- data.frame(
    metabolite = metabolite, 
    nr_variants = nrow(inp),
    nr_maf001_variants = nrow(inp_f),
    nr_chr = length(unique(inp$CHROM)),
    lambda = median(inp$CHISQ) / qchisq(0.5, df = 1),
    lambda_maf001 = lambda_001,
    max_beta = max(inp$BETA),
    min_beta = min(inp$BETA),
    max_se = max(inp$SE),
    min_se = min(inp$SE),
    most_sig_logP = max(inp$LOG10P), 
    nr_lead_variants = nrow(lead_var)
  ) %>%
    dplyr::bind_cols(inp_f %>% 
                       dplyr::select(CHISQ, BETA, SE, LOG10P) %>% 
                       dplyr::filter(!is.na(CHISQ)) %>%
                       dplyr::summarise(max_beta_maf001 = max(BETA),
                                        min_beta_maf001 = min(BETA),
                                        max_se_maf001 = max(SE),
                                        min_se_maf001 = min(SE),
                                        most_sig_logP_maf0012 = max(LOG10P)
                       ),
    ) %>%
    dplyr::relocate(nr_lead_variants, .after = last_col())
  message(paste(metabolite, "processed!"))
  
  
  readr::write_tsv(x = temp, 
                   file = paste0(output_path, "/chunk", chunk_idx, "_GWAS_QC_metrics.txt"),
                   append = TRUE,
                   col_names = (i == 1)
  )
  
  message(paste(metabolite, "data written!"))

  gc()
}

message("Done!")