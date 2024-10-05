library("readr")
library("dplyr")
library("purrr")
library("arrow")
library("tidyr")

args = commandArgs(trailingOnly=TRUE)

chunk_idx = as.numeric(args[1]) #index from 1 to n_metabolites / chunk_size
output_folder = args[2] #where to write results
chunk_size = as.numeric(args[3]) #how many metabolites should be analysed by one job
metabolite_list = args[4] #comma-separated names of metabolites
n_studies = as.numeric(args[5]) #number of studies included in meta-analysis

studies = c()
for (i in 1:n_studies){
  studies <- c(studies, args[5+i]) #directories where parquet files for metabolites are kept for different studies. Parquet files have to be kept in separate folders named after metabolites
}

studies <- as.list(studies)

print(paste("Chunk:", chunk_idx))
print(paste("Output folder:", output_folder))
print(paste("Chunk size:", chunk_size))
print(paste("Phenotype list:", metabolite_list))
print(paste("Number of studies:", n_studies))
print(paste("Included studies:", studies))

dir.create(file.path(output_folder), showWarnings = FALSE) #make output folder

all_metabolites <- colnames(readr::read_csv(metabolite_list)) #read metabolite names

# NB! Chunking copied from Urmo Võsa
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

#get metabolites for current job
metabolites <- unlist(chunks[[chunk_idx]])

for (i in 1:length(metabolites)){
  
  metabolite = metabolites[i]
  
  message(paste0("Working on metabolite ", i,": ", metabolite))
  
  chromosomes = 1:23
  
  #chromosome-wise solution
  merged_metabolite <- purrr::map_dfr(
    chromosomes,
    function(chromosome){
      message(chromosome)
      #find files
      metabolite_files <- lapply(studies, function(x) list.files(path = paste0(x, metabolite, "/"), pattern = ".parquet$", full.names = TRUE))
      
      chr_df_list <- lapply(metabolite_files, function(x){
        gwas_name <- basename(dirname(dirname(dirname(x)))) #only works on certain file tree
        
        open_dataset(x) %>% 
          filter(CHROM == chromosome) %>%
          mutate(AC = A1FREQ * N * 2) %>%
          mutate(w = 1 / (SE*SE), #weight
                 wES = w * BETA, #weighted effect size
                 sign = sign(BETA)) %>% #direction
          select(CHROM, GENPOS, ALLELE0, ALLELE1, N, INFO, AC, w, wES, sign) %>%
          collect() %>%
          rename_with(.fn = ~paste0(., "_", gwas_name), .cols = 5:10)
      }) 
      
      ci_multiplier = qnorm(0.975) #for 95% confidence interval
      
      chr_df <- purrr::reduce(chr_df_list, full_join) %>% #join by CHROM, GENPOS, ALLELE0, ALLELE1
        mutate(across(dplyr::starts_with("w_"), ~ifelse(is.na(.x),0,.x))) %>%
        mutate(across(dplyr::starts_with("wES_"), ~ifelse(is.na(.x),0,.x))) %>%
        mutate(across(dplyr::starts_with("AC_"), ~ifelse(is.na(.x),0,.x))) %>%
        mutate(across(dplyr::starts_with("N_"), ~ifelse(is.na(.x),0,.x))) %>%
        mutate(w_sum = rowSums(select(., starts_with("w_"))),
               wES_sum = rowSums(select(., starts_with("wES_"))),
               AC = rowSums(select(., starts_with("AC_"))),
               N = rowSums(select(., starts_with("N_"))),
               Effect = wES_sum/w_sum,
               StdErr = sqrt(1/w_sum),
               Z = abs(Effect) / StdErr,
               LOG10P = -pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE) * log10(exp(1)) - log10(2) #log10(exp(1)) converts log2 values to log10 values, subtracting log10(2) is the same as multiplying pvalue by 2 (necessary for two-tailed tests such as this)
        ) %>%
        mutate(across(dplyr::starts_with("sign_"), ~ifelse(is.na(.x),0,.x))) %>%
        unite("Direction", starts_with("sign_"), remove = TRUE, sep="") %>%
        mutate(Direction = gsub("-1", "-", Direction),
               Direction = gsub("0", "?", Direction),
               Direction = gsub("1", "+", Direction)) %>%
        select(CHROM, GENPOS, ALLELE0, ALLELE1, Effect, StdErr, Direction, LOG10P, N, AC, dplyr::starts_with("AC_"), dplyr::starts_with("INFO")) %>%
        arrange(CHROM, GENPOS)
    }
  )
  
  arrow::write_parquet(x = merged_metabolite, sink = paste0(output_folder, metabolite, "_EstBB_UKBB_full_metaanalysis.parquet"))
  
  rm(merged_metabolite)
  
  gc()
  
}

print("Done!")