library(dplyr)
library(readr)
library(GenomicRanges)

# Original created by Ida Rahu, recursively shrinks regions until they fit into win_size
find_regions <- function(data, win_size, win_shrink_ratio, max_region_width) {
  final <- data[-c(1:nrow(data)), ]
  data_copy <- data
  
  while (nrow(data) > 0) {
    lead_snp <- data[abs(data$log10P) == max(data$log10P), ]
    final <- rbind(final, lead_snp)
    data <- data[!(data$chromosome == lead_snp$chromosome &
                     data$position > lead_snp$position-win_size &
                     data$position < lead_snp$position+win_size), ]
  }
  final$start <- final$position - win_size
  final$end <- final$position + win_size
  
  final$region <-  ifelse (final$start < 0, paste0(final$chromosome, ':', 0, '-', final$end),
                           paste0(final$chromosome, ':', final$start, '-', final$end))
  
  final <- final %>% arrange(chromosome, position)
  
  gr_object <- GRanges(seqnames=final$chromosome, ranges=IRanges(start=final$start, end=final$end))
  reduced_gr <- reduce(gr_object)
  
  regions_to_keep <- reduced_gr[width(reduced_gr) <= max_region_width]
  regions_to_shrink <- reduced_gr[width(reduced_gr) > max_region_width]
  
  data_gr <- with(data_copy, GRanges(seqnames=chromosome, ranges=IRanges(start=position, end=position)))
  overlap_indexes <- findOverlaps(data_gr, regions_to_shrink)
  subset_data <- data_copy[queryHits(overlap_indexes), ]
  
  if (nrow(subset_data) > 0) {
    regions_to_keep_new <- find_regions(subset_data, win_size=win_size * win_shrink_ratio,
                           win_shrink_ratio=win_shrink_ratio, max_region_width=max_region_width)
    final_regions <- c(regions_to_keep, regions_to_keep_new)
  } else {
    final_regions <- regions_to_keep
  }
  
  return(sort(final_regions))
}

# clips away parts of regions that fall into MHC 
clip_MHC_regions <- function(regions, MHC_start = 28510120, MHC_end = 33480577){
  chr6_clipped_regions <- regions %>%
    filter(chromosome == 6) %>%
    mutate(is_start_in_MHC = start >= MHC_start & start <= MHC_end) %>%
    mutate(is_end_in_MHC = end >= MHC_start & end <= MHC_end) %>%
    mutate(is_completely_in_MHC = is_start_in_MHC & is_end_in_MHC) %>%
    filter(!is_completely_in_MHC) %>%
    mutate(start = ifelse(test = is_start_in_MHC, yes = MHC_end, no = start)) %>%
    mutate(end = ifelse(test = is_end_in_MHC, yes = MHC_start, no = end)) %>%
    mutate(width = end-start + 1) %>%
    select(1:4)
  
  result = regions %>% 
    filter(chromosome != 6) %>%
    rbind(chr6_clipped_regions) %>%
    arrange(chromosome, start, end) 
  
  return(result)
}

#if two regions overlap, clips the wider one
clip_overlapping_regions <- function(regions){
  overlapping <- regions %>% 
    group_by(chromosome) %>%
    mutate(next_start = start[c(2:n(), NA)]) %>%
    mutate(prev_end = end[c(NA, 1:(n()-1))]) %>%
    rowwise() %>%
    mutate(is_overlapping_with_next = !is.na(next_start) & next_start <= end) %>%
    mutate(is_overlapping_with_prev = !is.na(prev_end) & prev_end >= start) %>%
    mutate(is_overlapping = is_overlapping_with_prev || is_overlapping_with_next) %>%
    ungroup() %>%
    filter(is_overlapping) %>%
    mutate(overlap_size = ifelse(is_overlapping_with_prev, prev_end-start, end-next_start)) %>%
    mutate(next_width = width[c(2:n(), NA)]) %>%
    mutate(prev_width = width[c(NA, 1:(n()-1))]) %>%
    mutate(region_name=paste(chromosome, start, end, sep = "_")) %>%
    arrange(chromosome, start, end)
  
  clipped <- overlapping %>%
    mutate(new_end = ifelse(is_overlapping_with_next & (width > next_width), next_start, end)) %>%
    mutate(new_start = ifelse(is_overlapping_with_prev & (width > prev_width), prev_end, start)) %>%
    mutate(new_width = new_end-new_start + 1) %>%
    select(chromosome, start = new_start, end = new_end, width = new_width)
    
  
  result <- regions %>% 
    mutate(region_name=paste(chromosome, start, end, sep = "_")) %>%
    filter(!(region_name %in% overlapping$region_name)) %>%
    select(-region_name) %>%
    rbind(clipped) %>%
    arrange(chromosome, start, end)
  
  return(result)
}

# Wrapper to find regions and clip out MHC if needed
find_regions_wrapper <- function(data, win_size, win_shrink_ratio, max_region_width, MHC_start = 28510120, MHC_end = 33480577, remove_MHC = TRUE, clip_overlaps = FALSE) {
  regions <- find_regions(data, win_size, win_shrink_ratio, max_region_width) %>%
    as.data.frame() %>%
    select(1:4) %>%
    mutate(start = pmax(1, start))
  
  lvls <- levels(regions$seqnames)
  
  regions <- regions %>% 
    mutate(chromosome = as.numeric(lvls[seqnames])) %>%
    select(chromosome, start, end, width) %>%
    arrange(chromosome) %>%
    distinct()
  
  if (remove_MHC) {
    regions <- clip_MHC_regions(regions)
  }
  
  if (clip_overlaps){
    regions <- clip_overlapping_regions(regions)
  }
  
  return(regions)
}

# Keeps unique leads by chromosome and position. Retains highest log10P value across metabolites.
# Input - dataframe with columns c("metabolite", "rsid", "chromosome", "position", "ALL0", "ALL1")
# Output - dataframe with unique rows of c("rsid", "chromosome", "position", "log10P")
get_unique_leads <- function(df){
  return(df %>%
           group_by(chromosome, position, ALL0, ALL1) %>%
           filter(log10P == max(log10P)) %>%
           ungroup() %>%
           select(-c(metabolite, ALL0, ALL1)) %>%
           distinct() %>%
           arrange(chromosome, position))
}

window = 1*10^6 #in case of a single lead in the region, creates regions of lead +/- window bp (default +/- 1Mb)
max_width = 6*10^6 #maximum width of a region (default 6Mb)
win_shrink_ratio = 0.95 #how much to shrink window in every recursive step (default 5%) 

EST_leads <- read_tsv("~/Documents/tööasjad/biit/NMR/LD_2025/loci_2025_MAF_filtered/01_reformatted/EST_lead_variants_reformatted_MAF_filtered.tsv") %>%
  select(metabolite, rsid, chromosome = CHR, position = POS, ALL0, ALL1, log10P = LOG10P) 

meta_EUR_leads <- read_tsv("...") %>%
  select(metabolite, rsid, chromosome = CHR, position = POS, ALL0, ALL1, log10P = LOG10P) 

combined_EST <- get_unique_leads(EST_leads)

combined_EST_regions <- find_regions_wrapper(data = combined_EST, 
                                             win_size = window, 
                                             win_shrink_ratio = win_shrink_ratio, 
                                             max_region_width = max_width, 
                                             remove_MHC = T, 
                                             clip_overlaps = F) 

higher_RAM_limit <- 3*10^6 #from which region width to assign more RAM
higher_RAM_symbolically <- "3" #what to write in the filename :)

#these are submit_plink2_with_parameters.sh inputs (be sure to set correct RAM limits in .sh file)
write_tsv(x = combined_EST_regions %>% filter(width < higher_RAM_limit), 
          file = paste0("combined_EST_regions_under_", higher_RAM_symbolically, "Mb.tsv"),
          col_names = F)
write_tsv(x = combined_EST_regions %>% filter(width >= higher_RAM_limit), 
          file = paste0("combined_EST_regions_over_", higher_RAM_symbolically, "Mb.tsv"),
          col_names = F)

#finds all metabolites that have leads in each region
get_all_metabolites_in_regions <- function(leads_df_with_metabolites, regions_df) {
  result <- data.frame()
  
  for (chr in unique(leads_df_with_metabolites$chromosome)) {
    chr_leads <- leads_df_with_metabolites %>% filter(chromosome == chr) %>% 
      select(-log10P) %>%
      select(-rsid) %>%
      distinct()
    
    chr_regions <- regions_df %>% filter(chromosome == chr)
    
    chr_joined <- left_join(chr_leads, chr_regions, by = "chromosome", relationship = "many-to-many") %>%
      filter(between(x = position, left = start, right = end)) %>%
      distinct() %>%
      group_by(start, end) %>%
      mutate(n_metabolites = n_distinct(metabolite)) %>%
      select(chromosome, start, end, width, metabolite, n_metabolites) %>%
      distinct()
    
    result <- rbind(result, chr_joined)
  }
  
  return(result %>% arrange(chromosome, start, end, metabolite))
}

#divides each regions metabolites into chunks of a maximum of chunk_size
divide_into_chunks <- function(metabolites_in_EST_regions_df, chunk_size){
  result <- metabolites_in_EST_regions %>%
    group_by(chromosome, start, end) %>%
    mutate(n_metabolites = n_distinct(metabolite)) %>%
    group_by(chromosome, start, end, n_metabolites) %>%
    mutate(chunk_id = (row_number() - 1) %/% chunk_size + 1) %>%
    group_by(chromosome, start, end, width, chunk_id) %>%
    summarise(
      metabolites_combined = paste(metabolite, collapse = ";"),
      n_metabolites_total_region = dplyr::first(n_metabolites),
      n_metabolites_in_chunk = n_distinct(metabolite),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    select(1:4, metabolites_combined, chunk_id, n_metabolites_in_chunk, n_metabolites_total_region)
  
  return(result)
}

metabolites_in_EST_regions <- get_all_metabolites_in_regions(leads_df_with_metabolites = EST_leads, 
                                                             regions_df = combined_EST_regions)

chunk_size = 25 # maximum number of metabolites in one chunk

metabolites_chunked <- divide_into_chunks(metabolites_in_EST_regions, 
                                          chunk_size)

#these are the input files for submit_susie_chunked_with_parameters.sh (be sure to set correct RAM limits in .sh file)
write_tsv(x = metabolites_chunked %>% filter(width < higher_RAM_limit), 
          file = paste0("combined_EST_regions_for_finemapping_under_", higher_RAM_symbolically, "Mb.tsv"),
          col_names = F)
write_tsv(x = metabolites_chunked %>% filter(width >= higher_RAM_limit), 
          file = paste0("combined_EST_regions_for_finemapping_over_", higher_RAM_symbolically, "Mb.tsv"),
          col_names = F)
