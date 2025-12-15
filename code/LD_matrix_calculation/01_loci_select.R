library(dplyr)
library(readr)
library(GenomicRanges)
library(ggplot2)

#plotting aid, made by chatgpt
assign_tracks <- function(regions) {
  regions <- regions %>% arrange(start, end)
  tracks <- list()
  
  for (i in 1:nrow(regions)) {
    assigned <- FALSE
    for (track in seq_along(tracks)) {
      if (regions$start[i] > tracks[[track]]) {
        regions$track[i] <- track
        tracks[[track]] <- regions$end[i]
        assigned <- TRUE
        break
      }
    }
    if (!assigned) {
      # new track
      tracks[[length(tracks) + 1]] <- regions$end[i]
      regions$track[i] <- length(tracks)
    }
  }
  return(regions)
}

#plotting function, mostly made by chatgpt
plot_regions <- function(regions_df, leads_df, too_near_distance){
  regions_with_tracks <- regions_df %>%
    mutate(track = NA_integer_) %>%
    group_by(chromosome) %>%
    mutate(n_regions_in_chr = n()) %>%
    group_modify(~assign_tracks(.x)) %>%
    ungroup() %>%
    dplyr::rename(region_track = track) %>%
    mutate(label = paste0("Chromosome ", chromosome, " (", n_regions_in_chr, " regions)")) %>%
    arrange(chromosome, start, end) %>%
    mutate(label = factor(label, levels = unique(label[order(chromosome)])))
  
  leads_in_regions <- get_all_metabolites_in_regions(leads_df, regions_df) %>%
    select(-metabolite) %>%
    distinct() %>%
    left_join(regions_with_tracks) %>%
    mutate(distance_from_start = lead_position - start,
           distance_to_end = end - lead_position,
           distance_from_closest_edge = pmin(distance_from_start, distance_to_end),
           too_near = (start >= too_near_distance & distance_from_start < too_near_distance) | distance_to_end < too_near_distance,
           point_track = 1.5,
           point_color = ifelse(too_near, "red", "blue"),
           point_alpha = ifelse(too_near, 1, 0.4)) %>%
    group_by(chromosome, lead_position) %>%
    filter(distance_from_closest_edge == max(distance_from_closest_edge)) %>%
    ungroup() %>%
    arrange(chromosome, lead_position)
  
  
  p <- ggplot() +
    geom_rect(data = regions_with_tracks, mapping = aes(xmin = start, xmax = end, ymin = region_track - 0.4, ymax = region_track + 0.4), fill = "grey", color = "black") +
    geom_point(data = leads_in_regions, mapping = aes(x=lead_position, y=point_track, fill = point_color, colour = point_color, alpha = point_alpha)) +
    facet_wrap(~label, scales = "fixed", ncol = 1) +
    labs(x = "Genomic Position", y = paste0(nrow(regions_with_tracks), " total regions")) +
    theme_minimal() +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_alpha_identity() +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      panel.spacing = unit(0.5, "lines"),
      axis.text.y = element_blank(),        # Remove y-axis tick labels
      axis.ticks.y = element_blank(),       # Remove y-axis ticks
      panel.grid.major.y = element_blank(), # Remove horizontal grid lines
      panel.grid.minor.y = element_blank()  # Just in case minor grid lines exist
    ) +
    guides(alpha = "none")
  
  print(p)
}

#for each region, determine which metabolites have leads in it
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
      select(chromosome, start, end, width, metabolite, n_metabolites, lead_position = position) %>% #NB! lisasin positioni juurde
      distinct()
    
    result <- rbind(result, chr_joined)
  }
  
  return(result %>% arrange(chromosome, start, end, metabolite) %>% ungroup())
}

#create ranges of (lead-window...lead+window) for each lead and merge overlapping ones into a single region (can produce very wide regions)
get_grange_regions <- function(unique_leads_df, window = 1e6){
  df <- unique_leads_df %>%
    mutate(start = pmax(1, position-window)) %>%
    mutate(end = (position + window)) %>%
    arrange(chromosome, position)
  
  gr <- GRanges(seqnames = df$chromosome, ranges = IRanges(df$start, df$end))
  rd <- reduce(gr) 
  
  return(rd)
}

#chop regions that are larger than keep_width_limit into chunks that are split_width wide and overlap by split_width - step
split_regions <- function(gr, keep_width_limit = 3.5e6, split_width = 3e6, step = 1e6) {
  regions_to_keep <- gr[width(gr) <= keep_width_limit]
  regions_to_split <- gr[width(gr) > keep_width_limit]
  
  split <- slidingWindows(regions_to_split, width = split_width, step = 1e6) %>% unlist()
  
  return(c(regions_to_keep, split))
}

#convert granges object to dataframe
granges_to_df <- function(gr){
  lvls <- levels(gr@seqnames)
  
  result <- as.data.frame(gr) %>%
    select(1:4) %>%
    mutate(chromosome = as.numeric(lvls[seqnames])) %>%
    select(chromosome, start, end, width) %>%
    arrange(chromosome, start) %>%
    mutate(width = end-start + 1)
  
  return(result)
}

#if the last region in a group is under tail_size_limit wide, merge it with the second to last
merge_little_tails <- function(gr_df, tail_size_limit = 2.5e6){
  df_overlaps <- gr_df %>%
    arrange(chromosome, start) %>%
    group_by(chromosome) %>%
    mutate(next_end = lead(end, n = 1)) %>%
    mutate(next_width = lead(width, n = 1)) %>%
    mutate(next_start = lead(start, n = 1)) %>%
    ungroup() %>%
    mutate(overlaps_with_next = next_start < end) %>%
    mutate(merge_with_next = (next_width < tail_size_limit) & overlaps_with_next) %>%
    mutate(merge_with_next = ifelse(is.na(merge_with_next), FALSE, merge_with_next)) %>%
    mutate(end = ifelse(merge_with_next, next_end, end)) %>%
    mutate(remove_this = lag(merge_with_next, 1)) %>%
    mutate(remove_this = ifelse(is.na(remove_this), FALSE, remove_this)) %>%
    filter(!remove_this) %>%
    mutate(width = end-start+1) %>%
    select(chromosome, start, end, width) 
  
  return(df_overlaps)
}

remove_regions_with_no_leads <- function(gr_df, unique_leads_df, distance_allowed=1e6, MHC_start = 28510120, MHC_end = 33480577){
  #take out MHC
  unique_leads_df <- unique_leads_df %>%
    filter(!(chromosome == 6 & between(position, MHC_start, MHC_end))) 
  
  leads_in_regions <- unique_leads_df %>%
    left_join(gr_df, relationship = "many-to-many") %>%
    filter(between(position, start, end)) %>%
    mutate(distance_to_start = ifelse(start == 1, distance_allowed + 1, position-start)) %>% #if near beginning of chromosome, ignore
    mutate(distance_to_start = ifelse((chromosome == 6 & start == MHC_end + 1), distance_allowed + 1, distance_to_start)) %>% #if near end of MHC, ignore
    mutate(distance_from_end = end-position) %>%
    mutate(distance_from_end = ifelse((chromosome == 6 & end == MHC_start -1 ), distance_allowed + 1, distance_from_end)) %>% #if near beginning of MHC, ignore
    mutate(distance_to_edge = pmin(distance_from_end, distance_to_start)) %>%
    mutate(variant_fits_in_this_region = distance_to_edge >= distance_allowed) %>%
    filter(variant_fits_in_this_region) %>%
    group_by(chromosome, position) %>%
    mutate(n_suitable_regions = sum(variant_fits_in_this_region)) %>%
    ungroup() %>%
    group_by(chromosome, start, end) %>%
    mutate(all_variants_in_region_fit_in_multiple = all(n_suitable_regions > 1)) %>%
    ungroup() %>%
    mutate(region_name = paste(chromosome, start, end, sep = "_"))
  
  if(nrow(leads_in_regions) != nrow(unique_leads_df)){
    message("check rows")
    return(NULL)
  }
  if(!all(leads_in_regions$variant_fits_in_this_region)){
    message("check other")
    return(NULL)
  }
  
  df_empty_removed <- gr_df %>%
    mutate(region_name = paste(chromosome, start, end, sep = "_")) %>%
    filter(region_name %in% leads_in_regions$region_name) %>%
    select(-region_name)
  
  return(df_empty_removed)
}

#clip the regions overlapping MHC
clip_MHC <- function(gr, MHC_start = 28510120, MHC_end = 33480577){
  MHC_region <- GRanges(seqnames = "6", ranges = IRanges(MHC_start, MHC_end))
  
  return(setdiff(gr, MHC_region))
}

#check if all leads (except those in MHC) are in a region
are_all_leads_in_a_region <- function(leads_df, region_df, MHC_start = 28510120, MHC_end = 33480577) {
  leads_to_check <- leads_df %>%
    distinct(chromosome, position) %>%
    filter(!(chromosome == 6 & between(position, MHC_start, MHC_end)))
  
  checking_df <- leads_to_check %>%
    left_join(region_df, by = "chromosome", relationship = "many-to-many") %>%
    filter(between(position, start, end))
  
  leads_missing <- anti_join(leads_to_check, checking_df, by = c("chromosome", "position"))
  
  return(nrow(leads_missing) == 0)
}


meta_EUR_leads <- read_tsv("...") %>%
  select(metabolite, rsid, chromosome = CHR, position = POS, ALL0, ALL1, log10P = LOG10P) 

EUR_leads <- read_tsv("...") %>%
  select(metabolite, rsid, chromosome = CHR, position = POS, ALL0, ALL1, log10P = LOG10P) 

combined_EUR_leads <- rbind(meta_EUR_leads, EUR_leads) %>%
  arrange(chromosome, position, ALL0, ALL1)

combined_EUR_unique_leads <- rbind(meta_EUR_leads, EUR_leads) %>%
  distinct(chromosome, position) 

combined_EUR_regions <- get_grange_regions(combined_EUR_unique_leads, window = 1e6) %>% 
  clip_MHC() %>%
  split_regions(keep_width_limit = 3.5e6, split_width = 3e6, step = 2e6) %>% 
  granges_to_df() %>%
  merge_little_tails(tail_size_limit = 2.5e6) %>%
  remove_regions_with_no_leads(unique_leads_df = combined_EUR_unique_leads) 


#check if all variants are in a region somewhere
are_all_leads_in_a_region(leads_df = combined_EUR_leads, region_df = combined_EUR_regions)

plot_regions(regions_df = combined_EUR_regions,
             leads_df = combined_EUR_leads,
             too_near_distance = 1e6)

write_tsv(x = combined_EUR_regions, file = "meta_EUR_loci_DNAnexus.tsv", col_names = F)
