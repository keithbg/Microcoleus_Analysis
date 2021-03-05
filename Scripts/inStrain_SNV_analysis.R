## Investigate the SNVs.tsv file output
## AC= allele count
library(tidyverse)
library(vegan)


#### INPUT FILES ####
in_dir <- "Data/inStrain_data/inStrain_output"

species_lookup <- read_tsv("Data/inStrain_data/inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

snv_files <- list.files(in_dir, pattern= "species_1.pid96.SNVs.tsv")


## Read in all SNVs.tsv files
# AC = allele count
snv_list_AC1 <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                           mutate(pos_id= str_c(scaffold, position, sep= "-")) %>% 
                           filter(allele_count == 1) %>% 
                           select(pos_id, con_base)) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))

snv_list_AC2 <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                         mutate(pos_id= str_c(scaffold, position, sep= "-")) %>% 
                         filter(allele_count == 2) %>% 
                         select(pos_id, con_base)) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))


## Transform into data frame
snv_df_conBase_AC1 <- bind_rows(snv_list_AC1, .id= "sample") %>% 
  mutate(site= do.call(rbind, str_split(.$sample, "\\."))[, 1],
         species= do.call(rbind, str_split(.$sample, "\\."))[, 2]) %>% 
  inner_join(species_lookup, ., by= c("site", "species")) %>% # Filter to only include species recovered from each site
  mutate(con_base= ifelse(con_base == "A", 1,
                         ifelse(con_base == "T", 2,
                                ifelse(con_base == "C", 3, 4))))


snv_df_conBase_AC2 <- bind_rows(snv_list_AC2, .id= "sample") %>% 
  mutate(site= do.call(rbind, str_split(.$sample, "\\."))[, 1],
         species= do.call(rbind, str_split(.$sample, "\\."))[, 2]) %>% 
  inner_join(species_lookup, ., by= c("site", "species")) %>% # Filter to only include species recovered from each site
  mutate(con_base= ifelse(con_base == "A", 1,
                         ifelse(con_base == "T", 2,
                                ifelse(con_base == "C", 3, 4))))


## Make wide df of SNV position presence/absence 0 and 1 values
snv_binary_AC1 <- snv_df_conBase_AC1 %>% 
  select(site, pos_id, con_base) %>% 
  mutate(con_base= ifelse(con_base > 0, 1, 0)) %>% 
  spread(key= pos_id, value= con_base, fill= 0)

snv_binary_AC2 <- snv_df_conBase_AC2 %>% 
  select(site, pos_id, con_base) %>% 
  mutate(con_base= ifelse(con_base > 0, 1, 0)) %>% 
  spread(key= pos_id, value= con_base, fill= 0)


## SNV sharing (Jaccard distance, i.e. % shared SNV sites between pairwise comparisons)
  ## Calculate distance between sites
  dist.rows.AC1 <- vegan::vegdist(as.matrix(snv_binary_AC1[, -1]), method= "jaccard")
  dist.rows.AC2 <- vegan::vegdist(as.matrix(snv_binary_AC2[, -1]), method= "jaccard")

  snv_dist_mat.AC1 <- as.matrix(dist.rows.AC1)
  snv_dist_mat.AC2 <- as.matrix(dist.rows.AC2)
  
  colnames(snv_dist_mat.AC1) <- snv_binary_AC1$site
  rownames(snv_dist_mat.AC1) <- snv_binary_AC1$site

  colnames(snv_dist_mat.AC2) <- snv_binary_AC2$site
  rownames(snv_dist_mat.AC2) <- snv_binary_AC2$site

## Write files
write_tsv(as.data.frame(snv_dist_mat.AC1), path= "Data/inStrain_data/snv_pos_AC1_jaccard.tsv")
write_tsv(as.data.frame(snv_dist_mat.AC2), path= "Data/inStrain_data/snv_pos_AC2_jaccard.tsv")



