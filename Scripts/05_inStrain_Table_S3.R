## Script to generate Supplementary Table 3 (Table S3). 
## Each row is the output of the inStrain genome profile analysis for a single genome
# https://instrain.readthedocs.io/en/latest/user_manual.html#profile

library(tidyverse)
in_dir <- "Data/inStrain_data/inStrain_output"
genome_info_files <- list.files(in_dir, pattern= "pid96_genome_info")


gi_list <- map(genome_info_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                 mutate(species= ifelse(str_detect(genome, "12U"), "species_1",
                                        ifelse(str_detect(genome, "13D"), "species_2", "species_3")))) %>% 
                 # mutate(sample= str_replace(x, "_gene_info.tsv", ""),
                 #        site= str_split(sample, "\\.")[[1]][1],
                 #        species= str_split(sample, "\\.")[[1]][2])) %>% 
  setNames(str_replace(genome_info_files, ".species.*.tsv", ""))

gi_df <- bind_rows(gi_list, .id= "sample") %>% 
  select(sample, species, everything())

write_tsv(gi_df, path= file.path("Data/inStrain_data", "Table_S3.tsv"))
