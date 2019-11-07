## Combine the linkage data with mutation type
## reads in the XX_SNP_mutation_types.tsv and XX.linkage.tsv from inStrain
## exports a linkage file with the mutation type of the linkage file included


## Library
library(tidyverse)


#### INPUT FILES ####

link.files <- list.files("inStrain/inStrain_linkage_output/")
snv.files <- list.files("inStrain/inStrain_gene_profile_output/", pattern= "SNP_mutation_types")

link.list <- map(link.files, function (x) read_tsv(file.path("inStrain", "inStrain_linkage_output", x))) %>% 
  setNames(str_replace(link.files, ".linkage.tsv", ""))
snv.list <- map(snv.files, function (x) read_tsv(file.path("inStrain", "inStrain_gene_profile_output", x))) %>% 
  setNames(str_replace(snv.files, "_SNP_mutation_types.tsv", ""))

## Check to make sure lists are in the same order
table(names(link.list) == names(snv.list))


#### COMBINE MUTATION TYPE AND LINKAGE TOGETHER ###
linkage_type <- map2(snv.list, link.list, function(SNV, LINK){
  snv_red <- SNV %>% 
    select(scaffold, gene, position, mutation_type)
  
  # Join on positions A
  NS_positions_A <- inner_join(snv_red, LINK, by= c("scaffold", "position" = "position_A")) %>% 
    rename(mutation_type_A= mutation_type,
           position_A= position,
           gene_A= gene) %>% 
    select(scaffold, gene_A, position_A, mutation_type_A) %>% 
    distinct()
  
  # Join on position B
  NS_positions_B <- inner_join(snv_red, LINK, by= c("scaffold", "position" = "position_B")) %>% 
    rename(mutation_type_B= mutation_type,
           position_B= position,
           gene_B= gene) %>% 
    select(scaffold, gene_B, position_B, mutation_type_B) %>% 
    distinct()
  
  # Combine mutation types on positions A and B with the original data
  linkage.NS <- left_join(LINK, NS_positions_A, by= c("scaffold", "position_A")) %>% 
    left_join(., NS_positions_B, by= c("scaffold", "position_B")) %>% 
    mutate(link_type= str_c(mutation_type_A, mutation_type_B, sep= "-"))
  
  return(linkage.NS)
}
) %>% setNames(str_replace(snv.files, "_SNP_mutation_types.tsv", ""))

rm(link.list, snv.list)


## Transform into a dataframe
linkage_type_df <- bind_rows(linkage_type, .id= "sample") %>% 
  mutate(species= str_split(sample, "\\.")[[1]][2]) %>% 
  select(sample, species, scaffold, everything())

write_tsv(linkage_type_df, path= "inStrain/linkage_NS.tsv")

