library(tidyverse)
library(ggplot2)



atx_cluster <- read_tsv("ATX_gene_cluster/atx_gene_cluster_features.tsv")


species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

link.sp2 <- read_tsv("inStrain/linkage_NS.tsv") %>% 
  mutate(site= str_replace(sample, "\\.species.*", "")) %>% 
  left_join(species_lookup, ., by= c("site", "species")) %>%  # Filter to only include species recovered from each site
  filter(species == "species_2") %>% 
  filter((mutation_type_A == "N" | mutation_type_A == "S") & (mutation_type_B == "N" | mutation_type_B == "S"))

link.atx.scaf <- link.sp2 %>% 
  filter(str_detect(scaffold, "scaffold_594c")) %>% 
  left_join(., atx_cluster, c("gene_A" = "gene"))

table(link.atx.scaf$sample, link.atx.scaf$haplotype)



#### N:S RATIOS
in_dir <- "inStrain/inStrain_gene_profile_output"
snv_files <- list.files(in_dir, pattern= ".pid96_SNP")
snv_list <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                  mutate(sample= str_replace(x, "_SNP_mutation_types.tsv", "")) %>% 
                  mutate(site= str_split(sample, "\\.")[[1]][1],
                         species= str_split(sample, "\\.")[[1]][2])) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))

# Filter to only include species recovered from each site
snv_atx_scaf <- bind_rows(snv_list) %>% left_join(species_lookup, ., by= c("site", "species"))  %>% 
  filter(species == "species_2") %>% 
  filter(str_detect(scaffold, "scaffold_594c"))



snv_atx_NS_ratios <- snv_atx_scaf %>% 
  filter(mutation_type == "N" | mutation_type == "S") %>% 
  group_by(sample, gene) %>% 
  summarize(
    n= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N")),
    S_count= sum(str_detect(mutation_type, "S")),
    NS= N_count/S_count) %>% 
  ungroup() %>% 
  filter(n > 1) %>% 
  inner_join(., atx_cluster, by= "gene")


#### SNV per mBP ####

snvs_mbp_df <- snv_atx_scaf %>% 
  group_by(gene) %>% 
  summarize(SNVs= length(mutation_type)) %>%
  ungroup() %>% 
  mutate(SNV_mbp= SNVs / genome_mbp)









#### MAKE FIGURES ####
source("Scripts/ggplot_themes.R")
ggplot(snv_atx_NS_ratios) +
  geom_hline(yintercept = 1, size= 0.5, linetype= "dashed") +
  geom_point(aes(x= name, y= NS, size= n), pch= 21, fill= "gray70", position= "jitter") +
  labs(x= "Anatoxin biosynthesis genes", y= "N:S ratio") +
  scale_size_continuous(name= "SNV count") +
  scale_y_continuous(expand= c(0.02, 0)) +
  theme_strains
ggsave(last_plot(), filename= "atx_NS.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)







