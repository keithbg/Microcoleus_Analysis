## Investigate the ATX biosynthesis genes

## Librarires
library(tidyverse)
library(ggplot2)


## Format data
atx_cluster <- read_tsv("ATX_gene_cluster/atx_gene_cluster_features.tsv")

species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

#### SNV LINKAGE ####

link.sp2 <- read_tsv("inStrain/linkage_NS.tsv") %>% 
  mutate(site= str_replace(sample, "\\.species.*", "")) %>% 
  left_join(species_lookup, ., by= c("site", "species")) %>%  # Filter to only include species recovered from each site
  filter(species == "species_2") %>% 
  filter((mutation_type_A == "N" | mutation_type_A == "S") & (mutation_type_B == "N" | mutation_type_B == "S"))

link.atx.scaf <- link.sp2 %>% 
  filter(str_detect(scaffold, "scaffold_594c")) %>% 
  left_join(., atx_cluster, c("gene_A" = "gene"))

table(link.atx.scaf$sample, link.atx.scaf$haplotype)
table(link.atx.scaf$name, link.atx.scaf$haplotype)


#### SNV per mBP ####
atx_gi_df <- read_tsv("inStrain/output_tables/gene_info_df.tsv") %>% 
  inner_join(., atx_cluster, by= "gene") %>% 
  mutate(SNV_mbp= SNPs_per_bp*1000000)

summary(atx_gi_df$SNV_mbp)
table(atx_gi_df$SNV_mbp > 5000)



#### N:S RATIOS ####

# Read SNV data and filter to only include species recovered from each site
snv_atx_scaf <- read_tsv("inStrain/output_tables/snp_mutation_type_df.tsv") %>% 
  left_join(species_lookup, ., by= c("site", "species"))  %>% 
  filter(str_detect(scaffold, "scaffold_594c")) # scaffold with ATX biosynthesis cluster

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
  inner_join(., atx_cluster, by= "gene") %>% 
  left_join(., atx_gi_df)


#### MAKE FIGURES ####
source("Scripts/ggplot_themes.R")


ggplot(atx_gi_df, aes(x= name, y= SNV_mbp)) +
  geom_boxplot() +
  geom_point(position= "jitter", color= "black", alpha= 0.5, size= 2) +
  labs(x= "Anatoxin biosynthesis genes", y= "SNVs per mbp") +
  scale_y_continuous(expand= c(0.02, 0)) +
  theme_strains
ggsave(last_plot(), filename= "atx_SNVs_mbp.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)



ggplot(snv_atx_NS_ratios) +
  geom_hline(yintercept = 1, size= 0.5, linetype= "dashed") +
  geom_point(aes(x= name, y= NS, size= n), pch= 21, fill= "gray70", position= "jitter") +
  labs(x= "Anatoxin biosynthesis genes", y= "N:S ratio") +
  scale_size_continuous(name= "SNV count") +
  scale_y_continuous(expand= c(0.02, 0)) +
  theme_strains
ggsave(last_plot(), filename= "atx_NS.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)


## SNVs per mBP and N:S Ratio
ggplot(snv_atx_NS_ratios, aes(x= SNV_mbp, y= NS)) +
  #geom_vline(xintercept = 1500, size= 0.5, linetype= "dashed") +
  geom_point(aes(color= name), size= 3) +
  labs(x= "SNVs per mbp", y= "N:S ratio") +
  scale_color_discrete(name= "Anatoxin \ngene") +
  scale_x_continuous(limits= c(0, 35000),
                     breaks= seq(0, 35000, by= 5000), 
                     labels= c("0", "", "10000", "", "20000", "", "30000", ""),
                     expand= c(0.01, 0)) +
  scale_y_continuous(expand= c(0.02, 0)) +
  theme_strains
ggsave(last_plot(), filename= "atx_NS_SNPmbp.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)









