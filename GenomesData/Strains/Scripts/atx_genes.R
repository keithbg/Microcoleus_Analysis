## Investigate the ATX biosynthesis genes

## Librarires
library(tidyverse)
library(ggplot2)


## Format data
atx_cluster <- read_tsv("ATX_gene_cluster/atx_gene_cluster_features.tsv")
pks_domains <- read_tsv("ATX_gene_cluster/atx_pks_domains.tsv") %>% 
  mutate(domainID= str_c(domain, domain_start, sep= "-"))
  

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


test <- read_tsv("inStrain/output_tables/gene_info_df.tsv") %>% 
  filter(str_detect(scaffold, "594c"))
#### N:S RATIOS ####

# Read SNV data and filter to only include species recovered from each site
snv_atx_scaf <- read_tsv("inStrain/output_tables/snp_mutation_type_df.tsv") %>% 
  left_join(species_lookup, ., by= c("site", "species"))  %>% 
  filter(str_detect(scaffold, "scaffold_594c"))  # scaffold with ATX biosynthesis cluster

snv_atx <- read_tsv("inStrain/output_tables/snp_mutation_type_df.tsv") %>% 
  left_join(species_lookup, ., by= c("site", "species"))  %>% 
  #filter(str_detect(scaffold, "scaffold_594c"))  %>% 
  filter(str_detect(scaffold, "scaffold_594c")) %>% 
  inner_join(., atx_cluster, by= "gene")

snv_pks <- snv_atx %>% 
  filter(name == "anaE" | name == "anaF" | name == "anaG") %>% 
  left_join(., pks_domains)#, by= c("name", "gene"))

  
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
  left_join(., atx_gi_df) %>% 
  filter(name != "anaJ")


#### SNV sharing ####
snv_postions_counts <- snv_atx_scaf %>% 
  inner_join(., atx_cluster, by= "gene") %>% 
  count(name, gene, position) %>% 
  inner_join(., atx_cluster, by= c("name", "gene")) %>% 
  filter(name != "anaJ")
  

pks_snv_counts <- snv_atx_scaf %>% 
  inner_join(., pks_domains, by= "gene") %>% 
  filter(name == "anaE" | name == "anaF" | name == "anaG") %>% 
  count(name, gene, position) %>% 
  inner_join(., atx_cluster, by= c("name", "gene")) %>% 
  filter(name != "anaJ")


ggplot(data= snv_postions_counts, aes(x= name, y= n)) +
  geom_point(position= "jitter")





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
  geom_point(aes(color= name, size= coverage)) +
  labs(x= "SNVs per mbp", y= "N:S ratio") +
  scale_color_discrete(name= "Anatoxin \ngene") +
  scale_x_continuous(limits= c(0, 12500),
                     breaks= seq(0, 12500, by= 2500),
                     labels= c("0", "", "5000", "", "10000", ""),
                     expand= c(0.01, 0)) +
  scale_y_continuous(expand= c(0.02, 0)) +
  theme_strains
ggsave(last_plot(), filename= "atx_NS_SNPmbp.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)



## SNV sharing

ggplot(data= snv_postions_counts, aes(x= n)) +
  geom_histogram(binwidth = 1, fill= "black", color= "gray50") +
  labs(x= "Number of genomes sharing a SNV position", y= "Count") +
  scale_x_continuous(breaks = seq(1, 6),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  facet_grid(name ~.) +
  theme_strains
ggsave(last_plot(), filename= "atx_SNV_sharing.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)



line.end <- "round"
y.start= 0.1
y.end= 0.1
line.size= 4

anaB <- geom_segment(aes(x= 11267, xend= 12412, y= y.start, yend= y.end, color= "anaB"), lineend= line.end, size= line.size)
anaC <- geom_segment(aes(x= 12425, xend= 14044, y= y.start, yend= y.end, color= "anaC"), lineend= line.end, size= line.size)
anaD <- geom_segment(aes(x= 14059, xend= 14322, y= y.start, yend= y.end, color= "anaD"), lineend= line.end, size= line.size)
anaE <- geom_segment(aes(x= 14475, xend= 20360, y= y.start, yend= y.end, color= "anaE"), lineend= line.end, size= line.size)
anaF <- geom_segment(aes(x= 20882, xend= 26423, y= y.start, yend= y.end, color= "anaF"), lineend= line.end, size= line.size)
anaG <- geom_segment(aes(x= 26552, xend= 29771, y= y.start, yend= y.end, color= "anaG"), lineend= line.end, size= line.size)
anaA <- geom_segment(aes(x= 34957, xend= 35718, y= y.start, yend= y.end, color= "anaA"), lineend= line.end, size= line.size)
anaJ <- geom_segment(aes(x= 34086, xend= 34295, y= y.start, yend= y.end, color= "anaJ"), lineend= line.end, size= line.size)
anaI <- geom_segment(aes(x= 36094, xend= 37482, y= y.start, yend= y.end, color= "anaI"), lineend= line.end, size= line.size)

ggplot(data= snv_atx_scaf) +
  geom_histogram(aes(x= position),
                 binwidth = 10,
                 boundary= 1) +
  anaB +
  anaC +
  anaD +
  anaE +
  anaF +
  anaG +
  anaA +
  anaJ +
  anaI +
  scale_x_continuous(limits= c(11267, 37482)) +
  scale_y_continuous(expand= c(0.0, 0)) +
  scale_color_discrete(name= "Anatoxin \ngene") +
  theme_strains


ggplot(data= snv_atx_scaf) +
  geom_bar(aes(x= position)) +
  anaB +
  anaC +
  anaD +
  anaE +
  anaF +
  anaG +
  anaA +
  anaJ +
  anaI +
  scale_x_continuous(limits= c(11267, 37482)) +
  scale_color_discrete(name= "Anatoxin \ngene") +
  theme_strains

ggplot(data= snv_atx) +
  geom_histogram(aes(x= position),
                 binwidth = 5,
                 boundary= 1,
                 fill= "black") +
 # geom_bar(aes(x= position), fill= "black") +
  geom_vline(xintercept = c(34262, 34296)) +
  geom_point(aes(x= gene_start, y= 0), color= "transparent") +
  geom_point(aes(x= gene_end, y= 0), color= "transparent") +
  geom_line(data= pks_domains_plots, aes(x= domain_position, y= 0.5, group= domainID, color= domain), size= 3) +
  labs(x= "Contig position", y= "SNV count") +
  scale_color_discrete(name= "PKS \ndomain") +
  scale_y_continuous(breaks= seq(0, 10, by= 2), expand= c(0, 0)) +
  scale_x_continuous(expand= c(0, 0)) +
  facet_wrap(.~name, ncol= 2, scales= "free_x") +
  theme_strains +
  theme(legend.position = "top")


### PKS Domains
snv_pks_plots <- snv_pks %>% 
  mutate(domainID= str_c(domain, domain_start, sep= "-")) %>% 
  gather(key= domains_start_end, value= domain_position, domain_start:domain_end)

pks_domains_plots <- pks_domains %>% 
  gather(key= domains_start_end, value= domain_position, domain_start:domain_end)


ggplot(data= snv_pks) +
  geom_histogram(aes(x= position),
                 binwidth = 3,
                 boundary= 1,
                 fill= "black") +
  geom_line(data= pks_domains_plots, aes(x= domain_position, y= 0.5, group= domainID, color= domain), size= 3) +
  # geom_bar(aes(x= position), fill= "black") +
  labs(x= "Contig position", y= "SNV count") +
  scale_color_discrete(name= "PKS \ndomain") +
  scale_y_continuous(expand= c(0, 0)) +
  scale_x_continuous(expand= c(0, 0)) +
  facet_wrap(.~name, ncol= 1, scales= "free_x") +
  theme_strains 
  




