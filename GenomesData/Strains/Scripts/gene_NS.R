# Investigate nonsynonymous and synonymous SNV results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.SNP_mutation_types.tsv


## Libraries
library(tidyverse)
library(ggplot2)


#### INPUT FILES ####
in_dir <- "inStrain/inStrain_gene_profile_output"
species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)
gg_anno <- read_tsv("inStrain/ggkbase_anno.tsv") # ggkbase annotations


snv_files <- list.files(in_dir, pattern= ".pid96_SNP")
snv_list <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                  mutate(sample= str_replace(x, "_SNP_mutation_types.tsv", "")) %>% 
                  mutate(site= str_split(sample, "\\.")[[1]][1],
                         species= str_split(sample, "\\.")[[1]][2]))
names(snv_list) <- str_replace(snv_files, ".tsv", "")

# Filter to only include species recovered from each site
snv_df <- as_tibble(do.call(rbind, snv_list)) %>% left_join(species_lookup, ., by= c("site", "species"))  

# Get genome sizes
genome.size <- read_delim("genome_lengths.txt", delim= " ", col_names= FALSE) %>% 
  slice(c(13:16, 21:22))

genome_size_df <- tibble(species= c("species_1", "species_2", "species_3"),
                         genome_mbp= c(as.numeric(genome.size[2, ]) / 1000000, 
                                       as.numeric(genome.size[4, ]) / 1000000, 
                                       as.numeric(genome.size[6, ]) / 1000000))

#### CALCULATE SNVs per MBP ####
snvs_mbp_df <- snv_df %>% 
  left_join(., genome_size_df) %>% 
  group_by(sample, species, genome_mbp) %>% 
  summarize(SNVs= length(mutation_type)) %>%
  ungroup() %>% 
  mutate(SNV_mbp= SNVs / genome_mbp)


#### CALCULATE N:S RATIOS ####
# for both genome and individual genes
NS_genome_ratios <- snv_df %>% 
  group_by(sample, species) %>% 
  summarize(
    SNV_count= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N")),
    S_count= sum(str_detect(mutation_type, "S")),
    I_count= sum(str_detect(mutation_type, "I")),
    M_count= sum(str_detect(mutation_type, "M")),
    NS= N_count/S_count) 


NS_gene_ratios <- snv_df %>% 
  filter(mutation_type == "N" & mutation_type == "S") %>% 
  group_by(sample, gene) %>% 
  summarize(
    n= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N")),
    S_count= sum(str_detect(mutation_type, "S")),
    NS= N_count/S_count) 




## JOIN THE SNVS PER MBP AND NS RATIO DATA
snvs_genome_df <- left_join(snvs_mbp_df, NS_genome_ratios, by= c("sample", "species"))



#### MAKE FIGURES ####
source("Scripts/ggplot_themes.R")

ggplot(data= snvs_genome_df, aes(x= SNV_mbp, y= NS)) +
  geom_vline(xintercept= 1500, color= "gray50", size= 0.25) +
  geom_point(aes(fill= species, shape= species), size= 3) +
  labs(x= "SNVs per mbp", y= "N:S ratio") +
  scale_x_continuous(limits= c(0, 8500),
                     breaks= seq(0, 8000, by= 1000), 
                     labels= c("0", "", "2000", "", "4000", "", "6000", "", "8000"),
                     expand= c(0, 0)) +
  scale_fill_manual(values= species.colors,
                     labels= c("1", "2", "3"),
                     name= "Species") +
  
  scale_shape_manual(values= species.shapes,
                     labels= c("1", "2", "3"),
                     name= "Species") +
  theme_strains +
  theme(legend.position = c(0.92, 0.85))
ggsave(last_plot(), filename = "snv_mbp.pdf", height= 180*0.75, width= 180, units= "mm", device= cairo_pdf,
       path= "Output_figures")

