# Investigate nonsynonymous and synonymous SNV results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.SNP_mutation_types.tsv

### The SNP_mutation_types.tsv file filters SNVs from the SNVs.tsv file by morphia == 2 and cryptic == FALSE. 
## This is why there are fewer SNVs in the SNP_mutation_types.tsv file compared to the SNVs.tsv file

## Input df generated in format_inStrain_output.R

## Libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
source("Scripts_manuscript/ggplot_themes.R")

setwd("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis")

## SNV data (input table generated in: format_inStrain_output.R)
snv_df <- read_tsv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/Strains/inStrain/output_tables/snp_mutation_type_df.tsv")

# Get genome sizes
genome.size <- read_delim("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/Strains/genome_lengths.txt", delim= " ", col_names= FALSE) %>% 
  slice(c(13:16, 21:22))

genome_size_df <- tibble(species= c("species_1", "species_2", "species_3"),
                         genome_mbp= c(as.numeric(genome.size[2, ]) / 1000000, 
                                       as.numeric(genome.size[4, ]) / 1000000, 
                                       as.numeric(genome.size[6, ]) / 1000000))

## CALCULATE SNVs per MBP
snvs_mbp_df <- snv_df %>% 
  left_join(., genome_size_df) %>% 
  group_by(sample, species, genome_mbp) %>% 
  summarize(SNVs= length(mutation_type),
            cov_mean= mean(baseCoverage, na.rm= TRUE),
            cov_sd= sd(baseCoverage, na.rm= TRUE)) %>%
  ungroup() %>% 
  mutate(SNV_mbp= SNVs / genome_mbp)

## CALCULATE N:S RATIOS per genome
NS_genome_ratios <- snv_df %>% 
  group_by(sample, species) %>% 
  summarize(
    SNV_count= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N")),
    S_count= sum(str_detect(mutation_type, "S")),
    I_count= sum(str_detect(mutation_type, "I")),
    M_count= sum(str_detect(mutation_type, "M")),
    NS= N_count/S_count) 

snvs_genome_df <- left_join(snvs_mbp_df, NS_genome_ratios, by= c("sample", "species")) %>% 
  mutate(pop_age= ifelse(SNV_mbp < 1550, "Young", "Old"),
         pop_age_binary= ifelse(SNV_mbp < 1550, 1, 0),
         ggkbase_id= str_replace(sample, "\\.species.*$", ""))
  #left_join(., latlong, by= "ggkbase_id") %>% 
  #left_join(., watershed.area, by= "ggkbase_id") %>% 
 # rename(species= species.x)

## CALCULATE MINOR ALLELE FREQUENCIES FOR SPECIES 1
snv.freq.sp1 <- snv_df %>% 
  filter(species == "species_1") %>% 
  #filter(., site == "PH2015_03U") %>% 
  mutate(varFreq_r3 = round(varFreq, 3),
         varFreq_r2 = round(varFreq, 2)) %>% 
  group_by(site, varFreq_r2) %>% 
  summarize(n= length(varFreq)) %>% 
  ungroup() %>% 
  left_join(., select(snvs_genome_df, ggkbase_id, pop_age, NS, SNV_mbp), by= c("site" = "ggkbase_id")) %>% 
  mutate(facet_label= str_sub(site, start= 3, end= 13))

## Samples with secondary SNV peaks
secondary.snv.peaks <- c("2015_03D", "2015_03U", "2015_04D", "2015_04U", "2015_10S", "2017_02_FOX", 
                         "2017_03_ELD", "2017_04_SCI", "2017_05_CCC", "2017_06_SFM", "2017_07_MST")



#### FIGURES ####



secondary.peaks <- ggplot(filter(snv.freq.sp1, facet_label %in% secondary.snv.peaks),  aes(x= varFreq_r2, y= n)) +
  geom_point(color= species.colors[1], size= 1, alpha= 0.7) +
  #geom_point(aes(color= NS), size= 2) +
  #scale_color_gradientn(colors= colorRampPalette(c("snow2", species.colors[1]))(2),
  #                       name= "N:S") +
  geom_smooth(method= "gam", se= FALSE, size= 0.75, color= "black") +
  scale_x_continuous(limits= c(0.05, 0.5), 
                     breaks= seq(0.05, 0.5, by= 0.05),
                     labels= c("", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5"),
                     expand= c(0.02,0)) +
  #labs(x= "Minor allele frequency", y= "Count") +
  lemon::facet_rep_wrap(~facet_label, ncol= 6, scales= "free_y") +
  theme_strains +
  theme(text = element_text(size= 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t= 1, r=0.25, b= 0.5, l= 0.5, "cm"))


no.secondary.peaks <- ggplot(filter(snv.freq.sp1, !(facet_label %in% secondary.snv.peaks) & (facet_label != "2015_01D") & (facet_label != "2015_01U")),  aes(x= varFreq_r2, y= n)) +
  geom_point(color= species.colors[1], size= 1, alpha= 0.7) +
  #geom_point(aes(color= NS), size= 2) +
  #scale_color_gradientn(colors= colorRampPalette(c("snow2", species.colors[1]))(2),
  #                       name= "N:S") +
  geom_smooth(method= "gam", se= FALSE, size= 0.75, color= "black") +
  scale_x_continuous(limits= c(0.05, 0.5), 
                     breaks= seq(0.05, 0.5, by= 0.05),
                     labels= c("", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5"),
                     expand= c(0.02,0)) +
  #labs(x= "Minor allele frequency", y= "Count") +
  lemon::facet_rep_wrap(~facet_label, ncol= 6, scales= "free_y") +
  theme_strains +
  theme(text = element_text(size= 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t= 1, r=0.25, b= 0.5, l= 0.5, "cm"))

snv.freq.fig <- plot_grid(secondary.peaks, no.secondary.peaks,
                          labels= c("Secondary SNV peaks", "No secondary SNV peaks"),
                          # align= "v",
                          nrow= 2,
                          rel_heights = c(0.5, 1),
                          hjust= -0.2)

snv.freq.fig.anno <- annotate_figure(snv.freq.fig,
                                     left= text_grob("Number of SNV sites", rot= 90, vjust= 2),
                                     bottom= text_grob("Minor allele frequency", vjust= -1))
snv.freq.fig.anno
#ggsave(snv.freq.fig, filename = "snv_freq_sp1.pdf", height= 180, width= 180, units= "mm", device= cairo_pdf,
#      path= "Output_figures")
ggsave(snv.freq.fig.anno, filename = "Fig_S8.png", height= 180, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


