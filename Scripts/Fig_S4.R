## Figure S4

## Nucleotide diversity (pi) density distributions for Microcoleus species

## Investigate nucleotide diversity (pi) results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.gene_profile.tsv


#### Libraries #################################################################
library(tidyverse)
source("Scripts/ggplot_themes.R")
################################################################################


#### INPUT FILES
## Input df (gene_info_filt_df_v1.4.tsv) generated in format_inStrain_output.R
gi_filt_df <- read_tsv("Data/inStrain_data/gene_info_filt_df_v1.4.tsv")



#### SUMMARIZE PI VALUES ACROSS THE GENOME ####
gi_filt_summary <- gi_filt_df %>% 
  group_by(sample, site, species) %>% 
  summarize(n= length(pi),
            mean_pi= mean(pi, na.rm= TRUE),
            median_pi= median(pi, na.rm= TRUE),
            sd_pi= sd(pi, na.rm= TRUE),
            min_pi= min(pi, na.rm= TRUE),
            max_pi= max(pi, na.rm= TRUE),
            median_cov= median(coverage, na.rm= TRUE)) %>% 
  ungroup() 

## Summary values
summary(gi_filt_summary$mean_pi)
summary(gi_filt_summary$mean_pi2)

hist(gi_filt_summary$median_pi2 - gi_filt_summary$mean_pi2)
hist(gi_filt_summary$mean_pi2)
hist(gi_filt_summary$median_pi2)


#### MAKE FIGURE ####
gi_filt_df$species_facet <- factor(gi_filt_df$species,
                                   labels= c(expression(italic("Microcoleus")~"sp. 1"),
                                             expression(italic("Microcoleus")~"sp. 2"),
                                             expression(italic("Microcoleus")~"sp. 3")))


## Nucleotide diversity density curves
ggplot(data= gi_filt_df) +
  geom_density(aes(x= nucl_diversity, color= sample)) +
  labs(x= "Nucleotide diversity", y= "Density") +
  scale_x_log10(limits= c(0.00002, 0.52), # the minimum pi value (apart from 0) is 0.000026
                breaks= c(0.0001, 0.001, 0.01, 0.1),
                labels= c("0.0001", "0.001", "0.01", "0.1"),
                expand= c(0, 0)) +
  annotation_logticks(sides= "b") +
  scale_y_continuous(expand= c(0.01, 0)) +
  scale_color_discrete(guide= FALSE) +
  facet_rep_grid(species_facet~., scales= "free_y", labeller= label_parsed) +
  theme_strains

ggsave(last_plot(), filename = "Fig_S4_v1.4.png", dpi= 320, height= 180*0.75, width= 180, units= "mm",
       path= "Output_figures")



