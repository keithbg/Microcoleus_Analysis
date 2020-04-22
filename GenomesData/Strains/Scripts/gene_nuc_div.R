## Investigate nucleotide diversity (pi) results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.gene_profile.tsv

## Input df generated in format_inStrain_output.R


## Libraries
library(tidyverse)
library(ggplot2)


#### INPUT FILES
gi_filt_df <- read_tsv("inStrain/output_tables/gene_info_filt_df.tsv")


## Watershed area data
dir_input_watershed <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/EnvData"
watershed.area <-
  read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2) %>% 
  rename(site= ggkbase_id)




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
  ungroup() %>% 
  left_join(., watershed.area, by= "site") # COMBINE WITH WATERSHED AREA

#write_tsv(gi_filt_summary, "Output_tables/nuc_div_summary.txt")

summary(gi_filt_summary$mean_pi)
summary(gi_filt_summary$mean_pi2)

hist(gi_filt_summary$median_pi2 - gi_filt_summary$mean_pi2)
hist(gi_filt_summary$mean_pi2)
hist(gi_filt_summary$median_pi2)


#### MAKE FIGURES ####

## ggplot themes
source("Scripts/ggplot_themes.R")


## MEDIAN NUC. DIVERSITY X  WATERSHED AREA


nucDiv.watershed.plot <- ggplot(data= gi_filt_summary, aes(x= watershed_km2, y= median_pi)) +
  geom_point(aes(fill= species, shape= species), size= 4, color= "gray40") +
  labs(x= expression(paste("Watershed area (", km^{2}, ")")), y= "Median nucleotide diversity") +
  #  geom_smooth(method= "lm") +
  # scale_x_continuous(breaks= seq(0, 8000, by= 1000), 
  #                    labels= c("0", "", "2000", "", "4000", "", "6000", "", "8000"),
  #                    expand= c(0.02, 0)) +
  scale_x_log10(limits= c(1, 10000),
                expand= c(0, 0)) +
  annotation_logticks() +
  scale_y_continuous(limits= c(0, 0.0023), 
                     expand= c(0, 0)) +
  scale_fill_manual(values= species.colors,
                    labels= c("1", "2", "3"),
                    name= "Species") +
  scale_shape_manual(values= species.shapes,
                     labels= c("1", "2", "3"),
                     name= "Species") +
  theme_strains +
  theme(legend.position = c(0.92, 0.85))

nucDiv.watershed.plot
# ggsave(nucDiv.watershed.plot, filename = "nuc_div_watershed_gene.pdf", height= 180*0.75, width= 180, units= "mm", device= cairo_pdf,
#        path= "Output_figures")
# ggsave(nucDiv.watershed.plot, filename = "nuc_div_watershed_gene.png", dpi= 320, height= 180*0.75, width= 180, units= "mm",
#        path= "Output_figures")

# hist(log(gi_filt_summary$median_pi))
# 
# pi_fit1 <- lm(log(median_pi) ~ log10(watershed_km2), data= gi_filt_summary)
# pi_fit2 <- lm(log(median_pi) ~ log10(watershed_km2)*species, data= gi_filt_summary)
# summary(pi_fit2)
# anova(pi_fit1)
# anova(pi_fit1, pi_fit2)
# plot(pi_fit1)
# summary(lm(median_pi ~ watershed_km2, data= gi_sp1_filt_summary))



## Nucleotide diversity density curves
ggplot(data= gi_filt_df) +
  geom_density(aes(x= pi, color= sample)) +
  #geom_density(aes(x= pi)) +
  labs(x= "Nucleotide diversity", y= "Density") +
  scale_x_log10(limits= c(0.00002, 0.52), # the minimum pi value (apart from 0) is 0.000026
                breaks= c(0.0001, 0.001, 0.01, 0.1),
                labels= c("0.0001", "0.001", "0.01", "0.1"),
                expand= c(0, 0)) +
  annotation_logticks(sides= "b") +
  # scale_y_continuous(limits= c(0, 9),
  #                    breaks= c(0, 2, 4, 6, 8),
  #                    expand= c(0.02, 0)) +
  scale_color_discrete(guide= FALSE) +
  facet_grid(species~., scales= "free_y") +
  theme_strains
#ggsave(last_plot(), filename = "nuc_div_density.png", dpi= 320, height= 180*0.75, width= 180, units= "mm",
#       path= "Output_figures")


# 
# 
# ggplot(data= gi_filt_df) +
#   geom_histogram(aes(x= coverage, fill= species)) +
#   facet_wrap(~sample, nrow= 8, scales= "free_x")
# 
# ggplot(data= gi_filt_df) +
#   geom_point(aes(x= coverage, y= pi)) +
#   labs(x= "Coverage", y= "Nucleotide diversity") +
#   facet_grid(.~species) +
#   theme_strains
# 
# ggplot(data= gi_filt_df) +
#   geom_point(aes(x= coverage, y= SNPs_per_bp)) +
#   facet_grid(.~species) +
#   theme_bw()
# 
# ggplot(data= gi_filt_df) +
#   geom_point(aes(x= pi, y= SNPs_per_bp)) +
#   facet_grid(.~species) +
#   theme_bw()
# 
# ggplot(data= gi_filt_df) +
#   geom_boxplot(aes(x= sample, y= pi)) +
#   theme_bw()
# 
# ggplot(data= gi_filt_df) +
#   geom_point(aes(x= gene, y= pi, color= species)) +
#   theme(axis.text.x= element_blank()) +
#   facet_wrap(~sample, nrow= 8, scales= "free_x")
# 
# ggplot(data= gi_filt_df) +
#   geom_point(aes(x= gene, y= SNPs_per_bp, color= species)) +
#   theme(axis.text.x= element_blank()) +
#   facet_wrap(~sample, nrow= 8, scales= "free_x")
# 
# 
# ## SUMMARIZED PI VALUES
# 
# ggplot(data= gi_filt_df) +
#   geom_boxplot(aes(x= sample, y= pi)) +
#   theme_strains
# 
# 
# ggplot(data= gi_filt_summary) +
#   geom_point(aes(x= sample, y= sd_pi, color= species)) +
#   theme_bw()
# 
# ggplot(data= gi_filt_summary) +
#   geom_point(aes(x= sample, y= median_pi, shape= species)) +
#   geom_point(aes(x= sample, y= mean_pi, shape= species), color= "red") +
#   theme_bw()
# 
# 
# ggplot(data= gi_filt_summary) +
#   geom_point(aes(x= median_cov, y= median_pi, shape= species)) +
#   theme_strains
# 
# 
# 
# 
# 
# 
# ## Check for coverage and watershed area relationship
# ggplot(data= gi_filt_summary, aes(x= watershed_km2, y= median_cov)) +
#   geom_point(aes(color= species, shape= species), size= 3) +
#   labs(x= expression(paste("Watershed area (", km^{2}, ")")), y= "Median coverage") +
#   #scale_x_continuous(expand= c(0.02, 0)) +
#   scale_x_log10(expand= c(0.02, 0)) +
#   geom_smooth(method= "lm") +
#  # scale_y_continuous(limits= c(0, 0.0023), expand= c(0.02, 0)) +
#   theme_classic(base_size= 20)
# 
# 
# fit1 <- lm(median_pi ~ log10(watershed_km2), data= filter(gi_filt_summary, watershed_km2 < 2000))
# summary(fit1)
# plot(fit1)
# 
# getwd()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
