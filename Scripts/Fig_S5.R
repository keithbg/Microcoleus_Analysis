## Figure S5

## Histogram of SNV site sharing at fixed and bi-allelic sites

#### Libraries #################################################################
library(tidyverse)
source("Scripts/ggplot_themes.R")
################################################################################


## SNV SHARING X WATERSHED AREA
# jaccard distance files generated in: Scripts/inStrain_SNVS_analysis.R
# AC1: allele count = 1, fixed sites
# AC2: allele count = 2, biallelic sites
snv_dist_mat.AC1 <- read_tsv("Data/inStrain_data/snv_pos_AC1_jaccard.tsv") %>% 
  mutate(allele_count= "AC1")
snv_dist_mat.AC2 <- read_tsv("Data/inStrain_data/snv_pos_AC2_jaccard.tsv") %>% 
  mutate(allele_count= "AC2")

snv_dist_list <- list(AC1= snv_dist_mat.AC1, AC2= snv_dist_mat.AC2)


# Make SNV_distance long and merge with river distance  
snv_dist_df_list <- map(snv_dist_list, function(x) x %>% 
                          mutate(siteA= names(.)[-48]) %>% 
                          mutate_all(as.character) %>% 
                          pivot_longer(., names_to = "siteB", values_to= "snv_distance", 
                                       cols= starts_with("PH")) %>% 
                          mutate(snv_distance= as.numeric(snv_distance)))
                        

snv_dist_df <- bind_rows(snv_dist_df_list) %>% 
  filter(siteA != siteB)



#### MAKE FIGURE ####
allele_count.labels <- as_labeller(c(`AC1` = "Fixed sites", `AC2` = "Bi-allelic sites"))

snv_histogram <-  ggplot() +
  geom_histogram(data= filter(snv_dist_df, allele_count == "AC1", siteB != "PH2015_12D", siteA != "PH2015_12D", siteB != "PH2015_12U", siteA != "PH2015_12U"),
                 aes(x= (1-snv_distance)*100),
                 binwidth= 1, boundary=1,  fill= "black") +
  geom_histogram(data= filter(snv_dist_df, allele_count == "AC2"),
                 aes(x= (1-snv_distance)*100),
                 binwidth= 1, boundary=1,  fill= "black") +
  labs(x= "SNV site similarity (%)", y= "Paired sample comparisons") +
  scale_y_continuous(expand= c(0, 0)) +
  scale_x_continuous(limits= c(0, 100),
                     breaks= seq(0, 100, 10),
                     expand= c(0, 0)) +
  facet_rep_wrap(~allele_count, nrow=3,
                 labeller= labeller(allele_count= allele_count.labels)) +
  theme_strains


ggsave(snv_histogram, filename= "Fig_S5_v1.4.png", height= 180*.75, width= 180, units= "mm", dpi= 320, 
       path= "Output_figures")
