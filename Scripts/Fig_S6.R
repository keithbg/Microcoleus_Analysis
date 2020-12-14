## Histogram of SNV site s sharing at fixed and bi-allelic sites


library(tidyverse)
library(RColorBrewer)
library(lemon)
library(cowplot)
source("Scripts/ggplot_themes.R")


## SNV SHARING X WATERSHED AREA
#source("Scripts/SNVS_analysis.R")
snv_dist_mat.m1 <- read_tsv("Data/inStrain_data/snv_pos_M1_jaccard.tsv") %>% 
  mutate(morphia= "m1")
snv_dist_mat.m2 <- read_tsv("Data/inStrain_data/snv_pos_M2_jaccard.tsv") %>% 
  mutate(morphia= "m2")

snv_dist_list <- list(m1= snv_dist_mat.m1, m2= snv_dist_mat.m2)


# Make SNV_distance long and merge with river distance  
snv_dist_df_list <- map(snv_dist_list, function(x) x %>% 
                          mutate(siteA= names(.)[-48]) %>% 
                          mutate_all(as.character) %>% 
                          pivot_longer(., names_to = "siteB", values_to= "snv_distance", 
                                       cols= starts_with("PH")) %>% 
                          mutate(snv_distance= as.numeric(snv_distance))) #%>% 
                          #left_join(., river_dist_df) %>% 
                          #inner_join(., ani.dist.df) %>% 
                          #left_join(., select(watershed_area, siteA, watershed_km2_A)) %>% 
                          #left_join(., select(watershed_area, siteB, watershed_km2_B)) %>% 
                          #mutate(watershed_diff= abs(watershed_km2_A - watershed_km2_B)))

snv_dist_df <- bind_rows(snv_dist_df_list) %>% 
  filter(siteA != siteB)



#### MAKE FIGURE ####
morphia.labels <- as_labeller(c(`m1` = "Fixed sites", `m2` = "Bi-allelic sites"))

snv_histogram <-  ggplot() +
  geom_histogram(data= filter(snv_dist_df, morphia == "m1", siteB != "PH2015_12D", siteA != "PH2015_12D", siteB != "PH2015_12U", siteA != "PH2015_12U"),
                 aes(x= (1-snv_distance)*100),
                 binwidth= 1, boundary=1,  fill= "black") +
  geom_histogram(data= filter(snv_dist_df, morphia == "m2"),
                 aes(x= (1-snv_distance)*100),
                 binwidth= 1, boundary=1,  fill= "black") +
  labs(x= "SNV site similarity (%)", y= "Paired sample comparisons") +
  scale_y_continuous(expand= c(0, 0)) +
  scale_x_continuous(limits= c(0, 100), 
                     expand= c(0, 0)) +
  facet_rep_wrap(~morphia, nrow=3,
                 labeller= labeller(morphia= morphia.labels)) +
  theme_strains


ggsave(snv_histogram, filename= "Fig_S6.png", height= 180, width= 180, units= "mm", dpi= 320, 
       path= "Output_figures")
