## Figures 6 and S6

# Minor allele frequencies in Species 1

#### Libraries #################################################################
library(tidyverse)
source("Scripts/ggplot_themes.R")
################################################################################



#### SNV data (input table generated in: format_inStrain_output.R) ###
snv_df <- read_tsv(file.path("Data/inStrain_data", "snv_df_filt_v1.4.tsv"))
snvs_genome_df <- read_tsv(file.path("Data/inStrain_data", "snvs_genome_summary_v1.4.tsv"))


## CALCULATE MINOR ALLELE FREQUENCIES FOR SPECIES 1
## Samples with secondary SNV peaks
secondary.snv.peaks <- c("2015_03D", "2015_03U", "2015_04D", "2015_04U", "2015_10S", "2017_02_FOX", 
                         "2017_03_ELD", "2017_04_SCI", "2017_05_CCC", "2017_06_SFM", "2017_07_MST")

snv.freq.sp1 <- snv_df %>% 
  filter(species == "species_1") %>% 
  #filter(., site == "PH2015_03U") %>% 
  mutate(varFreq_r3 = round(var_freq, 3),
         varFreq_r2 = round(var_freq, 2)) %>% 
  group_by(site, varFreq_r2) %>% 
  summarize(n= length(var_freq)) %>% 
  ungroup() %>% 
  left_join(., select(snvs_genome_df, ggkbase_id, NS, SNV_mbp), by= c("site" = "ggkbase_id")) %>% 
  mutate(facet_label= str_sub(site, start= 3, end= 13)) %>% 
  mutate(sec.peak= ifelse(facet_label %in% secondary.snv.peaks, "Y", "N"))

snv.freq.NS <- snv.freq.sp1 %>% 
  select(site, NS, SNV_mbp, facet_label, sec.peak) %>% 
  distinct() 
  
#### STATISTICS ####
#summary(filter(snv.freq.NS, sec.peak == "Y")$NS)
#summary(filter(snv.freq.NS, sec.peak == "N")$NS)
#summary(lm(NS ~ sec.peak, data= snv.freq.NS))
#plot(lm(NS ~ sec.peak, data= snv.freq.NS))


## Samples 12U and 12D are large N:S outliers. These are also the reference genome, which introduces some bias into the results.
## Therefore I am removing these genomes for the analysis and report the values below in the manuscript
summary(filter(snv.freq.NS, sec.peak == "Y" & site != "PH2015_12U" & site != "PH2015_12D")$NS)
summary(filter(snv.freq.NS, sec.peak == "N" & site != "PH2015_12U" & site != "PH2015_12D")$NS)


boxplot(NS ~ sec.peak, data= filter(snv.freq.NS, site != "PH2015_12U" & site != "PH2015_12D"))
summary(lm(NS ~ sec.peak, data= filter(snv.freq.NS, site != "PH2015_12U" & site != "PH2015_12D")))


#### FIGURE 6 ####
ggplot(snv.freq.NS, aes(x= sec.peak, y= NS)) +
  geom_boxplot() +
  geom_point(position= "jitter") +
  theme_strains


sec.peaks.combined <- ggplot(filter(snv.freq.sp1, facet_label != "2015_01D" & sec.peak == "Y"), 
                             aes(x= varFreq_r2, y= n, group= facet_label)) +
  geom_point(aes(color= facet_label), size= 1, alpha= 0.7) +
  geom_smooth(aes(color= facet_label), method= "loess", span= 0.2, se= FALSE, size= 0.75) +
  labs(x= "Minor allele frequency", y= "Number of SNV sites") +
  scale_x_continuous(limits= c(0.04, 0.5), 
                     breaks= seq(0.05, 0.5, by= 0.05),
                     labels= c("", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5"),
                     expand= c(0,0)) +
  scale_y_continuous(limits= c(0, 3500)) +
  scale_color_manual(values= c("purple", pal_npg("nrc")(10)), guide= FALSE) +
  lemon::facet_rep_wrap(~facet_label, ncol= 4, scales= "free_x") +
  theme_strains
#sec.peaks.combined
ggsave(sec.peaks.combined, filename = "Fig_6_v1.4.png", height= 180*0.66, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")
  
#### FIGURE S6 ####
secondary.peaks <- ggplot(filter(snv.freq.sp1, facet_label %in% secondary.snv.peaks & varFreq_r2 >= 0.05 & varFreq_r2 <= 0.5),
                          aes(x= varFreq_r2, y= n)) +
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


no.secondary.peaks <- ggplot(filter(snv.freq.sp1, !(facet_label %in% secondary.snv.peaks) & varFreq_r2 >= 0.05 & varFreq_r2 <= 0.5 & (facet_label != "2015_01D") & (facet_label != "2015_01U")),  
                             aes(x= varFreq_r2, y= n)) +
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

ggsave(snv.freq.fig.anno, filename = "Fig_S6_v1.4.png", height= 180, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


