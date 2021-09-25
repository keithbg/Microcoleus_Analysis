## Figure S7

## Boxplots of inter-species Fst values

## PopGenome statistics on genes extracted from Roary output
## Plotting the data that was formatted in PopGenome_Roary_Analysis.R,


#### Libraries #################################################################
library(tidyverse)
source("Scripts/ggplot_themes.R")
################################################################################


#### READ DATA
pg.btw.filt <- read_tsv(file.path("Data/PopGenome_data", "PopGenome_btw_filt.tsv"))



#### CALC SUMMARY STATISTICS
pg.btw.filt %>% 
  filter(FST > 0) %>% 
  group_by(pops) %>% 
  summarize(
    min= min(FST, na.rm= TRUE),
    mean= mean(FST, na.rm= TRUE),
    sd= sd(FST, na.rm= TRUE),
    med= median(FST, na.rm= TRUE),
    max= max(FST, na.rm= TRUE)) %>% 
  ungroup()


#### MAKE FIGURE ####
ggplot(data= pg.btw.filt) +
  geom_boxplot(aes(x= pops, y= FST), outlier.size = 0.5, size= 0.3) +
  labs(x= NULL, y= expression("F"[st])) +
  scale_x_discrete(labels= c("Species\n1 & 2", "Species\n1 & 3", "Species\n2 & 3")) +
  scale_y_continuous(limits= c(0, 1), expand= c(0, 0), breaks= seq(0, 1, by= 0.1), 
                     labels= c("0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0")) +
  theme_strains +
  theme(text= element_text(size= 8))

ggsave(last_plot(), file= "Fig_S7.png", width= 90, height= 90, units= "mm", dpi= 320,
       path= "Output_figures") 


