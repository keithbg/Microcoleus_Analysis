## Histogram of populuation and consensus ANI values

library(tidyverse)
source("Scripts/ggplot_themes.R")

ani_sum <- read_tsv("Data/inStrain_data/ani_summary_v1.4.tsv") # ani_summary.tsv generateed in "Scripts/ANI_scaffold_data.R"
#ani_sum08 <- read_tsv("Data/inStrain_data/ani_summary_v0.8.tsv") 

## 99.6% conANI threshold

# 67 genomes with conANI > 99.6%
ani_sum %>% 
  count(mean_conANI >= 0.996)
67/1081 # 6.6% of genome pairs

ani_sum_996 <-  ani_sum %>% 
  filter(mean_conANI >= 0.996)



conANI_hist <- ggplot(ani_sum, aes(x= round(mean_conANI, 4))) +
  geom_histogram(binwidth= 0.0005, boundary= 1, fill= "black", color= "gray75") +
  labs(x= "Mean consensus ANI", y= "Count") +
  scale_x_continuous(breaks= seq(0.992, 1, by= 0.001),
                     labels= c("99.2", "99.3", "99.4", "99.5", "99.6", "99.7", "99.8", "99.9", "100"),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  theme_strains

popANI_hist <- ggplot(ani_sum, aes(x= round(mean_popANI, 4))) +
  geom_histogram(binwidth= 0.0005, boundary= 1, fill= "black", color= "gray75") +
  labs(x= "Mean population ANI", y= "Count") +
  scale_x_continuous(limits= c(0.9915, 1),
                     breaks= seq(0.992, 1, by= 0.001),
                     labels= c("99.2", "99.3", "99.4", "99.5", "99.6", "99.7", "99.8", "99.9", "100"),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  theme_strains

ANI_hist_combined <- plot_grid(popANI_hist + theme(axis.title.x = element_blank()), 
                               conANI_hist + labs(x= "Average nucleotide identity (%)"),
                               nrow= 2,
                               labels= c("A", "B")) +
  draw_label(label= "Population\nANI", x= 0.8, y= 0.91, hjust= 0) +
  draw_label(label= "Consensus\nANI", x= 0.8, y= 0.41, hjust= 0)

ggsave(ANI_hist_combined, filename = "Fig_S3_v1.4.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")
