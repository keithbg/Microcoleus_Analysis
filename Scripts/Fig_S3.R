## Histogram of populuation and consensus ANI values

library(tidyverse)
source("Scripts/ggplot_themes.R")

ani_sum <- read_tsv("Data/inStrain_data/ani_summary.tsv") # ani_summary.tsv generateed in "Scripts/ANI_scaffold_data.R"


conANI_hist <- ggplot(ani_sum, aes(x= mean_conANI)) +
  geom_histogram(binwidth= 0.0005, boundary= 1, fill= "black", color= "gray75") +
  labs(x= "Mean consensus ANI", y= "Count") +
  scale_x_continuous(breaks= seq(0.987, 1, by= 0.001),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  theme_strains

popANI_hist <- ggplot(ani_sum, aes(x= mean_popANI)) +
  geom_histogram(binwidth= 0.0005, boundary= 1, fill= "black", color= "gray75") +
  labs(x= "Mean population ANI", y= "Count") +
  scale_x_continuous(limits= c(0.987, 1),
                     breaks= seq(0.987, 1, by= 0.001),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  theme_strains

ANI_hist_combined <- plot_grid(conANI_hist + theme(axis.title.x = element_blank()), 
                               popANI_hist + labs(x= "Average nucleotide identity"),
                               nrow= 2,
                               labels= c("A", "B")) +
  draw_label(label= "Consensus\nANI", x= 0.12, y= 0.91, hjust= 0) +
  draw_label(label= "Population\nANI", x= 0.12, y= 0.41, hjust= 0)

ggsave(ANI_hist_combined, filename = "Fig_S3.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")
