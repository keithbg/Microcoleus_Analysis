

library(tidyverse)
library(ggplot2)
source("Scripts_manuscript/Fig_2_data.R")
source("Scripts_manuscript/ggplot_themes.R")



ggplot(data= ani_sum) +
  geom_vline(xintercept = 99.35,linetype= "dashed", color= "gray60", size= 0.5) +
  geom_point(aes(x= mean_conANI*100, y= 1 - frac_popSNPs, fill= pi_avg_mean, size= pi_diff_mean), 
             pch= 21,  color= "black") +
  labs(x= "Consensus ANI (%)", y= "Percent shared minor alleles (pMA; %)") +
  scale_x_continuous(limits= c(98.74, 100),
                     breaks= seq(98.80, 100, by= 0.1), 
                     labels= c("98.8", "", "99.0", "", "99.2", "", "99.4", "", "99.6", "", "99.8", "", "100"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.12, 1.025),
                     breaks= seq(0, 1, by= 0.1),
                     labels= c("0", "", "20", "", "40", "", "60", "", "80", "", "100"),
                     expand= c(0, 0)) +
  scale_fill_viridis_c(name= "Avg.\nnucleotide\ndiversity",
                       limits= c(0.0004, 0.0034)) +
  scale_size_binned(name= "Difference in\nnucleotide\ndiversity",
                    range= c(1, 7),
                    n.breaks= 3,
                    nice.breaks = TRUE) +
  guides(fill = guide_colourbar(order = 1),
         size_binned = guide_legend(order = 2)) +
  theme_strains +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.key.width= unit(10, "mm"),
        legend.box.just = "bottom",
        legend.box.spacing = unit(0, "mm"),
        legend.text= element_text(size= 8),
        legend.title= element_text(size= 10),
        legend.background = element_rect(color= "transparent", fill= "transparent"))

ggsave(last_plot(), filename = "Fig_5.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/Output_figures")


