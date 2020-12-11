setwd("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis")


library(tidyverse)
library(ggplot2)
library(cowplot)
source("Scripts/ANI_scaffold_data.R")
source("Scripts/ggplot_themes.R")


#### STATISTICS ####
fit.popANI <- lm((mean_popANI*100) ~ riv_dist, ani_sum)
fit.conANI <- lm((mean_conANI*100) ~ riv_dist, ani_sum)

conANI_25km_fit1 <- lm(log(mean_conANI*100) ~ riv_dist + watershed_diff, data= filter(ani_sum, riv_dist < 25000))
#summary(conANI_25km_fit1)
conANI_25km_fit2 <- lm(mean_conANI*100 ~ riv_dist * watershed_diff, data= filter(ani_sum, riv_dist < 25000))
#anova(conANI_25km_fit2)
anova(conANI_25km_fit1, conANI_25km_fit2)


#### FIGURES ####

## popANI
popANI_rivDist <- ggplot(data= ani_sum) +
  geom_point(aes(x= riv_dist, y= mean_popANI*100), pch= 21, fill= species.colors[1], color= "black", size= 3, alpha= 0.5) +
  geom_abline(intercept = fit.popANI$coefficients["(Intercept)"], slope= fit.popANI$coefficients["riv_dist"],
              color= "black", size= 1) +
  labs(x= "River network distance (km)", y= "Population ANI (%)") +
  scale_x_continuous(limits= c(0, 350000),
                     breaks= seq(0, 350000, by= 25000),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5000)) +
  theme_strains +
  theme(legend.position = "top")

## conANI
conANI_rivDist <- ggplot(data= ani_sum, aes(x= riv_dist, y= mean_conANI*100)) +
  geom_hline(yintercept = 99.35, linetype= "dashed", color= "gray60", size= 0.5) +
  annotate("text", x= 350*1000, y= 99.4, label= "99.35%", hjust= 1, vjust= 0, color= "gray30", size= 3) +
  geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_abline(intercept = fit.conANI$coefficients["(Intercept)"], slope= fit.conANI$coefficients["riv_dist"],
              color= "black", size= 1) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_x_continuous(limits= c(0, 350000),
                     breaks= seq(0, 350000, by= 25000),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5000)) +
  scale_y_continuous(breaks= seq(98.80, 100.00, by= 0.2),
                     labels= c("98.8", "99.0", "99.2", "99.4", "99.6", "99.8", "100.0")) +
  theme_strains +
  theme(legend.position = "top")
#conANI_rivDist

## conANI < 25 km
conANI_rivdist25km <- ggplot(data= filter(ani_sum, riv_dist < 25000), aes(x= riv_dist, y= mean_conANI*100)) +
  geom_hline(yintercept = 99.35, linetype= "dashed", color= "gray60", size= 0.5) +
  annotate("text", x= 25000, y= 99.4, label= "99.35%", hjust= 1, vjust= 0, color= "gray30",size= 3) +
  #geom_point(aes(fill= watershed_diff), size= 4, pch= 21, color= "black") +
  geom_point(aes(fill= watershed_diff, size= watershed_diff), shape= 21, color= "gray50") +
  geom_abline(intercept = conANI_25km_fit2$coefficients["(Intercept)"], slope= conANI_25km_fit2$coefficients["riv_dist"],
              color= "black", size= 1) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_fill_viridis_c(name= expression("Watershed difference (km"^2*")"),
                       option= "plasma") +
  scale_size_continuous(range= c(2.5, 6), guide= FALSE) +
  scale_x_continuous(limits= c(0, 25001),
                     breaks= seq(0, 25000, by= 5000),
                     expand= c(0.02, 0)) +
  scale_y_continuous(breaks= seq(98.80, 100.00, by= 0.2),
                     labels= c("98.8", "99.0", "99.2", "99.4", "99.6", "99.8", "100.0")) +
  theme_strains +
  theme(legend.position = c(0.7, 0.85),
        legend.direction = "horizontal",
        legend.background = element_rect(color= "transparent", fill= "transparent"))
conANI_rivdist25km


ANI_rivdist_combined <- plot_grid(popANI_rivDist+ theme(axis.title.x = element_blank()), 
                                  conANI_rivDist + theme(axis.title.x = element_blank()), 
                                  conANI_rivdist25km,
                                  nrow= 3,
                                  labels= c("A", "B", "C"))


ggsave(ANI_rivdist_combined, filename = "Fig_2.png", height= 180, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



