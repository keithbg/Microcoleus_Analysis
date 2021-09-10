## Figure 2

## Relationship between popANI, conANI, and river network distance

#### Libraries #################################################################
library(tidyverse)
source("Scripts/ggplot_themes.R")
################################################################################


#### READ DATA ####
## Read in ANI Summary data
ani_sum <- read_tsv("Data/inStrain_data/ani_summary_v1.4.tsv") # .tsv file generated in: Scripts/inStrain_ANI_summary.R

## Create data frame with RivDist and ANI_Summary
load("Data/flowDist_Vectors.Rdata") # from RiverDistances.R
## Combine data frame and flow-distance vectors together
site_pairs_xyVert_flow <- bind_cols(site_pairs_xyVert, 
                                    tibble(flowConnected, flowDistTotal, flowDistNet)) %>% 
  mutate(FlowConnection= ifelse(is.na(flowConnected), "No", "Yes")) 

ani_rivDist <- left_join(select(ani_sum, name1, name2, mean_conANI, mean_popANI, riv_dist, euc_dist, watershed_diff),
                         select(site_pairs_xyVert_flow, name1, name2, FlowConnection, flowDistTotal, flowDistNet)) %>% 
  mutate(riv_dist= riv_dist/1000,
         euc_dist= euc_dist/1000,
         flowDistTotal= abs(flowDistTotal/1000),
         dist_diff= riv_dist - flowDistTotal,
         mean_conANI= mean_conANI*100,
         mean_popANI= mean_popANI*100)


#### STATISTICS ####
fit.noflow2.conANI <- lm(mean_conANI ~ flowDistTotal*FlowConnection, data= ani_rivDist)
summary(fit.noflow2.conANI)
anova(fit.noflow2.conANI)

fit.noflow2.popANI <- lm(mean_popANI ~ flowDistTotal*FlowConnection, data= ani_rivDist)
summary(fit.noflow2.popANI)
anova(fit.noflow2.popANI)
plot(fit.noflow2.popANI)

conANI_25km_fit1 <- lm((mean_conANI) ~ flowDistTotal + watershed_diff, data= filter(ani_rivDist, flowDistTotal < 25 & FlowConnection == "Yes"))
summary(conANI_25km_fit1)

conANI_25km_fit2 <- lm((mean_conANI) ~ flowDistTotal * watershed_diff, data= filter(ani_rivDist, flowDistTotal < 25 & FlowConnection == "Yes"))
summary(conANI_25km_fit2)
anova(conANI_25km_fit2)

anova(conANI_25km_fit1, conANI_25km_fit2)



#### FIGURES ####
legend_title <- "Connected flow"

## Consensus ANI
conANI2 <- ggplot(data= ani_rivDist, aes(x= flowDistTotal, y= mean_conANI)) +
  geom_point(aes(fill= FlowConnection, shape= FlowConnection), size= 3, color= "black", alpha= 0.5) +
  geom_line(aes(x= flowDistTotal, y= predict(fit.noflow2.conANI), color= FlowConnection, group= FlowConnection), size= 2) +
  scale_fill_manual(values= c(species.colors[1], "wheat2"), name= legend_title) +
  scale_color_manual(values= c(species.colors[1], "wheat2"), name= legend_title) +
  scale_shape_manual(values= c(23, 21), name= legend_title) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_x_continuous(limits= c(0, 350),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains
conANI2

## Population ANI
popANI2 <- ggplot(data= ani_rivDist, aes(x= flowDistTotal, y= mean_popANI)) +
  geom_point(aes(fill= FlowConnection, shape= FlowConnection), size= 3, color= "black", alpha= 0.5) +
  geom_line(aes(x= flowDistTotal, y= predict(fit.noflow2.popANI), color= FlowConnection, group= FlowConnection), size= 2) +
  scale_fill_manual(values= c(species.colors[1], "wheat2"), name= legend_title) +
  scale_color_manual(values= c(species.colors[1], "wheat2"), name= legend_title) +
  scale_shape_manual(values= c(23, 21), name= legend_title) +
  labs(x= "River network distance (km)", y= "Population ANI (%)") +
  scale_x_continuous(limits= c(0, 350),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains
#popANI2


## conANI < 25 km
short_dist <- filter(ani_rivDist, flowDistTotal < 1)

sd2 <- short_dist %>% 
  filter((name1 == "PH2015_02D" & name2 == "PH2015_03U") | 
           (name1 == "PH2015_02D" & name2 == "PH2015_03D") |
           (name1 == "PH2015_02D" & name2 == "PH2015_04U") |
           (name1 == "PH2015_02D" & name2 == "PH2015_04D") |
           (name1 == "PH2015_02D" & name2 == "PH2017_05_CCC_O_A") |
           (name1 == "PH2015_02U" & name2 == "PH2015_03U") | 
           (name1 == "PH2015_02U" & name2 == "PH2015_03D") |
           (name1 == "PH2015_02U" & name2 == "PH2015_04D") |
           (name1 == "PH2015_02U" & name2 == "PH2015_04U") |
           (name1 == "PH2015_02U" & name2 == "PH2017_05_CCC_O_A") |
           (name1 == "PH2017_01_SCD_O_A" & name2 == "PH2015_03U") | 
           (name1 == "PH2017_01_SCD_O_A" & name2 == "PH2015_03D") |
           (name1 == "PH2017_01_SCD_O_A" & name2 == "PH2015_04D") |
           (name1 == "PH2017_01_SCD_O_A" & name2 == "PH2015_04U") |
           (name2 == "PH2017_01_SCD_O_A" & name1 == "PH2015_03U") | 
           (name2 == "PH2017_01_SCD_O_A" & name1 == "PH2015_03D") |
           (name2 == "PH2017_01_SCD_O_A" & name1 == "PH2015_04D") |
           (name2 == "PH2017_01_SCD_O_A" & name1 == "PH2015_04U") |
           (name1 == "PH2017_01_SCD_O_A" & name2 == "PH2017_05_CCC_O_A"))
           
sd2 %>% 
  select(name1, name2, mean_conANI) %>% 
  pivot_wider(names_from = name2, values_from= mean_conANI)

sd3 <- short_dist %>% 
  filter(name1 == "PH2015_02U" | name1 == "PH2015_02D" |
           name1 == "PH2015_03U" | name1 == "PH2015_03D" |
           name1 == "PH2015_04U" | name1 == "PH2015_04D" |
           name1 == "PH2017_01_SCD_O_A" | name1 == "PH2017_05_CCC_O_A") %>% 
  filter(name2 == "PH2015_02U" | name2 == "PH2015_02D" |
           name2 == "PH2015_03U" | name2 == "PH2015_03D" |
           name2 == "PH2015_04U" | name2 == "PH2015_04D" |
           name2 == "PH2017_01_SCD_O_A" | name2 == "PH2017_05_CCC_O_A") #%>% 
  filter((!(name1 == "PH2015_02D" & name2 == "PH2015_02U") &
           !(name1 == "PH2015_02D" & name2 == "PH2017_01_SCD_O_A") &
           !(name1 == "PH2015_02U" & name2 == "PH2017_01_SCD_O_A")) |
           (name1 == "PH2015_03D" & name2 == "PH2017_01_SCD_O_A") |
           (name1 == "PH2015_03U" & name2 == "PH2017_01_SCD_O_A") |
           (name1 == "PH2015_04D" & name2 == "PH2017_01_SCD_O_A") |
           (name1 == "PH2015_04U" & name2 == "PH2017_01_SCD_O_A") |
           (name1 == "PH2017_01_SCD_O_A" & name2 == "PH2017_05_CCC_O_A"))

"PH2017_01_SCD_O_A"
conANI_25km <- ggplot(data= filter(ani_rivDist, flowDistTotal < 25 & FlowConnection == "Yes"), aes(x= flowDistTotal, y= mean_conANI)) +
  geom_point(aes(fill= watershed_diff, size= watershed_diff), shape= 21, color= "gray50") +
  #geom_line(aes(x= flowDistTotal, y= predict(conANI_25km_fit1)), size= 2) +
  geom_abline(intercept = conANI_25km_fit2$coefficients["(Intercept)"], slope= conANI_25km_fit2$coefficients["flowDistTotal"],
              color= "black", size= 1) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_fill_viridis_c(name= expression("Watershed difference (km"^2*")"),
                       option= "plasma") +
  scale_size_continuous(range= c(2.5, 6), guide= FALSE) +
  scale_x_continuous(limits= c(0, 25.001),
                     breaks= seq(0, 25, by= 2.5),
                     labels= c("0", "", "5", "", "10", "", "15", "", "20", "", "25"),
                     expand= c(0.02, 0)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.00, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains +
  theme(legend.position = c(0.7, 0.85),
        legend.direction = "horizontal",
        legend.background = element_rect(color= "transparent", fill= "transparent"))
#onANI_rivdist25km




ANI_rivdist_combined2 <- plot_grid(conANI2 + theme(axis.title.x = element_blank(),
                                                   legend.direction = "horizontal",
                                                   legend.position = c(0.8, 0.95)), 
                                   popANI2 + theme(axis.title.x = element_blank(),
                                                   legend.position = "none"), 
                                   conANI_25km + theme(axis.title.x = element_blank()),
                                   nrow= 3,
                                   labels= c("A", "B", "C"))
#ANI_rivdist_combined2

ANI_rivdist_combined3 <-  annotate_figure(ANI_rivdist_combined2, 
                                          bottom= text_grob(label= "River network distance (km)", vjust= -0.5))

#ANI_rivdist_combined3

ggsave(ANI_rivdist_combined3, filename = "Fig_2B_v1.4.png", height= 180, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


