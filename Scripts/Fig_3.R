library(tidyverse)
library(ggplot)
library(RColorBrewer)
library(cowplot)
library(ggpubr) #ggarrange
source("Scripts/ggplot_themes.R")

#### NUCLEOTIDE DIVERSITY X WATERSHED AREA ####################################
## Investigate nucleotide diversity (pi) results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.gene_profile.tsv

#### FORMAT DATA ####

## Input df generated in format_inStrain_output.R
gi_filt_df <- read_tsv("Data/inStrain_data/gene_info_filt_df.tsv")

## Watershed area data
watershed.area <-
  read_tsv(file.path("Data", "WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2) %>% 
  rename(site= ggkbase_id)

## Summarize pi values across genome 
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

#### STATISTICS ####
fit.pi <- lm(median_pi ~ watershed_km2, data= gi_filt_summary)
summary(fit.pi)

#### FIGURE ####

nucDiv.watershed.plot <- ggplot(data= gi_filt_summary, aes(x= watershed_km2, y= median_pi)) +
  geom_point(aes(fill= species, shape= species), size= 3, color= "black") +
  labs(x= expression(paste("Watershed area (", km^{2}, ")")), y= "Median nucleotide diversity") +
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
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1))
#nucDiv.watershed.plot
###############################################################################

#### popANI X WATERSHED AREA ##################################################
ani_sum <- read_tsv("Data/inStrain_data/ani_summary.tsv") # ani_sum generated in ANI_scaffold_data.R

#### FORMAT DATA ####
## Watershed area data
watershed.area.popANI <- read_tsv(file.path("Data", "Spatial_data", "WatershedArea_Combined.tsv")) %>%
  select(ggkbase_id, watershed_km2)

popANI.wide <- pivot_wider(select(ani_sum, name1, name2, mean_popANI), 
                           names_from = name2, 
                           values_from = mean_popANI) %>% 
  ungroup() %>% 
  mutate(PH2015_01D= NA) %>% 
  select(name1, PH2015_01D, everything()) %>% 
  add_row(., name1= "PH2017_40_RAT_O_B")

gdata::lowerTriangle(popANI.wide[, -1]) <- gdata::upperTriangle(popANI.wide[, -1], byrow=TRUE)

popANI.long <- pivot_longer(popANI.wide, names_to= "name2", values_to= "mean_popANI", -name1) %>% 
  # WATERSHED AREA JOIN
  left_join(., watershed.area.popANI, by= c("name1" = "ggkbase_id")) %>% 
  left_join(., watershed.area.popANI, by= c("name2" = "ggkbase_id")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y)

popANI.sum <- popANI.long %>% 
  group_by(name1, watershed_1) %>% 
  summarise(meanPOPani= mean(mean_popANI, na.rm=TRUE),
            sdPOPani= sd(mean_popANI, na.rm=TRUE)) %>% 
  ungroup()
popANI.sum

#### STATISTICS ####
fit.popANI1 <- lm(meanPOPani ~ log10(watershed_1), data= popANI.sum)
summary(fit.popANI1)
plot(fit.popANI1)
hist(log(popANI.sum$meanPOPani))

#### FIGURE #####
popANI_watershed <- ggplot(data= popANI.sum, aes(x= watershed_1, y= meanPOPani*100)) +
  #geom_errorbar(aes(ymin= meanPOPani*100 - sdPOPani*100, ymax= meanPOPani*100 + sdPOPani*100), size= 0.2, color= "gray50") +
  geom_point(shape= 21, fill= species.colors[1], color= "black", size= 3) +
  labs(x= expression("Watershed area (km"^2*")"), y= "Mean population ANI (%)") +
  scale_x_log10(limits=c(1, 2200),
                expand= c(0, 0)) +
  annotation_logticks(sides= "b") +
  scale_y_continuous(breaks= seq(99.4, 99.8, by= 0.2),
                     labels= c("99.4", "99.6", "99.8")) +
  #labels= c(" 0.994", " 0.996", " 0.998")) +
  theme_strains
popANI_watershed
###############################################################################



#### SNV SHARING X WATERSHED AREA #############################################
#source("Scripts/SNVS_analysis.R")

#### FORMAT DATA ####
## Input Jaccard index of SNV sharing sites
snv_dist_list <- list(m1= read_tsv("Data/inStrain_data/snv_pos_M1_jaccard.tsv") %>% 
                        mutate(morphia= "m1"),
                      m2= read_tsv("Data/inStrain_data/snv_pos_M2_jaccard.tsv") %>% 
                        mutate(morphia= "m2"))


## River network distance between sites in meters
river_dist_mat <- read_tsv("Data/Spatial_data/Distance_RiverNetwork_meters.tsv") %>% 
  rename(site= Site)

# Make river distances long
river_dist_df <- river_dist_mat %>% 
  rename(siteA= site) %>% 
  gather(key= siteB, value= riv_dist, -siteA)

# Watershed areas
watershed_area <- read_tsv("Data/Spatial_data/WatershedArea_Combined.tsv") %>% 
  select(ggkbase_id, watershed_km2) %>% 
  mutate(siteA= ggkbase_id,
         watershed_km2_A= watershed_km2,
         siteB= ggkbase_id,
         watershed_km2_B= watershed_km2) %>% 
  select(-ggkbase_id, -watershed_km2)

# Make SNV_distance long and merge with river distance  
snv_dist_list_mutate <- map(snv_dist_list, function(x) x %>% 
                              mutate(siteA= names(.)[-48]) %>% 
                              mutate_all(as.character) %>% 
                              pivot_longer(., names_to = "siteB", values_to= "snv_distance", 
                                           cols= starts_with("PH")) %>% 
                              mutate(snv_distance= as.numeric(snv_distance)) %>% 
                              left_join(., river_dist_df) %>% 
                              #inner_join(., ani.dist.df) %>% 
                              left_join(., select(watershed_area, siteA, watershed_km2_A)) %>% 
                              left_join(., select(watershed_area, siteB, watershed_km2_B)) %>% 
                              mutate(watershed_diff= abs(watershed_km2_A - watershed_km2_B)))

snv_dist_df <- bind_rows(snv_dist_list_mutate) %>% 
  filter(siteA != siteB)

#### STATISTICS ####
fit.m1 <- lm(((1-snv_distance)) ~ riv_dist * log10(watershed_km2_A), filter(snv_dist_df, morphia == "m1", siteB != "PH2015_12D", siteA != "PH2015_12D", siteB != "PH2015_12U", siteA != "PH2015_12U"))
summary(fit.m1)
anova(fit.m1)

fit.m2 <- lm(((1-snv_distance)) ~ riv_dist * log10(watershed_km2_A), filter(snv_dist_df, morphia == "m2", siteB != "PH2015_12D", siteA != "PH2015_12D", siteB != "PH2015_12U", siteA != "PH2015_12U"))
summary(fit.m2)
anova(fit.m2)

#### FIGURE #####
morphia.labels <- as_labeller(c(`all` = "Allele count = 1-4", `m234` = "Allele count = 2,3,4", `m34` = "Allele count = 3,4", `m12` = "Allele count = 1,2", `m1` = "Allele count = 1", `m2` = "Allele count = 2"))
cols <- c("snow1", brewer.pal(9, "BuGn")[5], species.colors[1])

SNV_diss_watershed.m1 <- filter(snv_dist_df, morphia == "m1", siteB != "PH2015_12D", siteA != "PH2015_12D", siteB != "PH2015_12U", siteA != "PH2015_12U") %>% 
  ggplot(., aes(x= watershed_km2_A, y= 1 - snv_distance, group= morphia)) +
  geom_point(aes(fill= riv_dist/1000), shape= 21, color= "black", size= 3) +
  geom_abline(intercept = fit.m1$coefficients["(Intercept)"],
              slope = fit.m1$coefficients["log10(watershed_km2_A)"],
              color= "black", size= 1) +
  annotate("text", x= 1.05, y= 0.97, label= "Fixed SNV sites", hjust= 0, size= 4) +
  labs(x= expression('Watershed area (km'^"2"*")"), y= "SNV site similarity (%)") +
  scale_fill_gradientn(colors= cols,
                       name= "River network\ndistance (km)") +
  scale_x_log10(limits= c(1, 2200),
                expand= c(0,0)) +
  annotation_logticks(sides= "b") +
  scale_y_continuous(limits= c(0, 1.02), 
                     expand= c(0, 0),
                     breaks= seq(0, 1, by= 0.25),
                     labels= c("0", "25", "50", "75", "100")) +
  theme_strains +
  theme(legend.position = c(0, 0.9),
        legend.justification = c(0, 1),
        legend.direction= "horizontal")
SNV_diss_watershed.m1

SNV_diss_watershed.m2 <- ggplot(filter(snv_dist_df, morphia == "m2"), aes(x= watershed_km2_A, y= 1 - snv_distance, group= morphia)) +
  geom_point(aes(fill= riv_dist/1000), shape= 21, color= "black", size= 3) +
  geom_abline(intercept = fit.m2$coefficients["(Intercept)"],
              slope = fit.m2$coefficients["log10(watershed_km2_A)"],
              color= "black", size= 1) +
  annotate("text", x= 1.05, y= 0.97, label= "Bi-allelic SNV sites", hjust= 0, size= 4) +
  labs(x= expression('Watershed area (km'^"2"*")"), y= "SNV site similarity (%)") +
  scale_fill_gradientn(colors= cols,
                       name= "River network\ndistance (km)") +
  scale_x_log10(limits= c(1, 2200),
                expand= c(0,0)) +
  annotation_logticks(sides= "b") +
  scale_y_continuous(limits= c(0, 1.02), 
                     expand= c(0, 0),
                     breaks= seq(0, 1, by= 0.25),
                     labels= c("0", "25", "50", "75", "100")) +
  theme_strains +
  theme(legend.position = c(0, 0.9),
        legend.justification = c(0, 1),
        legend.direction= "horizontal")
SNV_diss_watershed.m2
###############################################################################







#### COMBINE PANELS ###########################################################
watershedArea_combined <- plot_grid(nucDiv.watershed.plot + labs(y= "Nucleotide diversity") +  theme(legend.position= "none",
                                                                                                     axis.title.x = element_blank(),
                                                                                                     text= element_text(size= 10)),
                                    popANI_watershed + labs(y= "Mean pop. ANI (%)") + theme(axis.title.x = element_blank(), text= element_text(size= 10)),
                                    SNV_diss_watershed.m2 + theme(#legend.position= "none",
                                      axis.title.x = element_blank(),
                                      text= element_text(size= 10)),
                                    SNV_diss_watershed.m1 + theme(legend.position= "none",
                                                                  axis.title.x = element_blank(),
                                                                  text= element_text(size= 10)),
                                    nrow= 4,
                                    labels= c("A", "B", "C", "D"),
                                    align= "v")

species_legend <- get_legend(nucDiv.watershed.plot + theme(legend.position = "top",
                                                           legend.justification = "left",
                                                           legend.box.margin = margin(0, 0, 0, 0)))

watershedArea_combined2 <- plot_grid(species_legend, watershedArea_combined,
                                     nrow= 2,
                                     rel_heights = c(0.05, 1))

watershedArea_combined3 <-  annotate_figure(watershedArea_combined2, 
                                            bottom= text_grob(label= expression('Watershed area (km'^"2"*")"), vjust= -0.5))



ggsave(watershedArea_combined3, filename = "Fig_3.png", height= 180*1.25, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")




