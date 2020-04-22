##


library(tidyverse)
library(RColorBrewer)
source("Scripts/ggplot_themes.R")

## NUCLEOTIDE DIVERSITY X WATERSHED AREA
source("Scripts/gene_nuc_div.R")
nucDiv.watershed.plot <- ggplot(data= gi_filt_summary, aes(x= watershed_km2, y= median_pi)) +
  geom_point(aes(fill= species, shape= species), size= 3, color= "black") +
  labs(x= expression(paste("Watershed area (", km^{2}, ")")), y= "Median nucleotide diversity") +
  #  geom_smooth(method= "lm") +
  # scale_x_continuous(breaks= seq(0, 8000, by= 1000), 
  #                    labels= c("0", "", "2000", "", "4000", "", "6000", "", "8000"),
  #                    expand= c(0.02, 0)) +
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
  theme(legend.position = c(0.92, 0.85))

nucDiv.watershed.plot



## SNV SHARING X WATERSHED AREA
#source("Scripts/SNVS_analysis.R")
snv_dist_mat <- read_tsv("inStrain/output_tables/snv_pos_M234_jaccard.tsv")

## River network distance between sites in meters
river_dist_mat <- read_tsv("../../EnvData/Distance_RiverNetwork_meters.tsv") %>% 
  rename(site= Site)

# Make river distances long
river_dist_df <- river_dist_mat %>% 
  rename(siteA= site) %>% 
  gather(key= siteB, value= riv_dist, -siteA)

# Watershed areas
watershed_area <- read_tsv("../../EnvData/PhormMeta17_WatershedArea_Combined.tsv") %>% 
  select(ggkbase_id, watershed_km2) %>% 
  mutate(siteA= ggkbase_id,
         watershed_km2_A= watershed_km2,
         siteB= ggkbase_id,
         watershed_km2_B= watershed_km2) %>% 
  select(-ggkbase_id, -watershed_km2)

# Make SNV_distance long and merge with river distance  
snv_dist_df <- snv_dist_mat %>% 
  as_tibble() %>%
  mutate(siteA= names(.)) %>% 
  gather(key= siteB, value= snv_distance, -siteA) %>% 
  left_join(., river_dist_df) %>% 
  #inner_join(., ani.dist.df) %>% 
  left_join(., select(watershed_area, siteA, watershed_km2_A)) %>% 
  left_join(., select(watershed_area, siteB, watershed_km2_B)) %>% 
  mutate(watershed_diff= abs(watershed_km2_A - watershed_km2_B))

cols <- c("snow1", brewer.pal(9, "BuGn")[5], species.colors[1])

SNV_diss_watershed <- snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  mutate(watershed_diff= abs(watershed_km2_A - watershed_km2_B)) %>% 
  ggplot(., aes(x= watershed_km2_A, y= 1 - snv_distance)) +
  geom_point(aes(fill= riv_dist/1000), shape= 21, color= "black", size= 3) +
  #geom_point(fill= species.colors[1], shape= 21, color= "black", size= 3, alpha= 0.5) +
  geom_smooth(se= FALSE, method= "lm", color= "black", size= 1) +
  labs(x= expression('Watershed area (km'^"2"*")"), y= "SNV position\nJaccard similarity") +
  scale_fill_gradientn(colors= cols,
                    name= "River network\ndistance (km)") +
  scale_x_log10(limits= c(1, 2200),
                expand= c(0,0)) +
  annotation_logticks(sides= "b") +
  scale_y_continuous(limits= c(0, 1.02), 
                     expand= c(0, 0),
                     breaks= seq(0, 1, by= 0.25),
                     labels= c("0.00", "0.25", "0.50", "0.75", "1.00")) +
                     #labels= c("  0.00", "  0.25", "  0.50", "  0.75", "  1.00")) +
  theme_strains +
  theme(legend.position = "top")#,
        #legend.direction = "horizontal",
        #legend.justification = "left")

SNV_diss_watershed

fit.SNVdiss <- lm(sqrt((1-snv_distance)) ~ watershed_km2_A, filter(snv_dist_df, siteA != siteB))
summary(fit.SNVdiss)
plot(fit.SNVdiss)

hist(sqrt(1 - filter(snv_dist_df, siteA != siteB)$snv_distance))


## popANI X WATERSHED AREA
source("Scripts/ANI_scaffold.R")
ani_sum <- read_tsv("Output_tables/ani_summary.tsv")
## Watershed area data
dir_input_watershed <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/EnvData"
watershed.area <- read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
  # rename(site= ggkbase_id) %>% 
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
  left_join(., watershed.area, by= c("name1" = "ggkbase_id")) %>% 
  left_join(., watershed.area, by= c("name2" = "ggkbase_id")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y)

popANI.sum <- popANI.long %>% 
  group_by(name1, watershed_1) %>% 
  summarise(meanPOPani= mean(mean_popANI, na.rm=TRUE)) %>% 
  ungroup()
popANI.sum                  

popANI_watershed <- ggplot(data= popANI.sum, aes(x= watershed_1, y= meanPOPani*100)) +
  geom_point(shape= 21, fill= species.colors[1], color= "black", size= 3) +
  labs(x= expression("Watershed area (km"^2*")"), y= "Mean population ANI (%)") +
  scale_x_log10(limits=c(1, 2200),
                expand= c(0, 0)) +
  annotation_logticks() +
  scale_y_continuous(breaks= seq(99.4, 99.8, by= 0.2),
                     labels= c("99.4", "99.6", "99.8")) +
                     #labels= c(" 0.994", " 0.996", " 0.998")) +
  theme_strains

popANI_watershed

ggsave(popANI_watershed, filename = "popANI_watershed.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



watershedArea_combined <- plot_grid(nucDiv.watershed.plot + labs(y= "Nucleotide diversity") +  theme(legend.position= "none",
                                                                                                     axis.title.x = element_blank(),
                                                                                                     text= element_text(size= 10)),
                               SNV_diss_watershed + labs(y= "SNV similarity") + theme(axis.title.x = element_blank(),
                                                                                         text= element_text(size= 10),
                                                                                      legend.position = "none"), 
                               popANI_watershed + labs(y= "Population ANI (%)") + theme(text= element_text(size= 10)),
                               nrow= 3,
                               labels= c("A", "B", "C"),
                               align= "v")
                               # align= "h", 
                               # axis= "l")


species_legend <- get_legend(nucDiv.watershed.plot + theme(legend.position = "top",
                                                           legend.justification = "center",
                                                           legend.box.margin = margin(0, 0, 0, 0)))

rivDist_legend <- get_legend(SNV_diss_watershed + theme(legend.position = "top",
                                                           legend.justification = "center",
                                                           legend.box.margin = margin(0, 0, 0, 0)))

legends <- plot_grid(species_legend, rivDist_legend,
                     ncol= 2)

watershedArea_combined2 <- plot_grid(legends, watershedArea_combined,
                                     nrow= 2,
                                     rel_heights = c(0.05, 1))

watershedArea_combined2



ggsave(watershedArea_combined2, filename = "watershedArea_combined.png", height= 180, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



# ggplot(data= popANI.long, aes(x= reorder(name1, watershed_1), y= mean_popANI)) +
#   geom_point(color="white") +
#   geom_boxplot(aes(fill= log10(watershed_1))) +
#   labs(x= "Sample", y= "Mean population ANI") +
#   scale_fill_viridis_c(breaks= c(1, 2, 3), 
#                        labels= c(10, 100, 1000),
#                        name= expression("Watershed area (km"^2*")")) +
#   theme_strains +
#   #scale_x_continuous(labels= as.character(sort(unique(popANI.long$watershed_1)))) +
#   theme(#axis.text.x = element_text(angle= 90, vjust= 0.5),
#     axis.text.x= element_blank(),
#     legend.position = "top")
# ggsave(last_plot(), filename = "popANI_watershed.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
#        path= "Output_figures")
# 
# 
# ggplot(data= popANI.long, aes(x= log10(watershed_1), y= mean_popANI)) +
#   geom_point() +
#   geom_smooth(method= "lm")
# #  scale_x_log10() +
# # annotation_logticks() +
# scale_fill_discrete(guide=FALSE) +
#   theme_strains
# 
# 

# 
# 
# fit.1 <- lm(mean_popANI ~ log10(watershed_1), popANI.long)
# plot(mean_popANI ~ log10(watershed_1), popANI.long)
# summary(fit.1)
# anova(fit.1)
# plot(fit.1)
# 
# hist(sqrt(popANI.long$mean_popANI))
# 
# fit.2 <- lm(meanPOPani ~ watershed_1, popANI.sum)
# summary(fit.2)
# anova(fit.2)
# plot(fit.2)
# 






