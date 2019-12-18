# Investigate nonsynonymous and synonymous SNV results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.SNP_mutation_types.tsv

### The SNP_mutation_types.tsv file filters SNVs from the SNVs.tsv file by morphia == 2 and cryptic == FALSE. 
## This is why there are fewer SNVs in the SNP_mutation_types.tsv file compared to the SNVs.tsv file

## Input df generated in format_inStrain_output.R

## Libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)


#### INPUT FILES ####
snv_df <- read_tsv("inStrain/output_tables/snp_mutation_type_df.tsv")
dir_input_map <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","EnvData")

## Read in Lat Longs
latlong <- read_csv(file.path(dir_input_map, "PhormMeta17_LatLong_combined.csv")) %>%
  mutate(year= as.character(year))

## Watershed area data
dir_input_watershed <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/EnvData"
watershed.area <-
  read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2)


# Get genome sizes
genome.size <- read_delim("genome_lengths.txt", delim= " ", col_names= FALSE) %>% 
  slice(c(13:16, 21:22))

genome_size_df <- tibble(species= c("species_1", "species_2", "species_3"),
                         genome_mbp= c(as.numeric(genome.size[2, ]) / 1000000, 
                                       as.numeric(genome.size[4, ]) / 1000000, 
                                       as.numeric(genome.size[6, ]) / 1000000))

#### CALCULATE SNVs per MBP ####
snvs_mbp_df <- snv_df %>% 
  left_join(., genome_size_df) %>% 
  group_by(sample, species, genome_mbp) %>% 
  summarize(SNVs= length(mutation_type),
            cov_mean= mean(baseCoverage, na.rm= TRUE),
            cov_sd= sd(baseCoverage, na.rm= TRUE)) %>%
  ungroup() %>% 
  mutate(SNV_mbp= SNVs / genome_mbp)



#### CALCULATE N:S RATIOS ####
# for both genome and individual genes
NS_genome_ratios <- snv_df %>% 
  group_by(sample, species) %>% 
  summarize(
    SNV_count= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N")),
    S_count= sum(str_detect(mutation_type, "S")),
    I_count= sum(str_detect(mutation_type, "I")),
    M_count= sum(str_detect(mutation_type, "M")),
    NS= N_count/S_count) 


NS_gene_ratios <- snv_df %>% 
  filter(mutation_type == "N" & mutation_type == "S") %>% 
  group_by(sample, gene) %>% 
  summarize(
    n= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N")),
    S_count= sum(str_detect(mutation_type, "S")),
    NS= N_count/S_count) 




## JOIN THE SNVS PER MBP AND NS RATIO AND LAT/LONGAND WATERSHED DATA
snvs_genome_df <- left_join(snvs_mbp_df, NS_genome_ratios, by= c("sample", "species")) %>% 
  mutate(pop_age= ifelse(SNV_mbp < 1550, "Young", "Old"),
         pop_age_binary= ifelse(SNV_mbp < 1550, 1, 0),
         ggkbase_id= str_replace(sample, "\\.species.*$", "")) %>% 
  left_join(., latlong, by= "ggkbase_id") %>% 
  left_join(., watershed.area, by= "ggkbase_id")
  
  
  



#### MAKE FIGURES ####
source("Scripts/ggplot_themes.R")


## SNVs_mbp X NS ratio
ggplot(data= snvs_genome_df, aes(x= SNV_mbp, y= NS)) +
  geom_vline(xintercept= 1550, color= "gray50", size= 0.5, linetype= "dashed") +
  #geom_point(aes(fill= species, shape= species), size= 4, color= "gray40") +
  geom_point(aes(fill= species, shape= species, size= cov_mean), color= "gray40") +
  geom_text(aes(x= 750, y= 2.4, label = "Young \npopulations", hjust= "middle")) +
  geom_text(aes(x= 2500, y= 2.4, label = "Old \npopulations", hjust= "middle")) +
  labs(x= "SNVs / mbp", y= "N:S ratio") +
  scale_x_continuous(limits= c(0, 8500),
                     breaks= seq(0, 8000, by= 1000), 
                     labels= c("0", "", "2000", "", "4000", "", "6000", "", "8000"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0, 2.5),
                     expand= c(0, 0))+
  scale_fill_manual(values= species.colors,
                     labels= c("1", "2", "3"),
                     name= "Species") +
  scale_shape_manual(values= species.shapes,
                     labels= c("1", "2", "3"),
                     name= "Species") +
  scale_size_continuous(name= "Mean site \ncoverage",
                        breaks= c(50, 100, 200, 300, 400)) +
  theme_strains
ggsave(last_plot(), filename = "snv_mbp.pdf", height= 180*0.75, width= 180, units= "mm", device= cairo_pdf,
       path= "Output_figures")
ggsave(last_plot(), filename = "snv_mbp.jpg", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

table(snvs_genome_df$pop_age, snvs_genome_df$watershed_km2 < 500, snvs_genome_df$species)
table(snvs_genome_df$species)

## SNVs_mbp X Watershed area
ggplot(data= snvs_genome_df, aes(x= watershed_km2, y= SNV_mbp)) +
  geom_hline(yintercept= 1550, color= "gray50", size= 0.5, linetype= "dashed") +
  geom_point(aes(fill= species, size= cov_mean, shape= species), color= "gray40") +
  geom_smooth(se= FALSE) +
  geom_text(aes(x= 1.1, y= 1000, label = "Young \npopulations", hjust= "left")) +
  geom_text(aes(x= 1.1, y= 2000, label = "Old \npopulations", hjust= "left")) +
  labs(x= expression('Watershed area (km'^"2"*")"), y= "SNVs mbp") +
  scale_x_log10(limits= c(1, 10000),
                expand= c(0, 0)) +
  annotation_logticks(sides= "b") +
  scale_fill_manual(values= species.colors,
                    labels= c("1", "2", "3"),
                    name= "Species") +
  scale_shape_manual(values= species.shapes,
                     labels= c("1", "2", "3"),
                     name= "Species") +
  scale_size_continuous(name= "Mean site \ncoverage",
                        breaks= c(50, 100, 200, 300, 400)) +
  theme_strains

ggsave(last_plot(), filename = "snv_mbp_km2.pdf", height= 180*0.75, width= 180, units= "mm", device= cairo_pdf,
       path= "Output_figures")
ggsave(last_plot(), filename = "snv_mbp_km2.jpg", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



ggplot(data= snvs_genome_df, aes(x= pop_age, y= watershed_km2)) +
geom_boxplot() +
  scale_y_log10() +
  theme_strains

ggplot(data= filter(snvs_genome_df, species == "species_1"), aes(x= watershed_km2, y= pop_age_binary)) +
  geom_point(aes(size= cov_mean))+
  geom_smooth(se= FALSE) +
  labs(x= expression('Watershed area (km'^"2"*")"), y= "Population age") +
  scale_y_continuous(limits= c(0, 1),
                     breaks= c(0, 1), 
                     labels= c("Old", "Young")) +
  scale_x_log10(limits= c(1, 2500),
                expand= c(0, 0)) +
  annotation_logticks(sides = "b") +
  theme_strains

ggplot(data= filter(snvs_genome_df, species == "species_1"), aes(x= watershed_km2)) +
  geom_histogram(binwidth= 10) +
  facet_wrap(~pop_age, nrow= 2) +
  theme_strains





#### MAKE MAP



## Make base map of Eel and Russian River watersheds
dir_input_script <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis", "Scripts_cyano_metagenomics_2017")
source(file.path(dir_input_script, "Map_eel_russian.R"))
# Returns R object: PH2017_eel_russian_base_map

## Remove rows with duplicate sites to reduce the number of points
latlong.map <- latlong %>%
  distinct(acronym, .keep_all = TRUE)

PH2017_map_theme <- theme(text= element_text(size= 14),
                          panel.background = element_rect(fill= "light blue"),
                          panel.border = element_rect(color= "black", fill= NA),
                          legend.key= element_rect(fill= "transparent"),
                          plot.background = element_rect(fill= "transparent", color= NA),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank())


## Add information to base map
PH2017_eel_russian_base_map +
  geom_point(data= snvs_genome_df, aes(x= long, y= lat, fill= pop_age), color= "black", size= 3, shape= 21, alpha= 0.5) +
  # geom_label_repel(data= latlong.map, aes(x= long, y= lat, label= acronym)) +
  #annotate("text", x= -123.3, y= 40.35, label= "Eel River", angle= 310, size= 4) +
  #annotate("text", x= -122.7, y= 38.9, label= "Russian River", angle= 310, size= 4) +
  #scale_shape_manual(values= c(21, 23, 24)) +
  scale_fill_manual(values= c("black", "red")) +
  #facet_wrap(~species) +
  PH2017_map_theme

#ggsave(last_plot(), filename= "Map_sample_locations.pdf", width= 6, height= 8, units= "in",
#       path= dir_output_fig)



PH2017_eel_russian_base_map +
  geom_point(data= filter(snvs_genome_df, species == "species_1"), aes(x= long, y= lat, fill= pop_age), color= "black", size= 3, shape= 21, alpha= 0.5) +
  scale_fill_manual(values= c("black", "red")) +
  #facet_wrap(~species) +
  PH2017_map_theme

PH2017_eel_russian_base_map +
  geom_point(data= filter(snvs_genome_df, species == "species_2"), aes(x= long, y= lat, fill= pop_age), color= "black", size= 3, shape= 21, alpha= 0.5) +
  scale_fill_manual(values= c("black", "red")) +
  #facet_wrap(~species) +
  PH2017_map_theme

PH2017_eel_russian_base_map +
  geom_point(data= filter(snvs_genome_df, species == "species_3"), aes(x= long, y= lat, fill= pop_age), color= "black", size= 3, shape= 21, alpha= 0.5) +
  scale_fill_manual(values= c("black", "red")) +
  #facet_wrap(~species) +
  PH2017_map_theme





