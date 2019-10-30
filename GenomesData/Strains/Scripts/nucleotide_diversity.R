## Investigate microdiverisity (population heterogeneity) in plots 
## microdiversity = 1 - clonality
## clonality given in strainrep .log files and 
## is the probability that a read has the same bp as another read at a given site

#### Libraries #################################################################
library(tidyverse)
library(ggplot2)
################################################################################

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "Strains")
dir_output_fig <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "Strains", "Figures")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "Strains", "Output_tables")
################################################################################

### SOURCE FUNCTIONS
#source(file.path(dir_input, "strains_analysis", "R_scripts", "snv_linkage_functions.R"))
source("R_scripts/diversity_functions.R")

nd.df <- get_microdiversity()


## Read in watershed area data (km^2)

dir_watershed <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")

ws.df <- read_tsv(file.path(dir_watershed, "PhormMeta17_WatershedArea_Combined.txt")) %>% 
  rename(sample= ggkbase_id) %>% 
  select(sample, year, fork, watershed_km2)

nd.ws.df <- left_join(nd.df, ws.df)


### PLOTTING PARAMETERS ######################################################
x_axis_format_distance <- scale_x_continuous(expand= c(0, 5))
x_axis_format_window <- scale_x_discrete(labels= NULL, expand= c(0.02, 0))
y_axis_format <- scale_y_continuous(limits= c(0, 1),
                                    breaks= seq(0, 1, by= 0.25),
                                    expand= c(0.02, 0))



pid.facet.labels <- as_labeller(c(`98` = "PID = 0.98", `99` = "PID = 0.99" ))
species_match.facet.labels <- as_labeller(c(`N` = "Species Mis-match", `Y` = "Species Match" ))

facet.by.pid <- facet_grid(.~pid, labeller= labeller(pid= pid.facet.labels, scales= "free_y"))
facet.by.pid.species_match <- facet_grid(species_match~pid, 
                                         scales= "free_y",
                                         labeller= labeller(pid= pid.facet.labels, 
                                                            species_match= species_match.facet.labels))


## ggplot theme for snv linkage
theme_snv <- theme(panel.grid = element_blank(),
                   plot.margin = unit(c(1, 1, 1, 1), "cm"),
                   text = element_text(size= 14),
                   plot.background = element_rect(fill = "transparent"), # bg of the plot
                   panel.background = element_rect(fill= "transparent", color="black"),
                   axis.text = element_text(colour="black"),
                   axis.title.x = element_text(vjust = -0.75),
                   axis.title.y = element_text(vjust = 1.5),
                   legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                   legend.key = element_blank(),
                   strip.background=element_rect(fill="transparent", color="transparent"),
                   legend.position = "top")
################################################################################

# OLD PLOT, NOW THE BREADTH FILTERING IS DONE IN THE GET_MICRODIVERSITY() FUNCTION

## Microdiversity filtered by breadth > 0.6 (filtered in get_microdiversity)
ggplot(data= nd.df, aes(x= species, y= microdiversity)) +
  geom_boxplot(aes(fill= species), alpha= 0.3) +
  geom_point(position = "jitter", size= 2) +
  labs(x= "", y= expression("Nucleotide diversity ("~pi~")")) +
  scale_y_continuous(limits= c(0, .0025), expand= c(0.01, 0)) +
  scale_fill_discrete(guide= FALSE) +
  #facet.by.pid +
  theme_snv
ggsave(last_plot(), filename= "nucleotide_diversity.jpg", width= 8, height= 6, units= "in", dpi= 320, path= "Output_figures")
       



ggplot(data= nd.ws.df, aes(x= watershed_km2, y= microdiversity)) +
  geom_point(aes(color= species), size= 2) +
  stat_smooth(aes(color= species), method= "lm", se= FALSE) +
  labs(x= expression("Watershed area ("~km^2~")"), y= expression("Nucleotide diversity ("~pi~")")) +
  scale_y_continuous(limits= c(0, .0025), expand= c(0.01, 0)) +
  #scale_x_log10(limits= c(1, 10000), expand= c(0, 0)) +
  #annotation_logticks() +
  #facet_grid(.~species, scales= "free_x") +
  scale_fill_discrete(guide= FALSE) +
  theme_snv
ggsave(last_plot(), filename= "nucleotide_diversity_watershed.jpg", width= 8, height= 6, units= "in", dpi= 320, path= "Output_figures")




ggplot(data= filter(nd.ws.df, species == "1"), aes(x= watershed_km2, y= microdiversity)) +
  #geom_point(size= 2) +
  geom_point(aes(size= cov_mean)) +
  #stat_smooth(aes(color= species), method= "lm", se= FALSE) +
  labs(x= expression("Watershed area ("~km^2~")"), y= expression("Nucleotide diversity ("~pi~")"), title= "Species 1") +
  scale_y_continuous(limits= c(0, .0025), expand= c(0.01, 0)) +
  scale_x_continuous(limits= c(0, 1750), expand= c(0.02, 0), breaks= seq(0, 1750, by= 250), labels= c("0", "", "500", "", "1000", "", "1500", "")) +
 # scale_x_log10(limits= c(1, 3000), expand= c(0, 0)) +
  #annotation_logticks() +
  scale_fill_discrete(guide= FALSE) +
  theme_snv
ggsave(last_plot(), filename= "nucleotide_diversity_watershed_sp1.jpg", width= 8, height= 6, units= "in", dpi= 320, path= "Output_figures")



 ab.samples <- nd.ws.df$sample[str_detect(nd.ws.df$sample, "SFW_U|ELK_O|MFU_D|RDM_U|LGB_O|RAT_O|RUC_O")]
 ab.sample.names <- do.call(rbind, str_split(ab.samples, "_"))[, 3:4]
 
 ab.samples.df <- nd.ws.df %>% 
   filter(.$sample %in% ab.samples) %>% 
   mutate(A_B= ifelse(str_detect(.$sample, "_B"), "B", "A"),
          site=  str_c(ab.sample.names[, 1], ab.sample.names[, 2], sep="_"))
 
 ggplot(data= ab.samples.df, aes(x= A_B, y= microdiversity, group= site)) +
   geom_point(size= 2) +
   geom_line() +
   labs(x= "Sample A or B from cobble", y= expression("Nucleotide diversity ("~pi~")")) +
   scale_x_discrete(expand= c(0.2, 0)) +
   scale_y_continuous(limits= c(0, 0.0012), expand= c(0.02, 0)) +
   theme_snv
 ggsave(last_plot(), filename= "nucleotide_diversity_sampA-B_sp1.jpg", width= 4, height= 6, units= "in", dpi= 320, path= "Output_figures")
 
 
 ab.samples %in% nd.ws.df$sample 
 