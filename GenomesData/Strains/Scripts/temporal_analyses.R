## Investigate sites that were sampled in both 2015 and 2017

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
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "Strains", "Output_tables")
################################################################################

### SOURCE FUNCTIONS
source(file.path(dir_input, "strains_analysis", "R_scripts", "snv_linkage_functions.R"))
  # Returns functions: analyse_linkage_data(), filter_snv_freq_window(), make_multi_panel_fig(), summarize_log_files()

## SAMPLE METADATA
samp.md <- read_tsv(file.path(dir_input, "sample.metadata.snv.linkage.tsv"))

## MULTI-SAMPLED SITES
ms.sites <- read_tsv(file.path(dir_input, "time_replicate_genomes.tsv")) %>% 
  mutate(sample= str_replace(ms.sites$genome, "_Osc.*$", ""),
         year= as.character(year))

#### Summarize data on .log files 
log_sp1 <- summarize_log_files(path= file.path(dir_input, "species_1", "log_files"))
log_sp2 <- summarize_log_files(path= file.path(dir_input, "species_2", "log_files"))
log_sp3 <- summarize_log_files(path= file.path(dir_input, "species_3", "log_files"))

# Combine into a master data frame
log_master <- do.call(rbind, list(log_sp1, log_sp2, log_sp3)) %>%
  left_join(., samp.md, by= "sample") %>% 
  filter(pid == "98")

# Filter by sites sampled multiple times and with breadth > 0.6
log_ms.sites <- inner_join(log_master, ms.sites) %>% 
  filter(breadth > 0.6)



## GGPLOT THEME
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


## Microdiversity changes in 2015 and 2017
ggplot(data= log_ms.sites, aes(x= year, y= (1 - clonality))) +
  geom_point(aes(color= ref_id, shape= year), size= 4) +
  labs(x= "", y= "Microdiversity") +
  scale_y_continuous(limits= c(0, .0015), expand= c(0.01, 0)) +
  #scale_color_discrete(guide=FALSE) +
  scale_shape_discrete(guide=FALSE) +
  facet_grid(.~site) +
  theme_snv
ggsave(last_plot(), filename= "microdiv_time.jpg", width= 8, height= 6, units= "in", dpi= 320, path= "Output_figures")
       
       
       

