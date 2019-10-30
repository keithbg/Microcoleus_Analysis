## Comparisons between linkage decay and other parameters
library(tidyverse)
library(ggplot2)

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



## REad in linkage decay parameters (ldp)
ldp <- read_tsv("Output_tables/nls_parameters.tsv") %>% 
  filter(sample != "PH2015_01D_0.98") %>% 
  mutate(sample= str_replace(sample, "_0.98$", ""),
         species= str_replace(species, "sp", "")) %>% 
  select(-std.error) %>% 
  spread(term, estimate)

## Get microdiveristy data
#source(file.path(dir_input, "strains_analysis", "R_scripts", "snv_linkage_functions.R"))
source("R_scripts/diversity_functions.R")

md.df <- get_microdiversity() 
md.df.trim <- md.df %>%
  select(sample, species, pid, species_present, multiple_species, species_match, microdiversity)


md.ldp <- full_join(md.df.trim, ldp)


ggplot(md.ldp) +
  geom_hline(yintercept= 0) +
  geom_point(aes(x= microdiversity, y= Asym, color= species)) +
  theme_snv

ggplot(md.ldp) +
  geom_point(aes(x= microdiversity, y= R0, color= species)) +
  theme_snv


ggplot(md.ldp) +
  geom_point(aes(x= microdiversity, y= lrc)) +
  theme_snv
