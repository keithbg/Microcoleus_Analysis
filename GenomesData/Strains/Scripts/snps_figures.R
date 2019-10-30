## Make figures for microdiversity and SNPS / Mbp
## Usess functions in: R_scripts/diversity_functions.R

source("R_scripts/diversity_functions.R")

## CALCULATE the SNP N:S RATIO

# Get files
NS_files_sp1 <- as.list(list.files("species_1", "*98_aa.freq$"))
NS_files_sp2 <- as.list(list.files("species_2", "*98_aa.freq$"))
NS_files_sp3 <- as.list(list.files("species_3", "*98_aa.freq$"))

# Calc NS ratio for each species
ns_sp1 <- calc_NS_ratio(file.list= NS_files_sp1, species= "species_1") %>% 
  mutate(species= "1")
ns_sp2 <- calc_NS_ratio(file.list= NS_files_sp2, species= "species_2") %>% 
  mutate(species= "2")
ns_sp3 <- calc_NS_ratio(file.list= NS_files_sp3, species= "species_3") %>% 
  mutate(species= "3")

# Merge files together
ns.df <- full_join(ns_sp1, ns_sp2) %>% 
  full_join(., ns_sp3) %>% 
  as_tibble()

## CALCULATE SNPS PER MEGA-BP
snps.df <- calc_snps_mbp()


## GET MICRODIVERSITY
microdiv <- get_microdiversity() %>% 
  mutate(sample= str_c(.$sample, "_0.98"))


# snps.ns.df <- left_join(snps.df, ns.df) %>% 
#   mutate(year= ifelse(str_detect(.$sample, "PH2015"), "2015", "2017"))

snps.ns.micro.df <- left_join(snps.df, ns.df) %>% 
  mutate(year= ifelse(str_detect(.$sample, "PH2015"), "2015", "2017")) %>% 
  inner_join(., select(microdiv, sample, species, microdiversity))





## ggplot themes
theme_strains <- theme(panel.grid = element_blank(),
                       plot.margin = unit(c(1, 1, 1, 1), "cm"),
                       text = element_text(size= 14),
                       plot.background = element_rect(fill = "transparent", color= "transparent"), # bg of the plot
                       panel.background = element_rect(fill= "transparent", color= "transparent"),
                       panel.border= element_rect(fill= NA, color= "black", linetype= "solid", size= 1),
                       panel.ontop = TRUE,
                       axis.text = element_text(colour="black"),
                       axis.title.x = element_text(vjust = -0.75),
                       axis.title.y = element_text(vjust = 1.5),
                       legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                       legend.key = element_blank(),
                       strip.background = element_rect(fill="transparent", color= "transparent"),
                       #axis.text.x = element_text(angle= 45, hjust= 1),
                       legend.position = "top")



ggplot(data= snps.ns.micro.df, aes(x= snps_mbp, y= microdiversity )) +
  geom_point(aes(color= year), size= 2) +
  labs(x= "SNPs / Mbp", y= expression("Nucleotide diversity ("~pi~")")) +
  facet_grid(.~species, scales= "free_x") +
  theme_strains
#ggsave(last_plot(), filename = "Output_figures/snps_mbp_microdiversity.pdf", height= 6, width= 8, units= "in", device= cairo_pdf)
ggsave(last_plot(), filename= "snps_mbp_microdiversity.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)



ggplot(data= snps.ns.micro.df, aes(x= snps_mbp, y= NS_ratio )) +
  geom_vline(xintercept= 1500, size= 0.5, color= "gray50") +
  geom_point(aes(color= year), size= 2) +
  labs(x= "SNPs / Mbp", y= "N:S ratio") +
  facet_grid(.~species) +
  theme_strains
#ggsave(last_plot(), filename = "Output_figures/snps_NS.pdf", height= 6, width= 8, units= "in", device= cairo_pdf)
ggsave(last_plot(), filename= "snps_NS.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)


ggplot(data= snps.ns.micro.df, aes(x= microdiversity, y= NS_ratio )) +
  geom_point(aes(color= year), size= 2) +
  labs(x= expression("Nucleotide diversity ("~pi~")"), y= "N:S ratio") +
  #scale_x_continuous(limits= c(0, 0.0035)) +
  facet_grid(.~species) +
  theme_strains
#ggsave(last_plot(), filename = "Output_figures/NS_microdiversity.pdf", height= 6, width= 8, units= "in", device= cairo_pdf)
ggsave(last_plot(), filename= "NS_microdiversity.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)




# source(file.path("strains_analysis", "R_scripts", "snv_linkage_functions.R"))
# # Returns functions: analyse_linkage_data(), filter_snv_freq_window(), make_multi_panel_fig(), summarize_log_files()
# 
# ## Sample metadata
# samp.md <- read_tsv(file.path("sample.metadata.snv.linkage.tsv"))
# 
# ## Summarize data on .log files 
# log_sp1 <- summarize_log_files(path= file.path(dir_input, "species_1", "log_files"))
# log_sp2 <- summarize_log_files(path= file.path(dir_input, "species_2", "log_files"))
# log_sp3 <- summarize_log_files(path= file.path(dir_input, "species_3", "log_files"))
# 
# ## Combine into a master data frame
# log_master <- do.call(rbind, list(log_sp1, log_sp2, log_sp3)) %>%
#   left_join(., samp.md, by= "sample") %>% 
#   filter(pid == "98") %>% 
#   mutate(sample= str_c(.$sample, "_0.98"),
#          species= str_replace(.$ref_id, "^[a-z]+_", ""),
#          microdiversity= 1-clonality)


# ## MERGE MICRODIVERSITY AND SNPS DATA
# snps.ns.micro.df <- left_join(snps.ns.df, select(log_master, sample, species, microdiversity)) %>% 
#   filter(sample != "PH2015_01U_0.98" & sample != "PH2015_01D_0.98") %>% 
#   filter(!(sample == "PH2015_07D_0.98" & species == "1")) %>% 
#   filter(!(sample == "PH2015_06S_0.98" & species == "1")) %>% 
#   filter(!(sample == "PH2015_08D_0.98" & species == "1")) %>% 
#   filter(!(sample == "PH2015_12U_0.98" & species == "1")) %>% 
#   filter(!(sample == "PH2017_06_SFM_O_A_0.98" & species == "3")) %>% 
#   filter(!(sample == "PH2017_07_MST_O_A_0.98" & species == "3"))


