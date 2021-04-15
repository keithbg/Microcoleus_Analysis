## Figure S4

## Relationship of environmental differences and nucleotide identities among samples

#### Libraries #################################################################
library(tidyverse)
source("Scripts/ggplot_themes.R")
################################################################################

## Function to reformat sample names for comparisons between ANI and environmental data
extract_envID <- function(x){
  if(str_detect(x, "PH2017")){
    #out <- str_split(x, "_")[[1]][3]
    out <- str_replace(str_split(x, "[0-9]_")[[1]][3], "_A|_B", "")
  } else {
    out <- str_split(x, "_")[[1]][2]
  }
  #if(str_detect(x, "PH2015")){
  #  out <- str_split(x, "_")[[1]][2]
  #}
}

#### IMPORT DATA ####
## ANI data
ani <- read_tsv("Data/inStrain_data/ani_summary_v1.4.tsv") %>% # ani_summary.tsv generated in ANI_scaffold_data.R
  mutate(envID1= do.call(rbind, map(name1, extract_envID))[, 1],
         envID2= do.call(rbind, map(name2, extract_envID))[, 1],
         compID= str_c(envID1, envID2, sep= "-"))

## PH2015 environmental data
env15 <- read_tsv(file.path("Data/Env_data", "PH2015_env_data.tsv")) %>% 
  mutate(NPOC_mgL= NPOC_mgL*1000,
         cond= cond/1000,
         canopy_avg= canopy_avg*100) %>% 
  rename(cond_ms= cond,do_mgL= DO_mgL, canopy_cover_percent= canopy_avg, DOC_ugL= NPOC_mgL, temp= temp_c, alk= Alk) %>% 
  dplyr::select(-original_sampleID)

## PH2017 environmental data
env17 <- read_tsv(file.path("Data/Env_data", "PH2017_env_data_formated.tsv"))

## NorWest Temperature data
norwest <- read_tsv(file.path("Data/Env_data", "temps_NorWest.tsv"))

## Combine all temperature data
env <- full_join(env15, env17)

tenv <- env %>% 
  left_join(., norwest, by= c("biotite_ID" = "site")) %>% 
  filter(!is.na(biotite_ID)) %>% 
  dplyr::select(biotite_ID, temp, alk, pH, do_mgL, mean_vX, TDP_ugL, NO3_ugL, NH4_ugL, DOC_ugL, TDN_ugL, canopy_cover_percent, cond_ms, temp_NorWest) %>% 
  arrange(biotite_ID) %>% 
  mutate(envID= do.call(rbind, map(biotite_ID, extract_envID))[, 1]) %>% 
  dplyr::select(biotite_ID, envID, everything())


## Calculate differences between environmental values
calc_diffs <- function(vec){
  df <- as_tibble(outer(vec, vec, `-`), .name_repair= "minimal")  
  colnames(df) <- tenv$envID
  df$envID1 <- tenv$envID
  
  # df.l <- pivot_longer(df, names_to= "name2", values_to = colname, contains("PH")) %>% 
  #   mutate(compID= str_c(name1, name2, sep= "-"))
  # 
  df.l <- pivot_longer(df, names_to= "envID2", values_to = "diff", -envID1) %>% 
    mutate(compID= str_c(envID1, envID2, sep= "-"),
           diff= abs(diff)) # absolue value of differences
}

comp_df <- map(tenv[, c(3:15)], calc_diffs) %>% 
  bind_rows(., .id="metric")  


## Join with ANI data
ani.env <- left_join(ani, comp_df) %>% 
  mutate(year= ifelse(str_detect(name1, "PH2015") & str_detect(name2, "PH2015"), "15-15",
                      ifelse(str_detect(name1, "PH2015") & str_detect(name2, "PH2017"), "15-17",
                             ifelse(str_detect(name1, "PH2017") & str_detect(name2, "PH2015"), "15-17", "17-17"))))
ani.env.plotting <- filter(ani.env, metric != "do_mgL" & metric != "mean_vX" & metric != "pH" & metric != "temp")


#### FIGURES ####
env.labels <- as_labeller(c(`alk` = "Alkalinity", `canopy_cover_percent` = "Canopy cover", `cond_ms` = "Conductivity", 
                            `DOC_ugL` = "DOC", `NH4_ugL` = "Ammonium", `NO3_ugL` = "Nitrate", 
                            `TDN_ugL` = "TDN", `TDP_ugL` = "TDP", `temp_NorWest` = "Mean temperature"))


env.watershed <- ggplot(ani.env.plotting, aes(x= watershed_diff + 0.1, y= diff)) +
  geom_point(aes(fill= year), pch= 21, color= "black", alpha= 0.3) +
  scale_fill_manual(values= year.fill.colors, name= "Year comparison") +
  scale_color_manual(values= year.colors, name= "Year comparison") +
  labs(x= expression("Watershed difference (km"^2*")"), y= "Pairwise parameter difference") +  
  geom_smooth(aes(color= year), size= 2, method= "lm", se= FALSE) +
  scale_x_log10() +
  annotation_logticks(side= "b") +
  facet_wrap(~metric, ncol= 2, scales= "free_y", labeller = labeller(metric= env.labels)) +
  theme_strains +
  theme(axis.text.x= element_text(angle= 45, vjust= 0.9, hjust= 0.9),
        legend.position = "top",
        legend.justification = c(0,0),
        legend.box.margin = margin(0, 0, 0, 0, unit= "pt"))


env.conANI <- ggplot(filter(ani.env.plotting, metric != "alk",  metric != "DOC_ugL", metric != "NH4_ugL", metric != "NO3_ugL"),
                     aes(x= diff, y= mean_conANI)) +
  geom_point(aes(fill= year), pch= 21, color= "black", alpha= 0.3) +
  scale_fill_manual(values= year.fill.colors, name= "Year comparison") +
  scale_color_manual(values= year.colors, name= "Year comparison") +
  #scale_y_continuous(breaks= seq(0.9875, 1, by=0.0025), labels= c("", "0.990", "", "0.995", "", "1.000")) +
  scale_y_continuous(breaks= seq(0.99, 1, by=0.001), labels= c("99.0", "", "99.2", "", "99.4", "", "99.6", "", "99.8", "", "100")) +
  scale_x_continuous(expand= c(0.04,0)) +
  labs(x= "Pairwise parameter difference", y= "Consensus ANI (%)") +  
  geom_smooth(aes(color= year),size= 2, method= "lm", se= FALSE) +
  facet_wrap(~metric, nrow= 5, scales= "free_x", labeller= labeller(metric= env.labels)) +
  theme_strains +
  theme(axis.text.x= element_text(angle= 45, vjust= 0.9, hjust= 0.9),
        legend.position = "top")


env.popANI <- ggplot(filter(ani.env.plotting, metric != "alk",  metric != "DOC_ugL", metric != "NH4_ugL", metric != "NO3_ugL"),
                     aes(x= diff, y= mean_popANI)) +
  geom_point(aes(fill= year), pch= 21, color= "black", alpha= 0.3) +
  scale_fill_manual(values= year.fill.colors, name= "Year comparison") +
  scale_color_manual(values= year.colors, name= "Year comparison") +
  labs(x= "Pairwise parameter difference", y= "Population ANI (%)") +  
  geom_smooth(aes(color= year), size= 2, method= "lm", se= FALSE) +
  scale_y_continuous(breaks= seq(0.993, 1, by=0.001), labels= c("", "99.4", "", "99.6", "", "99.8", "", "100")) +
  scale_x_continuous(expand= c(0.04,0)) +
  facet_wrap(~metric, nrow= 5, scales= "free_x", labeller= labeller(metric= env.labels)) +
  theme_strains +
  theme(axis.text.x= element_text(angle= 45, vjust= 0.9, hjust= 0.9),
        legend.position = "top")


year.legend <- get_legend(env.watershed)

text.size <- 8
env.plot.combined <- plot_grid(env.watershed + theme(legend.position = "none",
                                                     text= element_text(size= text.size)), 
                               env.conANI + theme(legend.position = "none",
                                                  text= element_text(size= text.size)), 
                               env.popANI + theme(legend.position = "none",
                                                  text= element_text(size= text.size)),
                               ncol= 3,
                               rel_widths = c(1.2, 0.8, 0.8),
                               labels= c("A", "B", "C"))


env.plot.combined2 <- plot_grid(year.legend, env.plot.combined, 
                                nrow=2,
                                rel_heights= c(0.1, 1))

ggsave(env.plot.combined2, filename = "Fig_S4_v1.4.png", dpi= 320, height= 180*1.2, width= 180, units= "mm",
       path= "Output_figures")





