



library(tidyverse)
library(ggplot2)

y=x
x=y[1]
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

#str_replace(str_split(x, "[0-9]_")[[1]][3], "_A", "")


ani <- read_tsv("Output_tables/ani_summary.tsv") %>% 
  #mutate(name1_acronym= do.call(rbind, map(ani$name1, extract_envID)),
  #       name2_acronym= do.call(rbind, map(ani$name2, extract_envID))) %>% 
  # mutate(envID1= str_replace(name1, "_A|_B$", ""),
  #        envID2= str_replace(name2, "_A|_B$", ""),
  #        compID= str_c(envID1, envID2, sep= "-"))
  mutate(envID1= do.call(rbind, map(name1, extract_envID))[, 1],
         envID2= do.call(rbind, map(name2, extract_envID))[, 1],
         compID= str_c(envID1, envID2, sep= "-"))

## PH2015 environmental data
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","EnvData")
env15 <- read_tsv(file.path(dir_input, "PH2015_env_data.tsv")) %>% 
  mutate(NPOC_mgL= NPOC_mgL*1000,
         cond= cond/1000,
         canopy_avg= canopy_avg*100) %>% 
  rename(cond_ms= cond,do_mgL= DO_mgL, canopy_cover_percent= canopy_avg, DOC_ugL= NPOC_mgL, temp= temp_c, alk= Alk) %>% 
  select(-original_sampleID)

## PH2017 Environmental data
env17 <- read_tsv(file.path(dir_input, "env_data_formated.tsv"))

env <- full_join(env15, env17)

tenv <- env %>% 
  filter(!is.na(biotite_ID)) %>% 
  select(biotite_ID, temp, alk, pH, do_mgL, mean_vX, TDP_ugL, NO3_ugL, NH4_ugL, DOC_ugL, TDN_ugL, canopy_cover_percent, cond_ms) %>% 
  arrange(biotite_ID) %>% 
  mutate(envID= do.call(rbind, map(biotite_ID, extract_envID))[, 1]) %>% 
  select(biotite_ID, envID, everything())


## Calculate differences between environmental values
calc_diffs <- function(vec){
  df <- as_tibble(outer(vec, vec, `-`))  
  colnames(df) <- tenv$envID
  df$envID1 <- tenv$envID

  # df.l <- pivot_longer(df, names_to= "name2", values_to = colname, contains("PH")) %>% 
  #   mutate(compID= str_c(name1, name2, sep= "-"))
  # 
  df.l <- pivot_longer(df, names_to= "envID2", values_to = "diff", -envID1) %>% 
    mutate(compID= str_c(envID1, envID2, sep= "-"),
           diff= abs(diff)) # absolue value of differences
}

comp_df <- map(tenv[, c(3:14)], calc_diffs) %>% 
  bind_rows(., .id="metric")  


## Join with ANI data
ani.env <- left_join(ani, comp_df) %>% 
  mutate(year= ifelse(str_detect(name1, "PH2015") & str_detect(name2, "PH2015"), "15-15",
                      ifelse(str_detect(name1, "PH2015") & str_detect(name2, "PH2017"), "15-17",
                             ifelse(str_detect(name1, "PH2017") & str_detect(name2, "PH2015"), "15-17", "17-17"))))

## Make plots
source("Scripts/ggplot_themes.R")


ggplot(filter(ani.env, metric != "do_mgL" & metric != "mean_vX" & metric != "pH"), aes(x= watershed_diff + 0.1, y= diff)) +
  geom_point(aes(color= year), alpha= 0.5) +
  scale_color_startrek() +
  labs(x= expression("Watershed difference (km"^2*")"), y= "Pairwise difference") +  
  geom_smooth(method= "lm", color= "black") +
  scale_x_log10() +
  annotation_logticks(side= "b") +
  facet_wrap(~metric, nrow= 4, scales= "free_y") +
  theme_strains +
  theme(axis.text.x= element_text(angle= 45, vjust= 0.9, hjust= 0.9),
        legend.position = "top")
ggsave(last_plot(), filename = "watershed_env_differences.png", dpi= 320, height= 180, width= 180, units= "mm",
       path= "Output_figures")

ggplot(filter(ani.env, diff > 0), aes(x= diff, y= mean_popANI)) +
  geom_point(aes(color= year), alpha= 0.5) +
  scale_color_startrek() +
  geom_smooth(method= "lm", color= "black") +
  facet_wrap(~metric, nrow= 4, scales= "free_x") +
  theme_strains

ggplot(filter(ani.env, diff > 0), aes(x= diff, y= mean_conANI)) +
  geom_point(aes(color= year), alpha= 0.5) +
  scale_color_startrek() +
  geom_smooth(method= "lm", color= "black") +
  facet_wrap(~metric, nrow= 4, scales= "free_x") +
  theme_strains



ggplot(ani.env, aes(x= diff, y= mean_conANI)) +
  geom_point() +
  geom_smooth(method= "lm") +
  facet_wrap(~metric, nrow= 4, scales= "free_x") +
  theme_bw()

ggplot(filter(ani.env, riv_dist < 25000), aes(x= diff, y= mean_conANI)) +
  geom_point() +
  geom_smooth(method= "lm") +
  facet_wrap(~metric, nrow= 4, scales= "free_x") +
  theme_bw()

ggplot(filter(ani.env, riv_dist < 25000 & metric == "cond_ms"), aes(x= riv_dist, y= mean_conANI, group= metric)) +
  geom_point(shape= 21, aes(fill= diff, size= diff)) +
  scale_fill_viridis_c() +
  facet_wrap(~metric, nrow= 4, scales= "free_x") +
  theme_strains

ggplot(filter(ani.env, riv_dist < 25000 & metric == "temp"), aes(x= riv_dist, y= mean_conANI, group= metric)) +
  geom_point(shape= 21, aes(fill= diff, size= diff)) +
  scale_fill_viridis_c() +
  facet_wrap(~metric, nrow= 4, scales= "free_x") +
  theme_strains

ggplot(filter(ani.env, riv_dist < 25000 & metric == "canopy_cover_percent"), aes(x= riv_dist, y= mean_conANI, group= metric)) +
  geom_point(shape= 21, aes(fill= diff, size= diff)) +
  scale_fill_viridis_c() +
  facet_wrap(~metric, nrow= 4, scales= "free_x") +
  theme_strains



#### STATISTICS
fitCon.cond <- lm(log(mean_conANI) ~ diff, filter(ani.env, metric == "cond_ms"))
summary(fitCon.cond)
anova(fitCon.cond)
plot(fitCon.cond)


fitKM2.ph <- lm(sqrt(diff) ~ watershed_diff, filter(ani.env, metric == "pH"))
summary(fitKM2.ph)
plot(fitCon.cond)

fitKM2.temp <- lm(sqrt(diff) ~ watershed_diff, filter(ani.env, metric == "temp"))
summary(fitKM2.temp)
plot(fitKM2.temp)

fitKM2.temp <- lm(sqrt(diff) ~ watershed_diff, filter(ani.env, metric == "cond_ms"))
summary(fitKM2.temp)
plot(fitKM2.temp)



hist(sqrt(filter(ani.env, metric == "temp" & diff > 0)$diff))



test <- ani.env %>% 
  filter(metric == "canopy_cover_percent") %>% 
  select(compID, watershed_diff, metric, diff)


