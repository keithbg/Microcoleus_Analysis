## Copied from Fig_S4.R script

source("Scripts/Fig_S4.R")

#### MODEL SELECTION ####
conANI.env.w <- ani.env %>% # ani.env generated in Fig_S4.R
  dplyr::select(mean_conANI, year, diff, metric, riv_dist, watershed_diff) %>% 
  pivot_wider(names_from = "metric", values_from = "diff") %>% 
  dplyr::select(-pH, -alk, -NH4_ugL, -mean_vX, -NO3_ugL, -DOC_ugL, -do_mgL, -temp)# %>% 
mutate_at(c("temp", "TDP_ugL",  "TDN_ugL", "canopy_cover_percent", "cond_ms"), 
          scale)
conANI.env.w

popANI.env.w <- ani.env %>% 
  dplyr::select(mean_popANI, year, diff, metric, riv_dist, watershed_diff) %>% 
  pivot_wider(names_from = "metric", values_from = "diff") %>% 
  dplyr::select(-pH, -alk, -NH4_ugL, -mean_vX, -NO3_ugL, -DOC_ugL, -do_mgL, -temp)



cor(conANI.env.w[, -c(1, 2)]) ## correlated variables removed


?cor
names(conANI.env.w)

lm.con.env <- lm(mean_conANI ~ ., data= dplyr::select(conANI.env.w, -riv_dist, -watershed_diff))
lm.con.env.dist <-  lm(mean_conANI ~ ., data= conANI.env.w)
anova(lm.con.env, lm.con.env.dist)

step1 <- step(lm.con.env.dist)
anova(step1)
?step
lm.pop.env <- lm(mean_popANI ~ ., data= dplyr::select(popANI.env.w, -riv_dist, -watershed_diff))
lm.pop.env.dist <-  lm(mean_popANI ~ ., data= popANI.env.w)

library(MASS)
modAIC.con.env <- stepAIC(lm.con.env)
summary(modAIC.con.env)
modAIC.con.env$anova

modAIC.con.env.dist <- stepAIC(lm.con.env.dist)
summary(modAIC.con.env.dist)
modAIC.con.env.dist$anova

anova(modAIC.con.env, modAIC.con.env.dist)

modAIC.pop.env <- stepAIC(lm.pop.env)
summary(modAIC.pop.env)
modAIC.pop.env$anova

modAIC.pop.env.dist <- stepAIC(lm.pop.env.dist)
summary(modAIC.pop.env.dist)
modAIC.pop.env.dist$anova
plot(modAIC.pop.env.dist)

anova(modAIC.pop.env, modAIC.pop.env.dist)




cAIC.con.env.dist <- stepcAIC(modAIC.con.env.dist, direction= "forward",
                              groupCandidates = "year",
                              trace= TRUE,
                              data= conANI.env.w)
## Variance of year as random factore ~0, so no need to add as Random effect. Keep it as a fixed effect.





modAIC.1$anova

car::vif(lm.1)

library(cAIC4)
lm.0 <- lm(mean_conANI ~ cond_ms + canopy_cover_percent, data= conANI.env.w)
lm.0 <- lm(mean_conANI ~ ., data= dplyr::select(conANI.env.w, -year))

lmer.0 <- lmer(log(mean_conANI) ~ 1 + (1|year), data= conANI.env.w)

lmer.1 <- lmer(mean_conANI ~ . + (1|year), data= conANI.env.w)

modcAIC.1 <- stepcAIC(lmer.1, trace= TRUE, direction= "backward", data= conANI.env.w)
modcAIC.1

cAIC.f.1 <- stepcAIC(lm.0, direction= "forward",
                     groupCandidates = "year",
                     trace= TRUE,
                     data= conANI.env.w,
                     numberOfSavedModels = 2)
summary(cAIC.f.1)
cAIC.f.1$additionalModels

str(cAIC.f.1)

anocAIC(c(lmer.0, lmer.1))
cAIC(lm.0)


(fm3 <- lmer(strength ~ 1 + (1|sample), Pastes))
fm3_step <- stepcAIC(fm3, direction = "backward", trace = TRUE, data = Pastes)




#### PRINCIPAL COMPONENTS ANALYSIS

watersheds <- ani.env %>% 
  select(name1, watershed_1) %>% 
  distinct() %>% 
  rename(ggkbase_id= name1) %>% 
  rbind(data.frame(ggkbase_id= "PH2017_40_RAT_O_B", watershed_1= 50.0))

## SNV data and young/old populations
snv_genomes <- read_tsv("Output_tables/snvs_genome_summary.tsv") %>% 
  filter(species == "species_1") %>% 
  select(ggkbase_id, SNV_mbp, pop_age)


## ORDINATE
tenv.ord <- tenv %>% 
  filter(biotite_ID %in% snv_genomes$ggkbase_id) %>% # include only Species 1
  left_join(., watersheds, by= c("biotite_ID" = "ggkbase_id")) %>% 
  select(-biotite_ID, -envID, -alk, -pH, -do_mgL, -mean_vX, -temp, -DOC_ugL, -NH4_ugL, -NO3_ugL) %>% 
  as.matrix(.) %>% 
  scale(.)
row.names(tenv.ord) <- filter(tenv, biotite_ID %in% snv_genomes$ggkbase_id)$biotite_ID

library(vegan)
env.pca <- rda(tenv.ord)
summary(env.pca)


## EXTRACT ORDINATION SCORES AND LOADINGS
pca.site.scores <- as.data.frame(scores(env.pca, choices= c(1, 2, 3,4))$sites) %>% 
  mutate(ggkbase_id= rownames(.)) %>% 
  as_tibble() %>% 
  left_join(., snv_genomes) %>% 
  mutate(year= ifelse(str_detect(.$ggkbase_id, "2015"), "2015", "2017"))
#row.names(pca.site.scores) <- str_replace(row.names(pca.site.scores), "PH20.*_", "")
pca.covariate.scores <- as.data.frame(scores(env.pca, choices= c(1, 2, 3, 4))$species) %>% 
  mutate(param= rownames(.)) %>% 
  as_tibble()



## PCA 1 & PCA 2
ggplot() +
  geom_hline(yintercept = 0, color= "gray") +
  geom_vline(xintercept= 0, color= "gray") +
  geom_point(data= pca.site.scores, aes(x= PC1, y= PC2, color= pop_age, shape= year), size= 2) +
  #geom_text_repel(data= pca.site.scores, aes(x= PC1, y= PC2, label= ord.sites), color= "black", size= 3) +
  stat_ellipse(data= pca.site.scores, aes(x= PC1, y= PC2, color= pop_age)) +
  geom_segment(data= pca.covariate.scores, aes(x= 0, xend= PC1, y=0, yend= PC2),
               arrow= arrow(length = unit(0.2, "cm")), size= 0.75,
               color= "black", alpha= 0.3) +
  ggrepel::geom_text_repel(data= pca.covariate.scores, aes(x= PC1, y= PC2, label= param), color= "Tomato", size= 5) +
  labs(x= "PC1 (34.7%)", y= "PC2 (28.4%)") +
  coord_equal() +
  theme_strains
#ggsave(last_plot(), file= "OldYoung_EnvData_PCA_1&2.pdf", height= 6, width= 6, units= "in", path= "Output_figures", device = cairo_pdf)
ggsave(last_plot(), file= "OldYoung_EnvData_PCA_1&2.png", height= 6, width= 6, units= "in", path= "Output_figures", dpi= 300)










