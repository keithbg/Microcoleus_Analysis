## Copied from Fig_S4.R script

source("Scripts/Fig_S4.R")

#### CONSENSUS ANI STATISTICS ####
filter(ani.env, metric == "canopy_cover_percent") %>% 
  plot(log(mean_conANI) ~ diff, data= .)

fitCon.canopy <- lm(log(mean_conANI) ~ diff, filter(ani.env, metric == "canopy_cover_percent"))
summary(fitCon.canopy)
plot(fitCon.canopy)
anova(fitCon.canopy)


filter(ani.env, metric == "cond_ms") %>% 
  plot(log(mean_conANI) ~ diff, data= .)

fitCon.cond <- lm(log(mean_conANI) ~ diff, filter(ani.env, metric == "cond_ms"))
summary(fitCon.cond)
anova(fitCon.cond)
plot(fitCon.cond)


filter(ani.env, metric == "temp_NorWest") %>% 
  plot(log(mean_conANI) ~ diff, data= .)

fitCon.NorWest <- lm(log(mean_conANI) ~ diff, filter(ani.env, metric == "temp_NorWest"))
summary(fitCon.NorWest)
anova(fitCon.NorWest)

## TDN and TDP variables no longer significant under inStrain v1.4
filter(ani.env, metric == "TDN_ugL") %>% 
  plot(log(mean_conANI) ~ diff, data= .)

fitCon.TDN <- lm(log(mean_conANI) ~ diff, filter(ani.env, metric == "TDN_ugL"))
summary(fitCon.TDN)
anova(fitCon.TDN)

filter(ani.env, metric == "TDP_ugL") %>% 
  plot(log(mean_conANI) ~ diff, data= .)

fitCon.TDP <- lm(log(mean_conANI) ~ diff, filter(ani.env, metric == "TDP_ugL"))
summary(fitCon.TDP)
anova(fitCon.TDP)









#### MODEL SELECTION ####
conANI.env.w <- ani.env %>% # ani.env generated in Fig_S4.R
  dplyr::select(mean_conANI, year, diff, metric, riv_dist, watershed_diff) %>% 
  pivot_wider(names_from = "metric", values_from = "diff") %>% 
  dplyr::select(-pH, -alk, -NH4_ugL, -mean_vX, -NO3_ugL, -DOC_ugL, -do_mgL, -temp) %>% 
  mutate_at(c("TDP_ugL",  "TDN_ugL", "canopy_cover_percent", "cond_ms"), 
            scale)
conANI.env.w

popANI.env.w <- ani.env %>% 
  dplyr::select(mean_popANI, year, diff, metric, riv_dist, watershed_diff) %>% 
  pivot_wider(names_from = "metric", values_from = "diff") %>% 
  dplyr::select(-pH, -alk, -NH4_ugL, -mean_vX, -NO3_ugL, -DOC_ugL, -do_mgL, -temp)


## Check that variables are not correlated
cor(conANI.env.w[, -c(1, 2)]) 

## Check for year as random effect
cAIC.con.env.dist <- cAIC4::stepcAIC(modAIC.con.env.dist, direction= "forward",
                              groupCandidates = "year",
                              trace= TRUE,
                              data= conANI.env.w)
## Variance of year as random factore ~0, so no need to add as random effect. 
##Keep it as a fixed effect.




## Linear models with all variables
lm.con.env <- lm(mean_conANI ~ ., data= dplyr::select(conANI.env.w, -riv_dist, -watershed_diff))
summary(lm.con.env)
lm.con.env.dist <-  lm(mean_conANI ~ ., data= conANI.env.w)
summary(lm.con.env.dist)
anova(lm.con.env, lm.con.env.dist)

lm.pop.env <- lm(mean_popANI ~ ., data= dplyr::select(popANI.env.w, -riv_dist, -watershed_diff))
lm.pop.env.dist <-  lm(mean_popANI ~ ., data= popANI.env.w)

## Select model based on AIC ##
# conANI
modAIC.con.env <- MASS::stepAIC(lm.con.env)
summary(modAIC.con.env)
modAIC.con.env$anova

modAIC.con.env.dist <- MASS::stepAIC(lm.con.env.dist)
summary(modAIC.con.env.dist)
modAIC.con.env.dist$anova

## Test improvement of adding distance parameters to the statistical fit
anova(modAIC.con.env, modAIC.con.env.dist)

# popANI
modAIC.pop.env <- stepAIC(lm.pop.env)
summary(modAIC.pop.env)
modAIC.pop.env$anova

modAIC.pop.env.dist <- stepAIC(lm.pop.env.dist)
summary(modAIC.pop.env.dist)
modAIC.pop.env.dist$anova
plot(modAIC.pop.env.dist)

anova(modAIC.pop.env, modAIC.pop.env.dist)







