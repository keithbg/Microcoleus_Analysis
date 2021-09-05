



#### Libraries #################################################################
#devtools::install_github("mbtyers/riverdist")
library(tidyverse)
library(riverdist) # https://rdrr.io/cran/riverdist/f/vignettes/riverdist_vignette.Rmd
library(ggplot2)
library(sp)
source("Scripts/ggplot_themes.R")
#source("/Users/kbg/R_functions/ggplot_scalebar_north_arrow.R")
################################################################################

# https://github.com/mbtyers/riverdist/blob/master/vignettes/riverdist_vignette.Rmd
# https://rdrr.io/cran/riverdist/f/vignettes/riverdist_vignette.Rmd
# https://cran.r-project.org/web/packages/riverdist/riverdist.pdf


## Read in shape files
eel.shp.path <- "Data/Spatial_data/Eel_river_network_shapefile"
eel_sp <- rgdal::readOGR(dsn = eel.shp.path, layer = "eel_rivernetwork_utm", verbose = TRUE)


## Subset by drainiage size (km2)
#sp2 <- subset(sp, CUMDRAINAG > 20) # 13 seconds
eel_sp2 <- subset(eel_sp, CUMDRAINAG > 2.5) # 100 seconds

## Create river network object
system.time(eel_sp2.l2n <- line2network(eel_sp2))

## Clean up river network
eel_sp2_clean <- cleanup(eel_sp2.l2n)
# Dissolve = Y
# Split segments = Y
# Insert vertices = Y
# Vertex length = 10 m
# River mouth = segment 77, vertex = 37 (91, 65)
# NO Remove additional segments = Y, 77

## IGNORE BRAIDING THERE ARE LOTS OF DETECTIONS, AND TAKING TOO LONG TO REMOVE THEM ALL
## DOESN'T IMPACT THE FINAL RESULTS

# Braiding detected = Yes, segment 81, remove 81
# 81, 91, 92, 53
# Recheck for braiding, still detected = Yes, remove 92
# Recheck for braiding, still detected = Yes, remove 91
# Recheck for braiding, still detected = Yes, remove 83 (probably should have been 53)
# New numbers: 82, 92, 93, 53, 83

## Set mouth
# https://rdrr.io/cran/riverdist/f/vignettes/riverdist_vignette.Rmd#Incorporating%20flow%20direction
zoomtoseg(seg= c(74, 76, 77, 78, 79, 90, 91, 95), rivers=eel_sp2_clean)
zoomtoseg(seg= c(91), rivers=eel_sp2_clean) # mouth is end of segment 91

showends(seg=90 , rivers=eel_sp2_clean)

eel_sp2_clean <- setmouth(seg=91, vert=65, rivers=eel_sp2_clean)

test2 <- checkbraidedTF(eel_sp2_clean, toreturn= "routes")

zoomtoseg(test2$route2, rivers= eel_sp2_clean)
zoomtoseg(c(72, 79, 80, 81), rivers= eel_sp2_clean)
zoomtoseg(c(52,53, 59, 82, 83), rivers= eel_sp2_clean)



str(eel_sp2_clean$lines[[93]])
summary(eel_sp2_clean)


checkbraidedTF(rivers=KilleyW, toreturn="routes")

# Build segment routes (can't do this because of braiding)
#system.time(eel_sp2_clean2 <- buildsegroutes(eel_sp2_clean, lookup= FALSE))



# Save river network objects
save(eel_sp2_clean, eel_sp2, eel_sp, file = "Data/Eel_RiverDistNetworks.RData")
# To load the data again
load("Data/Eel_RiverDistNetworks.RData")
#topologydots(eel_sp2_clean)

## Read site locations
sites <- read_csv("Data/Spatial_data/PhormMeta17_LatLong_combined.csv")
site_pairs <- read_tsv("Data/site_pairs.tsv") ## Get combination of site pairs from inStrain_ANI_summary.R

#site_pairs <- ani_sum %>% 
#  select(name1, name2)
#write_tsv(site_pairs, "Data/site_pairs.tsv")

## Convert sites to SP object
#https://cengel.github.io/R-spatial/intro.html#with-sp
sites.ll <- sites[c("long", "lat")]
coordinates(sites.ll) <- ~long+lat
#class(sites.ll)

proj4string(sites.ll) <- CRS("+init=epsg:4326") # this is WGS84

## Reproject to UTM to match the river network
sites.utm <- spTransform(sites.ll, CRS(proj4string(eel_sp)))
proj4string(sites.utm) 

## Snap the sites to the river network
eel_xyVert <- xy2segvert(sites.utm@coords[, "long"], 
                         sites.utm@coords[, "lat"], 
                         eel_sp2_clean) %>% 
  cbind(sites[, "ggkbase_id"], .) # add site name column
  

## Get vertices for every site pair
site_pairs_xyVert <- left_join(site_pairs, eel_xyVert, by= c("name1" = "ggkbase_id")) %>% 
  rename(seg_n1= seg, vert_n1= vert, snapdist_n1= snapdist) %>% #n1 = name 1
  left_join(., eel_xyVert, by= c("name2" = "ggkbase_id")) %>% 
  rename(seg_n2= seg, vert_n2= vert, snapdist_n2= snapdist) #n2 = name 2
  




## Calculate Upstream/Downstream location
## Returns NA if flow is unconnected

flowConnected <- map_chr(1:nrow(site_pairs_xyVert), function(row) 
{riverdirection(startseg = pull(site_pairs_xyVert[row, "seg_n1"]),
                startvert = pull(site_pairs_xyVert[row, "vert_n1"]),
                endseg = pull(site_pairs_xyVert[row, "seg_n2"]),
                endvert = pull(site_pairs_xyVert[row, "vert_n2"]),
                rivers= eel_sp2_clean,
                flowconnected = TRUE,
                algorithm= "Dijkstra")}
    )
table(flowconnected_vec, useNA= "always")

flowDistTotal <- map_dbl(1:nrow(site_pairs_xyVert), function(row) {
  upstream(startseg = pull(site_pairs_xyVert[row, "seg_n1"]),
           startvert = pull(site_pairs_xyVert[row, "vert_n1"]),
           endseg = pull(site_pairs_xyVert[row, "seg_n2"]),
           endvert = pull(site_pairs_xyVert[row, "vert_n2"]),
           rivers= eel_sp2_clean,
           flowconnected = FALSE,
           net= FALSE,
           algorithm= "Dijkstra")}
)

flowDistNet <- map_dbl(1:nrow(site_pairs_xyVert), function(row) {
  upstream(startseg = pull(site_pairs_xyVert[row, "seg_n1"]),
           startvert = pull(site_pairs_xyVert[row, "vert_n1"]),
           endseg = pull(site_pairs_xyVert[row, "seg_n2"]),
           endvert = pull(site_pairs_xyVert[row, "vert_n2"]),
           rivers= eel_sp2_clean,
           flowconnected = FALSE,
           net= TRUE,
           algorithm= "Dijkstra")}
)

save(site_pairs_xyVert, flowConnected, flowDistTotal, flowDistNet,  file = "Data/flowDist_Vectors.Rdata")
load("Data/flowDist_Vectors.Rdata")


## Calculate fraction of distance that is downstream
flowDownstreamFrac <- flowDistNet / flowDistTotal
flowUpstreamFrac <- (1 - ifelse(is.nan(flowDownstreamFrac), 1, flowDownstreamFrac))# NaN are sites at the same location or zero distance, 0/0 = NaN

## Combine data frame and flow-distance vectors together
site_pairs_xyVert_flow <- bind_cols(site_pairs_xyVert, 
                                    tibble(flowConnected, flowDistTotal, flowDistNet, flowUpstreamFrac)) %>% 
  mutate(FlowConnection= ifelse(is.na(flowConnected), "No", "Yes")) 


ani_sum <- read_tsv("Data/inStrain_data/ani_summary_v1.4.tsv") # .tsv file generated in: Scripts/inStrain_ANI_summary.R

ani_rivDist <- left_join(select(ani_sum, name1, name2, mean_conANI, mean_popANI, riv_dist, euc_dist),
                         select(site_pairs_xyVert_flow, name1, name2, FlowConnection, flowDistTotal, flowDistNet, flowUpstreamFrac)) %>% 
  mutate(riv_dist= riv_dist/1000,
         euc_dist= euc_dist/1000,
         flowDistTotal= abs(flowDistTotal/1000),
         dist_diff= riv_dist - flowDistTotal,
         mean_conANI= mean_conANI*100,
         mean_popANI= mean_popANI*100)

test <- filter(ani_rivDist, FlowConnection == "Yes")

#### STATISTICS ####
fit.flow.popANI <- lm(mean_popANI ~ flowDistTotal, filter(ani_rivDist, FlowConnection == "Yes"))
#summary(fit.flow.popANI)
#anova(fit.flow.popANI)

fit.flow.conANI <- lm(mean_conANI ~ flowDistTotal, filter(ani_rivDist, FlowConnection == "Yes"))
#summary(fit.flow.conANI)
#anova(fit.flow.conANI)


conANI_flow <- ggplot(data= filter(ani_rivDist, FlowConnection == "Yes"), aes(x= flowDistTotal, y= mean_conANI)) +
  geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_abline(intercept = fit.flow.conANI$coefficients["(Intercept)"], slope= fit.flow.conANI$coefficients["flowDistTotal"],
              color= "black", size= 1) +
 # annotate("text", x= 250, y= 99.92, label= "Sites connected \nby flow", hjust= 0) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_x_continuous(limits= c(0, 150),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains

popANI_flow <- ggplot(data= filter(ani_rivDist, FlowConnection == "Yes"), aes(x= flowDistTotal, y= mean_popANI)) +
  geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_abline(intercept = fit.flow.popANI$coefficients["(Intercept)"], slope= fit.flow.popANI$coefficients["flowDistTotal"],
              color= "black", size= 1) +
 # annotate("text", x= 110, y= 99.92, label= "Sites connected \nby flow", hjust= 0) +
  labs(x= "River network distance (km)", y= "Population ANI (%)") +
  scale_x_continuous(limits= c(0, 150),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains


#### STATISTICS ####
fit.noflow.popANI <- lm(mean_popANI ~ flowDistTotal, filter(ani_rivDist, FlowConnection == "No"))
#summary(fit.noflow.popANI)
#anova(fit.noflow.popANI)

fit.noflow.conANI <- lm(mean_conANI ~ flowDistTotal, filter(ani_rivDist, FlowConnection == "No"))
#summary(fit.noflow.conANI)
#anova(fit.noflow.conANI)


conANI_noflow <- ggplot(data= filter(ani_rivDist, FlowConnection == "No"), aes(x= riv_dist, y= mean_conANI)) +
  geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_abline(intercept = fit.noflow.conANI$coefficients["(Intercept)"], slope= fit.noflow.conANI$coefficients["flowDistTotal"],
              color= "black", size= 1) +
  #annotate("text", x= 250, y= 99.92, label= "Sites unconnected \nby flow", hjust= 0) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_shape_manual(values= c(21, 24)) +
  scale_x_continuous(limits= c(0, 350),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains 

popANI_noflow <- ggplot(data= filter(ani_rivDist, FlowConnection == "No"), aes(x= flowDistTotal, y= mean_popANI)) +
  geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_abline(intercept = fit.noflow.popANI$coefficients["(Intercept)"], slope= fit.noflow.popANI$coefficients["flowDistTotal"],
              color= "black", size= 1) +
  #annotate("text", x= 110, y= 99.92, label= "Sites unconnected \nby flow", hjust= 0) +
  labs(x= "River network distance (km)", y= "Population ANI (%)") +
  scale_x_continuous(limits= c(0, 350),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains


#### STATISTICS ####
fit.euc.popANI <- lm(mean_popANI ~ euc_dist, ani_rivDist)
#summary(fit.euc.popANI)
#anova(fit.euc.popANI)

fit.euc.conANI <- lm(mean_conANI ~ euc_dist, ani_rivDist)
#summary(fit.euc.conANI)
#anova(fit.euc.conANI)

conANI_euc <- ggplot(data= ani_rivDist, aes(x= euc_dist, y= mean_conANI)) +
  geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_abline(intercept = fit.euc.conANI$coefficients["(Intercept)"], slope= fit.euc.conANI$coefficients["euc_dist"],
              color= "black", size= 1) +
  #annotate("text", x= 100, y= 99.92, label= "Euclidean distance \nbetween sites", hjust= 0) +
  labs(x= "Euclidean distance (km)", y= "Consensus ANI (%)") +
  scale_shape_manual(values= c(21, 24)) +
  scale_x_continuous(limits= c(0, 150),
                     breaks= seq(0, 150, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains 

popANI_euc <- ggplot(data= ani_rivDist, aes(x= euc_dist, y= mean_popANI)) +
  geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_abline(intercept = fit.euc.popANI$coefficients["(Intercept)"], slope= fit.euc.popANI$coefficients["euc_dist"],
              color= "black", size= 1) +
  #annotate("text", x= 110, y= 99.92, label= "Sites connected \nby flow", hjust= 0) +
  labs(x= "Euclidean distance (km)", y= "Population ANI (%)") +
  scale_x_continuous(limits= c(0, 150),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains
t_marg <- 0.2
r_marg <- 0.2
l_marg <- 0.2
b_marg <- 0.2


ANI_rivdist_combined1 <- plot_grid(conANI_flow + theme(axis.title.x = element_blank(),
                                                       plot.margin = unit(c(t_marg, r_marg, b_marg, l_marg), "cm")), 
                                  conANI_noflow + theme(axis.title.x = element_blank(), 
                                                        axis.title.y = element_blank(),
                                                        plot.margin = unit(c(t_marg, r_marg, b_marg, l_marg), "cm")), 
                                  conANI_euc + theme(axis.title.x = element_blank(),
                                                     axis.title.y = element_blank(),
                                                     plot.margin = unit(c(t_marg, r_marg, b_marg, l_marg), "cm")),
                                  popANI_flow + theme(axis.title.x = element_blank(),
                                                      plot.margin = unit(c(t_marg, r_marg, b_marg, l_marg), "cm")), 
                                  popANI_noflow + theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        plot.margin = unit(c(t_marg, r_marg, b_marg, l_marg), "cm")), 
                                  popANI_euc + theme(axis.title.x = element_blank(),
                                                     axis.title.y = element_blank(),
                                                     plot.margin = unit(c(t_marg, r_marg, b_marg, l_marg), "cm")),
                                  nrow= 2,
                                  byrow= TRUE,
                                  #labels= c("A", "B", "C", "D", "E", "F"),
                                  labels= NULL,
                                  align= "hv",
                                  axis= "l")

ANI_rivdist_combined2 <-  annotate_figure(ANI_rivdist_combined1, 
                                          bottom= text_grob(label= "Distance (km)", vjust= -0.1),
                                          top= text_grob(label= c("Connected flow", "Unconnected flow", "Euclidean distance"),  
                                                         x= c(0.2, 0.55, 0.85),
                                                         hjust= 0.5,
                                                         vjust= 1,
                                                         face= "bold"))


ANI_rivdist_combined2


ggsave(ANI_rivdist_combined2, filename = "RivDist_test.png", height= 200*0.75, width= 200, units= "mm", dpi= 320,
       path= "Output_figures")





fit.noflow2.conANI <- lm(mean_conANI ~ flowDistTotal*FlowConnection, data= ani_rivDist)
summary(fit.noflow2.conANI)
anova(fit.noflow2.conANI)

fit.noflow2.popANI <- lm(mean_popANI ~ flowDistTotal*FlowConnection, data= ani_rivDist)
summary(fit.noflow2.popANI)
anova(fit.noflow2.popANI)


conANI2 <- ggplot(data= ani_rivDist, aes(x= flowDistTotal, y= mean_conANI)) +
  #geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_point(aes(fill= FlowConnection), size= 3, pch= 21, color= "black", alpha= 0.5) +
  geom_line(aes(x= flowDistTotal, y= predict(fit.noflow2.conANI), color= FlowConnection, group= FlowConnection), size= 2) +
  #geom_abline(intercept = fit.flow.conANI$coefficients["(Intercept)"], slope= fit.flow.conANI$coefficients["flowDistTotal"],
  #            color= "black", size= 1) +
  # annotate("text", x= 250, y= 99.92, label= "Sites connected \nby flow", hjust= 0) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_x_continuous(limits= c(0, 350),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains



popANI2 <- ggplot(data= ani_rivDist, aes(x= flowDistTotal, y= mean_popANI)) +
  #geom_point(size= 3, pch= 21, fill= species.colors[1], color= "black", alpha= 0.5) +
  geom_point(aes(fill= FlowConnection), size= 3, pch= 21, color= "black", alpha= 0.5) +
  geom_line(aes(x= flowDistTotal, y= predict(fit.noflow2.popANI), color= FlowConnection, group= FlowConnection), size= 2) +
  #geom_abline(intercept = fit.flow.conANI$coefficients["(Intercept)"], slope= fit.flow.conANI$coefficients["flowDistTotal"],
  #            color= "black", size= 1) +
  # annotate("text", x= 250, y= 99.92, label= "Sites connected \nby flow", hjust= 0) +
  labs(x= "River network distance (km)", y= "Population ANI (%)") +
  scale_x_continuous(limits= c(0, 350),
                     breaks= seq(0, 350, by= 25),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5)) +
  scale_y_continuous(limits= c(99.19, 100),
                     breaks= seq(99.2, 100.0, by= 0.1),
                     labels= c("99.2", "", "99.4", "", "99.6", "", "99.8", "", "100.0"),
                     expand= c(0.035,0)) +
  theme_strains


combined <- plot_grid(popANI2, conANI2, nrow= 2)


combined.euc <- plot_grid(popANI_euc, conANI_euc, nrow= 2)


ggsave(combined, filename = "RivDist2.png", height= 200, width= 200, units= "mm", dpi= 320,
       path= "Output_figures")

ggsave(combined.euc, filename = "EucDist2.png", height= 200, width= 200, units= "mm", dpi= 320,
       path= "Output_figures")


#### OLD CODE ##################################################




hist(ani_rivDist$flowDownstreamFrac)

fit1 <- lm((mean_popANI*100) ~ riv_dist, ani_rivDist)
summary(fit1)
anova(fit1)

fit2 <- lm((mean_popANI*100) ~ riv_dist + FlowConnection, ani_rivDist)
summary(fit2)
anova(fit2)


fit3 <- lm((mean_popANI*100) ~ riv_dist * FlowConnection, ani_rivDist)
summary(fit3)
anova(fit3)


fit4 <- lm((mean_popANI*100) ~ FlowConnection, ani_rivDist)
summary(fit4)
anova(fit4)
plot(fit4)

fit5 <- lm((mean_popANI*100) ~ riv_dist + FlowConnection, 
           filter(ani_rivDist, riv_dist > 30000))
summary(fit5)
anova(fit5)



boxplot((mean_popANI*100)  ~ FlowConnection, 
        filter(ani_rivDist, riv_dist > 40000))


boxplot(riv_dist ~ FlowConnection, filter(ani_rivDist, riv_dist > 40000))

plot(mean_popANI*100 ~ riv_dist + FlowConnection, ani_rivDist)
hist(filter(ani_rivDist, FlowConnection == "Yes")$riv_dist)
hist(filter(ani_rivDist, FlowConnection == "No")$riv_dist)
?hist

table(ani_rivDist$FlowConnection)

anova(fit1, fit2, fit3)

cor.test(ani_rivDist$riv_dist, ani_rivDist$flowUpstreamFrac, method= "kendall")

cor(ani_rivDist$riv_dist, ani_rivDist$flowUpstreamFrac, method= "spearman")


fitA <- lm((mean_conANI*100) ~ riv_dist, ani_rivDist)
summary(fitA)
anova(fitA)

fitB <- lm((mean_conANI*100) ~ riv_dist + flowUpstreamFrac, ani_rivDist)
summary(fitB)
anova(fitB)
plot(fitB)
predict(fitB)

ggplot()


ggplot(ani_rivDist, aes(x= riv_dist, y= mean_conANI*100)) +
  geom_point(aes(color= flowUpstreamFrac)) +
  geom_point(aes(x= ani_rivDist$riv_dist, y= predict(fitB)), color= "red") +
  geom_point(aes(x= ani_rivDist$riv_dist, y= predict(fitC)), color= "blue") +
  labs(x= "River distance (km)", y= "Mean consensus ANI (%)") +
  viridis::scale_color_viridis(discrete = FALSE,
                               name= "Upstream \nfraction") +
  ggtitle("Consensus ANI. Red= no interaction, Blue= interaction")
ggsave(last_plot(), filename= "flow_upstream_Interaction.png", height= 5, width= 8, 
       path= "/Users/kbg/Downloads")




fitB2 <- lm((mean_popANI*100) ~  flowUpstreamFrac + riv_dist, ani_rivDist)
summary(fitB2)
anova(fitB2)


fitC <- lm((mean_conANI*100) ~ riv_dist * flowUpstreamFrac, ani_rivDist)
summary(fitC)
anova(fitC)
plot(fitC)

fitD <- lm((mean_popANI*100) ~ flowUpstreamFrac * riv_dist, ani_rivDist)
summary(fitD)
anova(fitD)
plot(fitD)


lm(abs(flowDistTotal) ~ riv_dist, data= ani_rivDist)

ggplot(ani_rivDist, aes(x= riv_dist, y= abs(flowDistTotal))) +
  geom_point() +
  geom_abline(slope=1, intercept= 0, color= "blue")

ggplot(ani_rivDist, aes(x= riv_dist/1000, y= flowUpstreamFrac)) +
  geom_point(aes(color= mean_conANI), size=3) +
  labs(y= "Fraction of flow distance upstream", x= "River distance (km)") +
  viridis::scale_color_viridis(name= "Mean consensus ANI (%)") +
  theme(legend.position = "bottom")
ggsave(last_plot(), filename= "flow_upstream.png", height= 5, width= 8, path= "/Users/kbg/Downloads")


ggplot(ani_rivDist, aes(x= flowUpstreamFrac, y= mean_conANI*100)) +
  geom_point(aes(color= riv_dist/1000)) +
  labs(x= "Fraction of flow distance upstream", y= "Mean consensus ANI (%)") +
  viridis::scale_color_viridis(discrete = FALSE,
                               name= "Connected \nflow") +
  geom_smooth() +
  ggtitle("Consensus ANI")
ggsave(last_plot(), filename= "flow_conANI.png", height= 5, width= 8, path= "/Users/kbg/Downloads")

ggplot(ani_rivDist, aes(x= riv_dist/1000, y= mean_conANI*100)) +
  geom_point(aes(color= flowUpstreamFrac)) +
  labs(x= "River distance (km)", y= "Mean consensus ANI (%)") +
  viridis::scale_color_viridis(discrete = FALSE,
                               name= "Upstream \nfraction") +
  ggtitle("Consensus ANI")
ggsave(last_plot(), filename= "rivDist_conANI.png", height= 5, width= 8, path= "/Users/kbg/Downloads")


ggplot(ani_rivDist, aes(x= riv_dist/1000, y= mean_popANI*100)) +
  geom_point(aes(color= flowUpstreamFrac)) +
  labs(x= "River distance (km)", y= "Mean population ANI (%)") +
  viridis::scale_color_viridis(discrete = FALSE,
                               name= "Upstream \nfraction") +
  ggtitle("Consensus ANI")
ggsave(last_plot(), filename= "rivDist_popANI.png", height= 5, width= 8, path= "/Users/kbg/Downloads")



ggplot(ani_rivDist, aes(x= flowUpstreamFrac, y= mean_popANI*100)) +
  geom_point(aes(color= riv_dist/1000)) +
  labs(x= "Fraction of flow distance upstream", y= "Mean population ANI (%)") +
  viridis::scale_color_viridis(discrete = FALSE,
                               name= "Connected \nflow") +
  geom_smooth() +
  ggtitle("Population ANI")
ggsave(last_plot(), filename= "flow_popANI.png", height= 5, width= 8, path= "/Users/kbg/Downloads")

 
hist(ani_rivDist$flowDistTotal)


