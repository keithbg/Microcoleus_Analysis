



#### Libraries #################################################################
#devtools::install_github("mbtyers/riverdist")
library(tidyverse)
library(riverdist) # https://rdrr.io/cran/riverdist/f/vignettes/riverdist_vignette.Rmd
library(ggplot2)
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
ani_rivDist <- left_join(select(ani_sum, name1, name2, mean_conANI, mean_popANI, riv_dist),
                         select(site_pairs_xyVert_flow, name1, name2, FlowConnection, flowDistTotal, flowDistNet, flowUpstreamFrac))

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

cor.test(ani_rivDist$riv_dist, ani_rivDist$flowDownstreamFrac, method= "kendall")

cor(ani_rivDist$riv_dist, ani_rivDist$flowDownstreamFrac, method= "spearman")


fitA <- lm((mean_conANI*100) ~ riv_dist, ani_rivDist)
summary(fitA)
anova(fitA)

fitB <- lm((mean_conANI*100) ~  flowUpstreamFrac + riv_dist, ani_rivDist)
summary(fitB)
anova(fitB)


fitC <- lm((mean_conANI*100) ~ flowUpstreamFrac * riv_dist, ani_rivDist)
summary(fitC)
anova(fitC)

fitD <- lm((mean_popANI*100) ~ riv_dist * flowUpstreamFrac, ani_rivDist)
summary(fitD)
anova(fitD)
plot(fitD)


lm(abs(flowDistTotal) ~ riv_dist, data= ani_rivDist)

ggplot(ani_rivDist, aes(x= riv_dist, y= abs(flowDistTotal))) +
  geom_point() +
  geom_abline(slope=1, intercept= 0, color= "blue")

ggplot(ani_rivDist, aes(x= riv_dist, y= flowUpstreamFrac)) +
  geom_point(aes(color= mean_conANI), size=3) +
  viridis::scale_color_viridis()

ggplot(ani_rivDist, aes(x= flowUpstreamFrac, y= mean_conANI)) +
  geom_point(aes(color= FlowConnection))

hist(ani_rivDist$flowDistTotal)


