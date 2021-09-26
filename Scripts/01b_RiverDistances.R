### Script to calculate river network distances

## Uses R package riverdist
# Vignette: https://rdrr.io/cran/riverdist/f/vignettes/riverdist_vignette.Rmd, 
# Package pdf: https://cran.r-project.org/web/packages/riverdist/riverdist.pdf

# Output units = meters
# Positive distance values indicate second location is upstream of first location
# Negative distance values indicate second location is downstream of first location
# Total flow distance is the total distance in meters regardless of upstream or downstream direction (see the package vignette for more details)
# Net flow distance is the upstream distance minus downstream distance between 2 sites (see the package vignette for more details)


#### Libraries #################################################################
#devtools::install_github("mbtyers/riverdist")
library(tidyverse)
library(riverdist) # https://rdrr.io/cran/riverdist/f/vignettes/riverdist_vignette.Rmd, https://cran.r-project.org/web/packages/riverdist/riverdist.pdf
library(sp)
################################################################################

## Read site locations
sites <- read_csv("Data/Spatial_data/PhormMeta17_LatLong_combined.csv")
site_pairs <- read_tsv("Data/site_pairs.tsv") ## Get combination of site pairs from 02_inStrain_ANI_summary.R

# Load River network objects created by function, create_RiverDist_object
load("Data/Spatial_data/Eel_RiverDistNetworks.RData")

## Create RiverDist object from shapefile
create_RiverDist_object <- function(shapefile_path, layer_name){
## Read in shape files
eel_sp <- rgdal::readOGR(dsn = shapefile_path, layer = layer_name, verbose = TRUE)

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
# Remove additional segments = N
return(list(eel_sp2_clean, eel_sp2, eel_sp))
}

# eel_sp_list <- create_RiverDist_object(shapefile_path = "Data/Spatial_data/Eel_river_network_shapefile", 
#                      layer_name = "eel_rivernetwork_utm")

## Extract river network objects from output list
# eel_sp2_clean <- eel_sp_list[[1]]
# eel_sp2 <- eel_sp_list[[2]]
# eel_sp <- eel_sp_list[[3]]

## Save river network objects
# save(eel_sp2_clean, eel_sp2, eel_sp, file = "Data/Spatial_data/Eel_RiverDistNetworks.RData")



## Function to calculate river distances between pairs of sites
calculate_river_distances <- function(sites, site_pairs, rivdist_obj, Output_Rdata_filename){

## Convert sites to SP object
message("Convert sites to SP object")
#https://cengel.github.io/R-spatial/intro.html#with-sp
sites.ll <- sites[c("long", "lat")]
coordinates(sites.ll) <- ~long+lat
#class(sites.ll)

proj4string(sites.ll) <- CRS("+init=epsg:4326") # this is WGS84

## Reproject to UTM to match the river network
#eel_sp CRS proj4string = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs "
#sites.utm <- spTransform(sites.ll, CRS(proj4string(eel_sp)))
sites.utm <- spTransform(sites.ll, "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs")

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
message("Calculate upstream/downstream location")
flowConnected <- map_chr(1:nrow(site_pairs_xyVert), function(row) 
{riverdirection(startseg = pull(site_pairs_xyVert[row, "seg_n1"]),
                startvert = pull(site_pairs_xyVert[row, "vert_n1"]),
                endseg = pull(site_pairs_xyVert[row, "seg_n2"]),
                endvert = pull(site_pairs_xyVert[row, "vert_n2"]),
                rivers= eel_sp2_clean,
                flowconnected = TRUE,
                algorithm= "Dijkstra")}
    )

message("Calculate total flow distance")
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
message("Calculate net flow distance")
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
message("Save .Rdata file")
save(site_pairs_xyVert, flowConnected, flowDistTotal, flowDistNet,  file = file.path(Output_Rdata_filename))


## Combine data frame and flow-distance vectors together
site_pairs_xyVert_flow <- bind_cols(site_pairs_xyVert, 
                                    tibble(flowConnected, flowDistTotal, flowDistNet)) %>% 
  mutate(FlowConnection= ifelse(is.na(flowConnected), "No", "Yes")) 

save(site_pairs_xyVert_flow, site_pairs_xyVert, flowConnected, flowDistTotal, flowDistNet,  file = file.path(Output_Rdata_filename))


return(site_pairs_xyVert_flow)

}

# takes ~30  minutes to run function
system.time(site_pairs_xyVert_flow <- calculate_river_distances(sites, site_pairs, 
                                                                rivdist_obj= eel_sp2_clean, 
                                                                Output_Rdata_filename= "Data/Spatial_data/flowDist_Vectors.Rdata"))
#load("Data/Spatial_data/flowDist_Vectors.Rdata")

## Combine data frame and flow-distance vectors together
#site_pairs_xyVert_flow <- bind_cols(site_pairs_xyVert, 
                                    #tibble(flowConnected, flowDistTotal, flowDistNet)) %>% 
#  mutate(FlowConnection= ifelse(is.na(flowConnected), "No", "Yes")) 

