## Make a map of sampling locations in the Eel and Russian river network
## Samples collected 2015 and 2017

# HUC8 for Russian River: 18010110
# Drainage area column for Russian = TotDASqKM
# Drainage area column for Eel = CUMDRAINAG


#### Libraries #################################################################
library(tidyverse)
library(ggplot2)
library(sf)
library(ggspatial) # scale bar and north arrow
source("Scripts/ggplot_themes.R")
################################################################################

  
#### FILE PATHS ################################################################
shp_input_russian <- file.path("Data", "Spatial_data", "Russian_river_network")
shp_input_eel <- file.path("Data", "Spatial_data", "Eel_river_network_shapefile")
shp_input_CA <- file.path("Data", "Spatial_data", "CA_boundary_shapefile")
kml_input <- file.path("Data", "Spatial_data")
################################################################################

#### FORMAT DATA ################################################################

## River networks
river.network.format <- function(shp_pathway, drainage_column= NA, min_drainage_area= 0, crs_code= NA){
  require(sf)
  
  ## Read in shapefile
  river.network.sf <- st_read(shp_pathway)
  
  ## Filter by drainage area
  if(min_drainage_area > 0){
    print(paste("Removing drainage areas <", min_drainage_area))
    river.network.sf.filt <- filter(river.network.sf, get(drainage_column) > min_drainage_area)
  }
  
  ## Transform CRS
  if(is.na(crs_code) == FALSE)
  river.network.sf.filt.transform <- st_transform(river.network.sf.filt, crs= crs_code) 
  
  return(river.network.sf.filt.transform)
}



russian.network <- river.network.format(shp_pathway = shp_input_russian,
                                        drainage_column = "TotDASqKM",
                                        min_drainage_area = 20,
                                        crs_code= 32610) # WGS84 UTM 10

eel.network <- river.network.format(shp_pathway = shp_input_eel,
                                    drainage_column = "CUMDRAINAG",
                                    min_drainage_area = 20,
                                    crs_code= 32610) # WGS84 UTM 10

##  Watershed outlines
eel.watershed <- st_read(file.path(kml_input, "Eel_River_Drainage_Area.kml")) %>%
  st_transform(., crs= 32610) # WGS84 UTM 10

russian.watershed <- st_read(file.path(kml_input, "Russian_River_Drainage_Area.kml")) %>%
  st_transform(., crs= 32610) # WGS84 UTM 10


## California Outline
CA.outline <- st_read(shp_input_CA) %>%
  st_transform(., crs= 32610) # WGS84 UTM 10

## Sampling site location data
latlong.map <- read_csv(file.path("Data", "Spatial_data", "PhormMeta17_LatLong_combined.csv")) %>%
  mutate(year= as.character(year)) %>% 
  # Remove rows with duplicate sites to reduce the number of points
  distinct(acronym, .keep_all = TRUE) %>% 
  ## Convert to Simple Feature (sf) class
  st_as_sf(., coords = c("long", "lat"), 
           crs = 4326) %>%  # WGS84
  st_transform(., crs= 32610) #UTM 10

##### MAKE MAP ##################################################################

## Make base_map
species.location.map <- ggplot() +
  geom_sf(data= CA.outline, color= "black", fill= "snow1",  size= 0.25) +
  geom_sf(data= russian.watershed, size= 0.5, color= "black", fill= NA) +
  geom_sf(data= russian.network, size= 0.3, color= "#2A788EFF") +
  geom_sf(data= eel.watershed, size= 0.5, color= "black", fill= NA) +
  geom_sf(data= eel.network, size= 0.3, color= "#2A788EFF") +
  geom_sf(data= latlong.map, aes(fill= species, shape= species), color= "black", size= 4, alpha= 0.8) +
  scale_shape_manual(values= c(21, 23, 24, 25), name= "Species ID") +
  scale_fill_manual(values= c(species.colors[1], species.colors[2], "mediumpurple3", species.colors[3]), name= "Species ID") +
  annotation_scale(location = "bl", 
                   width_hint = 0.4,
                   height= unit(0.2, "cm")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         height= unit(0.75, "cm"),
                         width= unit(0.75*0.66, "cm"),
                         pad_y = unit(0.75, "cm"),
                         style = north_arrow_orienteering) +
  coord_sf(xlim= c(382000, 545000),
           ylim= c(4238984, 4508831)) +
  labs(x= "Longitude", y= "Latitude") +
  PH2017_map_theme# +
  #theme(panel.background = element_rect(fill= "#bbc2e0"))


## Make inset map
ca.inset <- ggplot() +
  geom_sf(data= CA.outline, size= 0.3, color= "black", fill= "white") +
  geom_sf(data= russian.watershed, size= 0.4, color= "black", fill= "#2A788EFF", alpha= 0.7) +
  geom_sf(data= eel.watershed, size= 0.4, color= "black", fill= "#2A788EFF", alpha= 0.7) +
  coord_sf() +
  PH2017_map_theme + 
  theme(panel.background = element_rect(fill= "transparent"),
        panel.border = element_rect(color= "transparent", fill= "transparent"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())

## Annotate map
species.location.map.anno <- species.location.map +
  annotate("text", x= 475000, y= 4470000, label= "Eel River", angle= 310, size= 4) +
  annotate("text", x= 523000, y= 4305000, label= "Russian River", angle= 310, size= 4) +
    PH2017_map_theme +
  theme(legend.position = c(0.15, 0.3))

## Place inset map
species.map.inset <- ggdraw() +
  draw_plot(species.location.map.anno) +
  draw_plot(ca.inset, x = 0.58, y = 0.7, width = 0.3, height = 0.3, hjust=0, vjust= 0)

## Write file
ggsave(species.map.inset, filename= "Fig_1.png", width= 120, height= 120*1.1, units= "mm", dpi= 300,
       bg= "transparent", path= "Output_figures")
ggsave(species.map.inset, filename= "Fig_1.tiff", width= 120, height= 120*1.1, units= "mm", dpi= 400,
       bg= "transparent", path= "Output_figures")
ggsave(species.map.inset, filename= "Fig_1.pdf", width= 120, height= 120*1.1, units= "mm",
       device= cairo_pdf, bg= "transparent", path= "Output_figures")


  


