## Make a map of the Eel and Russian river network

# HUC8 for Russian River: 18010110
# Drainage area column for Russian = TotDASqKM
# Drainage area column for Eel = CUMDRAINAG


#### Libraries #################################################################
library(tidyverse)
library(ggplot2)
library(ggmap)
source("Scripts/ggplot_scalebar_north_arrow.R")
source("Scripts/ggplot_themes.R")
################################################################################

#### IMPORT SHAPE FILES ########################################################
import_shape_files <- function(){

  require(tidyverse)
  require(ggmap)
  #library(ggsn)
  require(rgdal)
  require(broom)
#### FILE PATHS ################################################################
shp_input_russian <- file.path("Data", "Map_data", "Russian_river_network")
shp_input_eel <- file.path("Data", "Map_data", "Eel_river_network_shapefile")
shp_input_CA <- file.path("Data", "Map_data", "CA_boundary_shapefile")
kml_input <- file.path("Data", "Map_data")
dir_output <- file.path("Data", "Map_data")
################################################################################


## River networks
river.network.shp.import <- function(shp_pathway, layer_name, drainage_column= NA, min_drainage_area= 0){
  river.network.shp <- readOGR(dsn= shp_pathway, layer= layer_name)
  # Get column names in the attributes table
  print(ogrInfo(dsn= shp_pathway, layer= layer_name))
  
  # class(russian.network.shp) # SpatialLinesDataFrame
  # proj4string(russian.network.shp) # projection of shapefile
  if(min_drainage_area > 0){
    print(paste("Removing drainage areas <", min_drainage_area))
    col_num <- which(names(river.network.shp@data) %in% drainage_column) # get column number of drainage area data
    river.network.shp <- river.network.shp[river.network.shp@data[, col_num] >= min_drainage_area, ] # keep drainage areas >10 km^2
  }
  river.network.shp <- spTransform(river.network.shp, CRS("+proj=longlat +datum=WGS84")) # reproject to match ggmap
  river.network.shp <- tidy(river.network.shp)
  return(river.network.shp)
}

russian.network <- river.network.shp.import(shp_pathway = shp_input_russian,
                                            layer_name = "Russian_river_network",
                                            drainage_column = "TotDASqKM",
                                            min_drainage_area = 20)
eel.network <- river.network.shp.import(shp_pathway = shp_input_eel,
                                        layer_name = "eel_rivernetwork_utm",
                                        drainage_column = "CUMDRAINAG",
                                        min_drainage_area = 20)

##  Watershed outlines
eel.watershed <- readOGR(file.path(kml_input, "Eel_River_Drainage_Area.kml"), layer= "Eel_River_Drainage_Area") %>%
  spTransform(., CRS("+proj=longlat +datum=WGS84")) %>%  # reproject to match ggmap
  tidy(.)

russian.watershed <- readOGR(file.path(kml_input, "Russian_River_Drainage_Area.kml"), layer= "Russian_River_Drainage_Area") %>%
  spTransform(., CRS("+proj=longlat +datum=WGS84")) %>%  # reproject to match ggmap
  tidy(.)


## California Outline
CA.outline <- readOGR(dsn= shp_input_CA, layer= "CA_boundary_kbg") %>%
  spTransform(., CRS("+proj=longlat +datum=WGS84")) %>%  # reproject to match ggmap
  tidy(.)

return(list(eel_network= eel.network, eel_watershed= eel.watershed, 
            russian_network= russian.network, russian_watershed= russian.watershed,
            CA_outline= CA.outline))
}

shapefile.list <- import_shape_files()
eel.network <- shapefile.list[["eel_network"]]
eel.watershed <- shapefile.list[["eel_watershed"]]
russian.network <- shapefile.list[["russian_network"]]
russian.watershed <- shapefile.list[["russian_watershed"]]
CA.outline <- shapefile.list[["CA_outline"]]
rm(shapefile.list)

##### MAKE MAP ##################################################################

## Read in sampling site location data
latlong.map <- read_csv(file.path("Data", "Map_data", "PhormMeta17_LatLong_combined.csv")) %>%
  mutate(year= as.character(year)) %>% 
  # Remove rows with duplicate sites to reduce the number of points
  distinct(acronym, .keep_all = TRUE)

## Make base_map
PH2017_eel_russian_base_map <- ggplot() +
  geom_polygon(data= CA.outline, aes(x= long, y= lat, group=group), color= "black", fill= "snow1") +
  geom_polygon(data= russian.watershed, aes(x= long, y= lat, group=group), size= 0.5, color= "black", fill= NA) +
  geom_polygon(data= eel.watershed, aes(x= long, y= lat, group=group), size= 0.5, color= "black", fill= NA) +
  geom_path(data= russian.network, aes(x= long, y= lat, group=group), size= 0.3, color= "#2A788EFF") +
  geom_path(data= eel.network, aes(x= long, y= lat, group=group), size= 0.3, color= "#2A788EFF") +
  scale_bar(lon = -124.4, lat = 38.3,
            distance_lon = 20, distance_lat = 3, distance_legend = 6, dist_unit = "km",
            arrow_length= 10, arrow_distance = 8) +
  coord_map(xlim= c(-124.5, -122.4), ylim= c(38.27, 40.75)) +
  scale_y_continuous(breaks= c(38.5, 39, 39.5, 40, 40.5), labels= c("", "39", "", "40", "")) +
  scale_x_continuous(breaks= c(-124.5, -124, -123.5, -123.0, -122.5), labels= c("", "-124", "", "-123", "")) +
  labs(x= "Longitude", y= "Latitude") +
  PH2017_map_theme


## Make inset map
ca_inset <- ggplot() +
  geom_polygon(data= CA.outline, aes(x= long, y= lat, group=group), size= 0.3, color= "black", fill= "white") +
  geom_polygon(data= russian.watershed, aes(x= long, y= lat, group=group), size= 0.4, color= "black", fill= "#2A788EFF", alpha= 0.7) +
  geom_polygon(data= eel.watershed, aes(x= long, y= lat, group=group), size= 0.4, color= "black", fill= "#2A788EFF", alpha= 0.7) +
  coord_map() +
  PH2017_map_theme + 
  theme(panel.background = element_rect(fill= "transparent"),
        panel.border = element_rect(color= "transparent", fill= "transparent"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())


species.map <- PH2017_eel_russian_base_map +
  geom_point(data= latlong.map, aes(x= long, y= lat, fill= species, shape= species), color= "black", size= 4, alpha= 0.8) +
  #geom_label_repel(data= latlong.map, aes(x= long, y= lat, label= acronym)) +
  annotate("text", x= -123.3, y= 40.35, label= "Eel River", angle= 310, size= 4) +
  annotate("text", x= -122.7, y= 38.9, label= "Russian River", angle= 310, size= 4) +
  scale_shape_manual(values= c(21, 23, 24, 25), name= "Species ID") +
  scale_fill_manual(values= c(species.colors[1], species.colors[2], "mediumpurple3", species.colors[3]), name= "Species ID") +
  PH2017_map_theme +
  theme(legend.position = c(0.15, 0.3))

species.map.inset <- ggdraw() +
  draw_plot(species.map) +
  draw_plot(ca_inset, x = 0.59, y = 0.7, width = 0.3, height = 0.3, hjust=0, vjust= 0)

ggsave(species.map.inset, filename= "Fig_1.png", width= 120, height= 120*1.1, units= "mm", dpi= 300,
       path= dir_output_fig)

