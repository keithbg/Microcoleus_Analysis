library(tidyverse)
library(tidygraph)
library(ggplot2)
library(ggraph)

## Watershed area data
dir_input_watershed <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/EnvData"
watershed.area <-
  read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2) %>% 
  rename(site= ggkbase_id)

sp1.clusters <- read_tsv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/dRep/Output_tables/sp1_clusters.tsv")


# high_popANI <- read_tsv("Output_tables/popANI_90percentile.tsv") %>% 
#   filter(ani_dRep > 0.99)
high_dRep <- read_tsv("Output_tables/ani_summary.tsv") %>% 
  filter(ani_dRep > 0.99)

df= high_dRep
format_network <- function(df){

nodes <- unique(c(df$name1, df$name2)) %>% 
  tibble() %>% 
  #filter(!str_detect(., "_B$")) %>% 
  rowid_to_column(.) %>% 
  rename(id= rowid, label= ".")


edge.list <- df %>% 
  select(name1, name2)
 # rename(from= name1, to= name2)# %>% 
 # filter(!str_detect(from, "_B$")) %>% 
 # filter(!str_detect(to, "_B$"))


edges <- edge.list %>% 
  left_join(., nodes, by= c("name1" = "label")) %>% 
  rename(id_from= id) %>% 
  left_join(., nodes, by= c("name2" = "label")) %>%  
  rename(id_to= id) %>% 
  ungroup() %>%
  left_join(., select(df, name1, name2, riv_dist)) %>% 
  select(-name1, -name2)


network_object <- tbl_graph(nodes= nodes, edges= edges, directed= FALSE) %>% 
  activate(nodes) %>%
  mutate(degree = centrality_degree())
  
return(network_object)
}

high_dRep_network <- format_network(high_dRep)

node_degreeness_high <- high_dRep_network %>%
  activate(nodes) %>%
  mutate(importance = centrality_degree()) %>% 
  arrange(-importance) %>% 
  as_tibble() %>% 
  left_join(., watershed.area, by= c("label" = "site")) %>% 
  rename(importance_high = importance)


ggraph(high_dRep_network, layout= 'igraph', algorithm= "kk") + 
  #geom_edge_link(edge_width= 0.25, color= "gray50") + 
  geom_edge_link(edge_width= 1, aes(color= log10(riv_dist))) + 
  geom_node_point(size= 1) + 
  #geom_node_point(aes(size= as.character(degree))) + 
  geom_node_text(aes(label = label), repel = TRUE, size= 2) +
  scale_edge_color_viridis(option= "plasma",
                           label= c("10", "100", "1,000", "10,000", "100,000"),
                           name= "River distance (m)") +
  #scale_color_continuous() +
  #scale_color_manual(name="Type of Activity",guide=guide_legend(override.aes=list(label=c("w","s"))),
  #                   values=c("red","blue"),labels=c("Work","School"))
  theme_graph()

ggsave(last_plot(), filename = "network_dRep99.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "../dRep/Output_figures")


?ggraph

ggraph(net_high) + 
  geom_edge_link(edge_width= 0.25, color= "gray50") + 
  #geom_edge_arc(edge_width= 0.25)
  geom_node_point(size= 3) + 
  geom_node_text(aes(label = label), repel = TRUE, size= 3) +
  theme_graph()







node_degreeness_high <- net_high %>%
  activate(nodes) %>%
  mutate(importance = centrality_degree()) %>% 
  arrange(-importance) %>% 
  as_tibble() %>% 
  left_join(., watershed.area, by= c("label" = "site")) %>% 
  rename(importance_high = importance)



plot(importance_high ~ watershed_km2, data= node_degreeness_high)

high_popANI %>% 
  pivot_longer(names_to= "watershed", values_to= "km2", contains("watershed")) %>% 
  ggplot() +
  geom_histogram(aes(x= km2), binwidth= 50, boundary= 0)
plot(net)





nodes_low <- unique(c(low_popANI$name1, low_popANI$name2)) %>% 
  tibble() %>% 
  rowid_to_column(.) %>% 
  rename(id= rowid, label= ".")


edge.list.low <- low_popANI %>% 
  select(name1, name2) %>% 
  rename(from= name1, to= name2)

edges_low <- edge.list.low %>% 
  left_join(., nodes_low, by= c("from" = "label")) %>% 
  rename(id_from= id) %>% 
  left_join(., nodes_low, by= c("to" = "label")) %>%  
  rename(id_to= id) %>% 
  ungroup() %>% 
  select(-from, -to)


net_low <- tbl_graph(nodes= nodes_low, edges= edges_low, directed= FALSE)

node_degreeness_low <- net_low %>%
  activate(nodes) %>%
  mutate(importance = centrality_degree()) %>% 
  arrange(-importance) %>% 
  as_tibble() %>% 
  left_join(., watershed.area, by= c("label" = "site")) %>% 
  rename(importance_low = importance)

ggraph(net_low) + 
  geom_edge_link(edge_width= 0.25) + 
  geom_node_point(size= 3) + 
  geom_node_text(aes(label = label), repel = TRUE, size= 4) +
  theme_graph()

plot(importance ~ watershed_km2, data= node_degreeness_low)


node_degreeness <- full_join(select(node_degreeness_high, -id), select(node_degreeness_low, -id))

plot(importance_high ~ importance_low, data= node_degreeness)

sum(node_degreeness_high$label %in% node_degreeness_low$label)
sum(node_degreeness_low$label %in% node_degreeness_high$label)

length(node_degreeness_high$label)


## Make base map of Eel and Russian River watersheds
dir_input_script <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis", "Scripts_cyano_metagenomics_2017")
source(file.path(dir_input_script, "Map_eel_russian.R"))
# Returns R object: PH2017_eel_base_map

## Read in data
dir_latlong <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","EnvData")
latlong <- read_csv(file.path(dir_latlong, "PhormMeta17_LatLong_combined.csv")) %>%
  mutate(year= as.character(year))

high_dRep_map <- high_dRep %>% 
  select(name1, name2) %>% 
  mutate(nodes_id= str_c(name1, name2, sep= "-")) %>% 
  pivot_longer(names_to = "node", values_to = "ggkbase_id", name1:name2) %>% 
  left_join(., latlong) %>% 
  left_join(., node_degreeness_high, by= c("ggkbase_id" = "label"))

high_dRep_map_labels <- high_dRep_map %>% 
  select(long, lat, ggkbase_id) %>% 
  distinct()

PH2017_eel_base_map +
  geom_point(data= latlong, aes(x= long, y= lat), size= 3, pch=21, fill= "gray75", color= "black") +
  geom_line(data= high_dRep_map, aes(x= long, y= lat, group= nodes_id), color= "black", size= 0.3) +
  geom_point(data= high_dRep_map, aes(x= long, y= lat, fill= as.character(importance_high)), size= 4, pch= 21, color= "black") +
  #geom_label_repel(data= high_popANI_map_labels, aes(x= long, y= lat, label= ggkbase_id)) +
  scale_fill_viridis_d(name= "Node degreeness", direction= 1, begin= 0.2, option= "magma") +
  scale_size_continuous(guide = FALSE) +
  PH2017_map_theme
ggsave(last_plot(), filename = "network_map_dRep99.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "../dRep/Output_figures")




high_popANI_map <- high_popANI %>% 
  select(name1, name2) %>% 
  mutate(nodes_id= str_c(name1, name2, sep= "-")) %>% 
  pivot_longer(names_to = "node", values_to = "ggkbase_id", name1:name2) %>% 
  left_join(., latlong) %>% 
  left_join(., node_degreeness_high, by= c("ggkbase_id" = "label"))


high_popANI_map_labels <- high_popANI_map %>% 
  select(long, lat, ggkbase_id) %>% 
  distinct()


 hist(node_degreeness_high$importance_high)



sp1.clusters.map <- sp1.clusters %>% 
  left_join(., latlong, by= c("site" = "ggkbase_id")) %>% 
  mutate(sp1_clustID = as.character(sp1_clustID)) %>% 
  as_tibble()

PH2017_eel_base_map +
  geom_point(data= filter(sp1.clusters.map, sp1_clustID == "3"), aes(x= long, y= lat, fill= sp1_clustID), size=3, pch= 21, color= "black") +
  #geom_label_repel(data= high_popANI_map_labels, aes(x= long, y= lat, label= ggkbase_id)) +
  #scale_fill_viridis(name= "Node degreeness") +
  #scale_size_continuous(name= "Node degreeness") +
  #facet_wrap(~as.character(sp1.clusters.map$sp1_clustID), nrow= 3, scales= "free") +
  PH2017_map_theme






