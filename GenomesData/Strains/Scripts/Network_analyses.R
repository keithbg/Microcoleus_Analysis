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


high_popANI <- read_tsv("Output_tables/popANI_90percentile.tsv")



nodes <- unique(c(high_popANI$name1, high_popANI$name2)) %>% 
  tibble() %>% 
  filter(!str_detect(., "_B$")) %>% 
  rowid_to_column(.) %>% 
  rename(id= rowid, label= ".")


edge.list <- high_popANI %>% 
  select(name1, name2) %>% 
  rename(from= name1, to= name2) %>% 
  filter(!str_detect(from, "_B$")) %>% 
  filter(!str_detect(to, "_B$"))


edges <- edge.list %>% 
  left_join(., nodes, by= c("from" = "label")) %>% 
  rename(id_from= id) %>% 
  left_join(., nodes, by= c("to" = "label")) %>%  
  rename(id_to= id) %>% 
  ungroup() %>% 
  select(-from, -to)


net_high <- tbl_graph(nodes= nodes, edges= edges, directed= FALSE)

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

ggraph(net_high) + 
  geom_edge_link(edge_width= 0.25, color= "gray50") + 
  #geom_edge_arc(edge_width= 0.25)
  geom_node_point(size= 3) + 
  geom_node_text(aes(label = label), repel = TRUE, size= 3) +
  theme_graph()



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
