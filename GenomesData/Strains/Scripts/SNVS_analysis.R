## Investigate the SVNs.tsv file output

library(tidyverse)
library(vegan)
library(pheatmap)
library(ggplot2)

#### INPUT FILES ####
in_dir <- "inStrain/inStrain_gene_profile_output"

## River network distance between sites in meters
river_dist_mat <- read_tsv("../../EnvData/Distance_RiverNetwork_meters.tsv") %>% 
  rename(site= Site)

# Make river distances long
river_dist_df <- river_dist_mat %>% 
  rename(siteA= site) %>% 
  gather(key= siteB, value= riv_dist, -siteA)

## Distance metadata
source("../../Scripts_cyano_metagenomics_2017/ANI_distance_format.R")
rm(ani.dist1)

ani.dist.df <- ani.dist %>% 
  select(-querry, -reference, -euc_dist, -riv_dist, -primary_cluster, -ani) %>% 
  rename(siteA= sample.querry, forkA= fork.querry, speciesA= species.querry, yearA= year.querry, 
         forkB= fork.reference, siteB= sample.reference, speciesB= species.reference, yearB= year.reference)

# Watershed areas
watershed_area <- read_tsv("../../EnvData/PhormMeta17_WatershedArea_Combined.tsv") %>% 
  select(ggkbase_id, watershed_km2) %>% 
  mutate(siteA= ggkbase_id,
         watershed_km2_A= watershed_km2,
         siteB= ggkbase_id,
         watershed_km2_B= watershed_km2) %>% 
  select(-ggkbase_id, -watershed_km2)

# Table of species recovered at each sampling site
species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

snv_files <- list.files(in_dir, pattern= "species_1.pid96.SNVs.tsv")

## Read in all SNVs.tsv files
snv_list <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                  mutate(pos_id= str_c(scaffold, position, sep= "-")) %>% 
                  select(pos_id, refFreq)) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))

snv_list_morph234 <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                  mutate(pos_id= str_c(scaffold, position, sep= "-")) %>% 
                  filter(morphia != 1) %>% 
                  select(pos_id, conBase)) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))

snv_list_conBase <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                  mutate(pos_id= str_c(scaffold, position, sep= "-")) %>% 
                  select(pos_id, conBase)) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))


## Transform into data frame
snv_df <- bind_rows(snv_list, .id= "sample") %>% 
  mutate(site= do.call(rbind, str_split(.$sample, "\\."))[, 1],
         species= do.call(rbind, str_split(.$sample, "\\."))[, 2]) %>% 
  inner_join(species_lookup, ., by= c("site", "species")) # Filter to only include species recovered from each site

snv_df_conBase_morph234 <- bind_rows(snv_list_morph234, .id= "sample") %>% 
  mutate(site= do.call(rbind, str_split(.$sample, "\\."))[, 1],
         species= do.call(rbind, str_split(.$sample, "\\."))[, 2]) %>% 
  inner_join(species_lookup, ., by= c("site", "species")) %>% # Filter to only include species recovered from each site
  mutate(conBase= ifelse(conBase == "A", 1,
                         ifelse(conBase == "T", 2,
                                ifelse(conBase == "C", 3, 4))))

snv_df_conBase <- bind_rows(snv_list_conBase, .id= "sample") %>% 
  mutate(site= do.call(rbind, str_split(.$sample, "\\."))[, 1],
         species= do.call(rbind, str_split(.$sample, "\\."))[, 2]) %>% 
  inner_join(species_lookup, ., by= c("site", "species")) %>% # Filter to only include species recovered from each site
  mutate(conBase= ifelse(conBase == "A", 1,
                         ifelse(conBase == "T", 2,
                                ifelse(conBase == "C", 3, 4))))



## Make wide matrix of refFreq values
snv_matrix_refFreq <- snv_df %>% 
  select(site, pos_id, refFreq) %>% 
  spread(key= pos_id, value= refFreq, fill= 0)

snv_matrix_conBase <- snv_df_conBase %>% 
  select(site, pos_id, conBase) %>% 
  spread(key= pos_id, value= conBase, fill= 0)

snv_matrix_conBase_morph234 <- snv_df_conBase_morph234 %>% 
  select(site, pos_id, conBase) %>% 
  spread(key= pos_id, value= conBase, fill= 0)

## Make wide matrix of SNV position presence/absence 0 and 1 values
snv_matrix_binary <- snv_df %>% 
  select(site, pos_id, refFreq) %>% 
  mutate(refFreq= ifelse(refFreq > 0, 1, 0)) %>% 
  spread(key= pos_id, value= refFreq, fill= 0)

snv_matrix_binary_morph234 <- snv_df_conBase_morph234 %>% 
  select(site, pos_id, conBase) %>% 
  mutate(conBase= ifelse(conBase > 0, 1, 0)) %>% 
  spread(key= pos_id, value= conBase, fill= 0)

# Write files
write_tsv(snv_matrix_refFreq, path= "inStrain/output_tables/snv_pos_refFreq_sp1.tsv")
write_tsv(snv_matrix_binary, path= "inStrain/output_tables/snv_pos_binary_sp1.tsv")
write_tsv(snv_matrix_binary_morph234, path= "inStrain/output_tables/snv_pos_binaryM234_sp1.tsv")
write_tsv(snv_matrix_conBase, path= "inStrain/output_tables/snv_pos_conBase_sp1.tsv")
write_tsv(snv_matrix_conBase_morph234, path= "inStrain/output_tables/snv_pos_conBaseM234_sp1.tsv")



## SNV sharing
snv_sharing <- colSums(snv_matrix_binary_morph234[, -1])
table(snv_sharing)

## Calculate distances between sites
dist.rows <- vegdist(as.matrix(snv_matrix_binary_morph234[, -1]), method= "jaccard")

#dist.rows <- vegdist(as.matrix(snv_matrix_binary[, -1]), method= "jaccard")
snv_dist_mat <- as.matrix(dist.rows)

#colnames(snv_dist_mat) <- snv_matrix_binary$site
colnames(snv_dist_mat) <- snv_matrix_binary_morph234$site
rownames(snv_dist_mat) <- snv_matrix_binary_morph234$site
write_tsv(as.data.frame(snv_dist_mat), path= "inStrain/output_tables/snv_pos_M234_jaccard.tsv")
snv_dist_mat <- read_tsv("inStrain/output_tables/snv_pos_M234_jaccard.tsv")

head(as.data.frame(snv_dist_mat))
?write_csv


snv_nmds <- metaMDS(as.matrix(snv_matrix_binary_morph234[, -1]), distance= "jaccard", trymax= 100, k= 3)

## Set up data to plot NMDS in ggplot
data.scores <- as.data.frame(scores(snv_nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(snv_dist_mat)  # create a column of site names, from the rownames of data.scores

data.scores <- merge(data.scores, nonCyano_ward_groups_df)
data.scores$group <- as.factor(data.scores$group)
data.scores$site <- gsub(".*_", "", data.scores$site)
head(data.scores, 9)  #look at the data

species.scores <- as.data.frame(scores(snv_nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

# Without rpS3 gene information
ggplot(data= data.scores, aes(x=NMDS1,y=NMDS2)) +
  #geom_point(aes(color= atx_detect, shape= atx_detect), size= 3) +
  geom_point() +
 # ggrepel::geom_text_repel(aes(label=site), size= 3, show.legend= FALSE) +
  #annotate("text", label= "Stress= 0.085", x= -2.9, y= -1.2, size= 3) +
  # scale_color_manual(values= c("black", "#fc0000"), labels= c("No", "Yes"), name=
  #                      "ATX operon") +
  # scale_shape_manual(values= c(16, 15), labels= c("No", "Yes"), name=
  #                      "ATX operon") +
  coord_equal() +
  theme_strains + theme(plot.margin = margin(1, 1, 1, 1, "cm"),
                    legend.position = "top", panel.border = element_blank())





plot(snv_nmds, display= "sites")

# Make SNV_distance long and merge with river distance  
snv_dist_df <- snv_dist_mat %>% 
  as_tibble() %>%
  mutate(siteA= names(.)) %>% 
  gather(key= siteB, value= snv_distance, -siteA) %>% 
  left_join(., river_dist_df) %>% 
  #inner_join(., ani.dist.df) %>% 
  left_join(., select(watershed_area, siteA, watershed_km2_A)) %>% 
  left_join(., select(watershed_area, siteB, watershed_km2_B))




## Cluster sites
row_hclust <- hclust(dist.rows, method= "ward.D2")

row_hclust[["labels"]] <- snv_matrix_binary$site


pheatmap(snv_dist_mat,
                   #color= c("light gray", "black"),
                   #breaks= c(0, 0.5, 1),
                   #legend_breaks= c(0, 1),
                   #legend_labels= c("Absent", "Present"),
                   #scale= "none",
                   cluster_cols = TRUE,
                   cluster_rows = TRUE,
                   clustering_distance_rows= dist.rows,
                   clustering_distance_cols= dist.rows,
                   clustering_method= "ward.D2",
                   #cutree_cols = 10,
                   #cutree_rows= 5,
                   #annotation_col = col_clusters,
                   #annotation_colors = ann_colors,
                   #annotation_names_col = TRUE,
                   #margins= (c(10, 10)),
                   #display_numbers = TRUE,
                   show_rownames= TRUE,
                   show_colnames=  FALSE
                   #filename= fig.filename,
                   #height= 25,
                   #width= 40
)


# Cut clusters
# row_clusters <- cutree(row_hclust, k= 10) %>%
#   data.frame(family= names(.), cluster= ., stringsAsFactors = FALSE) %>%
#   mutate(cluster= as.character(cluster))
# row.names(col_clusters) <- col_clusters$family
# col_clusters <- col_clusters %>%
#   select(-family)


### PLOT SNV AND SPATIAL DISTANCES


  


#### MAKE FIGURES ####
source("Scripts/ggplot_themes.R")

snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  #filter(riv_dist < 1000) %>% 
ggplot(., aes(x= riv_dist, y= snv_distance)) +
  labs(x= "River network distance (meters)", y= "SNV position Jaccard dissimilarity") +
  #geom_point(aes(size= watershed_diff, fill= fork.match), pch= 21, color= "gray50") +
  geom_point(aes(size= watershed_diff), pch= 21, color= "gray50", fill= "black") +
  geom_smooth() +
  scale_y_continuous(limits= c(0, 1.01), 
                     expand= c(0, 0)) +
  # scale_x_log10(breaks= c(0.01, 0.1, 1, 100),
  #               labels= c("0.01", "0.1", "1", "100")) +
  scale_x_log10(limits= c(0.04, 400000),
                breaks= c(0.1, 1, 10, 100, 1000, 10000, 100000),
                labels= c("0.1", "1", "10", "100", "1,000", "10,000", "100,000"),
                expand= c(0, 0)) +
  annotation_logticks(side= "b") +
  #scale_color_discrete(name= "Region match") +
  #scale_fill_discrete(name= "Region match") +
  scale_size_continuous(name= expression('Watershed area difference (km'^"2"*")")) +
  theme_strains +
  theme(legend.position = "top",
        legend.title= element_text(size= 10))

ggsave(last_plot(), filename = "snv_diss_rivDist_M234.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

ggplot() +
  geom_histogram(aes(x= snv_dist_df$riv_dist), binwidth= 10000)


test <- snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  filter(watershed_diff < 400) %>% 
  filter(riv_dist < 5000 & fork.match == "Y" & (siteA != "PH2017_02_FOX_O_A" & siteB != "PH2017_02_FOX_O_A"))


facet_anno2 <- data.frame(x1= rep(1300, 2),
                         y1= rep(0.95, 2),
                         lab= c("Site A", "Site B"),
                         watershed_site= c("watershed_km2_A", "watershed_km2_B"),
                         stringsAsFactors = FALSE)


snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  filter(watershed_diff < 400) %>% 
  #filter(riv_dist < 5000 & fork.match == "Y") %>% 
  filter(riv_dist < 5000 & fork.match == "Y" & (siteA != "PH2017_02_FOX_O_A" & siteB != "PH2017_02_FOX_O_A")) %>% 
  gather(key= watershed_site, value= watershed_km2, watershed_km2_A:watershed_km2_B) %>% 
  ggplot(., aes(x= watershed_km2, y= snv_distance)) +
  labs(x= expression('Watershed area (km'^"2"*")"), y= "SNV position Jaccard dissimilarity") +
  geom_point(pch= 21, fill= "black", color= "gray50", size= 3) +
  geom_text(data= facet_anno2, aes(x= x1, y= y1, label= lab), size= 4) +
 # geom_smooth(method= "lm", se= FALSE) +
  scale_y_continuous(limits= c(0, 1), 
                     expand= c(0, 0)) +
  scale_x_log10(limits= c(10, 2000),
                expand= c(0,0)) +
  annotation_logticks(side= "b") +
  #scale_color_discrete(name= "Region match") +
  #scale_fill_discrete(name= "Region match") +
  scale_size_continuous(name= expression('Watershed area difference (km'^"2"*")")) +
  facet_grid(watershed_site~.) +
  theme_strains +
  theme(strip.text = element_blank(),
        legend.position = "top",
        legend.title= element_text(size= 10))

ggsave(last_plot(), filename = "snv_diss_km2_M234.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  mutate(watershed_diff= abs(watershed_km2_A - watershed_km2_B)) %>% 
  ggplot(., aes(x= watershed_km2_A, y= snv_distance)) +
  geom_point(aes(fill= riv_dist/1000), shape= 21, color= "gray50", size= 3) +
  geom_smooth(se= FALSE, method= "lm") +
  labs(x= expression('Watershed area (km'^"2"*")"), y= "SNV position Jaccard dissimilarity") +
  #scale_fill_viridis_c(name= expression('Watershed difference (km'^"2"*")")) +
  scale_fill_viridis_c(name= "River network distance (km)") +
  scale_x_log10(limits=c (1, 2200),
                expand= c(0,0)) +
  annotation_logticks(sides= "b") +
  scale_y_continuous(limits= c(0, 1), 
                     expand= c(0, 0)) +
  theme_strains +
  theme(legend.position = "top")

ggsave(last_plot(), filename = "snv_diss_km2_RivDist _M234.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
      path= "Output_figures")


fit <- lm(log(snv_distance) ~ watershed_km2_A*riv_dist,  data= filter(snv_dist_df, siteA != siteB))
summary(fit)
anova(fit)
plot(fit)


lowJaccard_highRivDist <- snv_dist_df %>% 
  filter(siteA != siteB, snv_distance < 0.5, riv_dist>10000)

snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  #mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  filter(riv_dist  < 1000) %>% 
  ggplot(., aes(x= watershed_diff, y= snv_distance)) +
  labs(x= "Watershed area difference (km2)", y= "SNV position Jaccard distance") +
  geom_point(aes(color= riv_dist)) +
  geom_smooth() +
  scale_y_continuous(limits= c(0, 1), 
                     expand= c(0.01, 0)) +
  # scale_x_log10(breaks= c(0.01, 0.1, 1, 100),
  #               labels= c("0.01", "0.1", "1", "100")) +
  # scale_x_log10(limits= c(0.04, 400000),
  #               breaks= c(0.1, 1, 10, 100, 1000, 10000, 100000),
  #               labels= c("0.1", "1", "10", "100", "1,000", "10,000", "100,000"),
  #               expand= c(0, 0)) +
  # annotation_logticks(side= "b") +
  theme_strains



snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  #mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  filter(riv_dist  < 1000) %>% 
  ggplot(., aes(x= watershed_km2_A, y= snv_distance)) +
  #labs(x= "Watershed area (km2)", y= "SNV position Jaccard distance") +
  geom_point(aes(size= watershed_diff)) +
  geom_smooth(se= FALSE) +
  scale_y_continuous(limits= c(0, 1), 
                     expand= c(0.01, 0)) +
  theme_strains




snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  #mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  filter(watershed_diff  < 500) %>% 
  ggplot(., aes(x= watershed_km2_B, y= snv_distance)) +
  #labs(x= "Watershed area (km2)", y= "SNV position Jaccard distance") +
  geom_point(aes(size= watershed_diff, color= fork.match)) +
  geom_smooth(aes(color= fork.match), se= FALSE) +
  scale_y_continuous(limits= c(0, 1), 
                     expand= c(0.01, 0)) +
  scale_x_log10() +
  annotation_logticks(side= "b") +
  theme_strains





snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  #mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  #filter(riv_dist  < 1000) %>% 
  ggplot(., aes(x= watershed_km2_B, y= watershed_diff, z= snv_distance)) +
  #labs(x= "Watershed area (km2)", y= "SNV position Jaccard distance") +
  geom_point(aes(color= snv_distance), size= 3) +
  scale_color_gradient2(midpoint= 0.5) +
  theme_strains




snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  filter(snv_distance  < 0.5) %>% 
  ggplot(., aes(x= watershed_diff, y= snv_distance)) +
  #labs(x= "Watershed area (km2)", y= "SNV position Jaccard distance") +
  geom_point(aes(color= fork.comparison), size= 3) +
  #scale_color_gradient2(midpoint= 0.5) +
  theme_strains




snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  filter(snv_distance  < 0.5) %>% 
  ggplot(., aes(x= watershed_km2_A, y= snv_distance)) +
  #labs(x= "Watershed area (km2)", y= "SNV position Jaccard distance") +
  geom_point(aes(color= fork.comparison), size= 3) +
  #scale_color_gradient2(midpoint= 0.5) +
  theme_strains


test <- snv_dist_df %>% 
  filter(siteA != siteB) %>% 
  #mutate(riv_dist= ifelse(riv_dist == 0, 0.05, riv_dist)) %>% 
  filter(watershed_diff  < 500) %>% 
  gather(key= watershed_site, value= watershed_km2, watershed_km2_A:watershed_km2_B)


