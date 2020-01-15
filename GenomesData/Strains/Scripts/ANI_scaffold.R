## Data from inStrain Comparer module

library(tidyverse)
library(ggplot2)
library(vegan)
library(pheatmap)


## Input files

## Watershed area data
dir_input_watershed <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/EnvData"
watershed.area <-
  read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2) %>% 
  rename(site= ggkbase_id)


sp1_sites <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present) %>% 
  filter(., species == "species_1") %>% 
  select(-multiple_species)

comp_sp1 <- read_tsv("inStrain/instrainComparer_NoLoc_comparisonsTable_sp1.tsv") %>% 
  mutate(name1= str_replace(name1, "^.*-vs-", ""),
         name2= str_replace(name2, "^.*-vs-", "")) %>% 
  mutate(name1= str_replace(name1, ".sorted.bam", ""),
         name2= str_replace(name2, ".sorted.bam", "")) %>% 
  left_join(sp1_sites, ., by= c("site" = "name1")) %>% # Filter only species one sites
  rename(name1= site) %>%
  select(-species) %>% 
  left_join(sp1_sites, ., by= c("site" = "name2"))  %>% 
  rename(name2= site) %>%
  select(-species) %>% 
  filter(!is.na(name1))


## Remove comparisons with percent_genome_compared < 0.75 
# this is 11.5% of all comparisons and includes many low ANI outliers
count(comp_sp1, percent_genome_compared > 0.75) %>% 
  mutate(frac= n/sum(n))

comp_sp1.F <- comp_sp1 %>% 
  filter(percent_genome_compared > 0.75) %>% 
  mutate(frac_popSNPs= population_SNPs / consensus_SNPs,
         frac_ANI= popANI  /conANI) # fraction of consensus_SNPs that are population_SNPs

# Make long for plotting
comp_sp1.FL <- comp_sp1.F %>% 
  gather(key= type, value= ANI, conANI:popANI) %>% 
  filter(!is.na(ANI))


#### DATA ANALYSIS

ani_sum.mat[28,]
test2 <- filter(ani_sum.wide, str_detect(name1, "BER"))

ani_sum <- comp_sp1.F %>% 
  group_by(name1, name2) %>% 
  summarize(scaf_num= length(scaffold),
            mean_popANI= mean(popANI),
            mean_conANI= mean(conANI),
            mean_frac_popSNPs= mean(frac_popSNPs, na.rm=TRUE),
            sd_popANI= sd(popANI),
            sd_conANI= sd(conANI),
            median_popANI= median(popANI),
            median_conANI= median(conANI)) %>% 
  left_join(., watershed.area, by= c("name1" = "site")) %>% 
  left_join(., watershed.area, by= c("name2" = "site")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y) %>% 
  left_join(., select(sp1.clusters, sp1_clustID, site), by= c("name1" = "site")) %>% 
  left_join(., select(sp1.clusters, sp1_clustID, site), by= c("name2" = "site")) %>% 
  mutate(clust_match = ifelse(sp1_clustID.x == sp1_clustID.y, "Yes", "No"))


ggplot(ani_sum, aes(x= mean_conANI, y= mean_popANI)) +
  geom_point(aes(color= clust_match)) +
  geom_abline(intercept= 0, slope= 1) +
  coord_equal() +
  theme_bw()

ggplot(ani_sum) +
  geom_boxplot(aes(x= "popANI", y= mean_popANI, fill= clust_match)) +
  geom_boxplot(aes(x= "conANI", y= mean_conANI))




# Make into a matrix
# First make wide, and add extra row and columns
ani_sum.wide <- pivot_wider(select(ani_sum, name1, name2, mean_popANI), 
            names_from = name2, 
            values_from = mean_popANI) %>% 
  ungroup() %>% 
  mutate(PH2015_01D= NA) %>% 
  select(name1, PH2015_01D, everything()) %>% 
  add_row(., name1= "PH2017_40_RAT_O_B")

# Transform into matrix
ani_sum.mat <- as.matrix(ani_sum.wide[, -1])
row.names(ani_sum.mat) <- ani_sum.wide$name1

# Add 1s for the diagonoal
diag(ani_sum.mat) <- 1

# Lower triangle all NAs. Fill in the lower triangle to make it a symmetrical matrix
ani_sum.mat[lower.tri(ani_sum.mat, diag= FALSE)] <- ani_sum.mat[upper.tri(ani_sum.mat, diag= FALSE)] 


# Calculate Euclidean distances
ani_sum.dist.row <- vegdist(ani_sum.mat, method= "bray")
ani_sum.dist.col <- vegdist(t(ani_sum.mat), method= "bray")

ani_sum.nmds <- metaMDS(ani_sum.mat, distance= "euclidean", trymax= 100, k= 3)
plot(ani_sum.nmds)

pheatmap(ani_sum.mat,
         #color= ani.colors,
         #breaks= ani.breaks,
         #legend_breaks= legend.breaks,
         #scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= ani_sum.dist.row,
         clustering_distance_cols= ani_sum.dist.col,
         clustering_method= "ward.D2",
         margins= (c(10, 10)),
         show_rownames= TRUE,
         show_colnames=  TRUE,
         #display_numbers = TRUE,
        # number_format = "%.3f",
         #fontsize_number = 4,
         #filename= fig.filename,
         #height= 10,
         #width= 15
)

pheatmap(ani_sum.mat,
         #color= ani.colors,
         #breaks= ani.breaks,
         #legend_breaks= legend.breaks,
         #scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= as.dist(1-t(ani_sum.mat)),
         clustering_distance_cols= as.dist(1-ani_sum.mat),
         clustering_method= "complete",
         margins= (c(10, 10)),
         show_rownames= TRUE,
         show_colnames=  TRUE,
         #display_numbers = TRUE,
         # number_format = "%.3f",
         #fontsize_number = 4,
         #filename= fig.filename,
         #height= 10,
         #width= 15
)




dist_conANI <- vegdist(ani_sum$mean_conANI, method= "euclidean")

dist_popANI <- vegdist(ani_sum$mean_popANI, method= "euclidean")

clust_popANI <- hclust(dist_popANI, method= "ward.D2")
clust_conANI <- hclust(dist_conANI, method= "ward.D2")

k_cut <- 4

plot(clust_popANI)
rect.hclust(clust_popANI, k= k_cut)

popANI.clusters <- data.frame(clusters= cutree(clust_popANI, k= k_cut))


#### Data Tables
high_popANI <- ani_sum %>% 
  filter(mean_popANI > 0.9982) %>%  #90 percentile
  left_join(., watershed.area, by= c("name1" = "site")) %>% 
  left_join(., watershed.area, by= c("name2" = "site")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y)

write_tsv(high_popANI, "Output_tables/popANI_90percentile.tsv")


low_popANI <- ani_sum %>% 
  filter(mean_popANI < 0.9941) %>%  #10th percentile
  left_join(., watershed.area, by= c("name1" = "site")) %>% 
  left_join(., watershed.area, by= c("name2" = "site")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y)

write_tsv(low_popANI, "Output_tables/popANI_10percentile.tsv")

table(high_popANI$watershed_1)
table(high_popANI$watershed_2)

summary(ani_sum$mean_popANI)
quantile(ani_sum$mean_popANI, probs= c(0.05, 0.1, 0.75, 0.8, 0.9, 0.95))


#### MAKE FIGURES

ggplot(comp_sp1.L, aes(x= type, y= ANI)) +
  geom_boxplot()

ggplot(comp_sp1.L, aes(x= coverage_overlap, y= ANI)) +
  geom_point(aes(color= type))

ggplot(comp_sp1.L, aes(x= percent_genome_compared, y= ANI)) +
  geom_point() +
  facet_grid(type~.)

ggplot(comp_sp1.L, aes(x= coverage_overlap, y= percent_genome_compared)) +
  geom_point() 


ggplot(comp_sp1.L, aes(x= coverage_overlap)) +
  geom_histogram() +
  facet_grid(type~.)

ggplot(comp_sp1.L, aes(x= percent_genome_compared)) +
  geom_histogram() +
  facet_grid(type~.)


ggplot(ani_sum) +
  geom_boxplot(aes(x= "popANI", y= mean_popANI, fill= clust_match)) +
  geom_boxplot(aes(x= "conANI", y= mean_conANI))


ggplot(ani_sum, aes(x= mean_conANI, y= mean_popANI)) +
  geom_point(aes(color= clust_match)) +
  geom_abline(intercept= 0, slope= 1) +
  coord_equal() +
  theme_bw()

ggplot(ani_sum, aes(x= watershed_1, y= mean_frac_popSNPs)) +
  geom_point() +
#  geom_abline(intercept= 0, slope= 1) +
 # coord_equal() +
  theme_bw()



ggplot(ani_sum, aes(x= mean_popANI)) +
  geom_histogram(binwidth = 0.0002,
                 boundary= 0,
                 color= "white")


ggplot(comp_sp1.F, aes(x= conANI, y= popANI)) +
  geom_point() 
  #geom_hex()

ggplot(comp_sp1.F, aes(x= frac_popSNPs)) +
  geom_density()

ggplot(comp_sp1.F, aes(x= frac_popSNPs)) +
  geom_histogram(binwidth = 0.01)

ggplot(comp_sp1.F, aes(x= frac_ANI)) +
  geom_histogram(binwidth= 0.0005)
