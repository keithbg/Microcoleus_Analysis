## Script to format and plot dRep ANI values


#### Libraries #################################################################
library(tidyverse)
library(vegan)
library(pheatmap)
library(grid)
################################################################################


## READ dRep OUTPUT
df <- read_csv(file.path("Data/dRep", "Ndb_20200422.csv")) %>%
  select(querry, reference, ani, alignment_coverage) %>% 
  mutate(querry= str_replace(querry, ".fa", ""),
         reference= str_replace(reference, ".fa", ""))

## Transform into a matrix
ani.mat <- df %>%
  select(querry, reference, ani) %>% 
  spread(key= querry, value= ani) %>%
  select(-reference) %>%
  as.matrix(.)
rownames(ani.mat) <- colnames(ani.mat)


## Remove NAs (these are generated because only samples within the same primary cluster have an NA value calculated)
ani.mat[is.na(ani.mat) == TRUE] <- 0.8

# Calculate Euclidean distance
ani.dist.row <- vegdist(ani.mat, method= "euclidean")

## Identify the 4 species
ani.clusters <- hclust(ani.dist.row, method = "ward.D2")
ani.species <- data.frame(ani_species= cutree(ani.clusters, k= 4))
ani.species$genome <- str_replace(row.names(ani.species), ".fa", "")

## Change the species numbers to match the PH2015 publication and write a TSV file
ani.species <- ani.species %>%
  mutate(ani_species= ifelse(ani_species == 1, 4,
                             ifelse(ani_species == 2, 1,
                                    ifelse(ani_species == 3, 2, 3)))) %>%
  select(genome, ani_species)
#write_tsv(ani.species, path= file.path("Data/dRep", "ani_species.tsv"))



#### SUMMARY STATISTICS BETWEEEN SPECIES ####
ani.df <- ani.mat %>% 
  as.data.frame() %>% 
  mutate(genome= str_replace_all(row.names(ani.mat), ".fa", "")) %>% 
  as_tibble() %>% 
  left_join(., ani.species) %>% 
  rename(genome_1= genome, species_1= ani_species)

## Make long
ani.df.l <- ani.df %>% 
  pivot_longer(names_to= "genome_2", values_to = "dRep_ANI", cols= contains("PH")) %>% 
  mutate(genome_2= str_replace_all(genome_2, ".fa", "")) %>% 
  left_join(., ani.species, by= c("genome_2" = "genome")) %>% 
  rename(species_2 = ani_species) %>% 
  select(genome_1, species_1, genome_2, species_2, dRep_ANI) %>% 
  mutate(species_comp= str_c(species_1, species_2, sep= "-")) %>% 
  left_join(., select(df, querry, reference, alignment_coverage), by= c("genome_1" = "querry", "genome_2" = "reference"))


## Summarize ANI values between species clusters
ani.species.comparison <- ani.df.l %>% 
  group_by(species_comp) %>% 
  summarize(n= length(dRep_ANI),
            mean_ani= mean(dRep_ANI),
            median_ani= median(dRep_ANI),
            sd_ani= sd(dRep_ANI),
            mean_cov= mean(alignment_coverage),
            median_cov= median(alignment_coverage),
            sd_cov= sd(alignment_coverage))

#### FIGURE #####

## Color ramp
ani.breaks <- c(0.85, seq(0.90, 1, by= 0.005))
legend.breaks <- c(0.85, seq(0.90, 1, by= 0.02))
legend.labels <- as.character(legend.breaks*100)
ani.colors <- c("gray70",  colorRampPalette(c("cornsilk", "snow", "goldenrod", "salmon", "tomato", "firebrick"))(21))


ani.heat <- pheatmap(ani.mat,
         color= ani.colors,
         breaks= ani.breaks,
         border_color = "gray60",
         legend_breaks= legend.breaks,
         legend_labels = legend.labels,
         scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= ani.dist.row,
         clustering_distance_cols= ani.dist.row,
         clustering_method= "ward.D2",
         margins= (c(10, 10)),
         show_rownames= TRUE,
         show_colnames=  FALSE,
         fontsize_row = 6,
)

png("Output_figures/Fig_S2b.png", height= 180, width= 180, units= "mm", res= 320, type= "cairo")
ani.heat

grid.text("sp. 1", x= 0.365, y= 0.475, just= "center", gp=gpar(fontsize=12, fontface= "bold", col= "white"))
grid.text("sp. 2", x= 0.148, y= 0.828, just= "center", gp=gpar(fontsize=12, fontface= "bold", col= "white"))
grid.text("sp. 3", x= 0.595, y= 0.09, just= "center", gp=gpar(fontsize=12, fontface= "bold", col= "white"))
grid.text("dRep ANI (%)", x= 0.99, y= 0.92, just= "right", gp=gpar(fontsize=10, col= "Black"))
dev.off()
