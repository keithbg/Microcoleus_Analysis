## Data from inStrain Comparer module

library(tidyverse)
library(ggplot2)
library(vegan)
library(pheatmap)
library(ggsci)
library(wesanderson)
source("Scripts/ggplot_themes.R")



## INPUT FILES
## River network distance data
source("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/Scripts_cyano_metagenomics_2017/ANI_distance_format.R")
rm(ani.dist)
ani.riv.dist <- ani.dist1 %>% 
  filter(primary_cluster == 1) %>% 
  select(-ani, -year.querry, -year.reference, -primary_cluster, -reference, -querry) %>% 
  rename(name1= sample.reference, name2= sample.querry)
rm(ani.dist1)

## Watershed area data
dir_input_watershed <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/EnvData"
watershed.area <-
  read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2) %>% 
  rename(site= ggkbase_id)

## Nucleotide diversity
nuc_div <- read_tsv("Output_tables/nuc_div_summary.txt") %>% 
  filter(species == "species_1") %>% 
  select(site, mean_pi, median_pi)

## dRep results
dir_input_dRep <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","GenomesData", "dRep","Output_tables")
sp1.clusters <- read_tsv(file.path(dir_input_dRep, "sp1_clusters.tsv"))
sp1.ani <- read_tsv(file.path(dir_input_dRep, "sp1_ani.tsv")) %>% 
  rename(ani_dRep= ani) %>% 
  mutate(name1= str_replace(row_name, "_s25.*$|_Oscill.*$", ""),
         name2= str_replace(col_name, "_s25.*$|_Oscill.*$", "")) %>% 
  mutate(ani99_dRep= ifelse(ani_dRep > 0.995, ">0.995", 
                            ifelse(ani_dRep <= 0.995 & ani_dRep >=0.99, "0.99-0.995", "<0.99")))
  
## Haplotype frequencies
haplos <- read_tsv(file= "inStrain/output_tables/haplotype_freqs.tsv") %>% 
  filter(species == "species_1")



## Species 1 sites
sp1_sites <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present) %>% 
  filter(., species == "species_1") %>% 
  select(-multiple_species)

## inStrain compare results
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


## Remove comparisons with percent_genome_compared < 0.25 
# this is 3.3% of all comparisons and includes many low ANI outliers
# 183 scaffolds in the reference genomes: length(unique(comp_sp1$scaffold))

# Had origiginally filtered by removing <0.75, when reducing to 0.25 it had negligible affect
# on the genome wide averages, so I lowered cutoff to 0.25. 

comp_sp1.F <- comp_sp1 %>% 
  filter(percent_genome_compared > 0.25) %>% 
  mutate(frac_popSNPs= population_SNPs / consensus_SNPs,
         frac_ANI= popANI  /conANI) # fraction of consensus_SNPs that are population_SNPs

# comp_sp1.L <- comp_sp1 %>%
#   gather(key= type, value= ANI, conANI:popANI) %>%
#   filter(!is.na(ANI))
# 
# ggplot(comp_sp1.L, aes(x= percent_genome_compared, y= ANI)) +
#   geom_point() +
#   facet_grid(type~.)
# 
# hist(comp_sp1$percent_genome_compared)
# 
# count(comp_sp1, percent_genome_compared > 0.75) %>% 
#   mutate(frac= n/sum(n))
# 
# count(comp_sp1, percent_genome_compared > 0.25) %>% 
#   mutate(frac= n/sum(n))
# 
# low_compare <- comp_sp1 %>% 
#   filter(percent_genome_compared < 0.25) %>% 
#   mutate(compareID= str_c(name2, name1, sep= "-"))
# 
# low_compare_count <- count(low_compare, compareID) %>% 
#   mutate(frac_low= n/183) %>% 
#   arrange(-n)
# 
# hist(low_compare_count$frac_low)
# 
# hist(count(low_compare, compareID)$n)

# Filter genome

# Make long for plotting
# comp_sp1.FL <- comp_sp1.F %>%
#   gather(key= type, value= ANI, conANI:popANI) %>%
#   filter(!is.na(ANI))


# no_snvs <- comp_sp1.F %>% 
#   filter(is.nan(frac_popSNPs))
# unique(no_snvs$scaffold)


#### DATA ANALYSIS

ani_sum <- comp_sp1.F %>% 
  group_by(name1, name2) %>% 
  summarize(scaf_num= length(scaffold),
            mean_popANI= mean(popANI),
            mean_conANI= mean(conANI),
            sum_pop_sites= sum(population_SNPs),
            sum_con_sites= sum(consensus_SNPs),
            frac_popSNPs= sum_pop_sites/sum_con_sites,
            #mean_frac_popSNPs= mean(frac_popSNPs, na.rm=TRUE),
            sd_popANI= sd(popANI),
            sd_conANI= sd(conANI),
            median_popANI= median(popANI),
            median_conANI= median(conANI)) %>% 
  ungroup() %>% 
  left_join(., watershed.area, by= c("name1" = "site")) %>% 
  left_join(., watershed.area, by= c("name2" = "site")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y) %>% 
  #left_join(., select(sp1.clusters, sp1_clustID, site), by= c("name1" = "site")) %>% 
  #left_join(., select(sp1.clusters, sp1_clustID, site), by= c("name2" = "site")) %>% 
  #rename(sp1_clustID.1= sp1_clustID.x, sp1_clustID.2= sp1_clustID.y) %>% 
  #mutate(clust_pair= str_c(sp1_clustID.1, sp1_clustID.2, sep= "-"),
  #       clust_match = ifelse(sp1_clustID.1 == sp1_clustID.2, "Yes", "No")) %>% 
  # DREP ANI JOIN
  left_join(., sp1.ani) %>% 
  # RIVER DISTANCE JOIN
  left_join(., ani.riv.dist, by= c("name1", "name2")) %>% 
  left_join(., ani.riv.dist, by= c("name1" = "name2", "name2" = "name1")) %>% 
  mutate(riv_dist= ifelse(is.na(riv_dist.x), riv_dist.y, riv_dist.x),
         euc_dist= ifelse(is.na(euc_dist.x), euc_dist.y, euc_dist.x)) %>% 
  # NUCLEOTIDE DIVERSITY JOIN
  left_join(., nuc_div, by= c("name1" = "site")) %>% 
  left_join(., nuc_div, by= c("name2" = "site")) %>% 
  rename(mean_pi.1= mean_pi.x, mean_pi.2= mean_pi.y,
         median_pi.1= median_pi.x, median_pi.2= median_pi.y) %>% 
  mutate(pi_diff_mean = abs(mean_pi.1 - mean_pi.2),
         pi_avg_mean = rowMeans(cbind(.$mean_pi.1, .$mean_pi.2)),
         pi_avg_median = rowMeans(cbind(.$median_pi.1, .$median_pi.2))) %>% 
  select(-contains(".x"), -contains(".y"))

#write_tsv(ani_sum, "Output_tables/ani_summary.tsv")

ani_sum_haplo <- ani_sum %>% 
  left_join(., select(haplos, site, haplotype, freq), by= c("name1" = "site")) %>% 
  left_join(., select(haplos, site, haplotype, freq), by= c("name2" = "site")) %>% 
  rename(haplotype.1= haplotype.x, haplotype.2= haplotype.y, haplo_freq.1= freq.x, haplo_freq.2= freq.y)
  

# Make into a matrix
# First make wide, and add extra row and columns
ani_sum.wide <- pivot_wider(select(ani_sum, name1, name2, mean_conANI), 
            names_from = name2, 
            values_from = mean_conANI) %>% 
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
ani_sum.dist.row <- vegdist(ani_sum.mat, method= "euclidean")
#ani_sum.dist.col <- vegdist(t(ani_sum.mat), method= "euclidean")

#ani_sum.nmds <- metaMDS(ani_sum.mat, distance= "euclidean", trymax= 100, k= 3)
#plot(ani_sum.nmds)

pheatmap(ani_sum.mat,
         #color= ani.colors,
         #breaks= ani.breaks,
         #legend_breaks= legend.breaks,
         #scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= ani_sum.dist.row,
         clustering_distance_cols= ani_sum.dist.row,
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

dist_conANI <- vegdist(ani_sum$mean_conANI, method= "euclidean")
dist_popANI <- vegdist(ani_sum$mean_popANI, method= "euclidean")

clust_popANI <- hclust(dist_popANI, method= "ward.D2")
clust_conANI <- hclust(dist_conANI, method= "ward.D2")

k_cut <- 4

plot(clust_popANI)
rect.hclust(clust_popANI, k= k_cut)

popANI.clusters <- data.frame(clusters= cutree(clust_popANI, k= k_cut))


#### Data Tables
summary(ani_sum$mean_popANI)
quantile(ani_sum$mean_popANI, probs= c(0.05, 0.1, 0.75, 0.8, 0.9, 0.95))
quantile(ani_sum$ani_dRep, probs= c(0.05, 0.1, 0.75, 0.8, 0.9, 0.95))


high_popANI <- ani_sum %>% 
  #filter(mean_popANI > 0.9982) #%>%  #90th percentile
  filter(mean_popANI > 0.995) #%>%  
  #left_join(., watershed.area, by= c("name1" = "site")) %>% 
  #left_join(., watershed.area, by= c("name2" = "site")) %>% 
  #rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y)

write_tsv(high_popANI, "Output_tables/popANI_90percentile.tsv")


low_popANI <- ani_sum %>% 
  filter(mean_popANI < 0.9941) #%>%  #10th percentile
  #left_join(., watershed.area, by= c("name1" = "site")) %>% 
  #left_join(., watershed.area, by= c("name2" = "site")) %>% 
  #rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y)

#write_tsv(low_popANI, "Output_tables/popANI_10percentile.tsv")

table(high_popANI$watershed_1)
table(high_popANI$watershed_2)



#### MAKE FIGURES

ggplot(comp_sp1.L, aes(x= type, y= ANI)) +
  geom_boxplot()

ggplot(comp_sp1.L, aes(x= coverage_overlap, y= ANI)) +
  geom_point(aes(color= type))



ggplot(comp_sp1.L, aes(x= coverage_overlap, y= percent_genome_compared)) +
  geom_point() 


ggplot(comp_sp1.L, aes(x= coverage_overlap)) +
  geom_histogram() +
  facet_grid(type~.)

ggplot(comp_sp1.L, aes(x= percent_genome_compared)) +
  geom_histogram() +
  facet_grid(type~.)

## Genome averages ANI.SUM
ggplot(ani_sum, aes(x= mean_conANI)) +
  geom_histogram(binwidth= 0.0005, boundary= 1, fill= "black", color= "gray75") +
  labs(x= "Mean consensus ANI", y= "Count") +
  scale_x_continuous(breaks= seq(0.987, 1, by= 0.001),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  theme_strains
ggsave(last_plot(), filename = "conANI_histogram.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

ggplot(ani_sum) +
  geom_boxplot(aes(x= "popANI", y= mean_popANI, fill= clust_match)) +
  geom_boxplot(aes(x= "conANI", y= mean_conANI))

#ggplot(data= filter(ani_sum, frac_popSNPs < 0.7), aes(x= ani_dRep, y= mean_conANI, fill= frac_popSNPs)) +
ggplot(data= ani_sum, aes(x= ani_dRep, y= mean_conANI, fill= frac_popSNPs)) +
  geom_abline(intercept= 0, slope= 1, size= 0.25) +
  geom_point(size= 3, pch= 21, color= "black") +
  labs(x= "dRep ANI", y= "Mean consensus ANI") +
  # scale_x_continuous(limits= c(0.9875, 1),
  #                    breaks= seq(.9875, 1, by= 0.0025), 
  #                    labels= c("", "0.9900", "0.9925", "0.995", "0.9975", "1"),
  #                    expand= c(0, 0)) +
  # scale_y_continuous(limits= c(0.9925, 1.0005),
  #                    breaks= seq(.9925, 1, by= 0.0025), 
  #                    labels= c("", "0.9950", "0.9975", "1"),
  #                    expand= c(0.00, 0)) +
  # scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
  #                   breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  scale_fill_viridis_c(name= "popSNV sites / conSNV sites") +
  coord_equal() +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1),
        legend.position = "top") 
ggsave(last_plot(), filename = "conANI_dRep.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



ggplot(data= filter(ani_sum, clust_match == "Yes")) +
  geom_abline(intercept= 0, slope= 1, size= 0.25) +
  geom_point(data= select(ani_sum, mean_conANI, mean_popANI), aes(x= mean_conANI, y= mean_popANI), color= "gray75", size= 1) +
  geom_point(aes(x= mean_conANI, y= mean_popANI, fill= clust_pair), size= 3, pch= 21, color= "black") +
  labs(x= "Mean consensus ANI", y= "Mean population ANI") +
  scale_x_continuous(limits= c(0.9875, 1),
                     breaks= seq(.9875, 1, by= 0.0025), 
                     labels= c("", "0.9900", "0.9925", "0.995", "0.9975", "1"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.9925, 1.0005),
                     breaks= seq(.9925, 1, by= 0.0025), 
                     labels= c("", "0.9950", "0.9975", "1"),
                     expand= c(0.00, 0)) +
  scale_fill_d3(name= "Cluster \npair") +
  coord_equal() +
  facet_wrap(~clust_pair, nrow= 3) +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1))

ggsave(last_plot(), filename = "popANI_conANI.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


ggplot(data= ani_sum) +
  geom_abline(intercept= 0, slope= 1, size= 0.25) +
  geom_point(aes(x= mean_conANI, y= mean_popANI, fill= ani99_dRep), size= 3, pch= 21, color= "black") +
  labs(x= "Mean consensus ANI", y= "Mean population ANI") +
  scale_x_continuous(limits= c(0.9875, 1),
                     breaks= seq(.9875, 1, by= 0.0025), 
                     labels= c("", "0.9900", "0.9925", "0.995", "0.9975", "1"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.9925, 1.0005),
                     breaks= seq(.9925, 1, by= 0.0025), 
                     labels= c("", "0.9950", "0.9975", "1"),
                     expand= c(0.00, 0)) +
  scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
                    breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  coord_equal() +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1),
        legend.position = "top") 

ggsave(last_plot(), filename = "popANI_conANI_dRep99.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


ggplot(data= ani_sum) +
  geom_abline(intercept= 0, slope= 1, size= 0.25, linetype= "dashed") +
  geom_point(aes(x= mean_conANI, y= mean_popANI, fill= frac_popSNPs), size= 3, pch= 21, color= "black") +
  labs(x= "Mean consensus ANI", y= "Mean population ANI") +
  # scale_x_continuous(limits= c(0.9875, 1),
  #                    breaks= seq(.9875, 1, by= 0.0025), 
  #                    labels= c("", "0.9900", "0.9925", "0.995", "0.9975", "1"),
  #                    expand= c(0, 0)) +
  scale_x_continuous(limits= c(0.9875, 1),
                     breaks= seq(.9880, 1, by= 0.001), 
                     labels= c("0.988", "", "0.990", "", "0.992", "", "0.994", "", "0.996", "", "0.998", "", "1"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.9925, 1.0002),
                     breaks= seq(.993, 1, by= 0.001), 
                     labels= c("", "0.994", "", "0.996", "", "0.998", "", "1"),
                     expand= c(0.00, 0)) +
 # scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
  #                  breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  scale_fill_viridis_c(name= "% shared clonal sites",
                       direction = -1,
                       limits= c(0, 1.001),
                       breaks= seq(0, 1, by= 0.1),
                       labels= c("0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1")) +
  coord_equal() +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1),
        legend.position = "top",
        legend.key.width= unit(10, "mm")) 

ggsave(last_plot(), filename = "popANI_conANI_perSubs.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


ggplot(data= ani_sum) +
  geom_abline(intercept= 0, slope= 1, size= 0.25, linetype= "dashed") +
  geom_point(aes(x= mean_conANI, y= mean_popANI, fill= pi_avg_mean, size= pi_diff_mean), pch= 21, color= "gray30") +
  labs(x= "Mean consensus ANI", y= "Mean population ANI") +
  scale_x_continuous(limits= c(0.9875, 1),
                     breaks= seq(.9880, 1, by= 0.001), 
                     labels= c("0.988", "", "0.990", "", "0.992", "", "0.994", "", "0.996", "", "0.998", "", "1"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.9925, 1.0002),
                     breaks= seq(.993, 1, by= 0.001), 
                     labels= c("", "0.994", "", "0.996", "", "0.998", "", "1"),
                     expand= c(0.00, 0)) +
  scale_fill_viridis_c(name= "Avg.\nnucleotide diversity") +
  scale_size_continuous(name= "Difference in\nnucleotide diversity") +
  coord_equal() +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1),
        legend.position = "top")
        #legend.key.width= unit(10, "mm"),
        #legend.box= "vertical") 

ggsave(last_plot(), filename = "popANI_conANI_NucDiv.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


ggplot(data= ani_sum) +
  #geom_hline(yintercept= 0) +
  #geom_point(aes(x= mean_conANI, y= 1 - frac_popSNPs, fill= mean_popANI), size= 3, pch= 21, color= "black") +
  geom_point(aes(x= mean_conANI, y= 1 - frac_popSNPs, fill= pi_avg_mean, size= pi_diff_mean), pch= 21,  color= "black") +
  #labs(x= "Mean consensus ANI", y= "popANI sites / conANI sites") +
  labs(x= "Mean consensus ANI", y= "% shared minor alleles") +
  scale_x_continuous(limits= c(0.9875, 1),
                     breaks= seq(.9880, 1, by= 0.001), 
                     labels= c("0.988", "", "0.990", "", "0.992", "", "0.994", "", "0.996", "", "0.998", "", "1"),
                     expand= c(0, 0)) +
  # scale_y_continuous(limits= c(0, 0.84),
  #                    breaks= seq(0, 1, by= 0.1), 
  #                    #labels= c("", "0.9950", "0.9975", "1"),
  #                    expand= c(0.02, 0)) +
  # scale_y_continuous(limits= c(0, 1.02),
  #                    breaks= seq(0, 1, by= 0.1),
  #                    labels= c("0", "", "20", "", "40", "", "60", "", "80", "", "100"),
  #                    expand= c(0, 0)) +
  #scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
  #                  breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  scale_fill_viridis_c(name= "Avg.\nnucleotide diversity") +
  scale_size_continuous(name= "Difference in\nnucleotide diversity") +
  theme_strains +
  theme(#axis.text.x = element_text(angle= 45, hjust= 1),
    legend.position = "top")
    #legend.key.width= unit(10, "mm")) 

# ggsave(last_plot(), filename = "popANI_frac_dRep99.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
#        path= "Output_figures")
ggsave(last_plot(), filename = "minor_alleles_NucDiv.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")





ggplot(data= ani_sum) +
  #geom_hline(yintercept= 0) +
  #geom_point(aes(x= mean_conANI, y= 1 - frac_popSNPs, fill= mean_popANI), size= 3, pch= 21, color= "black") +
  geom_point(aes(x= mean_conANI, y= frac_popSNPs, fill= mean_popANI), size= 3, pch= 21,  color= "black") +
  #labs(x= "Mean consensus ANI", y= "popANI sites / conANI sites") +
  labs(x= "Mean consensus ANI", y= "% shared clonal sites") +
   scale_x_continuous(limits= c(0.9875, 1),
                      breaks= seq(.9880, 1, by= 0.001), 
                      labels= c("0.988", "", "0.990", "", "0.992", "", "0.994", "", "0.996", "", "0.998", "", "1"),
                      expand= c(0, 0)) +
   # scale_y_continuous(limits= c(0, 0.84),
   #                    breaks= seq(0, 1, by= 0.1), 
   #                    #labels= c("", "0.9950", "0.9975", "1"),
   #                    expand= c(0.02, 0)) +
  # scale_y_continuous(limits= c(0, 1.02),
  #                    breaks= seq(0, 1, by= 0.1),
  #                    labels= c("0", "", "20", "", "40", "", "60", "", "80", "", "100"),
  #                    expand= c(0, 0)) +
  #scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
  #                  breaks= c(">0.995", "0.99-0.995", "<0.99")) +
   scale_fill_viridis_c(name= "Mean population ANI",
                        limits= c(0.9927, 1),
                        breaks= c(0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999, 1),
                        labels= c("", "0.994", "", "0.996", "", "0.998", "", "1"),
                        expand= c(0, 0)) +
  theme_strains +
  theme(#axis.text.x = element_text(angle= 45, hjust= 1),
        legend.position = "top",
        legend.key.width= unit(10, "mm")) 

# ggsave(last_plot(), filename = "popANI_frac_dRep99.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
#        path= "Output_figures")
ggsave(last_plot(), filename = "substitution_percent.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")




ggplot(data= ani_sum) +
  #geom_point(aes(x= riv_dist, y= mean_popANI, fill= ani99_dRep), size= 3, pch= 21, color= "black") +
  geom_point(aes(x= riv_dist, y= mean_popANI, fill= ani99_dRep, size= sum_con_sites), pch= 21, color= "black") +
  labs(x= "River network distance (km)", y= "Mean population ANI") +
  geom_smooth(aes(x= riv_dist, y= mean_popANI), method= "lm") +
  #scale_fill_npg(name= "dRep ANI") +
  scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
                    breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  scale_x_continuous(limits= c(0, 350000),
                     breaks= seq(0, 350000, by= 25000),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5000)) +
  theme_strains +
  theme(legend.position = "top")
ggsave(last_plot(), filename = "popANI_rivDist.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

ggplot(data= ani_sum) +
  geom_point(aes(x= riv_dist, y= frac_popSNPs, fill= ani99_dRep), size= 3, pch= 21, color= "black") +
  labs(x= "River network distance (km)", y= "Mean population ANI") +
  #scale_fill_npg(name= "dRep ANI") +
  scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
                    breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  # scale_x_continuous(limits= c(0, 350000),
  #                    breaks= seq(0, 350000, by= 25000),
  #                    labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
  #                    expand= c(0, 5000)) +
  scale_x_log10() +
  annotation_logticks() +
  theme_strains +
  theme(legend.position = "top")
ggsave(last_plot(), filename = "popANI_frac_rivDist.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

ggplot(data= ani_sum, aes(x= riv_dist, y= mean_conANI)) +
  #geom_point(aes(x= riv_dist, y= mean_popANI, fill= ani99_dRep), size= 3, pch= 21, color= "black") +
  geom_hline(yintercept = 0.9935, linetype= "dashed") +
  geom_point(size= 3, pch= 21, fill= "gray75", color= "black") +
  #geom_smooth(color= "black", fill= "dodgerblue") +
  geom_smooth(color= "black", fill= "dodgerblue", method= "lm") +
  labs(x= "River network distance (km)", y= "Mean consensus ANI") +
  #scale_fill_npg(name= "dRep ANI") +
  #scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
  #                  breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  scale_x_continuous(limits= c(0, 350000),
                     breaks= seq(0, 350000, by= 25000),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5000)) +
  theme_strains +
  theme(legend.position = "top")
ggsave(last_plot(), filename = "conANI_rivDist.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



ggplot(data= ani_sum) +
  geom_abline(intercept= 0, slope= 1, size= 0.25) +
  geom_point(aes(x= mean_conANI, y= mean_popANI, fill= riv_dist/1000), size= 3, pch= 21, color= "black") +
  labs(x= "Mean consensus ANI", y= "Mean population ANI") +
  scale_x_continuous(limits= c(0.9875, 1),
                     breaks= seq(.9875, 1, by= 0.0025), 
                     labels= c("", "0.9900", "0.9925", "0.995", "0.9975", "1"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.9925, 1.0005),
                     breaks= seq(.9925, 1, by= 0.0025), 
                     labels= c("", "0.9950", "0.9975", "1"),
                     expand= c(0.00, 0)) +
    scale_fill_viridis_c(name= "dRep ANI") +
  coord_equal() +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1))


ggplot(ani_sum_haplo, aes(x=haplo_freq.1, y= mean_popANI)) +
  geom_point() +
  facet_wrap(.~haplotype.1, ncol= 4, scales= "free_x") +
  geom_smooth(aes(color= haplotype.1)) +
  theme_strains



ggplot(ani_sum, aes(x= watershed_1, y= frac_popSNPs)) +
  geom_point() +
#  geom_abline(intercept= 0, slope= 1) +
 # coord_equal() +
  theme_bw()



ggplot(ani_sum, aes(x= mean_popANI)) +
  geom_histogram(binwidth = 0.0002,
                 boundary= 0,
                 color= "white")

