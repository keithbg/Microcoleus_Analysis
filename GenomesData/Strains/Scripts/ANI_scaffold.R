## Data from inStrain Comparer module

library(tidyverse)
library(ggplot2)
library(vegan)
library(pheatmap)
library(ggsci)
library(wesanderson)
library(cowplot)
library(ggpubr)
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
watershed.area <- read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
 # rename(site= ggkbase_id) %>% 
  select(ggkbase_id, watershed_km2)



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

## SNV data and young/old populations
snv_genomes <- read_tsv("Output_tables/snvs_genome_summary.tsv") %>% 
  filter(species == "species_1") %>% 
  select(ggkbase_id, SNV_mbp, pop_age)



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
  # WATERSHED AREA JOIN
  left_join(., watershed.area, by= c("name1" = "ggkbase_id")) %>% 
  left_join(., watershed.area, by= c("name2" = "ggkbase_id")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y) %>% 
  # DREP ANI JOIN
  left_join(., sp1.ani) %>% 
  # SNV_MBP JOIN
  left_join(., snv_genomes, by= c("name1" = "ggkbase_id")) %>% 
  left_join(., snv_genomes, by= c("name2" = "ggkbase_id")) %>% 
  rename(pop_age.1= pop_age.x, pop_age.2= pop_age.y,
         SNV_mbp.1= SNV_mbp.x, SNV_mbp.2= SNV_mbp.y) %>% 
  mutate(pop_age_pair= str_c(pop_age.1, pop_age.2, sep= "-")) %>%
  mutate(pop_age_pair= ifelse(pop_age_pair == "Young-Old", "Old-Young", pop_age_pair)) %>% 
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
  mutate(watershed_diff= abs(watershed_2 - watershed_1)) %>% 
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
 ani_sum %>% 
   count(mean_conANI > 0.9935)

conANI_hist <- ggplot(ani_sum, aes(x= mean_conANI)) +
  geom_histogram(binwidth= 0.0005, boundary= 1, fill= "black", color= "gray75") +
  labs(x= "Mean consensus ANI", y= "Count") +
  scale_x_continuous(breaks= seq(0.987, 1, by= 0.001),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  theme_strains
ggsave(last_plot(), filename = "conANI_histogram.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

popANI_hist <- ggplot(ani_sum, aes(x= mean_popANI)) +
  geom_histogram(binwidth= 0.0005, boundary= 1, fill= "black", color= "gray75") +
  labs(x= "Mean population ANI", y= "Count") +
  scale_x_continuous(limits= c(0.987, 1),
                     breaks= seq(0.987, 1, by= 0.001),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  theme_strains
ggsave(last_plot(), filename = "popANI_histogram.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

ANI_hist_combined <- plot_grid(conANI_hist + theme(axis.title.x = element_blank()), 
                               popANI_hist + labs(x= "Average nucleotide identity"),
                               nrow= 2,
                               labels= c("A", "B")) +
  draw_label(label= "Consensus\nANI", x= 0.12, y= 0.91, hjust= 0) +
  draw_label(label= "Population\nANI", x= 0.12, y= 0.41, hjust= 0)
ANI_hist_combined

ggsave(last_plot(), filename = "ANI_histogram_combined.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
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
  geom_point(aes(x= mean_conANI, y= mean_popANI, fill= 1- frac_popSNPs), size= 3, pch= 21, color= "black") +
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
  scale_fill_viridis_c(name= "% shared minor alleles",
                       direction = 1,
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
  scale_size_continuous(name= "Difference in\nnucleotide diversity", guide= FALSE) +
  coord_equal() +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1),
        legend.position = "top")
        #legend.key.width= unit(10, "mm"),
        #legend.box= "vertical") 

ggsave(last_plot(), filename = "popANI_conANI_NucDiv.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

ggplot(data= ani_sum) +
  geom_abline(intercept= 0, slope= 1, size= 0.25, linetype= "dashed") +
  geom_point(aes(x= mean_conANI, y= mean_popANI, fill= pop_age_pair), pch= 21, size= 3, color= "gray30") +
  labs(x= "Mean consensus ANI", y= "Mean population ANI") +
  scale_x_continuous(limits= c(0.9875, 1),
                     breaks= seq(.9880, 1, by= 0.001), 
                     labels= c("0.988", "", "0.990", "", "0.992", "", "0.994", "", "0.996", "", "0.998", "", "1"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.9925, 1.0002),
                     breaks= seq(.993, 1, by= 0.001), 
                     labels= c("", "0.994", "", "0.996", "", "0.998", "", "1"),
                     expand= c(0.00, 0)) +
  scale_fill_viridis_d(name= "Pop. age") +
  coord_equal() +
  theme_strains +
  theme(axis.text.x = element_text(angle= 45, hjust= 1),
        legend.position = "top")
#legend.key.width= unit(10, "mm"),
#legend.box= "vertical") 


ggplot(data= ani_sum) +
  geom_point(aes(x= mean_conANI*100, y= 1 - frac_popSNPs, fill= pi_avg_mean, size= pi_diff_mean), 
             pch= 21,  color= "black") +
  labs(x= "Consensus ANI (%)", y= "Shared minor alleles (%)") +
  scale_x_continuous(limits= c(98.74, 100),
                     breaks= seq(98.80, 100, by= 0.1), 
                     labels= c("98.8", "", "99.0", "", "99.2", "", "99.4", "", "99.6", "", "99.8", "", "100"),
                     expand= c(0, 0)) +
  scale_y_continuous(limits= c(0.12, 1.025),
                     breaks= seq(0, 1, by= 0.1),
                     labels= c("0", "", "20", "", "40", "", "60", "", "80", "", "100"),
                     expand= c(0, 0)) +
  scale_fill_viridis_c(name= "Avg.\nnucleotide\ndiversity",
                       limits= c(0.0004, 0.0034)) +
  scale_size_binned(name= "Difference in\nnucleotide\ndiversity",
                    range= c(1, 7),
                    n.breaks= 3,
                    nice.breaks = TRUE) +
  guides(fill = guide_colourbar(order = 1),
         size_binned = guide_legend(order = 2)) +
  theme_strains +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.key.width= unit(10, "mm"),
        legend.box.just = "bottom",
        legend.box.spacing = unit(0, "mm"),
        legend.text= element_text(size= 8),
        legend.title= element_text(size= 10),
        legend.background = element_rect(color= "transparent", fill= "transparent"))

ggsave(last_plot(), filename = "minor_alleles_NucDiv.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

ggplot(ani_sum, aes(x= pi_diff_mean)) +
  geom_histogram()

ggplot(data= ani_sum) +
  #geom_hline(yintercept= 0) +
  #geom_point(aes(x= mean_conANI, y= 1 - frac_popSNPs, fill= mean_popANI), size= 3, pch= 21, color= "black") +
  geom_point(aes(x= mean_conANI, y= 1 - frac_popSNPs, fill= pop_age_pair), pch= 21, size= 3,  color= "black") +
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
  scale_fill_viridis_d(name= "Pop. Age") +
  theme_strains +
  theme(#axis.text.x = element_text(angle= 45, hjust= 1),
    legend.position = "top")
#legend.key.width= unit(10, "mm")) 




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



fit.popANI <- lm((mean_popANI*100) ~ riv_dist, ani_sum)
popANI_rivDist <- ggplot(data= ani_sum) +
  geom_point(aes(x= riv_dist, y= mean_popANI*100), pch= 21, fill= "gray75", color= "black", size= 3) +
  geom_abline(intercept = fit.popANI$coefficients["(Intercept)"], slope= fit.popANI$coefficients["riv_dist"],
              color= "black", size= 1) +
  labs(x= "River network distance (km)", y= "Population ANI (%)") +
  #geom_smooth(aes(x= riv_dist, y= mean_popANI), color= "black", fill= "dodgerblue", alpha= 0.3) +
  scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
                    breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  scale_x_continuous(limits= c(0, 350000),
                     breaks= seq(0, 350000, by= 25000),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5000)) +
  theme_strains +
  theme(legend.position = "top")
popANI_rivDist


gam.fit <- mgcv::gam(mean_popANI ~ s(riv_dist, bs= "cr"), data= ani_sum)
summary(gam.fit)
plot(gam.fit)
gam.check(gam.fit) # check to make sure k value is not too low
anova(gam.fit)


ggsave(last_plot(), filename = "popANI_rivDist.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



fit.conANI <- lm((mean_conANI*100) ~ riv_dist, ani_sum)
conANI_rivDist <- ggplot(data= ani_sum, aes(x= riv_dist, y= mean_conANI*100)) +
  geom_hline(yintercept = 99.35, linetype= "dashed") +
  geom_point(size= 3, pch= 21, fill= "gray75", color= "black") +
  geom_abline(intercept = fit.conANI$coefficients["(Intercept)"], slope= fit.conANI$coefficients["riv_dist"],
              color= "black", size= 1) +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_x_continuous(limits= c(0, 350000),
                     breaks= seq(0, 350000, by= 25000),
                     labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
                     expand= c(0, 5000)) +
  scale_y_continuous(breaks= seq(98.80, 100.00, by= 0.2),
                     labels= c("98.8", "99.0", "99.2", "99.4", "99.6", "99.8", "100.0")) +
  theme_strains +
  theme(legend.position = "top")
conANI_rivDist

ggsave(last_plot(), filename = "conANI_rivDist.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")




conANI_rivdist25km <- ggplot(data= filter(ani_sum, riv_dist < 25000), aes(x= riv_dist/1000, y= mean_conANI*100)) +
  geom_hline(yintercept= 99.35, size= 0.5, linetype= "dashed") +
  #geom_point(aes(fill= watershed_diff), size= 4, pch= 21, color= "black") +
  geom_point(aes(fill= watershed_diff, size= watershed_diff), shape= 21, color= "gray50") +
  labs(x= "River network distance (km)", y= "Consensus ANI (%)") +
  scale_fill_viridis_c(name= expression("Watershed difference (km"^2*")"),
                       option= "plasma") +
  scale_size_continuous(range= c(2.5, 6), guide= FALSE) +
  scale_x_continuous(limits= c(0, 25.01),
                     breaks= seq(0, 25, by= 5),
                     expand= c(0.02, 0)) +
  scale_y_continuous(breaks= seq(98.80, 100.00, by= 0.2),
                     labels= c("98.8", "99.0", "99.2", "99.4", "99.6", "99.8", "100.0")) +
  theme_strains +
  theme(legend.position = c(0.7, 0.85),
        legend.direction = "horizontal",
        legend.background = element_rect(color= "transparent", fill= "transparent"))
conANI_rivdist25km

ggsave(conANI_rivdist25km, filename = "conANI_rivDist_25km.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

conANI_fit1 <- lm(mean_conANI ~ riv_dist + watershed_diff, data= filter(ani_sum, riv_dist < 25000))
summary(conANI_fit1)
plot(conANI_fit1)
anova(conANI_fit1)




ANI_rivdist_combined <- cowplot::plot_grid(popANI_rivDist+ theme(axis.title.x = element_blank()), 
                                           conANI_rivDist + theme(axis.title.x = element_blank()), 
                                           conANI_rivdist25km,
                                           nrow= 3,
                                           labels= c("A", "B", "C"))

ANI_rivdist_combined


ggsave(ANI_rivdist_combined, filename = "ANI_rivDist_combined.png", height= 180, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")
?plot_grid

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

ggplot(data= ani_sum, aes(x= watershed_diff, y= mean_conANI)) +
  #geom_point(aes(x= riv_dist, y= mean_popANI, fill= ani99_dRep), size= 3, pch= 21, color= "black") +
  geom_point(size= 3, pch= 21, color= "black") +
  #geom_smooth(color= "black", fill= "dodgerblue") +
 # geom_smooth(color= "black", fill= "dodgerblue", alpha= 0.3) +
#  labs(x= "River network distance (km)", y= "Mean consensus ANI") +
  #scale_fill_npg(name= "dRep ANI") +
  #scale_fill_manual(values= wes_palette("Royal1", n= 3), name= "dRep ANI",
  #                  breaks= c(">0.995", "0.99-0.995", "<0.99")) +
  # scale_x_continuous(limits= c(0, 350000),
  #                    breaks= seq(0, 350000, by= 25000),
  #                    labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
  #                    expand= c(0, 5000)) +
  theme_strains +
  theme(legend.position = "top")
ggsave(last_plot(), filename = "conANI_rivDist.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")




ggplot(data= ani_sum, aes(x= riv_dist, y= watershed_diff)) +
  # geom_point(size= 3, pch= 21, fill= "gray75", color= "black") +
  geom_point(aes(fill= mean_conANI, size= mean_conANI), pch= 21, color= "black") +
#  geom_smooth(color= "black", fill= "dodgerblue", alpha= 0.3) +
  #labs(x= "River network distance (km)", y= "Mean consensus ANI") +
  scale_fill_viridis_c(name= "Mean conANI") +

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


ggplot(data= ani_sum) +
  geom_point(aes(x= watershed_1, y= mean_popANI), pch= 21, fill= "gray75", color= "black", size= 3) +
  labs(x= "River network distance (km)", y= "Mean population ANI") +
  #geom_smooth(aes(x= riv_dist, y= mean_popANI), color= "black", fill= "dodgerblue", alpha= 0.3) +
  #geom_smooth(aes(x= riv_dist, y= mean_popANI), method= 'lm', color= "black", fill= "dodgerblue", alpha= 0.3) +
  #scale_fill_npg(name= "dRep ANI") +
  # scale_x_continuous(limits= c(0, 350000),
  #                    breaks= seq(0, 350000, by= 25000),
  #                    labels= c("0", "", "50", "", "100", "", "150", "", "200", "", "250", "", "300", "", "350"),
  #                    expand= c(0, 5000)) +
  theme_strains +
  theme(legend.position = "top")


  table(ani_sum$name2)
  
  
  popANI.wide <- pivot_wider(select(ani_sum, name1, name2, mean_popANI), 
                              names_from = name2, 
                              values_from = mean_popANI) %>% 
    ungroup() %>% 
    mutate(PH2015_01D= NA) %>% 
    select(name1, PH2015_01D, everything()) %>% 
    add_row(., name1= "PH2017_40_RAT_O_B")


  gdata::lowerTriangle(popANI.wide[, -1]) <- gdata::upperTriangle(popANI.wide[, -1], byrow=TRUE)
  
popANI.long <- pivot_longer(popANI.wide, names_to= "name2", values_to= "mean_popANI", -name1) %>% 
  # WATERSHED AREA JOIN
  left_join(., watershed.area, by= c("name1" = "ggkbase_id")) %>% 
  left_join(., watershed.area, by= c("name2" = "ggkbase_id")) %>% 
  rename(watershed_1= watershed_km2.x, watershed_2= watershed_km2.y)

popANI.sum <- popANI.long %>% 
  group_by(name1, watershed_1) %>% 
  summarise(meanPOPani= mean(mean_popANI, na.rm=TRUE)) %>% 
  ungroup()
popANI.sum                  


ggplot(data= popANI.long, aes(x= reorder(name1, watershed_1), y= mean_popANI)) +
  geom_point(color="white") +
  geom_boxplot(aes(fill= log10(watershed_1))) +
  labs(x= "Sample", y= "Mean population ANI") +
  scale_fill_viridis_c(breaks= c(1, 2, 3), 
                       labels= c(10, 100, 1000),
                       name= expression("Watershed area (km"^2*")")) +
  theme_strains +
  #scale_x_continuous(labels= as.character(sort(unique(popANI.long$watershed_1)))) +
  theme(#axis.text.x = element_text(angle= 90, vjust= 0.5),
        axis.text.x= element_blank(),
        legend.position = "top")
ggsave(last_plot(), filename = "popANI_watershed.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


ggplot(data= popANI.long, aes(x= log10(watershed_1), y= mean_popANI)) +
  geom_point() +
  geom_smooth(method= "lm")
  #  scale_x_log10() +
  # annotation_logticks() +
  scale_fill_discrete(guide=FALSE) +
  theme_strains

  ggplot(data= popANI.sum, aes(x= watershed_1, y= meanPOPani)) +
    geom_point(size= 2) +
    labs(x= expression("Watershed area (km"^2*")"), y= "Mean population ANI") +
    scale_x_log10() +
    annotation_logticks() +
    theme_strains
  
  
  

fit.1 <- lm(mean_popANI ~ log10(watershed_1), popANI.long)
plot(mean_popANI ~ log10(watershed_1), popANI.long)
summary(fit.1)
anova(fit.1)
plot(fit.1)

hist(sqrt(popANI.long$mean_popANI))

fit.2 <- lm(meanPOPani ~ watershed_1, popANI.sum)
summary(fit.2)
anova(fit.2)
plot(fit.2)




