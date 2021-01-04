## Format and join data from inStrain and spatial data

library(tidyverse)


#### INPUT DATA ####

## River network distance data OK
# River network distance calculate on ArcGIS and CSV exported

euclidean.distance <- read_tsv(file.path("Data/Spatial_data", "Distance_Euclidean_meters.tsv")) %>% 
  rename(sample.querry = sample.1, sample.reference = sample.2)

river.distance <- read_tsv(file.path("Data/Spatial_data", "Distance_RiverNetwork_meters.tsv")) %>% 
  gather(key= sample.querry, value= riv_dist, -Site) %>% 
  rename(sample.reference= Site)

ani.riv.dist <- left_join(river.distance, euclidean.distance) %>% 
  rename(name1= sample.reference, name2= sample.querry)
rm(euclidean.distance, river.distance)

## Watershed area data OK
watershed.area <- read_tsv(file.path("Data/Spatial_data", "WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2)

## Species 1 sites OK
sp1_sites <- read_tsv(file.path("Data/inStrain_data/", "inStrain_sample_species_lookup.tsv")) %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present) %>% 
  filter(., species == "species_1") %>% 
  select(-multiple_species)

## Nucleotide diversity (generated in inStrain_format_output.R) OK
nuc_div <- read_tsv(file.path("Data/inStrain_data/", "nuc_div_summary_TEST.txt")) %>% 
  filter(species == "species_1") %>% 
  select(site, mean_pi, median_pi)

## SNV Genome summary data (generated in inStrain_format_output.R) OK
snv_genomes <- read_tsv("Data/inStrain_data/snvs_genome_summary_TEST.tsv")

## inStrain compare results OK
comp_sp1 <- read_tsv("Data/inStrain_data/compare_module_output_sp1_comparisonsTable.tsv") %>% 
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


## Remove inStrain comparisons with percent_genome_compared < 0.25 
# this is 3.3% of all comparisons and includes many low ANI outliers
# 183 scaffolds in the reference genomes: length(unique(comp_sp1$scaffold))

# Had origiginally filtered by removing <0.75, when reducing to 0.25 it had negligible affect
# on the genome wide averages, so I lowered cutoff to 0.25. 

comp_sp1.F <- comp_sp1 %>% 
  filter(percent_genome_compared > 0.25) %>% 
  # Calculate pMA
  mutate(pMA= 1 - (population_SNPs / consensus_SNPs)) # shared minor alleles: 1 - fraction of consensus_SNVs that are also population_SNVs

#### JOIN DATA ####

ani_sum <- comp_sp1.F %>% 
  group_by(name1, name2) %>% 
  summarize(scaf_num= length(scaffold),
            mean_popANI= mean(popANI),
            mean_conANI= mean(conANI),
            sum_pop_sites= sum(population_SNPs),
            sum_con_sites= sum(consensus_SNPs),
            #frac_popSNVs= sum_pop_sites/sum_con_sites,
            pMA= 1 - (sum_pop_sites/sum_con_sites),
            #mean_frac_popSNVs= mean(frac_popSNVs, na.rm=TRUE),
            sd_popANI= sd(popANI),
            sd_conANI= sd(conANI),
            median_popANI= median(popANI),
            median_conANI= median(conANI)) %>% 
  ungroup() %>% 
  # WATERSHED AREA JOIN
  left_join(., watershed.area, by= c("name1" = "ggkbase_id")) %>% 
  left_join(., watershed.area, by= c("name2" = "ggkbase_id")) %>% 
  rename(watershed.1= watershed_km2.x, watershed.2= watershed_km2.y) %>% 
  mutate(watershed_diff= abs(watershed.2 - watershed.1)) %>% 
  # SNV_MBP JOIN
  left_join(., select(snv_genomes, ggkbase_id, SNV_mbp), by= c("name1" = "ggkbase_id")) %>% 
  left_join(., select(snv_genomes, ggkbase_id, SNV_mbp), by= c("name2" = "ggkbase_id")) %>% 
  rename(SNV_mbp.1= SNV_mbp.x, SNV_mbp.2= SNV_mbp.y) %>% 
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

#write_tsv(ani_sum, "Data/inStrain_data/ani_summary_TEST.tsv")


## dRep results REMOVEE FROM ANALYSES
# dir_input_dRep <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","GenomesData", "dRep","Output_tables")
# sp1.ani <- read_tsv(file.path(dir_input_dRep, "sp1_ani.tsv")) %>% 
#   rename(ani_dRep= ani) %>% 
#   mutate(name1= str_replace(row_name, "_s25.*$|_Oscill.*$", ""),
#          name2= str_replace(col_name, "_s25.*$|_Oscill.*$", "")) %>% 
#   mutate(ani99_dRep= ifelse(ani_dRep > 0.995, ">0.995", 
#                             ifelse(ani_dRep <= 0.995 & ani_dRep >=0.99, "0.99-0.995", "<0.99")))

# DREP ANI JOIN
#left_join(., sp1.ani) %>% 

## Haplotype frequencies REMOVE FROM ANALYSES
# haplos <- read_tsv(file= "inStrain/output_tables/haplotype_freqs.tsv") %>% 
#   filter(species == "species_1")

## SNV data and young/old populations REMOVE FROM ANALYSES
# snv_genomes <- read_tsv("Output_tables/snvs_genome_summary.tsv") %>% 
#   filter(species == "species_1") %>% 
#   select(ggkbase_id, SNV_mbp, pop_age)

#rename(#pop_age.1= pop_age.x, pop_age.2= pop_age.y,
#       SNV_mbp.1= SNV_mbp.x, SNV_mbp.2= SNV_mbp.y) %>% 
#mutate(pop_age_pair= str_c(pop_age.1, pop_age.2, sep= "-")) %>%
#mutate(pop_age_pair= ifelse(pop_age_pair == "Young-Old", "Old-Young", pop_age_pair)) %>% 



