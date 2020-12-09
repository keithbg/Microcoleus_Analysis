
setwd("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/Strains")


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

write_tsv(ani_sum, "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/Data/inStrain_data/ani_summary.tsv")
