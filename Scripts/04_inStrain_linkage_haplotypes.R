## Calculate haplotypes from linkage data

## Combine the linkage data with mutation type and haplotypes
## reads in 2 files from inStrain profile module:
##  XX.SNVs.tsv files (https://instrain.readthedocs.io/en/latest/example_output.html#snvs-tsv)
## and XX.linkage.tsv files (https://instrain.readthedocs.io/en/latest/example_output.html#linkage-tsv)

## exports a linkage file with the mutation type and haploytpe of the linkage file included
## Exported file = "Data/inStrain_data/linkage_haplo_df_v1.4.tsv"

library(tidyverse)

#### INPUT FILES ####
link.files <- list.files("Data/inStrain_Data/inStrain_output/", pattern= "linkage")
snv.files <- list.files("Data/inStrain_Data/inStrain_output/", pattern= "SNVs")

link.list <- map(link.files, function (x) read_tsv(file.path("Data", "inStrain_Data", "inStrain_output", x))) %>% 
  setNames(str_replace(link.files, "_linkage.tsv", ""))

snv.list <- map(snv.files, function (x) read_tsv(file.path("Data", "inStrain_Data", "inStrain_output", x))) %>% 
  setNames(str_replace(snv.files, "_SNVs.tsv", ""))

species_lookup <- read_tsv("Data/inStrain_data/inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

## Check to make sure lists are in the same order
table(names(link.list) == names(snv.list))


#### COMBINE MUTATION TYPE AND LINKAGE TOGETHER ####
linkage_type <- map2(snv.list, link.list, function(SNV, LINK){
  snv_red <- SNV %>% 
    select(scaffold, gene, position, mutation_type)
  
  # Join on positions A
  NS_positions_A <- inner_join(snv_red, LINK, by= c("scaffold", "position" = "position_A")) %>% 
    rename(mutation_type_A= mutation_type,
           position_A= position,
           gene_A= gene) %>% 
    select(scaffold, gene_A, position_A, mutation_type_A) %>% 
    distinct()
  
  # Join on position B
  NS_positions_B <- inner_join(snv_red, LINK, by= c("scaffold", "position" = "position_B")) %>% 
    rename(mutation_type_B= mutation_type,
           position_B= position,
           gene_B= gene) %>% 
    select(scaffold, gene_B, position_B, mutation_type_B) %>% 
    distinct()
  
  # Combine mutation types on positions A and B with the original data
  linkage.NS <- left_join(LINK, NS_positions_A, by= c("scaffold", "position_A")) %>% 
    left_join(., NS_positions_B, by= c("scaffold", "position_B")) %>% 
    mutate(link_type= str_c(mutation_type_A, mutation_type_B, sep= "-")) %>% 
    mutate(link_type= ifelse(link_type == "N-S" | link_type == "S-N", "N-S", link_type))
  
  return(linkage.NS)
}
) %>% setNames(str_replace(snv.files, "_SNVs.tsv", ""))

rm(link.list, snv.list)

## Transform into a dataframe
linkage_type_df <- bind_rows(linkage_type, .id= "sample") %>% 
  mutate(species= do.call(rbind, str_split(.$sample, "\\."))[, 2]) %>% 
  select(sample, species, scaffold, everything())



## CALCULATE HAPLOTYPES
calc_haplotypes <- function(df){
  counts <-  as.matrix(df[, c('countAB','countAb','countaB','countab')])
  counts[which(counts != 0)] <-  1
  
  df$haplotype <- ifelse(rowSums(counts) == 4, "h4", 
                         ifelse(rowSums(counts) == 3, "h3",
                                ifelse(rowSums(counts) == 2, "h2",
                                       ifelse(rowSums(counts) == 1, "h1", "error"))))
  return(df)
}

linkage_type_df_haplo <- calc_haplotypes(linkage_type_df)

## Write file
#write_tsv(linkage_type_df_haplo, path= "Data/inStrain_data/linkage_haplo_df_v1.4.tsv")


# Filter to only include species recovered from each site
linkage_type_df_haplo.species <- linkage_type_df_haplo %>% 
  mutate(site= str_replace(sample, "\\.species.*", "")) %>% 
  left_join(species_lookup, ., by= c("site", "species"))  

## Haplotype frequencies
haplo.freq <- linkage_type_df_haplo.species %>% 
  group_by(site, sample, species, multiple_species, haplotype) %>% 
  summarize(n= length(haplotype)) %>% 
  mutate(freq= n/sum(n)) %>% 
  ungroup()

#write_tsv(haplo.freq, path= "Data/inStrain_data/linkage_haplo_freqs_v1.4.tsv")







#### OLD CODE ####
# library(tidyverse)
# library(ggplot2)
# 
# link <- read_tsv("Data/inStrain_data/linkage_haplo_df_v1.4.tsv") %>% 
#   filter((mutation_type_A == "N" | mutation_type_A == "S") & (mutation_type_B == "N" | mutation_type_B == "S")) %>% 
#   mutate(site= str_replace(sample, "\\.species.*", ""))
# 
# link.all <- read_tsv("inStrain/linkage_haplo_df_v1.4.tsv") %>% 
#   mutate(site= str_replace(sample, "\\.species.*", ""))
# 
# 
# species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
#   mutate(species= str_c("species_", species_present)) %>% 
#   rename(site= sample) %>% 
#   select(-species_present)
# 
# link.species <- left_join(species_lookup, link, by= c("site", "species"))  # Filter to only include species recovered from each site
# link.all.species <- left_join(species_lookup, link.all, by= c("site", "species"))  # Filter to only include species recovered from each site
# 
# link_sum <- link %>% 
#   filter(haplotype == "h3" | haplotype == "h4") %>% 
#   group_by(sample, species, haplotype, distance, link_type) %>% 
#   summarize(r2_norm_mean= mean(r2_normalized, na.rm= TRUE),
#             Dprime_norm_mean= mean(d_prime_normalized, na.rm= TRUE))
# 
# #test <- link.species %>% 
# #  filter(str_detect(scaffold, "PH2015_13D"))
# 
# ## HAPLOTYPES
# link.species %>% 
#   group_by(sample, species) %>% 
#   summarize(n= length(haplotype)) %>% 
#   ggplot() +
#   geom_boxplot(aes(x= species, y= n)) +
#   scale_y_log10() +
#   theme_bw()
# 
# 
# haplo.freq <- link.all.species %>% 
#   group_by(site, sample, species, multiple_species, haplotype) %>% 
#   summarize(n= length(haplotype)) %>% 
#   mutate(freq= n/sum(n)) %>% 
#   ungroup()
# 
# #write_tsv(haplo.freq, path= "inStrain/output_tables/haplotype_freqs.tsv")
# 
# 
# #with(haplo.freq, table(filter(., species == "species_1"), freq < 0.01))
# 
# 
# haplo.freq.sp1 <- haplo.freq %>% 
#   filter(species == "species_1")
# 
# 
# # Summarize haplotype frequencies
# haplo.freq %>% 
#   group_by(species, haplotype) %>% 
#   summarize(count= length(freq),
#             min_freq= min(freq),
#             mean_freq= mean(freq),
#             med_freq= median(freq),
#             max_freq= max(freq))
#           
# 
# haplo.freq %>% 
#   group_by(haplotype) %>% 
#   summarize(count= length(freq),
#             min_freq= min(freq),
#             mean_freq= mean(freq),
#             med_freq= median(freq),
#             max_freq= max(freq))
# 
# ## Are H4 frequencies higher at sites where multiple species were recovered?
# haplo.freq %>% 
#   filter(species != "species_1") %>% 
#   group_by(haplotype, multiple_species, species) %>% 
#   summarize(count= length(freq),
#             min_freq= min(freq),
#             mean_freq= mean(freq),
#             med_freq= median(freq),
#             max_freq= max(freq))
# 
# 
# #### MAKE FIGURES ####
# source("Scripts/ggplot_themes.R")
# 
# 
# 
# ggplot(data= link_sum, aes(x= distance, y= r2_norm_mean)) +
#   geom_point() +
#   facet_grid(link_type~haplotype) +
#   theme_bw()
# 
# ggplot(data= filter(link_sum, haplotype == "h4"), aes(x= distance, y= Dprime_norm_mean)) +
#   geom_point(aes(color= link_type)) +
#   facet_wrap(~sample) +
#   theme_bw()
# 
# 
# ## HAPLOTYPE PLOTS ##
# 
# haplo.freq %>% 
#   group_by(sample, species, haplotype) %>% 
#   count(n>0)
# 
# 
# ggplot(data= haplo.freq, aes(x= site, y= freq)) +
#   geom_bar(aes(fill= haplotype), stat= "identity") +
#   #scale_x_discrete(labels= unique(str_replace(haplo.freq$sample, "\\.species.*", ""))) +
#   scale_y_continuous(expand= c(0, 0)) +
#   facet_grid(species~.) +
#   theme_strains +
#   theme(axis.text.x= element_text(angle= 90, vjust= 0.5))
# ggsave(last_plot(), filename= "haplotype_freqs.png", path= "Output_figures", height= 180*0.75, width= 180, units= "mm", dpi= 320,)
# 
# 
# ggplot(data= haplo.freq, aes(x= site, y= freq)) +
#   geom_bar(aes(fill= haplotype), stat= "identity") +
#   #scale_x_discrete(labels= unique(str_replace(haplo.freq$sample, "\\.species.*", ""))) +
#   scale_y_continuous(expand= c(0, 0)) +
#   facet_grid(species~.) +
#   theme_strains +
#   theme(axis.text.x= element_text(angle= 90, vjust= 0.5))
# 
# 
# ggplot(data= haplo.freq, aes(x= site, y= freq)) +
#   geom_bar(aes(fill= haplotype), stat= "identity") +
#   #scale_x_discrete(labels= unique(str_replace(haplo.freq$sample, "\\.species.*", ""))) +
#   scale_y_continuous(expand= c(0, 0)) +
#   coord_flip() +
#   facet_grid(species~., scales= "free_y") +
#   theme_strains +
#   theme(axis.text.x= element_text(angle= 90, vjust= 0.5))
# ggsave(last_plot(), filename= "haplotype_freqs_tall.png", path= "Output_figures", width= 8, height= 20, units= "in", dpi= 320)
# 
# 
# ggplot(data= haplo.freq, aes(x= freq)) +
#   geom_histogram(binwidth= 0.05, 
#                  boundary= 0, 
#                  fill= "black", color= "white") +
#   scale_y_continuous(expand= c(0, 0)) +
#   # scale_x_continuous(breaks= seq(0, 1, by= 0.05),
#   #                    labels= c("0", "", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5",
#   #                              "", "0.6", "", "0.7", "", "0.8", "", "0.9", "", "1.0")) +
#   scale_x_continuous(breaks= seq(0, 1, by= 0.1),
#                      labels= c("0", "0.1",  "0.2",  "0.3", "0.4",  "0.5",
#                                "0.6", "0.7",  "0.8",  "0.9",  "1.0"),
#                      expand= c(0.01, 0)) +
#   facet_grid(species~haplotype, scales= "free_y") +
#   theme_strains
# ggsave(last_plot(), filename= "haplotype_freqs_histogram.png", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)
# 
# ggplot(data= haplo.freq, aes(x= haplotype, y= freq)) +
#   geom_boxplot(aes(fill= haplotype)) +
#   labs(x= "Haplotype", y= "Frequency") +
#   scale_x_discrete(labels= c("H1", "H2", "H3", "H4")) +
#   scale_y_continuous(limits= c(0, 1), 
#                      breaks= seq(0, 1, by= 0.1),
#                      labels= c("0", "0.1",  "0.2",  "0.3", "0.4",  "0.5",
#                                "0.6", "0.7",  "0.8",  "0.9",  "1.0"),
#                      expand= c(0.01, 0)) +
#   scale_fill_manual(values= wes_palette("FantasticFox1")[2:5],
#                     guide= FALSE) +
#   facet_grid(.~species, scales= "free_y", labeller= labeller(species = as_labeller(c(`species_1` = "Species 1",
#                                                                                      `species_2` = "Species 2",
#                                                                                      `species_3` = "Species 3")))) +
#   theme_strains
# ggsave(last_plot(), filename= "haplotype_freqs_boxplot.png", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)
# 
# 
# ggplot(data= filter(link_sum, haplotype == "h4"), aes(x= distance, y= Dprime_norm_mean)) +
#   geom_point(aes(color= link_type)) +
#   facet_wrap(~sample) +
#   theme_bw()
# 
# 
# 
# table(link_sum$haplotype)
