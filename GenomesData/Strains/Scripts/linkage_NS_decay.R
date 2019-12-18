## Investigate linkage decay in genomes

## linkage_haplo_df.tsv created in linkage_NS_format.R script
## linkage_NS.tsv file created in linkage_NS_format.R script

library(tidyverse)
library(ggplot2)

species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

#filter((mutation_type_A == "N" | mutation_type_A == "S") & (mutation_type_B == "N" | mutation_type_B == "S")) %>% 

link <- read_tsv("inStrain/output_tables/linkage_haplo_df.tsv") %>% 
  mutate(site= str_replace(sample, "\\.species.*", "")) %>% 
  left_join(species_lookup, ., by= c("site", "species")) # Filter to only include species recovered from each site
  
link_NS <- link %>% 
  filter((mutation_type_A == "N" | mutation_type_A == "S") & (mutation_type_B == "N" | mutation_type_B == "S"))

link_NS_sp3 <- link_NS %>% 
  filter(species == "species_3") %>% 
  filter(haplotype == "h3" | haplotype == "h4")


link_NS_sp3_sum <- link_NS_sp3 %>% 
  filter(is.na(r2_normalized) == FALSE) %>% 
  group_by(species, sample, distance, link_type) %>% 
  summarize(n= length(total),
            total_sum= sum(total),
            mean_r2= mean(r2, na.rm= TRUE),
            mean_r2_norm= mean(r2_normalized, na.rm= TRUE),
            mean_d_prime= mean(d_prime, na.rm= TRUE)) %>% 
  ungroup() %>% 
  mutate(NN_link= ifelse(link_type == "N-N", "N-N", "Other"))


link_NS_sp3_sum_D <- link_NS_sp3 %>% 
  filter(haplotype == "h4") %>% 
  filter(is.na(d_prime_normalized) == FALSE) %>% 
  group_by(species, sample, distance, link_type) %>% 
  summarize(n= length(total),
            total_sum= sum(total),
            mean_d_prime= mean(d_prime, na.rm= TRUE),
            mean_d_prime_norm= mean(d_prime_normalized, na.rm= TRUE)) %>% 
  ungroup() %>% 
  mutate(NN_link= ifelse(link_type == "N-N", "N-N", "Other"))



# link_sum <- link %>% 
#   filter(haplotype == "h3" | haplotype == "h4") %>% 
#   group_by(sample, species, haplotype, distance, link_type) %>% 
#   summarize(r2_norm_mean= mean(r2_normalized, na.rm= TRUE),
#             Dprime_norm_mean= mean(d_prime_normalized, na.rm= TRUE))




#### MAKE FIGURES #####
source("Scripts/ggplot_themes.R")

link_NS_sp3 %>% 
  group_by(gene_A) %>% 
  count() %>% 
  ggplot() +
  geom_histogram(aes(x= n), binwidth= 10)


n(link_NS_sp3$gene_A)


ggplot(data= filter(link_NS_sp3, haplotype == "h4")) +
  geom_point(aes(x= distance, y= r2_normalized, color= link_type)) +
  geom_smooth(aes(x= distance, y= r2_normalized, color= link_type)) +
  facet_wrap(~sample, ncol=2) +
  theme_strains

ggplot() +
  geom_point(data= link_NS_sp3_sum, aes(x= distance, y= mean_r2_norm, size= total_sum)) +
  facet_grid(sample~NN_link) +
  theme_strains


### D prime
ggplot() +
  geom_point(data= link_NS_sp3_sum_D, aes(x= distance, y= mean_d_prime, size= total_sum)) +
  facet_grid(sample~NN_link) +
  theme_strains

ggplot() +
  geom_point(data= link_NS_sp3_sum_D, aes(x= distance, y= mean_d_prime_norm, size= total_sum)) +
  facet_grid(sample~NN_link) +
  theme_strains



ggplot(data= link_NS_sp3_sum, aes(x= distance, y= mean_r2, color= link_type)) +
  geom_point(aes(size= total_sum)) +
  #geom_smooth() +
  facet_grid(sample~link_type) +
  theme_strains


ggplot(data= link_NS_sp3_sum, aes(x= total_sum, y= mean_r2_norm, color= link_type)) +
  geom_point(aes(size= total_sum)) +
#  facet_grid(sample~link_type) +
  theme_strains



ggplot(data= link_NS_sp3_sum, aes(x= distance, y= n, color= link_type)) +
  geom_point(aes(size= total_sum)) +
  #  facet_grid(sample~link_type) +
  theme_strains










## HAPLOTYPES

# haplo.freq <- link.species %>% 
#   group_by(site, sample, species, multiple_species, haplotype) %>% 
#   summarize(n= length(haplotype)) %>% 
#   mutate(freq= n/sum(n))
# 
# 
# with(haplo.freq, table(filter(., species == "species_1"), freq < 0.01))
# 
# 
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
# yggplot(data= haplo.freq, aes(x= site, y= freq)) +
#   geom_bar(aes(fill= haplotype), stat= "identity") +
#   #scale_x_discrete(labels= unique(str_replace(haplo.freq$sample, "\\.species.*", ""))) +
#   scale_y_continuous(expand= c(0, 0)) +
#   facet_grid(species~.) +
#   theme_strains +
#   theme(axis.text.x= element_text(angle= 90, vjust= 0.5))
# ggsave(last_plot(), filename= "haplotype_freqs.jpg", path= "Output_figures", width= 10, height= 8, units= "in", dpi= 320)
# 
# ggplot(data= haplo.freq, aes(x= site, y= freq)) +
#   geom_bar(aes(fill= haplotype), stat= "identity") +
#   #scale_x_discrete(labels= unique(str_replace(haplo.freq$sample, "\\.species.*", ""))) +
#   scale_y_continuous(expand= c(0, 0)) +
#   coord_flip() +
#   facet_grid(species~., scales= "free_y") +
#   theme_strains +
#   theme(axis.text.x= element_text(angle= 90, vjust= 0.5))
# ggsave(last_plot(), filename= "haplotype_freqs_tall.jpg", path= "Output_figures", width= 8, height= 20, units= "in", dpi= 320)
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
# ggsave(last_plot(), filename= "haplotype_freqs_histogram.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)
# 
# 
# 
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
