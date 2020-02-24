## Phycobilisome gene analyses
## KEGG HMMs 

library(tidyverse)
library(ggplot2)


pbs_dir <- "../Photosystem_genes"
pbs_gene_names <- read_tsv(file.path(pbs_dir, "pbs_antennae_KEGG_hmms.txt"))

species_genomes <- tibble(genome= c("PH2015_12U_Oscillatoriales_45_315", "PH2015_13D_Oscillatoriales_45_19", "PH2017_22_RUC_O_B_Oscillatoriales_46_93"),
                  sample= c("PH2015_12U", "PH2015_13D", "PH2017_22_RUC_O_B"),
                  species= c("1", "2", "3"))


read_hmm_table <- function(hmm_file){
  require(tidyverse)
  read_delim(file= hmm_file, delim= " ", skip= 3, col_names = FALSE) %>% 
    select(X1, X3, X4, X5, X6, X8, X9) %>% 
    rename(target= X1, query= X3, query_accession= X4, eval_full= X5, score_full= X6, eval_domain= X8, score_domain= X9) %>% 
    mutate_all(funs(str_trim)) %>% 
    filter(target != "#") %>% 
    mutate(score_full= as.numeric(score_full),
           score_domain= as.numeric(score_domain),
           eval_full= as.numeric(eval_full),
           eval_domain= as.numeric(eval_domain),
           sample= str_replace(target, "_scaffold.*$", ""),
           scaffold= str_replace(target, "_[0-9]+$", "")) 
  
  
}


pbs_hmm <- read_hmm_table(file.path(pbs_dir, "pbs_hmm_results_table_20191219.txt")) %>% 
  mutate(kegg_hmm= str_replace(query, ".gCluster.*$", "")) %>% 
  left_join(., pbs_gene_names, by= "kegg_hmm") %>% 
  left_join(species_genomes, ., by= "sample") %>% 
  arrange(gene, species, eval_full)

gene_counts <- pbs_hmm %>% 
  group_by(gene, species, sample) %>% 
  count() %>% 
  arrange(-n)

min_eval <- pbs_hmm %>% 
  group_by(gene, species) %>% 
  filter(eval_full == min(eval_full)) %>% 
  arrange(gene, species)

ggplot(data= pbs_hmm, aes(x= reorder(target, -score_full), y= score_full)) +
  geom_point(aes(color= sample)) +
  facet_wrap(~gene, ncol= 6, scales= "free_x") +
  theme(axis.text.x= element_blank())

ggplot(data= pbs_hmm, aes(x= reorder(target, -score_full), y= score_full)) +
  geom_point(aes(color= sample)) +
  facet_grid(sample~gene, scales= "free_x") +
  theme(axis.text.x= element_blank())







ggplot(data= comEC_hmm, aes(x= reorder(target, -score_full), y= score_full)) +
  geom_point()

ggplot(data= dprA_hmm, aes(x= reorder(target, -score_full), y= score_full)) +
  geom_point()

ggplot(data= dprA_hmm, aes(x= reorder(target, -score_domain), y= score_domain)) +
  geom_point()


comEC_sample_counts <- comEC_hmm %>% 
  filter(score_full > 160) %>% 
  count(sample)

table(comEC_sample_counts$n > 1)

length(unique(comEC_sample_counts$sample))
length(unique(dprA_sample_counts$sample))


dprA_sample_counts <- dprA_hmm %>% 
  filter(score_full > 100) %>% 
  count(sample)
table(dprA_sample_counts$n > 1)



