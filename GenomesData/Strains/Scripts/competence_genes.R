
library(tidyverse)
library(ggplot2)


comp_dir <- "../competence"

comEC_hmm <- read_delim(file= file.path(comp_dir, "competence_hmm_results_table_20190517.txt"), delim= " ", skip= 3, col_names = FALSE) %>% 
  slice(-c(64:73)) %>% 
  mutate(gene= "comEC") %>% 
  select(X1, X3, gene, X4, X6, X9) %>% 
  rename(target= X1, query= X3, accession= X4, score_full= X6, score_domain= X9) %>% 
  mutate_all(funs(str_trim)) %>% 
  mutate(score_full= as.numeric(score_full),
         score_domain= as.numeric(score_domain),
         sample= str_replace(target, "_scaffold.*$", ""),
         scaffold= str_replace(target, "_[0-9]+$", ""))


dprA_hmm <- read_delim(file= file.path(comp_dir, "DNA_processg_A_hmm_results_table_20190517.txt"), delim= " ", skip= 3, col_names = FALSE) %>% 
  slice(-c(126:135)) %>% 
  mutate(gene= "dprA") %>% 
  select(X1, X3, gene, X4, X6, X9) %>% 
  rename(target= X1, query= X3, accession= X4, score_full= X6, score_domain= X9) %>% 
  mutate_all(funs(str_trim)) %>% 
  mutate(score_full= as.numeric(score_full),
         score_domain= as.numeric(score_domain),
         sample= str_replace(target, "_scaffold.*$", ""),
         scaffold= str_replace(target, "_[0-9]+$", ""))

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



