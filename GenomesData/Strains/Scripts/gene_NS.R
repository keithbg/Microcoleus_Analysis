library(tidyverse)
library(ggplot2)


#### INPUT FILES
sp1_samples <- read_tsv("sample_sp1_list.tsv")
in_dir <- "inStrain/inStrain_gene_profile_output"


snv_files <- list.files(in_dir, pattern= "species_1.pid96_SNP")
snv_list <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))))
names(snv_list) <- str_replace(snv_files, ".tsv", "")

snv.df <- snv_list[[1]]

## Remove all intergenic SNVs
snv_genes <- map(snv_list, function(x) filter(x, is.na(gene) == FALSE))
snv_genes.df <- snv_genes[[1]]

snv_genes_NS <- map(snv_genes, function(x) x %>% 
                      group_by(gene) %>% 
                      count(mutation_type) %>% 
                      spread(key=mutation_type, value= n))
test <- snv_genes.df %>% 
  group_by(gene) %>% 
  count(mutation_type) %>% 
  spread(key=mutation_type, value= n)
snv_genes_NS[[1]]

# make_gi_df <- function(dir_in, gi_file_names, samp_file){
#   require(tidyverse)
#   
#   # Read in gene info files
#   gi_list <- map(gi_file_names, function(x) suppressMessages(read_tsv(file.path(dir_in, x))) %>% 
#                    mutate(pi= 1- clonality,
#                           sample= str_replace(x, "_gene_info.tsv", "")))
#   names(gi_list) <- str_replace(gi_file_names, ".tsv", "")
#   
#   # Match samples in list with samples with the selected ANI species
#   sp_names <- do.call(rbind, map(names(gi_list), function(x) any(str_detect(x, samp_file$sample))))
#   gi_list_sp_subset <- gi_list[sp_names]
#   
#   # Transform to a data frame
#   gi_df <- as_tibble(do.call(rbind, gi_list_sp_subset))
#   return(gi_df)
# }
# 
# 
# gi_sp1_df <- make_gi_df(dir_in = in_dir,
#                       gi_file_names = gi_files,
#                       samp_file = sp1_samples)
# 
# gi_sp1_df_filt <- gi_sp1_df %>% 
#   filter(coverage > 2)
# 
# 
# 
# ggplot(data= gi_sp1_df_filt) +
#   geom_histogram(aes(x= coverage)) +
#   facet_wrap(~sample, nrow= 10, scales= "free_x")
# 
# ggplot(data= gi_sp1_df_filt) +
#   geom_point(aes(x= gene, y= pi, color= log10(coverage))) +
#   theme(axis.text.x= element_blank()) +
#   facet_wrap(~sample, nrow= 8)
# 
# 
# 
# ggplot(data= gi_sp1_df_filt) +
#   geom_point(aes(x= coverage, y= pi, color= log10(coverage)))
#   #scale_x_continuous(limits= c(2, 10))
# 
# 
# 
# 
# 
# # gi_pi_summary <- map(gi_list, function(x) unclass(summary(x$pi)))
# # 
# # gi_pi_summary_df <- as.data.frame(do.call(rbind, gi_summary))
# # gi_pi_summary_df$id <- row.names(gi_summary_df)
# # gi_pi_summary_df <- as_tibble(gi_summary_df)
# # 
# # ggplot(data= gi_pi_summary_df) +
# #   geom_jitter(aes(x= NA, y= `1st Qu.`))
# # 
# # ggplot(data= gi_pi_summary_df) +
# #   geom_jitter(aes(x= NA, y= Median))
# # ggplot(data= gi_pi_summary_df) +
# #   geom_jitter(aes(x= NA, y= Min.))
# # 
# # 
# # ggplot(data= filter(gi.df, coverage > 1 & breadth > 0.7)) +
# #   geom_point(aes(x= coverage, y= clonality))
# # 
# # 
# # filter(gi.df, coverage > 2 & breadth > 0.7)
# # table(gi.df$coverage > 1)
# # 
# # ggplot() +
# # geom_histogram(aes(x= gi.df$coverage), binwidth= 3)
