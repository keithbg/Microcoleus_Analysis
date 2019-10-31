library(tidyverse)
library(ggplot2)


#### INPUT FILES
sp1_samples <- read_tsv("sample_sp1_list.tsv")
in_dir <- "inStrain/inStrain_gene_profile_output"
gi_files <- list.files(in_dir, pattern= "species_1.pid96_gene_info")

#snp_files <- list.files(in_dir, pattern= "species_1.pid96_SNP")
# snp_list <- map(snp_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))))
# names(snp_list) <- str_replace(snp_files, ".tsv", "")

#gi.df <- read_tsv(file.path(in_dir, gi_files[3]))
#snp.df <- read_tsv(file.path(in_dir, snp_files[3]))


make_gi_df <- function(dir_in, gi_file_names, samp_file){
  require(tidyverse)
  
  # Read in gene info files
  gi_list <- map(gi_file_names, function(x) suppressMessages(read_tsv(file.path(dir_in, x))) %>% 
                   mutate(pi= 1- clonality,
                          sample= str_replace(x, "_gene_info.tsv", "")))
  names(gi_list) <- str_replace(gi_file_names, ".tsv", "")
  
  # Match samples in list with samples with the selected ANI species
  sp_names <- do.call(rbind, map(names(gi_list), function(x) any(str_detect(x, samp_file$sample))))
  gi_list_sp_subset <- gi_list[sp_names]
  
  # Transform to a data frame
  gi_df <- as_tibble(do.call(rbind, gi_list_sp_subset))
  return(gi_df)
}


gi_sp1_df <- make_gi_df(dir_in = in_dir,
                      gi_file_names = gi_files,
                      samp_file = sp1_samples)

#count(gi_sp1_df, coverage < 2)

cov_quantiles <- gi_sp1_df %>% 
  filter(coverage > 5) %>% # coverage threshold specified in inStrain
  group_by(sample) %>% 
  summarize(
    quant_0.1= quantile(coverage, probs= 0.1, na.rm= TRUE),
    quant_0.5= quantile(coverage, probs= 0.5, na.rm= TRUE),
    quant_0.9= quantile(coverage, probs= 0.9, na.rm= TRUE))


## FILTER OUT LOW COVERAGE GENES ##
# Deciding to filter out the lower 10% of the coverage distributions
# for samples there is a long low coverage tail of the distribution. 
# These low distribution samples may be a result of mis-mapping or poor binning
# These values also contain all the outlier pi values >0.1, which makes these high pi values suspect
# It will be more conservative to remove them from analyses
# The high end of the distribution is tighter, so will not be filtered

# Breadth >0.9
filter_coverage_breadth <- function(gene_df, quant_df, samp_name, breadth_thresh){
  quant_values <- quant_df %>% filter(sample == samp_name)
  
  gene_df_filtered <- gene_df %>% 
    filter(sample == samp_name) %>% 
    filter(breadth > breadth_thresh) %>% 
    filter(coverage > quant_values$quant_0.1)
  return(gene_df_filtered)
  
}

gi_sp1_filt <- map(cov_quantiles$sample, function(x) filter_coverage_breadth(gene_df= gi_sp1_df, 
                                                                             quant_df = cov_quantiles,
                                                                             breadth_thresh= 0.9,
                                                                             samp_name = x))
gi_sp1_filt_df <- do.call(rbind, gi_sp1_filt)




ggplot(data= gi_sp1_filt_df) +
  geom_histogram(aes(x= coverage)) +
  facet_wrap(~sample, nrow= 8, scales= "free_x")

ggplot(data= gi_sp1_filt_df) +
  geom_histogram(aes(x= breadth)) +
  facet_wrap(~sample, nrow= 8, scales= "free_x")

ggplot(data= gi_sp1_filt_df) +
  geom_point(aes(x= coverage, y= pi)) +
  theme_bw()

ggplot(data= gi_sp1_filt_df) +
  geom_point(aes(x= gene, y= pi)) +
  theme(axis.text.x= element_blank()) +
  facet_wrap(~sample, nrow= 8)





