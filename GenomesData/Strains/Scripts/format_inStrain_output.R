## Format the SNP_mutation_types.tsv and gene_info.tsv files from inStrain
## Reads in all the files for each sampe and combines into master data frames
## which are then exported to the inStrain/output_tables folder.


### The SNP_mutation_types.tsv file filters SNVs from the SNVs.tsv file by morphia == 2 and cryptic == FALSE. 
## This is why there are fewer SNVs in the SNP_mutation_types.tsv file compared to the SNVs.tsv file

## The gene_info.tsv and SNP_mutation_types.tsv are generated in the "gene_profile" command in inStrain



## Libraries
library(tidyverse)

#### FILE PATHS ##########################################################################################
in_dir <- "inStrain/inStrain_gene_profile_output"

species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

gg_anno <- read_tsv("inStrain/ggkbase_anno.tsv") # ggkbase annotations



#### SNP_MUTATION_TYPE.TSV FILE ##################################################################################

#### INPUT FILES ####
snv_files <- list.files(in_dir, pattern= ".pid96_SNP")
snv_list <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                  mutate(sample= str_replace(x, "_SNP_mutation_types.tsv", "")) %>% 
                  mutate(site= str_split(sample, "\\.")[[1]][1],
                         species= str_split(sample, "\\.")[[1]][2])) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))

# Filter to only include species recovered from each site
snv_df <- bind_rows(snv_list) %>% left_join(species_lookup, ., by= c("site", "species"))  

# Write files
write_tsv(snv_df, path= "inStrain/output_tables/snp_mutation_type_df.tsv")


#### GENE_INFO.TSV FILE ##################################################################################


## Investigate nucleotide diversity (pi) results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.gene_profile.tsv


#### INPUT FILES
gi_files <- list.files(in_dir, pattern= "pid96_gene_info")


## FORMAT GENE INFO FILES
make_gi_df <- function(dir_in, gi_file_names, samp_file){
  require(tidyverse)
  
  # Read in gene info files
  gi_list <- map(gi_file_names, function(x) suppressMessages(read_tsv(file.path(dir_in, x))) %>% 
                   mutate(pi= 1- clonality,
                          sample= str_replace(x, "_gene_info.tsv", "")) %>% 
                   mutate(site= str_split(sample, "\\.")[[1]][1],
                          species= str_split(sample, "\\.")[[1]][2])) %>% 
    setNames(str_replace(gi_file_names, ".tsv", ""))
  
  # Match samples in list with samples with the selected ANI species
  #sp_names <- do.call(rbind, map(names(gi_list), function(x) any(str_detect(x, samp_file$sample))))
  #gi_list_sp_subset <- gi_list[sp_names]
  
  # Transform to a data frame
  #gi_df <- as_tibble(do.call(rbind, gi_list_sp_subset))
  gi_df <- bind_rows(gi_list)
  
  return(gi_df)
}


gi_df <- make_gi_df(dir_in = in_dir,
                    gi_file_names = gi_files,
                    samp_file = species_lookup) %>% 
  left_join(species_lookup, ., by= c("site", "species"))  # Filter to only include species recovered from each site

# Write files
write_tsv(gi_df, path= "inStrain/output_tables/gene_info_df.tsv")



## FILTER OUT LOW COVERAGE GENES ##
# Deciding to filter out the lower 10% of the coverage distributions
# for samples there is a long low coverage tail of the distribution. 
# These low distribution samples may be a result of mis-mapping or poor binning
# These values also contain all the outlier pi values >0.1, which makes these high pi values suspect
# It will be more conservative to remove them from analyses
# The high end of the distribution is tighter, so will not be filtered

# Remove breadth <=0.9

## Calculate coverage quantiles
cov_quantiles <- gi_df %>% 
  filter(coverage > 5) %>% # coverage threshold specified in inStrain
  group_by(sample) %>% 
  summarize(
    quant_0.05= quantile(coverage, probs= 0.05, na.rm= TRUE),
    quant_0.1= quantile(coverage, probs= 0.1, na.rm= TRUE),
    quant_0.5= quantile(coverage, probs= 0.5, na.rm= TRUE),
    quant_0.9= quantile(coverage, probs= 0.9, na.rm= TRUE),
    quant_0.95= quantile(coverage, probs= 0.95, na.rm= TRUE))


## Filter out low coverage and low breadth genes
filter_coverage_breadth <- function(gene_df, quant_df, samp_name, breadth_thresh){
  quant_values <- quant_df %>% filter(sample == samp_name)
  
  gene_df_filtered <- gene_df %>% 
    filter(sample == samp_name) %>% 
    filter(breadth > breadth_thresh) %>% 
    filter(coverage > quant_values$quant_0.1)
  return(gene_df_filtered)
  
}


gi_filt <- map(cov_quantiles$sample, function(x) filter_coverage_breadth(gene_df= gi_df, 
                                                                         quant_df = cov_quantiles,
                                                                         breadth_thresh= 0.9,
                                                                         samp_name = x))
## COMBINE WITH GGKBASE ANNOTATIONS
gi_filt_df <- bind_rows(gi_filt) %>%  # transform list into a data frame
  left_join(., select(gg_anno, c(gene, uniref_anno, uniprot_anno, kegg_anno)), by= "gene") # combine with ggkbase annotations

# Write files
write_tsv(gi_filt_df, path= "inStrain/output_tables/gene_info_filt_df.tsv")


