## Format the SNP_mutation_types.tsv and gene_info.tsv files from inStrain
## Reads in all the files for each sampe and combines into master data frames
## which are then exported to the Data/inStrain_data folder.


### The SNP_mutation_types.tsv file filters SNVs from the SNVs.tsv file by morphia == 2 and cryptic == FALSE. 
## This is why there are fewer SNVs in the SNP_mutation_types.tsv file compared to the SNVs.tsv file

## The gene_info.tsv and SNP_mutation_types.tsv are generated in the "profile" command in inStrain



## Libraries
library(tidyverse)

#### FILE PATHS ##########################################################################################
in_dir <- "Data/inStrain_data/inStrain_output"

species_lookup <- read_tsv("Data/inStrain_data/inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)

#gg_anno <- read_tsv("inStrain/ggkbase_anno.tsv") # ggkbase annotations




#### GENE_INFO.TSV FILE ##################################################################################
## Investigate nucleotide diversity (pi) results from inStrain
## The data are generated in the "gene_profile" command in inStrain
## The files are SampleName.gene_profile.tsv


#### INPUT FILES
gi_files <- list.files(in_dir, pattern= "pid96_gene_info")


test <- read_tsv(file.path(in_dir, gi_files[5])) %>% 
  mutate(sample= str_replace(gi_files[5], "_gene_info.tsv", "")) %>% 
  mutate(site= str_split(sample, "\\.")[[1]][1],
         species= str_split(sample, "\\.")[[1]][2])


## FORMAT GENE INFO FILES
make_gi_df <- function(dir_in, gi_file_names){
  require(tidyverse)
  
  # Read in gene info files
  gi_list <- map(gi_file_names, function(x) suppressMessages(read_tsv(file.path(dir_in, x))) %>% 
                   mutate(sample= str_replace(x, "_gene_info.tsv", ""),
                          site= str_split(sample, "\\.")[[1]][1],
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

# samp_file = species_lookup
gi_df <- make_gi_df(dir_in = in_dir,
                    gi_file_names = gi_files) %>% 
  left_join(species_lookup, ., by= c("site", "species"))  # Filter to only include species recovered from each site

# Write files
#write_tsv(gi_df, path= "inStrain/output_tables/gene_info_df.tsv")



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
    filter(coverage > quant_values$quant_0.1) #lower 10% removed
  return(gene_df_filtered)
  
}


gi_filt_df <- map(cov_quantiles$sample, function(x) filter_coverage_breadth(gene_df= gi_df, 
                                                                            quant_df = cov_quantiles, #lower 10% removed
                                                                            breadth_thresh= 0.9,
                                                                            samp_name = x)) %>% 
  bind_rows(gi_filt)
## COMBINE WITH GGKBASE ANNOTATIONS
#gi_filt_df <- bind_rows(gi_filt) %>%  # transform list into a data frame
  #left_join(., select(gg_anno, c(gene, uniref_anno, uniprot_anno, kegg_anno)), by= "gene") # combine with ggkbase annotations

# Write files
write_tsv(gi_filt_df, path= "Data/inStrain_data/gene_info_filt_df_TEST.tsv")


#### SUMMARIZE PI VALUES ACROSS THE GENOME ####
gi_filt_summary <- gi_filt_df %>% 
  group_by(sample, site, species) %>% 
  summarize(n= length(nucl_diversity),
            mean_pi= mean(nucl_diversity, na.rm= TRUE),
            median_pi= median(nucl_diversity, na.rm= TRUE),
            sd_pi= sd(nucl_diversity, na.rm= TRUE),
            min_pi= min(nucl_diversity, na.rm= TRUE),
            max_pi= max(nucl_diversity, na.rm= TRUE),
            median_cov= median(coverage, na.rm= TRUE)) %>% 
  ungroup() #%>% 
  #left_join(., watershed.area, by= "site") # COMBINE WITH WATERSHED AREA

write_tsv(gi_filt_summary, "Data/inStrain_data/nuc_div_summary_TEST.txt")




#### SNVs.TSV FILE ##################################################################################

## Input files
snv_files <- list.files(in_dir, pattern= ".pid96_SNVs")
snv_list <- map(snv_files, function(x) suppressMessages(read_tsv(file.path(in_dir, x))) %>% 
                  mutate(sample= str_replace(x, "_SNVs.tsv", "")) %>% 
                  mutate(site= str_split(sample, "\\.")[[1]][1],
                         species= str_split(sample, "\\.")[[1]][2])) %>% 
  setNames(str_replace(snv_files, ".tsv", ""))

# Filter to only include species recovered from each site
### The SNP_mutation_types.tsv file filters SNVs from the SNVs.tsv file by morphia == 2 and cryptic == FALSE. 
## This is why there are fewer SNVs in the SNP_mutation_types.tsv file compared to the SNVs.tsv file

snv_df <- bind_rows(snv_list) %>% left_join(species_lookup, ., by= c("site", "species"))  
snv_df_filt <- snv_df %>% 
  filter(allele_count <= 2 & cryptic == FALSE)

# Write files
#write_tsv(snv_df_filt, path= "Data/inStrain_data/snv_df_filt_TEST.tsv")


# Get genome sizes
# genome_lengths.txt generated from bash command
genome.size <- read_delim("Data/inStrain_data/genome_lengths.txt", delim= " ", col_names= FALSE) %>% 
  slice(c(13:16, 21:22))

genome_size_df <- tibble(species= c("species_1", "species_2", "species_3"),
                         genome_mbp= c(as.numeric(genome.size[2, ]) / 1000000, 
                                       as.numeric(genome.size[4, ]) / 1000000, 
                                       as.numeric(genome.size[6, ]) / 1000000))

#### SNVs/Mbp and N:S ####
snvs_mbp_df <- snv_df_filt %>% 
  left_join(., genome_size_df) %>% 
  group_by(sample, species, genome_mbp) %>% 
  summarize(SNVs= length(mutation_type),
            cov_mean= mean(position_coverage, na.rm= TRUE),
            cov_sd= sd(position_coverage, na.rm= TRUE)) %>%
  ungroup() %>% 
  mutate(SNV_mbp= SNVs / genome_mbp)



## CALCULATE N:S RATIOS
# for both genome and individual genes
NS_genome_ratios <- snv_df_filt %>% 
  group_by(sample, species) %>% 
  summarize(
    SNV_count= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N"), na.rm = TRUE),
    S_count= sum(str_detect(mutation_type, "S"), na.rm = TRUE),
    I_count= sum(str_detect(mutation_type, "I"), na.rm = TRUE),
    M_count= sum(str_detect(mutation_type, "M"), na.rm = TRUE),
    NS= N_count/S_count) 


## JOIN THE SNVS PER MBP AND NS RATIO  DATA
snvs_genome_df <- left_join(snvs_mbp_df, NS_genome_ratios, by= c("sample", "species")) %>% 
  mutate(ggkbase_id= str_replace(sample, "\\.species.*$", ""))# %>% 
  
write_tsv(snvs_genome_df, "Data/inStrain_data/snvs_genome_summary_TEST.tsv")




