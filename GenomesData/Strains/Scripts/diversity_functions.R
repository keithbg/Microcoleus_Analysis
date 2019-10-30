## Reference genomes
# ANI Species 1: dRep cluster 1_17
  # PH2015_12U_Oscillatoriales_45_315.fa
# ANI Species 2: dRep cluster 2_2
  # PH2015_13D_Oscillatoriales_45_19.fa
# ANI Species 3: dRep cluster 1_1
# PH2017_22_RUC_O_B_Oscillatoriales_46_93.fa


## Calculate SNPs per mega base pair
## SNPs= number of rows in .freq file from strainRep.py
## Rows calculated with a BASH loop March 22, 2019 on Biotite

calc_snps_mbp <- function(){
  require(tidyverse)
  
  genome.size <- read_delim("genome_lengths.txt", delim= " ", col_names= FALSE) %>% 
    slice(c(13:16, 21:22))
  
  genome.sp1 <- as.numeric(genome.size[2, ]) / 1000000
  genome.sp2 <- as.numeric(genome.size[4, ]) / 1000000
  genome.sp3 <- as.numeric(genome.size[6, ]) / 1000000
  
  snps.sp1 <- read_delim("SNPs_species_1.txt", delim= " ", col_names= FALSE) %>% 
    rename(snps= X1, file= X2) %>% 
    mutate(species= "1",
           snps_mbp= snps / genome.sp1) %>% 
    filter(snps > 200) # filter to remove samples that did not have that species binned in that sample
  
  snps.sp2 <- read_delim("SNPs_species_2.txt", delim= " ", col_names= FALSE) %>% 
    rename(snps= X1, file= X2) %>% 
    mutate(species= "2",
           snps_mbp= snps / genome.sp2) %>% 
    filter(snps > 300)
  
  snps.sp3 <- read_delim("SNPs_species_3.txt", delim= " ", col_names= FALSE) %>% 
    rename(snps= X1, file= X2) %>% 
    mutate(species= "3",
           snps_mbp= snps / genome.sp3) %>% 
    filter(snps > 1300)
  
  snps.df <- full_join(snps.sp1, snps.sp2) %>% 
    full_join(., snps.sp3) %>% 
    mutate(file= str_replace(.$file, ".freq", "")) %>% 
    rename(sample= file)
    
  return(snps.df)
}




## Calculate SNPs per mega base pair
## SNPs= number of rows in .freq file from strainRep.py
## Rows calculated with a BASH loop March 22, 2019 on Biotite

## List _aa.freq files
calc_NS_ratio <- function(file.list, species){
  require(tidyverse)
  
  if(class(file.list) != "list"){
    stop("file.list must be a LIST of filenames")
  }
  
  
  count_Ns_Ss <- function(file){
    df <- suppressMessages(read_tsv(file.path(species, file)))
    Ns <- sum(str_detect(df$mutation, "^N:"))
    Ss <- sum(str_detect(df$mutation, "^S:"))
    data.frame(Ns, Ss)
  }
  
  NS.list <- lapply(file.list, count_Ns_Ss)
  
  NS.df <-  do.call(rbind, NS.list) %>% 
    mutate(NS_ratio= Ns / Ss,
           sample= str_replace(do.call(rbind, file.list), "_aa.freq", ""))
  
  return(NS.df)
}





## Get Microdiversity
get_microdiversity <- function(){
  require(tidyverse)
  
### SOURCE FUNCTIONS
source(file.path("strains_analysis", "R_scripts", "snv_linkage_functions.R"))
# Returns functions: analyse_linkage_data(), filter_snv_freq_window(), make_multi_panel_fig(), summarize_log_files()

## SAMPLE METADATA
samp.md <- read_tsv("sample.metadata.snv.linkage.tsv")

#### Summarize data on .log files 
#log.files <- list.files(file.path(dir_input, "species_1", "log_files"), "*.pid98.log")

log_sp1 <- summarize_log_files(path= file.path("species_1", "log_files"))
log_sp2 <- summarize_log_files(path= file.path("species_2", "log_files"))
log_sp3 <- summarize_log_files(path= file.path("species_3", "log_files"))

# Combine into a master data frame
log_master <- do.call(rbind, list(log_sp1, log_sp2, log_sp3)) %>%
  left_join(., samp.md, by= "sample") %>% 
  mutate(microdiversity= 1-clonality,
         ref_id= str_replace(.$ref_id, "[a-z]+_", "")) %>% 
  rename(species = ref_id)

# Check if species returned from read mapping matches the species binned from the sample
temp_species <- str_extract_all(log_master$species_present, "[0-9]", simplify = TRUE)
log_master <- log_master %>%
  mutate(species_match= ifelse(str_replace_all(log_master$species, "[^0-9]", "") == temp_species[, 1] | str_replace_all(log_master$species, "[^0-9]", "") == temp_species[, 2], "Y", "N")) %>% 
  filter(species_match == "Y" & pid == "98") %>% 
  filter(breadth > 0.1)


return(log_master)
}
