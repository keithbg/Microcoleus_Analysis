## # Generalized Dissimilarity Matrices of Microcoleus species 1

# Generalized Dissimilarity Matrices
# Fitzpatrick and Keller 2015, Ecology Letters

# R package, gdm: # https://cran.r-project.org/web/packages/gdm/index.html

## Average nucleotide identity (ANI) calculated with dRep

## ani.species.tsv generated in dRep_format.R


#### Libraries #################################################################
library(tidyverse)
# library(gdm) DON"T LOAD THIS PACKAGE AS IT CONFLICTS WITH TIDYVERSE
# https://cran.r-project.org/web/packages/gdm/index.html
################################################################################


#### FILE PATHS ################################################################
dir_input_scripts <- file.path("Scripts")
dir_input_ani <- file.path("Data", "dRep")
dir_input_species <- file.path("Data")
dir_input_env <- file.path("Data", "Env_data")
dir_input_spatial <- file.path("Data", "Spatial_data")
dir_output <- file.path("GenomesData", "GDM")
################################################################################

#### READ AND FORMAT DATA ######################################################

#### PH2015 ENVIRONMENTAL DATA ####

env.2015.df <- read_tsv(file.path(dir_input_env, "PH2015_env_data_V2.tsv")) %>% 
  select(-Alk) %>% 
  rename(site= biotite_ID)


env.2015.df <- env.2015.df %>%
  rename(temp= temp_c, do_mgL= DO_mgL, cond_ms= cond, canopy_cover_percent= canopy_avg, depth_cm= depth_cont, DOC_ugL= NPOC_mgL) %>%
  select(-watershed_km2) %>%
  mutate(DOC_ugL= DOC_ugL*1000, # convert from mg/L to ug/L
         flow_cont= flow_cont/100, # convert from cm/s to m/s
         NP= (TDN_ugL/14.0067)/(TDP_ugL/30.9737),
         canopy_cover_percent= canopy_cover_percent*100,
         cond_ms= cond_ms/1000)

#### PH2017 ENVIRONMENTAL DATA ####
##(variables have already been normalized by mean and variance)

env.2017.df <- read_tsv("Data/Env_data/PH2017_env_data_formated.tsv") %>%
  filter(uniqueID != "bear_us") %>%
  select(-site, -date, -uniqueID, -fork, -us_ds, -time, -mbars, -do_perL, -dist_cob,  -alk, -alk_conv, -flow_depth_dft, -biotite_acronym) %>%
  rename(site= biotite_ID, flow_cont= mean_vX)

#### MERGE 2015 and 2017 DATA TOGETHER ####
env.all.df <- full_join(env.2015.df, env.2017.df) %>%
                #mutate_at(vars(temp:NP), list(~scale(.))) %>% # scale variables by mean and std. deviation
                select(-pH, -DOC_ugL, -depth_cm) # Remove variables with correlations > 0.5 (see below)

## Test for correlations between variables
   # pairs(env.all.df[, -1])
   # cor(env.all.df[, -1]) > 0.5
   # Removed: pH, DOC, -depth_cm


#### LAT/LONG DATA
lat.long.info <-  read_csv(file.path(dir_input_spatial, "PhormMeta17_LatLong_combined.csv")) %>%
  select(ggkbase_id, year, fork, lat, long) %>%
  rename(site= ggkbase_id)


## Add lat/long data to environmental data
env.all.df <- lat.long.info %>%
               select(site, lat, long) %>%
               right_join(., env.all.df)


#### ADD GENOMES TO ENVIRONMENTAL DATA ####

## List of ANI species 1 genomes
ani.sp1.genomes <- read_tsv(file.path("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/Data/dRep/", "ani_species.tsv")) %>%
  filter(ani_species == 1) %>%
  select(genome)

genome.site <- read_tsv(file.path(dir_input_env, "site_genome_table_2.tsv")) %>%
                mutate(site2= str_replace(.$site, "[A|B]$", "")) %>% # Adjust PH2017 site id to account for A and B samples with the same env. data
                mutate(site3= str_replace(.$site2, "PH2017_[0-9].", "")) %>%
                select(-site, -site2)

env.all.df2 <- env.all.df %>%
                mutate(site2= str_replace(.$site, "[A|B]$", "")) %>%
                mutate(site3= str_replace(.$site2, "PH2017_[0-9].", "")) %>%
                select(-site, -site2) %>%
                full_join(genome.site, .) %>%
                select(-site3)

# Filter by only ANI species 1
env.sp1.df <- env.all.df2 %>%
  filter(env.all.df2$genome %in% ani.sp1.genomes$genome) %>%
  select(-original_sampleID) %>% 
  as.data.frame()

env.CC.df <- env.sp1.df[c(3, 5, 7, 15, 19),]# %>% 
  mutate(across(temp:NP, scale))

#### ANI DATA ####
ani.df <- read_csv(file.path(dir_input_ani, "Ndb_20200422.csv")) %>%
  select(querry, reference, ani) %>%
  spread(key= querry, value= ani) %>%
  select(-reference) %>%
  mutate(genome= str_replace(colnames(.), ".fa", "")) #%>% # Clean up genome names
  #mutate(genome= str_replace(.$genome, "_s25", ""))
colnames(ani.df) <- str_replace(colnames(ani.df), ".fa", "")
#colnames(ani.df) <- c(ani.df$genome, "genome")

ani.df[is.na(ani.df)] <- 0.5 # Add a low ANI value to all NAs

## Subset only species 1
ani.sp1 <- left_join(ani.sp1.genomes, ani.df)
ani.sp1 <- ani.sp1[, colnames(ani.sp1) %in% c(ani.sp1$genome, "genome")] # select only species 1 columns
ani.sp1.dissmat <- cbind(ani.sp1[, 1], (1 - ani.sp1[, -1]))

## Subset only Cedar Creek genomes
ani.CC.dissim <- left_join(ani.sp1.genomes[c(3, 5, 7, 15, 19), 1], ani.df) %>% 
  select(1, c(6, 8, 10, 29, 33)) %>% 
  mutate(across(where(is.numeric), ~ 1 - .x))



## Calculate ANI distances (method Euclidean)
ani.sp1.dist <- as.data.frame(as.matrix(dist(ani.sp1[, -1], method= "euclidean")))
row.names(ani.sp1.dist) <- ani.sp1$genome
colnames(ani.sp1.dist) <- row.names(ani.sp1.dist)

## Add genome as first column
ani.sp1.dist.mat <- as.matrix(data.frame(genome= row.names(ani.sp1.dist), ani.sp1.dist))
ani.sp1.dist.DF <- data.frame(genome= row.names(ani.sp1.dist), ani.sp1.dist, stringsAsFactors = FALSE)


#### RUN GDM MODEL ############################################################

ani.sp1$genome %in% env.sp1.df$genome

## Format data for model
gdm.format <- gdm::formatsitepair(bioData = ani.sp1.dissmat, #ani.sp1.dist.DF,
                             bioFormat = 3,
                             predData=   env.sp1.df,
                             siteColumn = "genome",
                             XColumn = "long",
                             YColumn = "lat")
# DOES NOT CONVERGE

## Cedar Creek only species
gdm.cc.format <- gdm::formatsitepair(bioData = ani.CC.dissim, #ani.sp1.dist.DF,
                                     bioFormat = 3,
                                     predData=   env.CC.df,
                                     siteColumn = "genome",
                                     XColumn = "long",
                                     YColumn = "lat")
# DOES NOT CONVERGE


