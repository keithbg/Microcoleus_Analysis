library(tidyverse)
library(ggplot2)


#### INPUT FILES
species_lookup <- read_tsv("inStrain_sample_species_lookup.tsv") %>% 
  mutate(species= str_c("species_", species_present)) %>% 
  rename(site= sample) %>% 
  select(-species_present)
in_dir <- "inStrain/inStrain_gene_profile_output"
gi_files <- list.files(in_dir, pattern= "pid96_gene_info")


## FORMAT gene info data
make_gi_df <- function(dir_in, gi_file_names, samp_file){
  require(tidyverse)
  
  # Read in gene info files
  gi_list <- map(gi_file_names, function(x) suppressMessages(read_tsv(file.path(dir_in, x))) %>% 
                   mutate(pi= 1- clonality,
                          sample= str_replace(x, "_gene_info.tsv", "")) %>% 
                   mutate(site= str_split(sample, "\\.")[[1]][1],
                          species= str_split(sample, "\\.")[[1]][2]))
  names(gi_list) <- str_replace(gi_file_names, ".tsv", "")
  
  # Match samples in list with samples with the selected ANI species
  #sp_names <- do.call(rbind, map(names(gi_list), function(x) any(str_detect(x, samp_file$sample))))
  #gi_list_sp_subset <- gi_list[sp_names]
  
  # Transform to a data frame
  #gi_df <- as_tibble(do.call(rbind, gi_list_sp_subset))
  gi_df <- as_tibble(do.call(rbind, gi_list))
  
  return(gi_df)
}


gi_df <- make_gi_df(dir_in = in_dir,
                        gi_file_names = gi_files,
                        samp_file = species_lookup) %>% 
  left_join(species_lookup, ., by= c("site", "species"))  # Filter to only include species recovered from each site
  
  

cov_quantiles <- gi_sp1_df %>% 
  filter(coverage > 5) %>% # coverage threshold specified in inStrain
  group_by(sample) %>% 
  summarize(
    quant_0.1= quantile(coverage, probs= 0.1, na.rm= TRUE),
    quant_0.5= quantile(coverage, probs= 0.5, na.rm= TRUE),
    quant_0.9= quantile(coverage, probs= 0.9, na.rm= TRUE))

cov_quantiles <- gi_df %>% 
  filter(coverage > 5) %>% # coverage threshold specified in inStrain
  group_by(sample) %>% 
  summarize(
    quant_0.05= quantile(coverage, probs= 0.05, na.rm= TRUE),
    quant_0.1= quantile(coverage, probs= 0.1, na.rm= TRUE),
    quant_0.5= quantile(coverage, probs= 0.5, na.rm= TRUE),
    quant_0.9= quantile(coverage, probs= 0.9, na.rm= TRUE),
    quant_0.95= quantile(coverage, probs= 0.95, na.rm= TRUE))

## FILTER OUT LOW COVERAGE GENES ##
# Deciding to filter out the lower 10% of the coverage distributions
# for samples there is a long low coverage tail of the distribution. 
# These low distribution samples may be a result of mis-mapping or poor binning
# These values also contain all the outlier pi values >0.1, which makes these high pi values suspect
# It will be more conservative to remove them from analyses
# The high end of the distribution is tighter, so will not be filtered

# Remove breadth <=0.9
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


gi_filt <- map(cov_quantiles$sample, function(x) filter_coverage_breadth(gene_df= gi_df, 
                                                                             quant_df = cov_quantiles,
                                                                             breadth_thresh= 0.9,
                                                                             samp_name = x))
gi_filt_df <- do.call(rbind, gi_filt)




ggplot(data= gi_filt_df) +
  geom_histogram(aes(x= coverage, fill= species)) +
  facet_wrap(~sample, nrow= 8, scales= "free_x")

ggplot(data= gi_filt_df) +
  geom_point(aes(x= coverage, y= pi)) +
  facet_grid(.~species) +
  theme_bw()

ggplot(data= gi_filt_df) +
  geom_boxplot(aes(x= sample, y= pi)) +
  theme_bw()

ggplot(data= gi_filt_df) +
  geom_point(aes(x= gene, y= pi, color= species)) +
  theme(axis.text.x= element_blank()) +
  facet_wrap(~sample, nrow= 8, scales= "free_x")



#### SUMMARIZE PI VALUES ACROSS THE GENOME ####
gi_sp1_filt_summary <- gi_sp1_filt_df %>% 
  group_by(sample) %>% 
  summarize(n= length(pi),
            mean_pi= mean(pi, na.rm= TRUE),
            median_pi= median(pi, na.rm= TRUE),
            sd_pi= sd(pi, na.rm= TRUE),
            min_pi= min(pi, na.rm= TRUE),
            max_pi= max(pi, na.rm= TRUE),
            median_cov= median(coverage, na.rm= TRUE)) %>% 
  ungroup() %>% 
  mutate(ggkbase_id= str_replace(sample, ".species_1.pid96", ""))

gi_filt_summary <- gi_filt_df %>% 
  group_by(sample, site, species) %>% 
  summarize(n= length(pi),
            mean_pi= mean(pi, na.rm= TRUE),
            median_pi= median(pi, na.rm= TRUE),
            sd_pi= sd(pi, na.rm= TRUE),
            min_pi= min(pi, na.rm= TRUE),
            max_pi= max(pi, na.rm= TRUE),
            median_cov= median(coverage, na.rm= TRUE)) %>% 
  ungroup() #%>% 
  #mutate(ggkbase_id= str_replace(sample, ".species_1.pid96", ""))



## Add watershed area to the data frame
dir_input_watershed <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/EnvData"
watershed.area <-
  read_tsv(file.path(dir_input_watershed, "PhormMeta17_WatershedArea_Combined.tsv")) %>% 
  select(ggkbase_id, watershed_km2) %>% 
  rename(site= ggkbase_id)
  
gi_filt_summary <- left_join(gi_filt_summary, watershed.area, by= "site")



ggplot(data= gi_filt_summary) +
  geom_point(aes(x= sample, y= sd_pi, color= species)) +
  theme_bw()

ggplot(data= gi_filt_summary) +
  geom_point(aes(x= sample, y= median_pi, shape= species)) +
  geom_point(aes(x= sample, y= mean_pi, shape= species), color= "red") +
  theme_bw()


ggplot(data= gi_filt_summary, aes(x= watershed_km2, y= median_pi)) +
  geom_point(aes(color= species), size= 3) +
  labs(x= expression(paste("Watershed area (", km^{2}, ")")), y= expression(paste("Median ", pi))) +
   scale_x_continuous(#limits= c(0, 2000),
                      expand= c(0.02, 0)) +
  scale_y_continuous(limits= c(0, 0.0023), 
                     expand= c(0.02, 0)) +
  scale_color_manual(values= c("black", "tomato", "dodgerblue"),
                     labels= c("1", "2", "3"),
                     name= "Species")+
  theme_bw(base_size= 20)

ggplot(data= gi_filt_summary, aes(x= watershed_km2, y= median_pi)) +
  geom_point(aes(color= species), size= 3) +
  labs(x= expression(paste("Watershed area (", km^{2}, ")")), y= expression(paste("Median ", pi))) +
  scale_x_continuous(#limits= c(0, 2000),
                     expand= c(0.02, 0)) +
  #scale_x_log10(#limits= c(0, 2000),
  #  expand= c(0.02, 0)) +
  scale_y_continuous(limits= c(0, 0.0023), 
                     expand= c(0.02, 0)) +
  facet_grid(.~species) +
  theme_bw(base_size= 20)



ggplot(data= gi_filt_summary, aes(x= watershed_km2, y= median_cov)) +
  geom_point(aes(color= species, shape= species), size= 3) +
  labs(x= expression(paste("Watershed area (", km^{2}, ")")), y= "Median coverage") +
  #scale_x_continuous(expand= c(0.02, 0)) +
  scale_x_log10(expand= c(0.02, 0)) +
  geom_smooth(method= "lm") +
 # scale_y_continuous(limits= c(0, 0.0023), expand= c(0.02, 0)) +
  theme_classic(base_size= 20)


fit1 <- lm(median_pi ~ log10(watershed_km2), data= filter(gi_filt_summary, watershed_km2 < 2000))
summary(fit1)
plot(fit1)

getwd()



summary(lm(median_pi ~ watershed_km2, data= gi_sp1_filt_summary))






