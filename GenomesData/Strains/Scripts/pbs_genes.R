## Phycobilisome gene analyses
## KEGG HMMs 

## Minimum coverage for nucleotide diversity = 15

library(tidyverse)
library(ggplot2)
library(vegan)

## Input data
pbs_dir <- "../Photosystem_genes"
pbs_gene_names <- read_tsv(file.path(pbs_dir, "pbs_antennae_KEGG_hmms.txt"))

species_genomes <- tibble(genome= c("PH2015_12U_Oscillatoriales_45_315", "PH2015_13D_Oscillatoriales_45_19", "PH2017_22_RUC_O_B_Oscillatoriales_46_93"),
                  sample= c("PH2015_12U", "PH2015_13D", "PH2017_22_RUC_O_B"),
                  species= c("species_1", "species_2", "species_3"))

## Function to parse output HMM table
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


pbs_hmm_raw <- read_hmm_table(file.path(pbs_dir, "pbs_hmm_results_table_20191219.txt")) %>% 
  mutate(kegg_hmm= str_replace(query, ".gCluster.*$", "")) %>% 
  left_join(., pbs_gene_names, by= "kegg_hmm") %>% 
  left_join(species_genomes, ., by= "sample") %>% 
  filter(species == "species_1") %>% 
  select(-genome, -sample) %>% 
  arrange(gene_name, eval_full) %>% 
  rename(gene= target)

## Select minimum eval_full value for each contig feature 
## and assign it that gene ID
pbs_hmm <- pbs_hmm_raw %>% 
  # group_by(gene) %>% 
  group_by(gene_name) %>% 
  filter(eval_full == min(eval_full)) %>% 
  ungroup()

test <- pbs_hmm %>% 
  count(gene_name)



#### SNV per mBP ####
gi_df <- read_tsv("inStrain/output_tables/gene_info_df.tsv") %>% 
  filter(species == "species_1") %>% 
  select(-sample)

pbs_gi_df <- gi_df %>% 
  inner_join(., pbs_hmm) %>% 
  filter(breadth > 0.95) %>% # Remove genes with breadth < 0.95
  filter(coverage > 15) %>%  # Remove all genes with coverage <15 (61 genes)
  mutate(SNV_mbp= SNPs_per_bp*1000000)
  


hist(pbs_gi_df$pi)
quantile(pbs_gi_df$pi, 0.9)

high_pi <- pbs_gi_df %>% 
  group_by(gene_name) %>% 
  filter(pi > quantile(pi, 0.8)) %>% 
  ungroup()


high_sites <- high_pi %>% 
  count(site, gene_name) %>% 
  count(site)

unique(high_pi$site)
  
summary(pbs_gi_df$SNV_mbp)


summary(pbs_gi_df$pi)
table(pbs_gi_df$SNV_mbp > 5000)
table(pbs_gi_df$gene_name, pbs_gi_df$SNV_mbp == 0)

#### N:S RATIOS ####

# Read SNV data
snv_pbs <- read_tsv("inStrain/output_tables/snp_mutation_type_df.tsv") %>% 
  inner_join(., pbs_hmm)

snv_pbs %>% 
  count(gene_name, site) %>% 
  count(gene_name)

pbs_gi_df %>% 
  count(gene_name)


snv_pbs_NS_ratios <- snv_pbs %>% 
  filter(mutation_type == "N" | mutation_type == "S") %>% 
  group_by(sample, gene) %>% 
  summarize(
    n= length(mutation_type),
    N_count= sum(str_detect(mutation_type, "N")),
    S_count= sum(str_detect(mutation_type, "S")),
    NS= N_count/S_count) %>% 
  ungroup() %>% 
  filter(n > 1) %>% 
  mutate(site= str_replace_all(sample, "\\.species.*$", "")) %>% 
  left_join(., pbs_gi_df)

  
#### SNV sharing ####
snv_postions_counts <- snv_pbs %>% 
  count(gene_name, gene, position) %>% 
  inner_join(., pbs_gene_names)

unique(snv_pbs$site)

snv_pbs_wide <- snv_pbs %>% 
  select(site, gene_name, position, refFreq) %>% 
  mutate(refFreq= ifelse(refFreq > 0, 1, 0),
         position= str_c(gene_name, position, sep= "-")) %>% 
  pivot_wider(-gene_name, names_from= position, values_from= refFreq, values_fill= list(refFreq = 0))

snv_pbs_dist <- as.matrix(vegdist(as.matrix(snv_pbs_wide[, -1]), 
                              method= "jaccard"))
row.names(snv_pbs_dist) <- snv_pbs_wide$site
colnames(snv_pbs_dist) <- row.names(snv_pbs_dist) 

pbs_nmds <- metaMDS(snv_pbs_dist, trymax= 200, k= 3)

## Set up data to plt NMDS in ggplot
nmds.data <- try(as_tibble(scores(pbs_nmds)))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
nmds.data$site <- try(rownames(snv_pbs_dist))  # create a column of site names, from the rownames of data.scores

ggplot(data= nmds.data, aes(x= NMDS1, y= NMDS2)) +
  geom_point(size= 2) +
  ggrepel::geom_text_repel(aes(label= site), size= 2, show.legend= FALSE) +
  coord_equal() +
 # facet_wrap(~gene_name, ncol= 5) +
  theme_strains


snv_pbs_list <- snv_pbs %>% 
  group_split(gene_name)
names(snv_pbs_list) <- sort(unique(snv_pbs$gene_name))


snv_matrix_list <- map(snv_pbs_list, function(x) {
  x %>% 
    select(site, position, refFreq) %>% 
    mutate(refFreq= ifelse(refFreq > 0, 1, 0),
           position= str_c("p", position)) %>% 
    spread(key= position, value= refFreq, fill= 0)
  })
snv_matrix_list[[4]]

## Calculate Jaccard Distances
snv_distance_list <- map(snv_matrix_list, function(x) {
  dist.mat <- as.matrix(vegdist(as.matrix(x[, -1]), 
                           method= "jaccard"))
  row.names(dist.mat) <- x$site
  return(dist.mat)})

## Run NMDS

nmds_list <- map(snv_distance_list[-4], function(x) {
  nmds <- try(metaMDS(x, trymax= 200, k= 3))
  ## Set up data to plot NMDS in ggplot
  data.scores <- try(as_tibble(scores(nmds)))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
  data.scores$site <- try(rownames(x))  # create a column of site names, from the rownames of data.scores
  return(data.scores)
})

# Remove list elements where the nmds failed (is_tibble == FALSE)
nmds_list2 <- nmds_list[unlist(map(nmds_list, is_tibble))]

# Make into dataframe
nmds_df <- bind_rows(nmds_list2, .id= "gene_name")



#### MAKE FIGURES
source("Scripts/ggplot_themes.R")


nmds_plot <- ggplot(data= nmds_df, aes(x= NMDS1, y= NMDS2)) +
  geom_point(size= 1) +
  ggrepel::geom_text_repel(aes(label= site), size= 1.5, show.legend= FALSE) +
  coord_equal() +
  facet_wrap(~gene_name, ncol= 5) +
  theme_strains
nmds_plot

## SNVs per mBP and N:S Ratio
ggplot(filter(snv_pbs_NS_ratios, is.infinite(snv_pbs_NS_ratios$NS) == FALSE), aes(x= SNV_mbp, y= NS)) +
  geom_vline(xintercept = 1500, size= 0.5, linetype= "dashed") +
  geom_point(aes(color= gene_name)) +
  labs(x= "SNVs per mbp", y= "N:S ratio") +
  scale_color_discrete(name= "PBS \ngene") +
  facet_wrap(~gene_name, ncol= 4) +
  # scale_x_continuous(limits= c(0, 16000),
  #                    breaks= seq(0, 15000, by= 2500),
  #                    labels= c("0", "", "5000", "", "10000", "", "15000"),
  #                    expand= c(0.01, 0)) +
  # scale_y_continuous(expand= c(0.02, 0)) +
  theme_strains
#ggsave(last_plot(), filename= "atx_NS_SNPmbp.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)


ggplot(filter(snv_pbs_NS_ratios, is.infinite(snv_pbs_NS_ratios$NS) == FALSE, is.na(gene_name) == FALSE), aes(x= reorder(gene_name, -NS), y= NS)) +
  #geom_vline(xintercept = 1500, size= 0.5, linetype= "dashed") +
  geom_boxplot() +
  geom_point(color= "gray50", position= "jitter", size= 2, alpha= 0.5) +
  labs(x= "Phycobilisome gene", y= "N:S ratio") +
  scale_y_continuous(expand= c(0.02, 0)) +
  theme_strains +
  theme(axis.text.x= element_text(angle= 45, vjust= 0.9, hjust= 0.9))

ggsave(last_plot(), filename = "pbs_NS.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



## Nucleotide diversity
ggplot(pbs_gi_df, aes(x= pi)) +
  geom_histogram(boundary= 1, binwidth= 0.0001) +
  theme_strains


ggplot(pbs_gi_df, aes(x= reorder(gene_name, -pi), y= pi)) +
  geom_boxplot() +
  #geom_point(alpha= 0.3) +
  scale_y_log10(breaks= c(0.0001, 0.001, 0.01),
    labels= c("0.0001", "0.001", "0.01")) +
  annotation_logticks(sides= "l") +
  theme_strains

avg_pi_order <- pbs_gi_df %>% 
  group_by(gene_name) %>% 
  summarize(median_pi= mean(pi)) %>% 
  arrange(-median_pi) %>% 
  pull(gene_name)


test <- pbs_gi_df %>% 
  select(gene_name, site, pi) %>% 
  group_by(gene_name) %>% 
  arrange(-pi) %>% 
  mutate(rank= row_number()) %>% 
  ungroup()


ggplot(test, aes(x= rank)) +
  geom_histogram(binwidth= 1) +
  facet_wrap(~site, ncol= 8) +
  theme_strains

ggplot(test, aes(x= rank)) +
  geom_density(aes(color= site)) +
  #facet_wrap(~site, ncol= 3) +
  theme_strains
  


ggplot(test, aes(x= reorder(site, -rank), y= pi)) +
  geom_point()+
  scale_y_log10(breaks= c(0.0001, 0.001, 0.01),
                labels= c("0.0001", "0.001", "0.01")) +
  annotation_logticks(sides= "l") +
  facet_wrap(~gene_name, nrow= 1, scales= "free_x") +
  theme_strains +
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5))

pbs_gi_df %>% 
  group_by(gene_name, site) %>% 
  summarize(mean_pi= mean(pi)) %>% 
  arrange(-median_pi) %>% 
  pull(gene_name)



ggplot(pbs_gi_df, aes(x= factor(gene_name, level= avg_pi_order), y= pi)) +
  geom_boxplot() +
  labs(x= "Phycobilisome gene", y= "Nucleotide diversity") +
  scale_y_log10(breaks= c(0.0001, 0.001, 0.01),
                labels= c("0.0001", "0.001", "0.01")) +
  annotation_logticks(sides= "l") +
  scale_x_discrete(expand= c(0, 1)) +
  theme_strains +
  theme(axis.text.x= element_text(angle= 45, vjust= 0.9, hjust= 0.9))
ggsave(last_plot(), filename = "pbs_nuc_div.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")


ggplot(pbs_gi_df, aes(x= gene_name, y= pi)) +
  geom_boxplot() +
  theme_str

ggplot(pbs_gi_df, aes(x= reorder(site, -pi), y= pi)) +
  geom_boxplot() +
  #geom_point(aes(color= gene_name)) +
  #geom_point(alpha= 0.3) +
  scale_y_log10(breaks= c(0.0001, 0.001, 0.01),
                labels= c("0.0001", "0.001", "0.01")) +
  annotation_logticks(sides= "l") +
  theme_strains +
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5))

ggplot(pbs_gi_df, aes(x= reorder(site, -pi), y= pi)) +
  #geom_boxplot() +
  geom_point() +
  scale_y_log10(breaks= c(0.0001, 0.001, 0.01),
                labels= c("0.0001", "0.001", "0.01")) +
  annotation_logticks(sides= "l") +
  facet_wrap(~gene_name, nrow= 1, scales= "free_x") +
  theme_strains +
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5))

ggplot(pbs_gi_df, aes(x= coverage, y= pi)) +
  #geom_boxplot() +
  geom_point(aes(color= gene_name)) +
 theme_strains

ggplot(pbs_gi_df, aes(x= coverage, y= pi)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_strains

ggplot(filter(pbs_gi_df, gene_name == "cpeU"), aes(x= coverage, y= pi)) +
  geom_point() +
  facet_wrap(.~gene, ncol= 4) +
  theme_strains


filter(pbs_gi_df, gene_name == "cpeU") %>% 
  count(gene)

## SNV sharing

ggplot(data= snv_postions_counts, aes(x= n)) +
  geom_histogram(binwidth = 1, fill= "black", color= "gray50") +
  labs(x= "Number of genomes sharing a SNV position", y= "Number of SNV positions shared") +
  scale_x_continuous(breaks= seq(1, 12, by= 2),
                     expand= c(0, 0)) +
  scale_y_continuous(expand= c(0, 0)) +
  #facet_grid(gene_name ~.) +
  theme_strains
#ggsave(last_plot(), filename= "pbs_SNV_sharing_hist.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)





gene_counts <- pbs_hmm %>% 
  group_by(gene, species, sample) %>% 
  count() %>% 
  arrange(-n)

min_eval <- pbs_hmm %>% 
  group_by(gene, species) %>% 
  filter(eval_full == min(eval_full)) %>% 
  arrange(gene, species)

ggplot(data= pbs_hmm_raw, aes(x= reorder(gene, -score_full), y= score_full)) +
  geom_point() +
  facet_wrap(~gene_name, ncol= 6, scales= "free_x") +
  theme(axis.text.x= element_blank())

ggplot(data= pbs_hmm, aes(x= reorder(target, -score_full), y= score_full)) +
  geom_point(aes(color= sample)) +
  facet_grid(sample~gene, scales= "free_x") +
  theme(axis.text.x= element_blank())

### SNV network
library(tidygraph)
library(ggraph)

create_network_df <- function(matrx){
sim.mat <- (1 - matrx)
colnames(sim.mat) <- rownames(matrx)
dist.df <- as_tibble(sim.mat)
dist.df$name1 <- rownames(sim.mat)

dist.df.l <- dist.df %>% 
  pivot_longer(names_to= "name2", values_to= "jaccard_sim", cols= -name1)
}

format_network <- function(df){
  
  nodes <- unique(c(df$name1, df$name2)) %>% 
    tibble() %>% 
    #filter(!str_detect(., "_B$")) %>% 
    rowid_to_column(.) %>% 
    rename(id= rowid, label= ".")
  
  
  edge.list <- df %>% 
    select(name1, name2)
  # rename(from= name1, to= name2)# %>% 
  # filter(!str_detect(from, "_B$")) %>% 
  # filter(!str_detect(to, "_B$"))
  
  
  edges <- edge.list %>% 
    left_join(., nodes, by= c("name1" = "label")) %>% 
    rename(id_from= id) %>% 
    left_join(., nodes, by= c("name2" = "label")) %>%  
    rename(id_to= id) %>% 
    ungroup() %>%
    left_join(., select(df, name1, name2, jaccard_sim, gene_name)) %>% 
    select(-name1, -name2)
  
  
  network_object <- tbl_graph(nodes= nodes, edges= edges, directed= FALSE) #%>% 
    # activate(nodes) %>%
    # mutate(degree = centrality_degree())
    # 
  return(network_object)
}

#cpcF.net <- format_network(cpcF.dist.df.l)

network_dist_list <- map(snv_distance_list, function(x)
  network_df <- create_network_df(x))

high_sim <- network_dist_list %>% 
  bind_rows(.id= "gene_name") %>% 
  filter(name1 != name2) %>% 
  filter(jaccard_sim > 0.95)


network_objects_list <- map2(snv_distance_list, names(snv_distance_list), function(x, x_name){
  network_df <- create_network_df(x)
  network_object <- format_network(network_df)

  return(network_object)
})

network_plots <- map2(network_objects_list, names(network_objects_list), function(x, x_name){
  network_graph <- ggraph(x, layout= "linear", circular= TRUE) + 
    geom_node_point(size= 3) + 
    geom_edge_link(aes(width= jaccard_sim, color= jaccard_sim, alpha= jaccard_sim)) + 
    scale_edge_color_viridis(option= "viridis",
                             limits= c(0, 1)) +
    scale_edge_width(guide= FALSE) +
    scale_edge_alpha(guide= FALSE) +
    geom_node_text(aes(label = label), repel = TRUE, size= 3) +
    ggtitle(x_name) +
    theme_graph()
  
  ggsave(network_graph, filename = str_c(x_name,".png"), height= 180*0.75, width= 180, units= "mm", dpi= 320,
         path= "Output_figures/pbs_network")
}
)


high_sim_network <- format_network(high_sim)

ggraph(high_sim_network, layout= "igraph", algorithm= "kk") + 
  #geom_edge_link(aes(width= jaccard_sim, color= jaccard_sim, alpha= jaccard_sim)) + 
  geom_edge_fan(aes(color= gene_name)) + 
  geom_node_point(size= 3) + 
  scale_edge_color_discrete(limits= c("cpeT", "cpeU", "cpcF", "cpcE", "cpcD", "cpcA")) +
  # scale_edge_width(guide= FALSE) +
  # scale_edge_alpha(guide= FALSE) +
  geom_node_text(aes(label = label), repel = TRUE, size= 3) +
  theme_graph()
ggsave(last_plot(), filename = "pbs_network.png", height= 180*0.75, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")

unique(high_sim$gene_name)

count(high_sim, gene_name) %>% 
  arrange(-n)


ggraph(cpcF.net, layout= "linear", circular= TRUE) + 
  geom_node_point(size= 3) + 
  geom_edge_link(aes(width= jaccard_sim, color= jaccard_sim, alpha= jaccard_sim)) + 
  scale_edge_color_viridis(option= "viridis") +
  scale_edge_width(guide= FALSE) +
  scale_edge_alpha(guide= FALSE) +
  geom_node_text(aes(label = label), repel = TRUE, size= 3) +
  theme_graph()











map(seq(1:11), function(x) map(seq(2:12), function(y) {
    pair <- str_c(cpcF.df[x, ]$site, cpcF.df[y, ]$site, sep= "-")
    snv_sum <- sum(colSums(cpcF.df[c(x,y), -1]) == 2)
    vec.out <- c(cpcF.df[x, ]$site, cpcF.df[y, ]$site, pair, snv_sum)
}
))
colSums(cpcF.df[c(11, 6), -1]) 

cpcF.df[11, ]$site
cpcF.df[6, ]$site

cpcF.df[1, ]$site


