## Recombination: Haplotype frequencies and popCOGenT results
# popCOGenT clonal divergence cutoff from paper: 0.000355362 
# (https://github.com/philarevalo/PopCOGenT/blob/master/src/PopCOGenT/cluster.py)

#setwd("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis")

## Results from popCOGenT

library(tidyverse)
library(tidygraph)
library(ggraph)
source("Scripts/ggplot_themes.R")





#### IMPORT DATA ####
## Haplotype frequencies
#haplo.freq <- read_tsv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/Strains/inStrain/output_tables/haplotype_freqs.tsv")
haplo.freq <- read_tsv("Data/inStrain_data/linkage_haplo_freqs_v1.4.tsv")

##  PopCOGentT results
#pgt <- read_csv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/popCOGenT/popCOGenT_RESULTS.csv") %>% 
pgt <- read_csv("Data/popCOGenT_data/popCOGenT_RESULTS.csv") %>% 
  rename(id= X1, name1= `Strain 1`, name2= `Strain 2`, divergence= `Initial divergence`, alignment_length= `Alignment size`, g1_length= `Genome 1 size`, g2_length= `Genome 2 size`,
         ssd_obs= `Observed SSD`, ssd_95CI_low= `SSD 95 CI low`, ssd_95CI_high= `SSD 95 CI high`) %>% 
  mutate(ssd_95CI_int= ssd_95CI_high - ssd_95CI_low)

## Identify high length bias nodes and edges
lb.high.cutoff <- quantile(pgt$ssd_obs, probs= 0.98) # Identify high length bias 98th percentils

pgt.high <- pgt %>% 
  filter(ssd_obs > lb.high.cutoff)
pgt.high.names <- unique(c(pgt.high$name1, pgt.high$name2))

## Infograph network
# infomap <- igraph::read_graph("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/popCOGenT/infomap/sp1_0.000355362.txt.unclust.graphml",
#                               format= "graphml") %>% 
infomap <- igraph::read_graph("Data/popCOGenT_data/infomap_data/sp1_0.000355362.txt.unclust.graphml",
                                format= "graphml") %>% 
  as_tbl_graph() %>% 
  # apply the high length bias cutoffs
  activate(nodes) %>% 
  mutate(lb.high.node= ifelse(id %in% pgt.high.names, id, NA)) %>% 
  mutate(lb.high.node= str_replace(lb.high.node, "_s25.*$|_Oscill.*$", "")) %>%
  activate(edges) %>% 
  mutate(lb.high.edge= ifelse(weight > lb.high.cutoff, "Y", "N"))



#### FIGURES ####
## Haplotypes
haplo.boxplot <- ggplot(data= haplo.freq, aes(x= haplotype, y= freq)) +
  geom_boxplot(aes(fill= haplotype)) +
  labs(x= "Haplotype", y= "Frequency") +
  scale_x_discrete(labels= c("H1", "H2", "H3", "H4")) +
  scale_y_continuous(limits= c(0, 1), 
                     breaks= seq(0, 1, by= 0.1),
                     labels= c("0", "0.1",  "0.2",  "0.3", "0.4",  "0.5",
                               "0.6", "0.7",  "0.8",  "0.9",  "1.0"),
                     expand= c(0.01, 0)) +
  scale_fill_manual(values= wes_palette("FantasticFox1")[2:5],
                    guide= FALSE) +
  facet_rep_grid(.~species, scales= "free_y", labeller= labeller(species = as_labeller(c(`species_1` = "Species 1",
                                                                                         `species_2` = "Species 2",
                                                                                         `species_3` = "Species 3")))) +
  theme_strains
#haplo.boxplot

## Length bias
length.bias.boxplot <- ggplot(pgt, aes(x= "Species 1 genomes", y= ssd_obs)) +
  geom_hline(yintercept = lb.high.cutoff, size= 0.2, linetype= "dotted") +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=27.9026), fill= "gray80", alpha= 50) + #max negative selection cutoff from popcogent_NegSelCutoff.R
  geom_boxplot(outlier.size = 0.2, outlier.color= "black", lwd= 0.2) +
  labs(x= "", y= "Length bias") +
  scale_y_log10(limits= c(7, 6000),
                breaks= c(10, 50, 100, 500, 1000, 5000),
                expand= c(0, 0)) +
  annotation_logticks(side= "l",
                      size= 0.3) +
  theme_strains +
  theme(text= element_text(size= 10),
        axis.text.y = element_text(size= 7),
        axis.text.x = element_text(size= 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 10),
        axis.ticks.y = element_line(size= 0.3))
#length.bias.boxplot


##  PopCOGenT network graph
infomap.graph <- #ggraph(infomap) + 
  ggraph(infomap, layout= 'igraph', algorithm= "fr") + 
  geom_edge_link(color= "gray85", width= 0.2) + 
  geom_node_point(color= "black", size= 2) + 
  scale_edge_width_discrete(range= c(0.2, 1.5), guide= FALSE) +
  theme_graph() +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit= "cm"))
#infomap.graph


## Combine figures
popCOGenT.panels <- plot_grid(length.bias.boxplot, infomap.graph,
                              nrow= 1,
                              rel_widths = c(0.4, 1),
                              labels= c("B", "C"))
#popCOGenT.panels


recombination.figure <- plot_grid(haplo.boxplot, popCOGenT.panels,
                                  nrow= 2,
                                  rel_heights= c(0.75, 1),
                                  labels= c("A"))
ggsave(recombination.figure, filename= "Fig_6_v1.4.png" , height= 180, width= 180, units= "mm", dpi= 320,
       path= "Output_figures")



