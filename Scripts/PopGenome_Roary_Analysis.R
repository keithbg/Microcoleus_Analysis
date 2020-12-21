## PopGenome statistics on genes extracted Roary output


#### POP GENOME INFORMATION ####################################################

## Valid sites in the alignment= sites without a gap or N, these sites are excluded from calculations

## Regions are each of the nucleotide files included in the folder input into readData
## Populations are organisms within each data file defined to be in differen populations

## Run a command such as diversity.stats(mydata)
## See output list contents get.diversity(mydata)
## See contents of each list get.diversity(mydata)[[1]]

#### Libraries #################################################################
library(tidyverse)
library(ggplot2)
library(ggExtra)
source("Scripts/PopGenome_Roary_Functions.R")

#### FILE PATHS ################################################################
dir_roary <- "Data/Roary_data/Output_PHall_bp90_c90"
dir_core_genes <- "Data/Roary_data/Output_PHall_bp90_c90/core_genome_sequences"
dir_out_table <- "Output_tables"
dir_out_figs <- "Output_figures"


## Extract core genes and put in "core_genomes_sequences" folder
#core_genes <- extract_core_aln(roary_path= dir_roary, core_prop= 0.2, move_files = TRUE)
core_genes <- extract_core_aln(roary_path= dir_roary, core_prop= 0.2, move_files = FALSE)


## Run Bash script using terminal from "core_genes_sequences" folder to format fasta headers for PopGenome
"cd Data/Roary_data/Output_PHall_bp90_c90/core_genome_sequences"
"bash Scripts/Roary_format_core_fasta_header.sh"

## Trim sequences that are not divisible by 3 (i.e. 3 nucleotides in a codon)
long_seqs <- trim_seq_length(core_genes_path= dir_core_genes, trim_seq= TRUE)


## Run PopGenome
# prop_valid_sites= 0.5, export_tables= FALSE
pg.list <- run_PopGenome(gene_path = dir_core_genes, pops= c(1, 2, 3), prop_valid_sites= 0.5, export_tables= dir_out_table)
#save(pg.list, file= "pg.list.RData")
#load("pg.list.RData")

pg.btw <- pg.list[["pop.stats.btw.df"]]
pg.wtn <- pg.list[["pop.stats.wtn.df"]]
#pg.object <- pg.list[["mydata"]]


#pg.object@region.names[1]
#pg.object@region.data@biallelic.sites[[1]]
# pg.object@region.stats@biallelic.structure[[1]]["POP"]
# pg.object@region.stats@biallelic.structure[[1]]["POP"][[1]][2, ]

## Identify genes not present in all genomes
gene_pa <- get_genes_by_species(gene_presence_absence = core_genes,
                                species_file= "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/genome_ani_species.tsv")




## Filter by proportion of genomes within a species the gene is within
## Gene must be in at least 50% of genomes in a species to be considered for within species statistics
pg.wtn.filt <- inner_join(pg.wtn, filter(gene_pa, prop > 0.5), by= c("geneID", "pops"))
length(unique(pg.wtn.filt$geneID)) # get number of genes for wtn analyses

## Genes must be in at least 50% of each species for all between species comparisons
genes.btw <- pg.wtn.filt %>% 
  group_by(geneID) %>% 
  summarize(pop1= any(str_detect(pops, "pop 1")),
            pop2= any(str_detect(pops, "pop 2")),
            pop3= any(str_detect(pops, "pop 3"))) %>% 
  ungroup() %>% 
  group_by(geneID) %>% 
  mutate(`pop1/pop2`= all(pop1, pop2),
         `pop1/pop3`= all(pop1, pop3),
         `pop2/pop3`= all(pop2, pop3)) %>% 
  select(geneID, `pop1/pop2`, `pop1/pop3`, `pop2/pop3`) %>% 
  gather(key= pops, value= present, `pop1/pop2`:`pop2/pop3`) %>% 
  filter(present == TRUE) %>% 
  select(geneID, pops)

pg.btw.filt <- inner_join(pg.btw, genes.btw, by= c("geneID", "pops"))
length(unique(pg.btw.filt$geneID)) # get number of genes for btw analyses


## Write output tables of filtered genes (including gene annotations from Roary (i.e. Prokka))
gene_annotations <- read_csv(file.path(dir_roary, "gene_presence_absence.csv")) %>%
  select(Gene, Annotation) %>% 
  rename(geneID= Gene, annotation= Annotation)

inner_join(pg.wtn.filt, gene_annotations, by= "geneID") %>% 
  write_tsv(path= file.path(dir_out_table, "PopGenome_wtn_filt.tsv"))

left_join(pg.btw.filt, gene_annotations, by= "geneID") %>% 
  write_tsv(path= file.path(dir_out_table, "PopGenome_btw_filt.tsv"))
