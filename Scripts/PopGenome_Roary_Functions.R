## PopGenome statistics on genes extracted from Roary output



#### POP GENOME INFORMATION 

## Valid sites in the alignment= sites without a gap or N, these sites are excluded from calculations

## Regions are each of the nucleotide files included in the folder input into readData
## Populations are organisms within each data file defined to be in differen populations

## Run a command such as diversity.stats(mydata)
## See output list contents get.diversity(mydata)
## See contents of each list get.diversity(mydata)[[1]]


#### EXTRACT CORE GENOME ALIGNMENTS FROM ROARY OUTPUT

## Extract the gene alignments of the core genes from the Roary output
## Must run Roary with the -z flag: "don't delete intermediate files"

## Gene alignments are in "pan_genome_sequences" folder

## roary_path= folder of Roary Output
## core_prop= Proportion of genomes with gene to be considered core (0-1)
## move_files= move core genes to their own folder

extract_core_aln <- function(roary_path, core_prop, move_files= FALSE){
  suppressMessages(require(tidyverse))
  
  if(core_prop > 1){
    stop("core_prop must be between 0 and 1")
  }
  
  ## Read in gene presence/absence table
  gene.df <- suppressMessages(read_tsv(file.path(roary_path, "gene_presence_absence.Rtab"))) %>%
    mutate(gene_counts= rowSums(.[, -1]))
  
  ## Count genomes in table
  genome.count <- sum(str_detect(names(gene.df), "PH"))
  
  ## Filter genes in the core
  gene.core <- gene.df %>%
    filter(gene_counts >= ceiling(genome.count * core_prop)) ## Proportion of genomes considered to be core
  
  if(move_files== TRUE){
    core.genes <- gene.core %>%
      select(Gene) %>%
      pull()
    
    ## List gene alignment sequences output by Roary
    gene.sequences <- list.files(file.path(roary_path, "pan_genome_sequences"))
    
    ## Subset the core sequences
    core.sequences <- gene.sequences[str_replace(gene.sequences, ".fa.aln", "") %in% core.genes]
    
    ## Move the core sequences to a new folder, "core_genome_sequences"
    if(dir.exists(file.path(roary_path, "core_genome_sequences")) == TRUE){
      unlink(file.path(roary_path, "core_genome_sequences"), recursive = TRUE)
    }
    
    dir.create(file.path(roary_path, "core_genome_sequences"))
    file.copy(from= file.path(roary_path, "pan_genome_sequences", core.sequences),
              to= file.path(roary_path, "core_genome_sequences"))
  }
  
  return(gene.core)
}




#### CHECK SEQUENCES TO MAKE SURE DIVISIBLE BY 3
## PopGenome throws an error if sequence lengths are not divisible by 3
## This script will trim the last 1 or 2 nucleotides from a sequence
## to make the length divisible by 3

## core_genes_path= path to gene fasta files
## trim_seq= trim the ends of the sequence to be divisible by 3

trim_seq_length <- function(core_genes_path, trim_seq= FALSE){
  suppressMessages(require(seqinr))
  gene_list <- list.files(core_genes_path)
  
  gene_length_df <- data.frame(geneID= rep(NA, length(gene_list)), 
                               gene_count= rep(NA, length(gene_list)), 
                               gene_length= rep(NA, length(gene_list)), 
                               div3= rep(NA, length(gene_list)))
  
  for(i in seq(1, length(gene_list))){
    
    seq <- read.fasta(file.path(core_genes_path, gene_list[i]))
    seq_length <- getLength(seq)
    seq_num <- length(seq_length)
    div3 <- sum(seq_length %% 3) / seq_num # The %% operator returns the remainder of the division
    
    gene_length_df[i, ] <- c(gene_list[i], seq_num, seq_length[1], div3)
  }
  
  gene_length_df <- gene_length_df %>% 
    as_tibble() %>% 
    mutate_at(vars(gene_count:div3), funs(as.numeric))
  
  ## Trim sequences and write new fasta file
  if(trim_seq== TRUE){
    message("Moving original long fasta files to ", file.path(str_replace("../Roary/Output_PHall_bp90_c90/core_genome_sequences", "\\/[^\\/]*$", ""), "core_long_seqs"))
    message("Writing new fasta files")
    
    dir.create(file.path(str_replace("../Roary/Output_PHall_bp90_c90/core_genome_sequences", "\\/[^\\/]*$", ""), "core_long_seqs"))
    
    trim_seqs <- gene_length_df %>% 
      filter(div3 != 0)
    
    for(i in seq(1, nrow(trim_seqs))){
      gene_name <- pull(trim_seqs[i, "geneID"])
      seq <- read.fasta(file.path(dir_core_genes, gene_name))
      trim_length <- pull(trim_seqs[i, "div3"]) # either 1 or 2
      gene_length <- pull(trim_seqs[i, "gene_length"])
      
      # Trim sequence by 1 or 2 nucleotides
      seq_frag <- getFrag(seq, begin= 1, end= gene_length - trim_length)
      
      ## Move original file to new folder
      file.rename(from= file.path(core_genes_path, gene_name),
                  to= file.path(str_replace("../Roary/Output_PHall_bp90_c90/core_genome_sequences", "\\/[^\\/]*$", ""), "core_long_seqs", paste0("orig_", gene_name)))
      
      # Write fasta file of trimmed sequence to original core_genomes folder
      write.fasta(sequences= seq_frag, names= names(seq), file.out = file.path(core_genes_path, gene_name))
    }
  }

  return(gene_length_df)
}


#### RUN POPGENOME AND CALCULATE POPULATION STATISTICS
## gene_path= folder where .aln files are located from Roary otput
## pops= which of the caynobacterial populations are included in the analysis, should be a vector
## prop_valid_sites= proportion of non-gaps n
## export_tables= either FALSE, or a pathway for where files should be written

run_PopGenome <- function(gene_path, pops= c(1, 2, 3), prop_valid_sites= 0.5, export_tables= FALSE){
  suppressMessages(require(tidyverse))
  suppressMessages(require(PopGenome))

   
  # Whole alignment FASTA
## Read in alignment of the genes
  message("Reading in fasta files to PopGenome")
  mydata <- readData(path= gene_path, format="FASTA", SNP.DATA = FALSE, big.data= TRUE)


#### SET POPULATIONS (i.e. ANI defined species 1 to 4)
#individual.names <- sort(get.individuals(mydata)[[2]])
individual.names <- suppressMessages(read_tsv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Microcoleus_Analysis/GenomesData/genome_ani_species.tsv")) %>% 
  select(genome) %>% 
  #filter(genome != "PH2015_01D_Oscillatoriales_44_513" | genome != "PH2015_02D_Oscillatoriales_45_1038") %>% 
  pull()


message("Assigning populations")

if(all(pops == c(1, 2, 3))){

populations <- list(c(individual.names[c(-1, -3, -c(12:19), -c(25:27), -c(37:39), -48, -49, -53)]),
                    c(individual.names[c(12, 14, 17, 25, 26, 27)]),
                    c(individual.names[c(13, 15:16, 18:19, 37:39, 48:49, 53)]))
}

if(all(pops == c(1, 2))){
  
  populations <- list(c(individual.names[c(-1, -3, -c(12:19), -c(25:27), -c(37:39), -48, -49, -53)]),
                      c(individual.names[c(12, 14, 17, 25, 26, 27)]))
}

if(all(pops == c(1, 3))){
  
  populations <- list(c(individual.names[c(-1, -3, -c(12:19), -c(25:27), -c(37:39), -48, -49, -53)]),
                      c(NA),
                      c(individual.names[c(13, 15:16, 18:19, 37:39, 48:49, 53)]))
}

if(all(pops == c(2, 3))){
  
  populations <- list(c(NA),
                      c(individual.names[c(12, 14, 17, 25, 26, 27)]),
                      c(individual.names[c(13, 15:16, 18:19, 37:39, 48:49, 53)]))
}




mydata <- set.populations(mydata, populations, diploid= FALSE)

mydata.sum <- as.data.frame(get.sum.data(mydata)) %>%
                mutate(geneID= str_replace(row.names(.), ".fa.aln", ""),
                       prop.valid.sites= round(n.valid.sites / n.sites, 4)) %>%
                select(geneID, everything()) %>%
                as_tibble()

#### CALCULATE STATISTICS 
# show.slots(mydata) # slots are the different type of analyses that can be conducted
message("Calculating population genetic statistics")

## Calculate Neutrality statistics
mydata <- neutrality.stats(mydata, detail= TRUE)

## Calculate F_ST and Diversity statistics
mydata <- F_ST.stats(mydata, mode= "nucleotide") # this also calculates diversity statistics

## Calculate MK Test
mydata <- MKT(mydata, new.populations= populations,
              do.fisher.test= TRUE,
              fixed.threshold.fst= FALSE)

## Calculate four-gamete test
# mydata2 <- recomb.stats(mydata)
# mydata3 <- recomb.stats(mydata, new.populations= populations)
# mydata4 <- recomb.stats(mydata, new.populations= TRUE)
# 
# mydata2@region.stats@
# recomb.values <- get.recomb(mydata2)
# recomb.values3 <- get.recomb(mydata3)
# recomb.values4 <- get.recomb(mydata4)
# head(recomb.values4[[1]]) # first population !
# mydata@region.stats@Hudson.RM
# mydata3@RM
# 
# table(recomb.values3[[3]])
# 
# mydata3@Pop_Recomb
# 
# 
# ?recomb.stats


## Extract FST between populations
pairwise.FST <- t(mydata@nuc.F_ST.pairwise) %>%
  as_tibble() %>%
  #select(-`pop1/pop4`, -`pop2/pop4`, -`pop3/pop4`) %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  #gather(key= pops, value= FST, `pop1/pop2`:`pop2/pop3`) %>%
  gather(key= pops, value= FST, contains("pop")) %>%
  arrange(geneID)


## Extract nucleotide diversity between and within populations
nucdiv.btw <- t(mydata@nuc.diversity.between / mydata.sum$n.sites) %>%
  as_tibble() %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  #gather(key= pops, value= nuc.btw.div, `pop1/pop2`:`pop2/pop3`) %>%
  gather(key= pops, value= nuc.btw.div, contains("pop")) %>%
  arrange(geneID)

nucdiv.within <- (mydata@nuc.diversity.within / mydata.sum$n.sites) %>%
  as_tibble() %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  #gather(key= pops, value= nuc.wtn.div, `pop 1`:`pop 3`) %>%
  gather(key= pops, value= nuc.wtn.div, contains("pop")) %>%
  arrange(geneID)

## Extract Tajima D
tajima <- mydata@Tajima.D %>%
  as_tibble() %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  gather(key= pops, value= tajimaD, contains("pop")) %>%
  arrange(geneID)

## Extract MKT data
format_mkt_data <- function(){
  ## Extract results as list and add gene names
  mkt.list <- mapply(cbind, get.MKT(mydata), mydata@region.names)
  
  ## Fill NULL results (where biallialec sites = 0) with a blank DF
  null.elements <- which(do.call(rbind, lapply(get.MKT(mydata), is.null)))
  mkt.list[null.elements] <- lapply(null.elements,
                                    function(x)
                                      mkt.list[x] <- cbind(matrix(NA, nrow= 3, ncol= 9), rep(mydata@region.names[x], 3))) # blank matrix of NAs, with geneID included
  
  ## Transform list into a dataframe
  mkt.df <- do.call(rbind, mkt.list) %>%
    as.data.frame() %>%
    mutate(pops = str_replace(row.names(.), "\\.[0-9].*$", "")) %>%
    mutate(pops = str_replace(pops, "\\.", "/")) %>% 
    rename("geneID" = V10) %>%
    as_tibble() %>%
    mutate_at(vars(P1_nonsyn:fisher.P.value), funs(as.numeric(as.character(.)))) %>% 
    mutate(geneID= as.character(.$geneID)) %>%   
    mutate(geneID = str_replace(as.character(geneID), ".fa.aln", ""))
  #left_join(., feature.annotation, by= "geneID") # combine with annotation information
  return(mkt.df)
}
mkt.df <- format_mkt_data()


#### Combine statistics into data frames
## Remove genes with valid sites < 50% of gene length
pop.stats.btw.df <- left_join(mydata.sum, pairwise.FST) %>% # BETWEEN STATISTICS
  left_join(., nucdiv.btw) %>%
  mutate(geneID= str_replace(.$geneID, ".fa", "")) %>%
  left_join(., mkt.df) %>%
  filter(prop.valid.sites > prop_valid_sites)


pop.stats.wtn.df <- left_join(mydata.sum, nucdiv.within) %>% # WITHIN STATISTICS
  left_join(., tajima) %>%
  mutate(geneID= str_replace(.$geneID, ".fa", "")) %>%
  filter(prop.valid.sites > prop_valid_sites)



#### EXPORT TABLES
if(export_tables != FALSE){
write_tsv(pop.stats.btw.df, path= file.path(export_tables, "PopGenome_btw.tsv"))
write_tsv(pop.stats.wtn.df, path= file.path(export_tables, "PopGenome_wtn.tsv"))
}


#### RETURN DATA FRAMES 

output_list <- list(mydata.sum, pop.stats.wtn.df, pop.stats.btw.df, mydata)
names(output_list) <- c("mydata.sum", "pop.stats.wtn.df", "pop.stats.btw.df", "mydata")
return(output_list)
}




#### CALCULATE PROPORTION OF GENOMES IN A SPECIES THAT POSSESS EACH GENE

## gene_presence_absence= the gene_presence_absence.Rtab matrix from Roary or the output from extract_core_aln function
get_genes_by_species <- function(gene_presence_absence, species_file){
  message("function takes the gene_presence_absence.Rtab matrix from Roary\n or the output from extract_core_aln")
  
  ## Read in table of species ID for each genome
  genome_species <- suppressMessages(read_tsv(species_file))
  
  ## Initialize empty tibble
  gene_occurrence <- tibble(Gene= core_genes$Gene)
  
  ## Loop through each species and count number of genomes in the species that possess each gene
  for(species in 1:3){
    
    ## Extract genomes in a species
    sp_list <- genome_species %>%
      filter(ani_species == species)
    sp_genomes <- names(core_genes)[names(core_genes) %in% sp_list$genome]
    
    ## Initialize empty column names
    var_count <- paste0("gene_count_sp", species)
    var_prop <- paste0("gene_prop_sp", species)
    
    ## Count number of genomes that possess each gene, then calculate proportion
    sp_genes <- core_genes %>%
      select(Gene, sp_genomes) %>%
      mutate(!!var_count := rowSums(.[, -1]),
             !!var_prop := rowSums(.[, -1]) / length(sp_genomes)) %>% 
      select(-contains("PH"))
    
    gene_occurrence <- left_join(gene_occurrence, sp_genes, by= "Gene")
    
  }
  
  ## Reformat pops column to match PopGenome output
  gene_occurrence <- gene_occurrence %>% 
    gather(key= pops, value= prop, contains("prop")) %>% 
    arrange(Gene) %>% 
    mutate(pops= ifelse(pops == "gene_prop_sp1", "pop 1",
                        ifelse(pops == "gene_prop_sp2", "pop 2", "pop 3"))) %>% 
    select(-contains("count")) %>% 
    rename(geneID= Gene)
  
  
  return(gene_occurrence)
}


#### PROCESS TABLES



  # message("Exporting tables")
  # # Fst
  # fst.outlier <- pop.stats.btw.df %>% 
  #   filter(FST < 0.5) %>% 
  #   mutate(outlier= "fst")
  # write_tsv(fst.outlier, file.path(export_tables, "Fst_low.txt"))
  # 
  # # Intra-specific
  # intra.outlier <- pop.stats.wtn.df %>% 
  #   filter(nuc.wtn.div > 0.01) %>% 
  #   mutate(outlier= "intra")
  # write_tsv(intra.outlier, file.path(export_tables, "Intra_specific_high.txt"))
  # 
  # # Inter-specific
  # inter.outlier <- pop.stats.btw.df %>% 
  #   filter(nuc.btw.div > 0.25) %>% 
  #   mutate(outlier= "inter")
  # write_tsv(inter.outlier, file.path(export_tables, "Inter_specific_high.txt"))
  # 
  # # Tajima's D
  # tajima.outlier <- pop.stats.wtn.df %>% 
  #   filter(tajimaD > 2 | tajimaD < -2) %>% 
  #   mutate(outlier= "tajima")
  # write_tsv(tajima.outlier, file.path(export_tables, "Tajima_2.txt"))
  # 
  # 
  # # Read in MK Test results
  # # mkt <- read_tsv(file.path(dir_output_table, "MK_test_sig.txt")) %>% 
  # #   select(geneID, feature, pops, neutrality.index, alpha, fisher.P.value, anno) %>% 
  # #   rename(attribute= feature) %>% 
  # #   mutate(outlier= "MKT",
  # #          pops= str_replace(.$pops, "\\.", "\\/")) %>% 
  # #   separate(anno, into= c("uniref_anno", "uniprot_anno", "kegg_anno"), sep= " -- ")
  # 
  # 
  # 
  # 
  # outlier.master <- full_join(fst.outlier, intra.outlier) %>% 
  #   full_join(., inter.outlier) %>% 
  #   full_join(., tajima.outlier) %>% 
  #   full_join(., mkt) %>% 
  #   select(geneID, attribute, outlier, pops, everything()) %>% 
  #   arrange(geneID)
  # write_tsv(outlier.master, file.path(export_tables, "PopGenome_outlier_genes.txt"))
  # 