## Concatenate ggkbase features tables (.ql files) from the three species


read_ql <- function(file_name, dir_input_ql, speciesID){
require(tidyverse)

## Read and format .ql file (files formated slightly differently from ggkbase for PH2015 and PH2017 samples)
  if(str_detect(file_name, "PH2015")){
    ql.df <- suppressMessages(read_tsv(file.path(dir_input_ql, file_name), col_names = FALSE))
    if(length(names(ql.df)) == 14){
      names(ql.df) <-  c("genome", "gene", "contig_ID", "feature_num",  "contig_length", "GC", "coverage", "codon_table", "winning_taxonomy", "begin_position", "end_position", "complement" , "annotation", "orf_taxonomy")
    }
    
    if(length(names(ql.df)) == 15){
      names(ql.df) <-  c("genome", "gene", "contig_ID", "feature_num",  "contig_length", "GC", "coverage", "codon_table", "winning_taxonomy", "begin_position", "end_position", "complement" , "annotation", "orf_taxonomy", "db_info")
    }
  }
  
  if(str_detect(file_name, "PH2017")){
    ql.df <- suppressMessages(read_tsv(file.path(dir_input_ql, file_name), col_names = FALSE))
    
    if(length(names(ql.df)) == 14){
      names(ql.df) <-  c("gene", "contig_ID", "feature_num",  "contig_length", "GC", "coverage", "codon_table", "winning_taxonomy", "begin_position", "end_position", "complement" , "annotation", "orf_taxonomy", "db_info")
    }
    
    # if(length(names(ql.df)) == 15){
    #   names(ql.df) <-  c("orf", "contig_ID", "gene", "feature_num",  "contig_length", "GC", "coverage", "codon_table", "winning_taxonomy", "begin_position", "end_position", "complement" , "annotation", "db_info")
    # }
  }
  
  # Separate annotation column into the 3 different databases
  clean.up.annotation.regex <- " Tax.*$| evalue.*$| \\{.*$| \\(.*$|,.*$| bin.*$"
  ql.df.clean.anno <- ql.df %>%
    separate(annotation, into= c("uniref_all", "uniprot_all", "kegg_all"), sep= "__ ") %>%
    mutate(uniref_anno= Hmisc::capitalize(str_replace(.$uniref_all, clean.up.annotation.regex, "")),
           uniprot_anno= Hmisc::capitalize(str_replace(.$uniprot_all, clean.up.annotation.regex, "")),
           kegg_anno= Hmisc::capitalize(str_replace(.$kegg_all, clean.up.annotation.regex, "")),
           gene_length= end_position - begin_position + 1) %>%
    mutate(species= speciesID) %>% 
    select(species, gene, contig_ID, feature_num, contig_length, gene_length, begin_position, end_position, contains("anno"), contains("all"))
  
  return(ql.df.clean.anno) 
}

sp1_ql <- read_ql(file_name = "PH2015_12U_Oscillatoriales_45_315.ql", 
                  speciesID= "species_1", 
                  dir_input_ql = "../ggkbase_features_tables")
sp2_ql <- read_ql(file_name = "PH2015_13D_Oscillatoriales_45_19.ql", 
                  speciesID= "species_2", 
                  dir_input_ql = "../ggkbase_features_tables")
sp3_ql <- read_ql(file_name = "PH2017_22_RUC_O_B_Oscillatoriales_46_93.ql", 
                  speciesID= "species_3", 
                  dir_input_ql = "../ggkbase_features_tables")

anno_df <- rbind(sp1_ql, sp2_ql, sp3_ql)

write_tsv(anno_df, path= "inStrain/ggkbase_anno.tsv")


