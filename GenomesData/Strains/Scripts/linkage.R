## inStrain linkage output

library(tidyverse)
library(ggplot2)

link <- read_tsv("inStrain/inStrain_linkage_output/PH2015_03D.species_1.pid96.linkage.tsv")

calc_haplotypes <- function(df){
  counts <-  as.matrix(df[, c('countAB','countAb','countaB','countab')])
  counts[which(counts != 0)] <-  1
  
  df$haplotype <- ifelse(rowSums(counts) == 4, "h4", 
                         ifelse(rowSums(counts) == 3, "h3",
                                ifelse(rowSums(counts) == 2, "h2",
                                       ifelse(rowSums(counts) == 1, "h1", "error"))))
  return(df)
}


link_haplo <- calc_haplotypes(link)

link_mean <- link_haplo %>% 
  group_by(distance, haplotype) %>% 
  summarize(r2_mean= mean(r2, na.rm= TRUE))



ggplot(data= filter(link_haplo, haplotype == "h3" | haplotype == "h4"), aes(x= distance, y= r2)) +
  geom_hex() +
  facet_grid(.~haplotype) +
  theme_bw()



ggplot(data= filter(link_mean, haplotype == "h3" | haplotype == "h4"), aes(x= distance, y= r2_mean)) +
  geom_point() +
  facet_grid(.~haplotype) +
  theme_bw()



ggplot(data= link_haplo, aes(x= distance, y= r2)) +
  geom_hex() +
  facet_grid(.~haplotype) +
  theme_bw()


table(is.na(link_haplo$r2))

## Why do some 
r2_na <- link_haplo %>% 
  filter(is.na(r2))

h2 <- link_haplo %>% 
  filter(haplotype == "h2")



table(link_haplo$haplotype
      )
