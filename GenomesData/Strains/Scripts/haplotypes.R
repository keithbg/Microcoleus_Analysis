# Analyze linkage decay from strainRep2.py
## strainRep2.py created by Alex Crits-Christoph: https://github.com/alexcritschristoph/strains_analysis


#### Libraries #################################################################
library(tidyverse)
library(ggplot2)
library(purrr)
################################################################################

dir_out_fig <- "Output_figures"

## ggplot theme
theme_snv <- theme(panel.grid = element_blank(),
                   plot.margin = unit(c(1, 1, 1, 1), "cm"),
                   text = element_text(size= 14),
                   plot.background = element_rect(fill = "transparent"), # bg of the plot
                   panel.background = element_rect(fill= "transparent", color="black"),
                   axis.text = element_text(colour="black"),
                   axis.title.x = element_text(vjust = -0.75),
                   axis.title.y = element_text(vjust = 1.5),
                   legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                   legend.key = element_blank(),
                   strip.background=element_rect(fill="transparent", color="transparent"),
                   legend.position = "top")




#### SOURCE FUNCTIONS ##########################################################
#source(file.path(dir_input, "strains_analysis", "R_scripts", "snv_linkage_functions.R"))
# Returns functions: analyse_linkage_data(), linkage_avg_plots(), filter_snv_freq_window(), make_multi_panel_fig(), summarize_log_files()



## Functions to calculate
calc_dprime <- function(df){
  
    mat <- as.matrix(df[ ,c('total', 'count_AB', 'count_Ab', 'count_aB', 'count_ab', 'linkD')])
    
    get_dprime <- function(row){
    p1 = as.numeric((row['count_AB'] + row['count_Ab']) / row['total'])
    p2 = as.numeric(1 - p1)
    q1 = as.numeric((row['count_AB'] + row['count_aB']) / row['total'])
    q2 = as.numeric(1 - q1)

    if(row['linkD'] < 0) {
      denom = max(c((-p1*q1),(-p2*q2)))
    } else if (row['linkD'] > 0) {
      denom = min(c((p1*q2), (p2*q1)))
    } else {
      return(0)
    }
    return(round(as.numeric(row['linkD']) / denom,3))
    }
    
    dprimes <- apply(mat, 1, FUN= get_dprime)
    return(dprimes)
}

# This is used in make_linkage_df()
calc_haplotypes <- function(df){
  counts <-  as.matrix(df[, c('count_AB','count_Ab','count_aB','count_ab')])
  counts[which(counts != 0)] <-  1
  
  df$haplotype <- ifelse(rowSums(counts) == 4, "h4", 
                    ifelse(rowSums(counts) == 3, "h3",
                           ifelse(rowSums(counts) == 2, "h2",
                                  ifelse(rowSums(counts) == 1, "h1", "error"))))
  return(df)
  }

## Extract linkage files and make a dataframe
make_linkage_df <- function(file, species, depth_min= 20){
  require(tidyverse)
  ## Get all linkage files into a list
  message("Reading in .linkage files")
  
  ## sp.linkage = dataframe of .linkage files in a folder
  sp.linkage <- suppressMessages(read_tsv(file) %>% 
    mutate(filepath= str_c(species, "/", filename)))
  
  sp.list <- apply(sp.linkage, 1, FUN= function(x) read_tsv(x['filepath']))
  names(sp.list) <- str_replace(sp.linkage$filename, ".linkage", "")
  
  ## Add haplotype column to each sample
  message("Calculating haplotype frequencies")
  sp.list <- lapply(sp.list, calc_haplotypes)
  
  sp.df <- map2(sp.list, names(sp.list), function(df, df_names) mutate(df, sample= df_names)) %>% 
    do.call(rbind, .) %>% 
    filter(total >= depth_min)
  return(sp.df)

}
sp1.df <- make_linkage_df(file= "species_1_linkage_files.tsv", species= "species_1")
sp2.df <- make_linkage_df(file= "species_2_linkage_files.tsv", species= "species_2", depth_min = 10)
sp3.df <- make_linkage_df(file= "species_3_linkage_files.tsv", species= "species_3")

## Calculate average r2 at each distance for each haplotype
sp1.df.means <- sp1.df %>% 
  group_by(sample, Distance, haplotype) %>% 
  summarize(n= length(r2),
            r2_mean = round(mean(r2, na.rm= TRUE), 4))

sp2.df.means <- sp2.df %>% 
  group_by(sample, Distance, haplotype) %>% 
  summarize(n= length(r2),
            r2_mean = round(mean(r2, na.rm= TRUE), 4))

sp3.df.means <- sp3.df %>% 
  group_by(sample, Distance, haplotype) %>% 
  summarize(n= length(r2),
            r2_mean = round(mean(r2, na.rm= TRUE), 4))


## Get linkage data frames into a list
get_linkage_list <- function(df, haplo){
  samps<- unique(df$sample)
  samp.list <- list()
  
  for(i in samps){
    samp.list[[i]] <- df %>% 
      filter(sample == i & haplotype == haplo) %>% 
      filter(Distance > 0)
  }
  return(samp.list)
}

sp1.list <- get_linkage_list(df= sp1.df.means, haplo= "h4")
sp3.list <- get_linkage_list(df= sp3.df.means, haplo= "h4")
  
# Function to fit NLS    
fit.nls <- function(df){
  require(nls.multstart)
  
  message(paste("Try NLS fit: ", df$sample[1]))
  # try({fit <- nls(r2_mean ~ SSasymp(log(Distance), Asym, R0, lrc), data = df, trace= TRUE)
  # fit.df <- broom::tidy(fit) %>%
  #   select(term, estimate, std.error)
  # })
  # 
  
  ## nls.multstart: https://cran.r-project.org/web/packages/nls.multstart/README.html
    ## used this function because regular nls was giving me errors and failing to find a regression
  # Define the asymptotic regression function (same as used in SSasym() function)
    # Asym= lower asymptote of the curve
    # R0= x intercept value
    # lrc= "curviness parameter"
  
  asym_regression <- function(input, Asym, R0, lrc){
    Y <- Asym + (R0 - Asym) * exp(-exp(lrc)* input)
    return(Y)  
  }
  
  ## Fit the model with nls_multstart
     # To determine starting value ranges, I looked at all the models that fit correctly with regular nls()
     # I used the min and max of these values to bound the starting range
     # nls_multstart is run with gridstart approach (see link above for more info)
  try({fit <- nls_multstart(r2_mean ~ asym_regression(input= Distance, Asym, R0, lrc),
                       data = df,
                       iter = rep(5, 3),
                       start_lower = c(Asym= -15, R0= 0.5, lrc= -5),
                       start_upper = c(Asym= 1, R0= 1, lrc= -0.1),
                       lower= c(Asym= 0, R0= -Inf, lrc= -Inf),
                       upper= c(Asym= 1, R0= 1, lrc= Inf),
                       supp_errors = 'Y',
                       na.action = na.omit)
      fit.df <- broom::tidy(fit) %>%
        select(term, estimate, std.error)
  })
  try(print(fit.df))
  
  if(exists("fit") == TRUE){
    predicted.values <- data.frame(Distance= df$Distance, n= df$n, r2_mean = df$r2_mean, r2_predict= predict(fit))
    fit.list <- list(nls.params= fit.df, nls.predicts= predicted.values )
    rm(fit)
    return(fit.list)
  } else {
    return(list(params= NA, predicts= NA))
  }
}

nls.sp1 <- lapply(sp1.list, fit.nls)
nls.sp3 <- lapply(sp3.list, fit.nls)

# Check to see how many samples failed to fit an nls curve
  #lapply(nls.sp1, is.na)
  #lapply(nls.sp3, is.na)

# Function to plot the NLS curves for each sample
nls.plots <- function(list.element, samp.name, species){
  if(is.na(list.element[[2]]) == FALSE){
    
    ## Extract parameter values
    ests <- as.character(round(list.element[[1]]$estimate, 3))
    anno.params <- str_c("Asym= ", ests[1], "  R0= ", ests[2], "  lrc= ", ests[3], sep="")
    
    ## Make plot
    p1 <-  ggplot(data= list.element[[2]]) +
      geom_point(aes(x= Distance, y= r2_mean, color= n), size= 2) +
      geom_line(aes(x= Distance, y= r2_predict), color= "red", size= 1) +
      annotate("text", label= anno.params, x= 1.5, y= 0.99, hjust= 0, size= 4) +
      scale_y_continuous(limits= c(0, 1), expand= c(0.02, 0)) +
      scale_x_continuous(limits= c(0, 450), expand= c(0.01, 0)) +
      labs(title= str_c(samp.name, species, sep="-")) +
      theme_bw()
    ggsave(p1, filename= str_c(samp.name, ".pdf"), path= file.path("Output_figures", str_c("linkage_", species)), height= 4, width= 6, units= "in", device= NULL)
    #ggsave(p1, filename= str_c(samp.name, ".pdf"), path= file.path("Output_figures", "linkage_sp1"), height= 4, width= 6, units= "in", device= NULL)
    
      }
}

## CHECK nls.plots as I added some species information in the ggsave and title command that needs to be made generalizable
#pmap(list(list.element= nls.sp1, samp.name= names(nls.sp1), species= "sp1"), nls.plots) # pmap is the purrr version of lapply, and lets you iterate over multiple arguments in a function
#pmap(list(list.element= nls.sp3, samp.name= names(nls.sp3), species= "sp3"), nls.plots) 


# Function to extract the predictor values from the NLS list
nls.preds <- function(list){
  params <- lapply(list, '[[', 2)
  # Add sample column
  preds <- map2(params, names(params), function(df_data, df_names) try(mutate(df_data, sample= df_names), silent= TRUE))
  
  # Get vector of TRUE/FALSE for nls models that worked
  worked <- do.call(rbind,lapply(preds, function(x) class(x) != "try-error"))[, 1]
  
  # Extract worked=TRUE from list and return a data frame
  preds.df <- do.call(rbind, preds[worked]) 
  return(preds.df)
}
nls.preds.sp1 <- nls.preds(nls.sp1) %>% mutate(species= "sp1")
nls.preds.sp3 <- nls.preds(nls.sp3) %>% mutate(species= "sp3")


# Function to extract the parameter values from the NLS list
nls.params <- function(list){
  params <- lapply(list, '[[', 1)
  # Add sample column
  params <- map2(params, names(params), function(df_data, df_names) try(mutate(df_data, sample= df_names), silent= TRUE))
  
  # Get vector of TRUE/FALSE for nls models that worked
  worked <- do.call(rbind,lapply(params, function(x) class(x) != "try-error"))[, 1]
  
  # Extract worked=TRUE from list and return a data frame
  params.df<-do.call(rbind, params[worked]) 
  return(params.df)
}
nls.params.sp1 <- nls.params(nls.sp1) %>% mutate(species= "sp1")
nls.params.sp3 <- nls.params(nls.sp3)  %>% mutate(species= "sp3")

## Combine species 1 and 3 parameters
nls.params.comb <- rbind(nls.params.sp1, nls.params.sp3)

nls.params.comb.w <- nls.params.comb %>% 
  select(-std.error) %>% 
  spread(term, estimate)

## Calculate summary statistics on the parameters
nls.params.sum <- nls.params.comb %>% 
  group_by(term) %>% 
  summarize(
    min= min(estimate),
    mean= mean(estimate),
    sd= sd(estimate),
    med= median(estimate),
    max= max(estimate)
  )
#write_tsv(nls.params.comb, path= "Output_tables/nls_parameters.tsv")

## Calc haplotypes frequency
  # map2 is the purrr version of lapply, and lets you iterate over two variables
  # This function takes each element in sp1.list calculates the frequency of each haplotype, and adds a column for sample name
  # map_dfr takes  list and returns a dataframe

sp1.haplo.freq <- sp1.df.means %>% 
  group_by(sample, haplotype) %>% 
  summarize(n= length(haplotype)) %>% 
  mutate(freq= n/sum(n))


sp2.haplo.freq <- sp2.df.means %>% 
  group_by(sample, haplotype) %>% 
  summarize(n= length(haplotype)) %>% 
  mutate(freq= n/sum(n))

sp3.haplo.freq <- sp3.df.means %>% 
  group_by(sample, haplotype) %>% 
  summarize(n= length(haplotype)) %>% 
  mutate(freq= n/sum(n))



#### PLOTS #########################################################################



## ggplot theme for snv linkage
theme_snv <- theme(panel.grid = element_blank(),
                   plot.margin = unit(c(1, 1, 1, 1), "cm"),
                   text = element_text(size= 14),
                   plot.background = element_rect(fill = "transparent"), # bg of the plot
                   panel.background = element_rect(fill= "transparent", color="black"),
                   axis.text = element_text(colour="black"),
                   axis.title.x = element_text(vjust = -0.75),
                   axis.title.y = element_text(vjust = 1.5),
                   legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                   legend.key = element_blank(),
                   strip.background=element_rect(fill="transparent", color="transparent"),
                   legend.position = "top")


## Species 1 h4
ggplot(data= nls.preds.sp1) +
  geom_point(aes(x= Distance, y= r2_mean), size= 2) +
  geom_line(aes(x= Distance, y= r2_predict), color= "red", size= 1) +
  #annotate("text", label= anno.params, x= 1.5, y= 0.99, hjust= 0, size= 2) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.02, 0)) +
  scale_x_continuous(limits= c(0, 440), expand= c(0.01, 0)) +
  labs(x= "Distance (bp)", y= expression("mean"~r^2), title= "Species 1 - h4") +
  facet_wrap(sample ~ ., nrow= 7) +
  theme_snv
#ggsave(last_plot(), filename= "linkage_sp1_h4.pdf", path= "Output_figures", height= 8, width= 10.5, units= "in", device= cairo_pdf)
ggsave(last_plot(), filename= "linkage_sp1_h4.jpg", path= "Output_figures", width= 6, height= 4, units= "in", dpi= 320)



## Species 3 h4
ggplot(data= nls.preds.sp3) +
  geom_point(aes(x= Distance, y= r2_mean), size= 2) +
  geom_line(aes(x= Distance, y= r2_predict), color= "red", size= 1) +
  #annotate("text", label= anno.params, x= 1.5, y= 0.99, hjust= 0, size= 2) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.02, 0)) +
  scale_x_continuous(limits= c(0, 440), expand= c(0.01, 0)) +
  labs(x= "Distance (bp)", y= expression("mean"~r^2), title= "Species 3 - h4") +
  facet_wrap(sample ~ ., nrow= 7) +
  theme_snv
#ggsave(last_plot(), filename= "linkage_sp3_h4.pdf", path= "Output_figures", height= 8, width= 10.5, units= "in", device= cairo_pdf)
ggsave(last_plot(), filename= "linkage_sp3_h4.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)



## Histograms of NLS parameters
ggplot(data= filter(nls.params.comb, term== "Asym")) +
  geom_histogram(aes(x= estimate), fill= "black", binwidth = 0.01) +
  labs(x= "Asym", y= "Count") +
  scale_x_continuous(limits= c(-0.005, 0.4),  expand= c(0.01, 0)) +
  scale_y_continuous(limits= c(0, 10), breaks= seq(0, 10, by= 1), labels= c("0", "", "2", "", "4", "", "6", "", "8", "", "10"),  expand= c(0, 0)) +
  theme_snv

ggplot(data= filter(nls.params.comb, term== "lrc")) +
  geom_histogram(aes(x= estimate), fill= "black", binwidth = 0.1) +
  labs(x= "lrc", y= "Count") +
  #scale_x_continuous(limits= c(-0.005, 0.4),  expand= c(0.01, 0)) +
  #scale_y_continuous(limits= c(0, 10), breaks= seq(0, 10, by= 1), labels= c("0", "", "2", "", "4", "", "6", "", "8", "", "10"),  expand= c(0, 0)) +
  theme_snv


ggplot(data= nls.params.comb.w) +
  geom_point(aes(x= lrc, y= Asym)) +
  labs(x= "lrc", y= "Asym") +
  #scale_x_continuous(limits= c(-0.005, 0.4),  expand= c(0.01, 0)) +
  #scale_y_continuous(limits= c(0, 10), breaks= seq(0, 10, by= 1), labels= c("0", "", "2", "", "4", "", "6", "", "8", "", "10"),  expand= c(0, 0)) +
  theme_snv



# sp1.df %>% 
#   filter(haplotype == "h4") %>% 
#   ggplot(., aes(x= Distance, y= r2)) +
#   geom_point() +
#   labs(title= "Haplotype= 4")
#   facet_wrap(~sample, nrow= 9) +
#   theme_snv
# 
#   
#   sp1.df.means %>% 
#     filter(haplotype == "h4") %>% 
#     ggplot(., aes(x= Distance, y= r2_mean)) +
#     geom_point() +
#     stat_smooth(method= "lm", se= FALSE, color= "tomato") +
#     labs(x= "Distance (bp)", y= expression(mean~r^2),  title= "Species_1 h4") +
#     scale_y_continuous(limits= c(0, 1)) +
#    # scale_y_log10(limits= c(0.001, 1)) +
#     #scale_x_log10() +
#     #annotation_logticks() +
#     facet_wrap(~sample, nrow= 9) +
#     theme_snv
#   #ggsave(last_plot(), file= "linkage_sp1_h4.pdf", width= 10, height= 10, units= "in", path= dir_out_fig, device= cairo_pdf)
#   
#   
#   sp1.df.means %>% 
#     filter(haplotype == "h3" | haplotype == "h4") %>% 
#     ggplot(., aes(x= Distance, y= r2_mean)) +
#     geom_point() +
#     stat_smooth(method= "lm", se= FALSE, color= "tomato") +
#    # scale_x_log10() +
#   #  scale_y_log10() +
#    # annotation_logticks() +
#     scale_y_continuous(limits= c(0, 1)) +
#     labs(x= "Distance (bp)", y= expression(mean~r^2),  title= "Species_1 h3 & h4") +
#     facet_wrap(~sample, nrow= 9) +
#     theme_snv
#   ggsave(last_plot(), file= "linkage_sp1_h3_h4.pdf", width= 10, height= 10, units= "in", path= dir_out_fig, device= cairo_pdf)
#   
#   sp3.df.means %>% 
#     filter(haplotype == "h3" | haplotype == "h4") %>% 
#     ggplot(., aes(x= Distance, y= r2_mean)) +
#     geom_point() +
#     stat_smooth(method= "lm", se= FALSE, color= "tomato") +
#     # scale_x_log10() +
#     #  scale_y_log10() +
#     # annotation_logticks() +
#     scale_y_continuous(limits= c(0, 1)) +
#     labs(x= "Distance (bp)", y= expression(mean~r^2),  title= "Species_3 h3 & h4") +
#     facet_wrap(~sample, nrow= 9) +
#     theme_snv
#   ggsave(last_plot(), file= "linkage_sp3_h3_h4.pdf", width= 10, height= 10, units= "in", path= dir_out_fig, device= cairo_pdf)
#   




### Make Haploytype frequency plots
ggplot(data= sp1.haplo.freq) +
  geom_bar( aes(x=sample, y=freq, fill= haplotype), stat= "identity") +
  coord_flip() +
  labs(title= "Species 1") +
  theme_snv
#ggsave(last_plot(), file= "sp1_haplotype_freqs.pdf", width= 10, height= 10, units= "in", path= dir_out_fig, device= cairo_pdf)
ggsave(last_plot(), filename= "sp1_haplotype_freqs.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)


ggplot(data= sp2.haplo.freq) +
  geom_bar( aes(x=sample, y=freq, fill= haplotype), stat= "identity") +
  coord_flip() +
  labs(title= "Species 2") +
  theme_snv
#ggsave(last_plot(), file= "sp2_haplotype_freqs.pdf", width= 10, height= 10, units= "in", path= dir_out_fig, device= cairo_pdf)
ggsave(last_plot(), filename= "sp2_haplotype_freqs.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)


ggplot(data= sp3.haplo.freq) +
  geom_bar( aes(x=sample, y=freq, fill= haplotype), stat= "identity") +
  coord_flip() +
  labs(title= "Species 3") +
  theme_snv
#ggsave(last_plot(), file= "sp3_haplotype_freqs.pdf", width= 10, height= 10, units= "in", path= dir_out_fig, device= cairo_pdf)
ggsave(last_plot(), filename= "sp3_haplotype_freqs.jpg", path= "Output_figures", width= 8, height= 6, units= "in", dpi= 320)

