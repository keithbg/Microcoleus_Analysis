## ggplot themes and plotting parameters for figures

library(ggplot2)
library(cowplot)
library(ggpubr) # ggarrange
library(lemon) #facet_rep_wrap
library(RColorBrewer)
library(wesanderson)
library(ggsci)



theme_strains <- theme(panel.grid = element_blank(),
                       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                       text = element_text(size= 12),
                       plot.background = element_rect(fill = "transparent", color= "transparent"), # bg of the plot
                       panel.background = element_rect(fill= "transparent", color= "transparent"),
                       #panel.border= element_rect(fill= "transparent", color= "black", linetype= "solid", size= 0.5),
                       panel.border= element_rect(fill= NA, color= NA, linetype= "solid", size= 1),
                       panel.ontop = TRUE,
                       panel.spacing.y = unit(1, "lines"),
                       panel.spacing.x = unit(0.5, "lines"),
                       axis.line = element_line(color= "black", size= 0.25),
                       axis.text = element_text(colour="black"),
                       axis.title.x = element_text(vjust = -0.75),
                       axis.title.y = element_text(vjust = 1.5),
                       #legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                       legend.box.background  = element_rect(color= "transparent", fill= "transparent"),
                       legend.background = element_rect(color= "transparent", fill= "transparent"),
                       legend.key = element_blank(),
                       strip.background = element_rect(fill="transparent", color= "transparent")
                       #axis.text.x = element_text(angle= 45, hjust= 1))
)


PH2017_map_theme <- theme(plot.margin= margin(0.2, 0.2, 0.2, 0.2, unit= "cm"),
                          text= element_text(size= 10),
                          #panel.background = element_rect(fill= "#c7cee7"),
                          #panel.border = element_rect(color= "black", fill= NA),
                          plot.background = element_rect(fill = "transparent", color= "transparent"), # bg of the plot
                          #panel.background = element_rect(fill= "transparent", color= "transparent"),
                          panel.background = element_rect(fill= "#bbc2e0", color= "transparent"),
                          #panel.border= element_rect(fill= "transparent", color= "black", linetype= "solid", size= 0.5),
                          panel.border= element_rect(fill= NA, color= NA, linetype= "solid", size= 1),
                          legend.key= element_rect(fill= "transparent"),
                          legend.background = element_rect(fill= "transparent"),
                          #plot.background = element_rect(fill= "transparent", color= NA),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank()
)


theme_freshSci <- theme(panel.grid = element_blank(),
                        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                        text = element_text(size= 10),
                        plot.background = element_rect(fill = "transparent", color= "transparent"), # bg of the plot
                        panel.background = element_rect(fill= "transparent", color= "transparent"),
                        panel.border= element_rect(fill= NA, color= NA, linetype= "solid", size= 1),
                        panel.ontop = TRUE,
                        axis.line = element_line(color= "black", size= 0.25),
                        axis.text = element_text(colour="black"),
                        axis.title.x = element_text(vjust = -0.75),
                        axis.title.y = element_text(vjust = 1.5),
                        #strip.background = element_rect(fill="transparent", color= "transparent"),
                        strip.background = element_blank(),
                        legend.background = element_rect(color=NULL, fill= "transparent"),
                        legend.key = element_blank(),
                        legend.justification = c(0, 1),
                        legend.box.margin = margin(0, 0, 0, 0, unit= "cm"),
                        legend.text = element_text(size= 8),
                        legend.title = element_text(size= 8)
)


## Species colors
library(wesanderson)


species.colors <- ggsci::pal_npg("nrc")(3)[c(3, 1, 2)]
#species.colors <- ggsci::pal_jco("default")(3)#[c(3, 1, 2)]
#species.colors <- wes_palette("Royal1")[1:3]
#species.colors <- c("black", "tomato", "dodgerblue")
species.shapes <- c(21, 22, 24)

## Year colors
year.fill.colors <- pal_npg("nrc")(9)[c(1, 4, 7)]
year.colors <- pal_npg("nrc")(9)[c(1, 4, 7)]


# scale_color_manual(values= c("black", "tomato", "dodgerblue"),
#                    labels= c("1", "2", "3"),
#                    name= "Species") +


## Make alternating blank axis tick labels

axis.labels.tick <- function(breaks.vec, first_blank=FALSE){ 
  labels.vec <- as.character(breaks.vec)
  
  if(first_blank == FALSE){
    labels.tick <- rep(c(FALSE, TRUE), length(labels.vec)/2) 
    labels.vec[labels.tick] <- ""
  } else {
    labels.tick <- rep(c(TRUE, FALSE), length(labels.vec)/2) 
    labels.vec[labels.tick] <- ""
  }
  return(labels.vec)
}
