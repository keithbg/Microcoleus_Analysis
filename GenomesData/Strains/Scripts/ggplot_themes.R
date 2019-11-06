## ggplot themes and plotting parameters for figures


theme_strains <- theme(panel.grid = element_blank(),
                       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                       text = element_text(size= 12),
                       plot.background = element_rect(fill = "transparent", color= "transparent"), # bg of the plot
                       panel.background = element_rect(fill= "transparent", color= "black"),
                       #panel.border= element_rect(fill= "transparent", color= "black", linetype= "solid", size= 0.5),
                       panel.ontop = TRUE,
                       axis.text = element_text(colour="black"),
                       axis.title.x = element_text(vjust = -0.75),
                       axis.title.y = element_text(vjust = 1.5),
                       legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                       legend.key = element_blank(),
                       strip.background = element_rect(fill="transparent", color= "transparent")
                       #axis.text.x = element_text(angle= 45, hjust= 1))
)

## Species colors
species.colors <- c("black", "tomato", "dodgerblue")
species.shapes <- c(21, 22, 24)
# scale_color_manual(values= c("black", "tomato", "dodgerblue"),
#                    labels= c("1", "2", "3"),
#                    name= "Species") +