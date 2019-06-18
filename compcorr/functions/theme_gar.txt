theme_gar <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, colour="black"),
        legend.key.width = unit(1.5,"cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text = element_text(size=20, face="bold", colour="white"),
        strip.background = element_rect(colour="black", fill="grey40"))