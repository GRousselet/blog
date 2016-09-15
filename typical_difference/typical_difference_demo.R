## typical_difference_demo

# Code to reproduce the figures in my blog post on the shift function
# Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

## ============================
## dependencies
# The text files should be in working directory of this script,
# or you can specify the path in each `source()` command:
source("Rallfun-v30.txt")
source("wilcox_modified.txt")
source("rgar_visualisation.txt")
source("rgar_utils.txt")
source("rgar_stats.txt")
# packages:
library("ggplot2")
library("cowplot")
library("tidyr")
## ============================

# --------------------------------
# Independent groups 
# --------------------------------

# make data: two skewed distributions
set.seed(1)
df <- 4 # Chi2 degrees of freedom
n <- 50 # sample size for both groups
g1 <- rnorm(n) + df
g1 <- g1 - hd(g1) + df + 1 # median centre + shift
g2 <- rchisq(n, df)
g2 <- g2 - hd(g2)  + df + 1 # median centre + shift

g2plot(g1,g2) # kernel density plot

t.test(g1,g2) # regular Welsh t-test

wilcox.test(g1,g2) # Mann-Whitney-Wilcoxon test (this test has nothing to do with Rand Wilcox!)

cidv2(g1,g2) # Cliff's delta test

ks(g1,g2) # Kolmogorov-Smirnov test
# ks(cont,expt,w=T) # uses a weighted version more sensitive to differences occuring in the tails

# compute shift function
set.seed(7)
out <- shifthd( g1, g2, nboot=200)

# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
sf

# make data frame
df <- mkdf2(g1,g2,group.label="Group")

# kernel density plots
kde <- plot.kde_rug_dec2(df)
kde

# strip charts (1D scatterplots)
set.seed(7)
scat <- plot.scat_dec(df,flip=TRUE) + 
        theme(axis.title.y = element_blank())  # Remove y-axis label
scat

# compute all pairwise differences + save data.frame
apd <- mkdf1(allpdiff(g1,g2))

# make kernel density plot
kde.diff <- plot.kde_rug_dec1(apd,fill.colour="#ffb347",fill.alpha=.9)
kde.diff

# strip chart
set.seed(7)
scat.diff <- plot.scat_dec(apd,flip=TRUE) + 
  theme(axis.title.y = element_blank(), # Remove y-axis label
        axis.text.y = element_blank())  # Remove name
scat.diff

# Wilcox's original difference asymmetry plot: q + (1-q) for several quantiles
# slow because of the bootstrap confidence intervals
set.seed(7)
dasym <- qwmwhd(g1,g2,q=seq(5,40,5)/100,plotit=TRUE,alpha=.05,nboot=1000)

# ggplot2 version
diff_asym <- plot.diff_asym(data=dasym$output)
diff_asym

# combine plots
plot_grid(scat, kde.diff, sf, diff_asym, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2, rel_heights = c(1, 1),label_size = 18,hjust = -1,scale=.95)

# save figure 1
ggsave(filename='typ_diff_fig1_ind.jpeg') #path=pathname

# ---------------------------------
# Example with linear shift effect:
# make data
set.seed(6)
g1 <- rnorm(n)
g2 <- rnorm(n)+1

ks(g1,g2) # Kolmogorov-Smirnov test

df <- mkdf2(g1,g2,group.label="Group") # make data frame

# strip charts (1D scatterplots)
set.seed(7)
scat <- plot.scat_dec(df,flip=TRUE) + 
  theme(axis.title.y = element_blank())  # Remove y-axis label

# shift function
set.seed(7)
out <- shifthd( g1, g2, nboot=200)
sf <- plot.sf(data=out) # function from rgar_visualisation.txt

# compute all pairwise differences + save data.frame
apd <- mkdf1(allpdiff(g1,g2))

# make kernel density plot
kde.diff <- plot.kde_rug_dec1(apd,fill.colour="#ffb347",fill.alpha=.9)

# difference distribution asymmetry figure
set.seed(7)
dasym <- qwmwhd(g1,g2,q=seq(5,40,5)/100,plotit=TRUE,alpha=.05,nboot=1000)
diff_asym <- plot.diff_asym(data=dasym$output)

# combine plots
plot_grid(scat, kde.diff, sf, diff_asym, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2, rel_heights = c(1, 1),label_size = 18,hjust = -1,scale=.95)

# save figure 2
ggsave(filename='typ_diff_fig2_ind_linear_effect.jpeg') #path=pathname

# ---------------------------------
# Example with no effect:
# make data
set.seed(7)
g1 <- rnorm(n)
g2 <- rnorm(n)

ks(g1,g2) # Kolmogorov-Smirnov test

df <- mkdf2(g1,g2,group.label="Group") # make data frame

# strip charts (1D scatterplots)
set.seed(7)
scat <- plot.scat_dec(df,flip=TRUE) + 
  theme(axis.title.y = element_blank())  # Remove y-axis label

# shift function
set.seed(7)
out <- shifthd( g1, g2, nboot=200)
sf <- plot.sf(data=out) # function from rgar_visualisation.txt

# compute all pairwise differences + save data.frame
apd <- mkdf1(allpdiff(g1,g2))

# make kernel density plot
kde.diff <- plot.kde_rug_dec1(apd,fill.colour="#ffb347",fill.alpha=.9)

# difference distribution asymmetry figure
set.seed(7)
dasym <- qwmwhd(g1,g2,q=seq(5,40,5)/100,plotit=TRUE,alpha=.05,nboot=1000)
diff_asym <- plot.diff_asym(data=dasym$output)

# combine plots
plot_grid(scat, kde.diff, sf, diff_asym, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2, rel_heights = c(1, 1),label_size = 18,hjust = -1,scale=.95)

# save figure 3
ggsave(filename='typ_diff_fig3_ind_no_effect.jpeg') #path=pathname

# ---------------------------------
# Another example with no effect:
# make data
set.seed(3) #set.seed(7)
g1 <- rnorm(n*2)
g2 <- rnorm(n*2)

ks(g1,g2) # Kolmogorov-Smirnov test

df <- mkdf2(g1,g2,group.label="Group") # make data frame

# kernel density plots
kde <- plot.kde_rug_dec2(df)

# shift function
set.seed(7)
out <- shifthd( g1, g2, nboot=200)
sf <- plot.sf(data=out) # function from rgar_visualisation.txt

# compute all pairwise differences + save data.frame
apd <- mkdf1(allpdiff(g1,g2))

# make kernel density plot
kde.diff <- plot.kde_rug_dec1(apd,fill.colour="#ffb347",fill.alpha=.9)

# difference distribution asymmetry figure
set.seed(7)
dasym <- qwmwhd(g1,g2,q=seq(5,40,5)/100,plotit=TRUE,alpha=.05,nboot=1000)
diff_asym <- plot.diff_asym(data=dasym$output)

# combine plots
plot_grid(kde, kde.diff, sf, diff_asym, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)

# save figure 4
ggsave(filename='typ_diff_fig4_ind_no_effect2.jpeg') #path=pathname

# --------------------------------
# Dependent groups 
# --------------------------------

# load data
load("paired_example.RData") # pdata
df <- pdata

# --------------------------------
# Make figure 5: data descriptions

# independent stripcharts (1D scatterplots)
# make long format
df.long <- gather(df,gr,data,condition1,condition2)
df.long$participant <- as.factor(df.long$participant)
df.long$gr <- as.factor(df.long$gr)
set.seed(7)
strip <- plot.scat_dec(df.long,flip=FALSE) + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16))   # Remove y-axis label
strip

t.test(df$condition1,df$condition2,paired=TRUE)

# linked stripcharts -----------------
linkedstrip <- ggplot(df.long, aes(x=gr, y=data, group=participant)) +
  geom_line(aes(colour=participant),size=1, alpha=.5,
            position=position_dodge(width=0.4)) +
  geom_point(aes(fill=participant,colour=gr),shape = 21,size = 4, stroke = 1,
             position=position_dodge(width=0.4),alpha=.5,colour="grey5") +
  theme_grey() +
  theme(axis.text.x = element_text(colour="grey20",size=16),
        axis.text.y = element_text(colour="grey20",size=16),  
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(colour="grey20",size=20),
        legend.position="none") +
  labs(title="Paired observations") +
  scale_y_continuous(limits=c(0, 27),breaks=seq(0,25,5)) 
linkedstrip

# scatterplot of paired observations -----------------
scatterdiff <- scat_paired_obs(df=df,xname="condition1",yname="condition2",
                               min.x=5,min.y=5,max.x=27,max.y=27,axis.steps=2)
scatterdiff

# 1D scatterplots = stripcharts of differences -------
paired_differences <- df.long$data[df.long$gr=="condition2"]-df.long$data[df.long$gr=="condition1"]
diff <- mkdf1(paired_differences)
set.seed(8)
cc <- "grey80" # colour to plot deciles
p <- ggplot(diff, aes(x=gr,y=data,fill=gr,colour=gr,shape=gr)) +
  geom_abline(intercept = 0, slope = 0)
p <- stat_summary_deciles(p,cc="grey80",width=0.5,size.d=1,size.md=2)
diffstrip <- p + geom_jitter(position=position_jitter(0.3), size=4, stroke=1, alpha=.5,shape=21,
              colour="black",fill="#ffb347") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="grey20",size=16),  
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(colour="grey20",size=20)) +
  labs(title="Differences: condition1 - condition2") +
  scale_y_continuous(limits=c(-2, 7),breaks=seq(-2,7,1)) 
diffstrip

# combine plots
plot_grid(strip, linkedstrip, scatterdiff, diffstrip, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2, rel_heights = c(1, 1),label_size = 18,hjust = -1,scale=.95)

# save figure 5
ggsave(filename='typ_diff_fig5_dep1.jpeg') #path=pathname

# ----------------------------------------------------------
# Make figure 6: descriptions + full distribution inferences

# compute shift function
set.seed(7)
g2 <- df.long$data[df.long$gr=="condition2"]
g1 <- df.long$data[df.long$gr=="condition1"]
out <- shiftdhd( g1, g2, nboot=200)

# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
sf

# asymmetry plot
set.seed(7)
dasym <- difQpci(paired_differences,q=seq(5,40,5)/100,plotit=TRUE,alpha=.05,nboot=1000,SEED=TRUE,LINE=TRUE)

# ggplot2 version
diff_asym <- plot.diff_asym(data=dasym$output)
diff_asym

# combine plots
strip.s <- strip + coord_flip()
plot_grid(strip.s, diffstrip, sf, diff_asym, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2, rel_heights = c(1, 1),label_size = 18,hjust = -1,scale=.95)

# save figure 6
ggsave(filename='typ_diff_fig6_dep2.jpeg') #path=pathname

# ----------------------------------------
# Figure 7: deciles + confidence intervals

# deciles + confidence intervals
set.seed(3)
out <- quantiles_pbci(paired_differences,q=seq(1,9)/10,nboot=2000,alpha=0.05)

# decile plot
decile_plot <- plot.deciles(out=out,plotzero=TRUE,xtitle="Differences",hjust=-.05,vjust=.2,size=5)
decile_plot

# save figure 7
ggsave(filename='typ_diff_fig7_dep3_decile_plot.jpeg',width=5,height=5) #path=pathname

# -----------------------------------
# figure 8: asymmetric differences

# increase magnitude of 4 largest differences
spd <- sort(paired_differences)
spd[seq(length(spd)-3,length(spd))] <- spd[seq(length(spd)-3,length(spd))]*1.5

# stripchart of differences
set.seed(8)
cc <- "grey80" # colour to plot deciles
p <- ggplot(mkdf1(spd), aes(x=gr,y=data,fill=gr,colour=gr,shape=gr)) +
  geom_abline(intercept = 0, slope = 0)
p <- stat_summary_deciles(p,cc="grey80",width=0.5,size.d=1,size.md=2)
diffstrip <- p + geom_jitter(position=position_jitter(0.3), size=4, stroke=1, alpha=.5,shape=21,
                             colour="black",fill="#ffb347") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="grey20",size=16),  
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(colour="grey20",size=20)) +
  labs(title="Differences: condition1 - condition2") +
  scale_y_continuous(limits=c(-2, 11),breaks=seq(-2,11,1)) 
diffstrip

# deciles + confidence intervals
set.seed(3)
out <- quantiles_pbci(spd,q=seq(1,9)/10,nboot=2000,alpha=0.05)

# decile plot
set.seed(3)
decile_plot <- plot.deciles(out=out,plotzero=TRUE,xtitle="Differences",hjust=-.1,vjust=.2,size=5)
decile_plot

# asymmetry plot
set.seed(7)
dasym <- difQpci(spd,q=seq(5,40,5)/100,plotit=FALSE,alpha=.05,nboot=1000)
diff_asym <- plot.diff_asym(data=dasym$output)
diff_asym

# combine plots
plot_grid(diffstrip, decile_plot, diff_asym, labels=c("A", "B", "C"), ncol = 3, nrow = 1, rel_widths = c(1,1,1),label_size = 18,hjust = -1,scale=.95)

# save figure 8
ggsave(filename='typ_diff_fig8_dep4_larger_diff.jpeg',width=15,height=5) #path=pathname
