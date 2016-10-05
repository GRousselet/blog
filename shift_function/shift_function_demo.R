## shift_function_demo

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
# packages:
library("ggplot2")
library("cowplot")
## ============================

# --------------------------------
# example 1: difference in spread 
# --------------------------------
set.seed(21)
g1 <- rnorm(1000) + 6
g2 <- rnorm(1000)*1.5 + 6

ks(g1,g2) # Kolmogorov-Smirnov test
# ks(cont,expt,w=T) # uses a weighted version more sensitive to differences occuring in the tails
t.test(g1,g2) # regular Welsh t-test

# make data frame
df <- mkdf2(g1,g2)

# kernel density estimate + rug plot + superimposed deciles
kde <- plot.kde_rug_dec2(df)
kde

# compute shift function
out <- shifthd( g1, g2, nboot=200)

# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
sf

# combine KDE + SF
plot_grid(kde, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)

# save figure
ggsave(filename='shift_function_ex1_spread.jpeg') #path=pathname

# --------------------------------
# example 2: complete shift 
# --------------------------------
set.seed(21)
g1 <- rnorm(1000) + 6
g2 <- rnorm(1000) + 6.5

ks(g1,g2) # Kolmogorov-Smirnov test
t.test(g1,g2) # regular Welsh t-test

# make data frame
df <- mkdf2(g1,g2)
# kernel density estimate + rug plot + superimposed deciles
kde <- plot.kde_rug_dec2(df)
# compute shift function
out <- shifthd( g1, g2, nboot=200)
# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
# combine KDE + SF
plot_grid(kde, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex2_complete.jpeg') #path=pathname

# --------------------------------
# example 3: one sided shift (I) 
# --------------------------------
set.seed(21)
g1 <- rnorm(1000)
g1[g1>0] <- g1[g1>0]*2
g2 <- rnorm(1000)

ks(g1,g2) # Kolmogorov-Smirnov test
t.test(g1,g2) # regular Welsh t-test

# make data frame
df <- mkdf2(g1,g2)
# kernel density estimate + rug plot + superimposed deciles
kde <- plot.kde_rug_dec2(df)
# compute shift function
out <- shifthd( g1, g2, nboot=200)
# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
# combine KDE + SF
plot_grid(kde, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex3_onesided1.jpeg') #path=pathname

# -----------------------------------
# example 4: one sided shift (II)
# -----------------------------------
set.seed(4)
n = 300
g1 <- rgamma(n, shape = 7, scale = 2)
g2 <- rgamma(n, shape = 7, scale = 2.1)
md.g2 <- hd(g2) # median
# g2plot(g1,g2,op=4);abline(v=md.g2)
g2[g2>md.g2] <- sort(g2[g2>md.g2]) * seq(1,1.3,length.out=sum(g2>md.g2))
g1 <- g1*4 + 300
g2 <- g2*4 + 300

ks(g1,g2) # Kolmogorov-Smirnov test
t.test(g1,g2) # regular Welsh t-test

# make data frame
df <- mkdf2(g1,g2)
# kernel density estimate + rug plot + superimposed deciles
kde <- plot.kde_rug_dec2(df)
# compute shift function
out <- shifthd( g1, g2, nboot=200)
# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
# combine KDE + SF
plot_grid(kde, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex4_onesided2.jpeg') #path=pathname

# --------------------------------
# example 5: no change
# --------------------------------
set.seed(7)
g1 <- rnorm(1000)
g2 <- rnorm(1000)

ks(g1,g2) # Kolmogorov-Smirnov test
t.test(g1,g2) # regular Welsh t-test

# make data frame
df <- mkdf2(g1,g2)
# kernel density estimate + rug plot + superimposed deciles
kde <- plot.kde_rug_dec2(df) + 
        scale_x_continuous(breaks = seq(-3, 3, 1))
# compute shift function
out <- shifthd( g1, g2, nboot=200)
# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
# combine KDE + SF
plot_grid(kde, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex5_nochange.jpeg') #path=pathname

# -----------------------------------
# example 6: ozone data
# -----------------------------------
# data used in:
# Doksum, K.A. & Sievers, G.L. (1976)
# Plotting with Confidence - Graphical Comparisons of 2 Populations
# Biometrika, 63, 421-434
# control group:
g1 <- c(41,38.4,24.4,25.9,21.9,18.3,13.1,27.3,28.5,-16.9,26,17.4,21.8,15.4,27.4,19.2,22.4,17.7,26,29.4,21.4,26.6,22.7)
# ozone group:
g2 <- c(10.1,6.1,20.4,7.3,14.3,15.5,-9.9,6.8,28.2,17.9,-9,-12.9,14,6.6,12.1,15.7,39.9,-15.9,54.6,-14.7,44.1,-9)

ks(g1,g2) # Kolmogorov-Smirnov test
t.test(g1,g2) # regular Welsh t-test

# make data frame
df <- mkdf2(g1,g2)
# 1D scatterplots + superimposed deciles
scat <- plot.scat_dec(df,size=4,stroke=1) + 
      labs(list(x = "Rat groups", y = "Weight gains in grams")) +
  scale_x_discrete(breaks=c("Group 1", "Group 2"),
                   labels=c("Control group","Ozone group")) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
# scat  

# compute shift function
set.seed(4)
out <- shifthd( g1, g2, nboot=200)
# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
# combine scatterplots + SF
plot_grid(scat, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex6_ozone.jpeg') #path=pathname

# -----------------------------------
# example 7: Guinea pig survival data
# -----------------------------------
# data from:
# Bjerkedal, T. (1960)
# Acquisition of resistance in guinea pigs infected with different doses of virulent tubercle bacilli
# American Journal of Hygiene, 72(1):130?48
# Table 6, columns 1 & 2
# and used in:
# Doksum, K. (1974) 
# Empirical Probability Plots and Statistical Inference for Nonlinear Models in the two-Sample Case
# Annals of Statistics, 2, 267-277
g1 <- scan("bjerkedal_control_group.txt")
g2 <- scan("bjerkedal_treatment_group.txt")

ks(g1,g2) # Kolmogorov-Smirnov test
t.test(g1,g2) # regular Welsh t-test

# make data frame
df <- mkdf2(g1,g2)
# 1D scatterplots + superimposed deciles
scat <- plot.scat_dec(df,size=4,stroke=1) + 
  labs(list(x = "Rat groups", y = "Survival time in days")) +
  scale_x_discrete(breaks=c("Group 1", "Group 2"),
                   labels=c("Control group","Treatment group")) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
scat
# compute shift function
set.seed(4)
out <- shifthd( g1, g2, nboot=200)
# plot shift function
sf <- plot.sf(data=out) # function from rgar_visualisation.txt
# combine KDE + SF
plot_grid(scat, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex7_guineapigs.jpeg') #path=pathname

# ----------------------------------------
# Other options to compute shift functions
# ----------------------------------------

# The `shifthd` function uses a percentile bootstrap estimation of the 
# standard error of the quantiles to compute the difference 
# confidence intervals (Wilcox 1995).
# An alternative strategy is to use straightforward percentile bootstrap 
# confidence intervals (Wilcox et al. 2012, 2014).
# I have modified these newer functions from Wilcox, so they can be used
# with the ggplot2 function plot.pbsf() in `rgar_visualisation.txt`. 
# The modified functions qcomhd() and Dqcomhd() are in `wilcox_modified.txt`
# here <https://github.com/GRousselet/rstats>.

# `qcomhd` can be used to compare the deciles using the Harrell-Davis estimator,
# similarly to shifhd:
qcomhd(g1,g2,q=seq(.1,.9,.1),xlab="Control",ylab="Control-treatment");abline(0,0)
# `qcomhd` can be also be used to:
# - compare quantiles other the deciles,
qcomhd(g1,g2,q=seq(.1,.9,.05),xlab="Control",ylab="Control-treatment");abline(0,0)
# - consider alpha other than 0.05
qcomhd(g1,g2,q=c(.1,.25,.5,.75,.9),xlab="Control",ylab="Control-treatment",alpha=0.1);abline(0,0)
# Also, unlike shifthd, qcomhd performs well even when there are tied values.

## ==========================================================================  
## Dependent groups
## ==========================================================================

### Example with repeated measures and systematic shift between conditions:
set.seed(4)
g1 <- rnorm(100)
g2 <- x+.1*runif(100,0,1)

### to look at deciles only, with alpha = 0.05 only, and no tied values
shiftdhd(g1,g2)

### Dqcomhd - equivalent to qcomhd, but for dependent groups
Dqcomhd(g1,g2)

### compare deciles using single order statistics + percentile bootstrap
qdec2ci(g1,g2)

# -----------------------------------
# example 8: test-retest
# -----------------------------------
# ses1 contains ERP onsets from 120 participants
# we only use participants tested twice
x <- read.table("onset.txt") 
twoses <- scan("onset_twoses.txt")
ses1 <- x[twoses==1,3]
# ses2 contains onsets from 74 participants who also provided ses1 onsets
x <- read.table("onset2.txt")
ses2 <- x[,3]

# make data frame
df <- mkdf2(ses1,ses2)
# kernel density estimate + rug plot + superimposed deciles
kde <- plot.kde_rug_dec2(df) + scale_x_continuous(breaks = seq(50,200,25))
# compute shift function
set.seed(4)
out <- shiftdhd( ses1, ses2, nboot=200) 
# plot shift function
sf <- plot.sf(data=out) +
  scale_x_continuous(breaks = seq(70,120,10))
# combine KDE + SF
plot_grid(kde, sf, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1.5, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex8_onsets.jpeg') #path=pathname

# with KDE of pairwise differences ----------------------------------
diff <- ses1-ses2 # pairwise differences
df1 <- mkdf1(diff) # make data frame
kde1 <- plot.kde_rug_dec1(df1) + xlab("Paired differences") +
        scale_x_continuous(breaks = seq(-50,60,10))

# combine KDE + KDE(difference) + SF
plot_grid(kde, kde1, sf, labels=c("A", "B", "C"), ncol = 1, nrow = 3, rel_heights = c(1.5, 1, 1),label_size = 18,hjust = -1,scale=.95)
# save figure
ggsave(filename='shift_function_ex9_onsets_diff.jpeg') #path=pathname

# with violinplot of pairwise differences ----------------------------------
p <- ggplot(df1, aes(x=1, y=data))
p1 <- p + geom_violin(fill="grey90",colour="black",adjust=.5,size=.5) +
  geom_abline(intercept = 0, slope = 0,colour="black",size=1,linetype=2) + 
  stat_summary(fun.data="plot.q123",geom="pointrange", color = "grey20",size=1) +
  scale_y_continuous(breaks = seq(-50, 60, 10)) +
  labs(y = "Pairwise differences") +
  theme(axis.text.y = element_text(colour="grey20",size=16),
        axis.title.y = element_text(colour="grey20",size=18),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_continuous(limits=c(0.2, 1.8))
# p1

toprow <- plot_grid(kde, p1, labels=c("A", "B"), ncol = 2, nrow = 1, rel_widths = c(3, 1),label_size = 18,hjust = -1)
plot_grid(toprow, sf, labels = c('', 'C'), ncol = 1, nrow = 2, rel_heights = c(1.5, 1), scale = 1,label_size = 18,
          align='h',hjust = -1)
# save figure
ggsave(filename='shift_function_ex10_onsets_diff_violin.jpeg') #path=pathname

# with scatterplot of pairwise differences ----------------------------------
p <- ggplot(df1, aes(x=gr, y=data))
p1 <- p + geom_jitter(position=position_jitter(0.3),fill="grey90",colour="grey50", 
                      shape=21,size=3,stroke=1,alpha=0.5) +
  geom_abline(intercept = 0, slope = 0,colour="black",size=1,linetype=2) +
  stat_summary(fun.data="plot.q123",geom="pointrange", color = "black",size=1) +
  scale_y_continuous(breaks = seq(-50, 60, 10)) +
  labs(y = "Pairwise differences") +
  theme(axis.text.y = element_text(colour="grey20",size=16),
        axis.title.y = element_text(colour="grey20",size=18),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# add horizintal lines markign the quartiles instead:
# p1 + stat_summary(fun.data="plot.q1", geom="errorbar", colour="grey20", width=.4, size=1) +
#   stat_summary(fun.data="plot.q2", geom="errorbar", colour="grey20", width=.4, size=2) +
#   stat_summary(fun.data="plot.q3", geom="errorbar", colour="grey20", width=.4, size=1)
# p1

toprow <- plot_grid(kde, p1, labels=c("A", "B"), ncol = 2, nrow = 1, rel_widths = c(3, 1),label_size = 18,hjust = -1)
plot_grid(toprow, sf, labels = c('', 'C'), ncol = 1, nrow = 2, rel_heights = c(1.5, 1), scale = 1,label_size = 18,
          align='h',hjust = -1)
# save figure
ggsave(filename='shift_function_ex11_onsets_diff_scatter.jpeg') #path=pathname
# ggsave(filename='shift_function_ex11_onsets_diff_scatter.jpeg',width=10,height=10) #path=pathname
