# define quartile functions to use with ggplot2 stat_summary
# GAR, University of Glasgow, 2016-07-07
plot.q1 <- function(x){
  m <- hd(x,.25)
  c(y = m, ymin = m, ymax = m)
}
plot.q2 <- function(x){
  m <- hd(x,.5)
  c(y = m, ymin = m, ymax = m)
}
plot.q3 <- function(x){
  m <- hd(x,.75)
  c(y = m, ymin = m, ymax = m)
}

plot.q123 <- function(x){
c(y = hd(x,.5), ymin = hd(x,.25), ymax = hd(x,.75))
}

plot.quartiles <- function(p,colour="grey70", width=.5){
p <- p + stat_summary(fun.data="plot.q1", geom="errorbar", colour=colour, width=width, size=1) +
  stat_summary(fun.data="plot.q2", geom="errorbar", colour=colour, width=width, size=2) +
  stat_summary(fun.data="plot.q3", geom="errorbar", colour=colour, width=width, size=1)
p
}

# define decile functions to use with ggplot2 stat_summary
# GAR, University of Glasgow, 2016-07-09
plot.d1 <- function(x){
  m <- hd(x,.1)
  c(y = m, ymin = m, ymax = m)
}
plot.d2 <- function(x){
  m <- hd(x,.2)
  c(y = m, ymin = m, ymax = m)
}
plot.d3 <- function(x){
  m <- hd(x,.3)
  c(y = m, ymin = m, ymax = m)
}
plot.d4 <- function(x){
  m <- hd(x,.4)
  c(y = m, ymin = m, ymax = m)
}
plot.d5 <- function(x){
  m <- hd(x,.5)
  c(y = m, ymin = m, ymax = m)
}
plot.d6 <- function(x){
  m <- hd(x,.6)
  c(y = m, ymin = m, ymax = m)
}
plot.d7 <- function(x){
  m <- hd(x,.7)
  c(y = m, ymin = m, ymax = m)
}
plot.d8 <- function(x){
  m <- hd(x,.8)
  c(y = m, ymin = m, ymax = m)
}
plot.d9 <- function(x){
  m <- hd(x,.9)
  c(y = m, ymin = m, ymax = m)
}

annotate.quartiles <- function(p,data,x=0,hjust=0,vjust=0,size=10){
# add text annotation of quartiles using ggplot2 annotate
# GAR, University of Glasgow, 2016-07-07
hdq = c(hd(out,.25),hd(out,.5),hd(out,.75)) # compute quartiles
caption <- as.character(round(hdq, digits=2)) # turn into characters
p <- p + annotate("text", x = x, y = hdq[1], label = caption[1],
                  hjust = hjust, vjust = vjust, size = size) +
         annotate("text", x = x, y = hdq[2], label = caption[2],
                  hjust = hjust, vjust = vjust, size = size) +
         annotate("text", x = x, y = hdq[3], label = caption[3],
                  hjust = hjust, vjust = vjust, size = size)
p
}

plot.sf <- function(data=sf){
# Plot shift function using output from shifthd or shiftdhd
# GAR, University of Glasgow, 2016-07-09
ylim <- max(max(abs(data$ci_upper)),max(abs(data$ci_lower)))
ylim <- c(-ylim,ylim)
if("group1" %in% colnames(data)){
xintercept <- data$group1[5]
xplot = "group1"
xlabel = "Group 1 quantiles"
}
if("condition1" %in% colnames(data)){
xintercept <- data$condition1[5]
xplot = "condition1"
xlabel = "Condition 1 quantiles"
}
  p <- ggplot(data, aes_string(x=xplot, y="difference")) +
  geom_hline(yintercept=0,linetype=2,alpha=0.5) + # x=0 reference line
  geom_vline(xintercept=xintercept,linetype=2,alpha=0.5) + # y=median reference line
  geom_linerange(aes(ymin=ci_lower, ymax=ci_upper), colour="#009E73") +
  geom_line(colour="#009E73", linetype="solid", size=1.5) +
  geom_point(colour="#009E73", size=4, shape=21, fill="white") + #999999
  xlab(xlabel) +
  ylab("Quantile differences") +
  theme_grey() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold")) +
  scale_y_continuous(limits = ylim)
p
}

plot.kde_rug_dec2 <- function(data=df){
# Plot kernel density estimates + rug plots + superimposed deciles
# for 2 groups stored in a data frame
# GAR, University of Glasgow, 2016-07-09
require(plyr)
cdat <- ddply(data, "gr", summarise, deciles=deciles(data))
p <- ggplot(data, aes(x=data, fill=gr)) + geom_density(alpha=.3) +
  facet_grid(gr ~ .) +
  geom_vline(data=cdat, aes(xintercept=deciles,  colour=gr),
             linetype="solid", size=1, alpha=.5) +
  geom_rug() +
  theme(legend.position="none",
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold"),
        strip.text.y = element_text(size = 20, colour = "white"),
        strip.background = element_rect(colour="darkgrey", fill="darkgrey")) +
        ylab("Density")
p
}

plot.kde_rug_dec1 <- function(data=df){
# Plot kernel density estimate + rug plot + superimposed deciles
# for 1 group stored in a data frame
# GAR, University of Glasgow, 2016-07-11
require(plyr)
cdat <- ddply(data, "gr", summarise, deciles=deciles(data))
p <- ggplot(data, aes(x=data)) +
  geom_density(alpha=.3,fill="grey50",colour="black") +
  geom_vline(data=cdat, aes(xintercept=deciles),
             linetype="solid", size=1, alpha=.5, colour="black") +
  geom_rug() +
  theme(legend.position="none",
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold"),
        strip.text.y = element_text(size = 20, colour = "white"),
        strip.background = element_rect(colour="darkgrey", fill="darkgrey")) +
        ylab("Density")
p
}

plot.scat_dec <- function(data=df,size=4, stroke=1){
# Plot stripchart (1D scatterplot) + superimposed deciles
# for 2 groups stored in a data frame
# GAR, University of Glasgow, 2016-07-09
cc <- "grey80" # colour to plot deciles
p <- ggplot(data, aes(x=gr, y=data, fill=gr, colour=gr, shape=gr)) +
  stat_summary(fun.data="plot.d1", geom="errorbar", colour=cc, width=0.5, size=1) +
  stat_summary(fun.data="plot.d2", geom="errorbar", colour=cc, width=0.5, size=1) +
  stat_summary(fun.data="plot.d3", geom="errorbar", colour=cc, width=0.5, size=1) +
  stat_summary(fun.data="plot.d4", geom="errorbar", colour=cc, width=0.5, size=1) +
  stat_summary(fun.data="plot.d5", geom="errorbar", colour=cc, width=0.5, size=2) +
  stat_summary(fun.data="plot.d6", geom="errorbar", colour=cc, width=0.5, size=1) +
  stat_summary(fun.data="plot.d7", geom="errorbar", colour=cc, width=0.5, size=1) +
  stat_summary(fun.data="plot.d8", geom="errorbar", colour=cc, width=0.5, size=1) +
  stat_summary(fun.data="plot.d9", geom="errorbar", colour=cc, width=0.5, size=1) +
  #geom_dotplot(binaxis='y', stackdir='center', method="histodot", dotsize=dotsize,binwidth=binwidth, alpha=.5) +
  geom_jitter(position=position_jitter(0.3), size=size, stroke=stroke, alpha=.5) +
  theme_bw() +
  scale_colour_manual(values = c("grey5","grey5")) +
  scale_shape_manual(values=c(21,21)) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold"),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14)) +
  coord_flip()
p
}
