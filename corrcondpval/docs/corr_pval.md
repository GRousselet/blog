Correlation estimation: reporting conditional on p values
================
Guillaume A. Rousselet
2018-06-13

-   [Dependencies](#dependencies)
-   [Rho = 0](#rho-0)
    -   [Parameters](#parameters)
    -   [Simulation](#simulation)
    -   [Make KDE](#make-kde)
    -   [Plot kernel density estimates](#plot-kernel-density-estimates)
-   [Rho = 0.4](#rho-0.4)
    -   [Parameters](#parameters-1)
    -   [Simulation](#simulation-1)
    -   [Make KDE](#make-kde-1)
    -   [Plot kernel density estimates](#plot-kernel-density-estimates-1)

Dependencies
============

``` r
source('./functions/akerd.txt')
library(ggplot2)
library(tibble)
library(viridis)
```

    ## Loading required package: viridisLite

``` r
# https://www.r-bloggers.com/computing-the-mode-in-r/
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
```

Reporting correlation results conditional on statistical significance (p&lt;=0.05) leads to inflated effect size estimation in the literature. This bias estimation increases with lower sample sizes.

Rho = 0
=======

Parameters
----------

``` r
nseq <- c(seq(10, 100, 10), 150, 200, 300) # sample sizes
Nn <- length(nseq)
pts <- seq(-1, 1, 0.025) # KDE points
Np <- length(pts)
p <- 2
mu <- 0
nsim <- 50000
rho <- 0
sigma <- diag(p)
sigma[sigma==0] <- rho
```

Simulation
----------

``` r
set.seed(21)

res.r <- matrix(data = 0, nrow = nsim, ncol = Nn) # R values
res.p <- matrix(data = 0, nrow = nsim, ncol = Nn) # p values

for(iter.n in 1:Nn){
  
  print(paste0("Sample size = ", nseq[iter.n]))
  
  for(iter in 1:nsim){
    
    data <- MASS::mvrnorm(n = nseq[iter.n], mu = rep(mu, p), Sigma = sigma)
    out <- cor.test(data[,1],data[,2], method = "pearson")
    res.r[iter,iter.n] <- out$estimate
    res.p[iter,iter.n] <- out$p.value

  }
}

hist(res.r)
mean(res.p<=0.05)

save(res.r, res.p,
     file = "./data/rpval.RData")
```

Make KDE
--------

``` r
# get data
load("./data/rpval.RData")

kde.r.ori <- matrix(data = 0, nrow = Np, ncol = Nn) 
kde.r.cond <- matrix(data = 0, nrow = Np, ncol = Nn) 

for(iter.n in 1:Nn){
  
  print(paste0("Sample size = ", nseq[iter.n]))
  
  kde.r.ori[,iter.n] <- akerd(res.r[,iter.n], pyhat=TRUE, pts=pts, plotit=FALSE)
  kde.r.cond[,iter.n] <- akerd(res.r[res.p[,iter.n]<=0.05,iter.n], pyhat=TRUE, pts=pts, plotit=FALSE)
  
}

save(kde.r.ori, kde.r.cond,
     file = "./data/rpval.kde.RData")
```

Plot kernel density estimates
-----------------------------

### All correlation estimates

``` r
# get data
load("./data/rpval.kde.RData")

# make data frame
fm <- array(0, dim = c(Np, Nn+1)) # make full matrix
fm[,1] <- pts
fm[,2:(Nn+1)] <- kde.r.ori
colnames(fm) <- c("x",nseq)
df <- as_tibble(fm)
df <- tidyr::gather(df, SS, Density,2:(Nn+1))
df[[2]] <- as.character(df[[2]])
df[[2]] <- factor(df[[2]], levels=unique(df[[2]]))

# make plot
p <- ggplot(df, aes(x, Density)) + theme_classic() +
          geom_line(aes(colour = SS), size = 1)  + 
          scale_color_viridis(discrete = TRUE) + 
          theme(axis.title.x = element_text(size = 18),
                axis.text = element_text(size = 14, colour = "black"),
                axis.title.y = element_text(size = 18),
                legend.key.width = unit(1.5,"cm"),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                legend.position = "right",
                plot.title = element_text(size = 20, colour = "black"),
                panel.background = element_rect(fill="grey90")) +
          scale_x_continuous(limits = c(-1, 1), 
                             breaks = seq(-1, 1, 0.2)) +
  labs(x = "Correlation estimates", y = "Density") +
  ggtitle("Correlation sampling distributions") +
  guides(colour = guide_legend(override.aes = list(size=3), # make thicker legend lines
        title="Sample size")) # change legend title
p
```

![](corr_pval_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# save figure
ggsave(filename='./figures/figure_rpval_ori.png',width=9,height=5) 
```

### Correlation estimates p&lt;=0.05

``` r
# get data
load("./data/rpval.kde.RData")

# make data frame
fm <- array(0, dim = c(Np, Nn+1)) # make full matrix
fm[,1] <- pts
fm[,2:(Nn+1)] <- kde.r.cond
colnames(fm) <- c("x",nseq)
df <- as_tibble(fm)
df <- tidyr::gather(df, SS, Density,2:(Nn+1))
df[[2]] <- as.character(df[[2]])
df[[2]] <- factor(df[[2]], levels=unique(df[[2]]))

# make plot
p <- ggplot(df, aes(x, Density)) + theme_classic() +
          geom_line(aes(colour = SS), size = 1)  + 
          scale_color_viridis(discrete = TRUE) + 
          theme(axis.title.x = element_text(size = 18),
                axis.text = element_text(size = 14, colour = "black"),
                axis.title.y = element_text(size = 18),
                legend.key.width = unit(1.5,"cm"),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                legend.position = "right",
                plot.title = element_text(size = 20, colour = "black"),
                panel.background = element_rect(fill="grey90")) +
          scale_x_continuous(limits = c(-1, 1), 
                             breaks = seq(-1, 1, 0.2)) +
  labs(x = "Correlation estimates", y = "Density") +
  ggtitle("Correlation sampling distributions") +
  guides(colour = guide_legend(override.aes = list(size=3), # make thicker legend lines
        title="Sample size")) # change legend title
p
```

![](corr_pval_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
# save figure
ggsave(filename='./figures/figure_rpval_cond.png',width=9,height=5) 
```

Rho = 0.4
=========

Parameters
----------

``` r
nseq <- c(seq(10, 100, 10), 150, 200, 300) # sample sizes
Nn <- length(nseq)
pts <- seq(-1, 1, 0.025) # KDE points
Np <- length(pts)
p <- 2
mu <- 0
nsim <- 50000
rho <- 0.4
sigma <- diag(p)
sigma[sigma==0] <- rho
```

Simulation
----------

``` r
set.seed(21)

res.r <- matrix(data = 0, nrow = nsim, ncol = Nn) # R values
res.p <- matrix(data = 0, nrow = nsim, ncol = Nn) # p values

for(iter.n in 1:Nn){
  
  print(paste0("Sample size = ", nseq[iter.n]))
  
  for(iter in 1:nsim){
    
    data <- MASS::mvrnorm(n = nseq[iter.n], mu = rep(mu, p), Sigma = sigma)
    out <- cor.test(data[,1],data[,2], method = "pearson")
    res.r[iter,iter.n] <- out$estimate
    res.p[iter,iter.n] <- out$p.value

  }
}

hist(res.r)
mean(res.p<=0.05)

save(res.r, res.p,
     file = "./data/rpval_04.RData")
```

Make KDE
--------

``` r
# get data
load("./data/rpval_04.RData")

kde.r.ori <- matrix(data = 0, nrow = Np, ncol = Nn) 
kde.r.cond <- matrix(data = 0, nrow = Np, ncol = Nn) 

for(iter.n in 1:Nn){
  
  print(paste0("Sample size = ", nseq[iter.n]))
  
  kde.r.ori[,iter.n] <- akerd(res.r[,iter.n], pyhat=TRUE, pts=pts, plotit=FALSE)
  kde.r.cond[,iter.n] <- akerd(res.r[res.p[,iter.n]<=0.05,iter.n], pyhat=TRUE, pts=pts, plotit=FALSE)
  
}

save(kde.r.ori, kde.r.cond,
     file = "./data/rpval.kde_04.RData")
```

Plot kernel density estimates
-----------------------------

### All correlation estimates

``` r
# get data
load("./data/rpval.kde_04.RData")

# make data frame
fm <- array(0, dim = c(Np, Nn+1)) # make full matrix
fm[,1] <- pts
fm[,2:(Nn+1)] <- kde.r.ori
colnames(fm) <- c("x",nseq)
df <- as_tibble(fm)
df <- tidyr::gather(df, SS, Density,2:(Nn+1))
df[[2]] <- as.character(df[[2]])
df[[2]] <- factor(df[[2]], levels=unique(df[[2]]))

# make plot
p <- ggplot(df, aes(x, Density)) + theme_classic() +
          geom_line(aes(colour = SS), size = 1)  + 
          scale_color_viridis(discrete = TRUE) + 
          theme(axis.title.x = element_text(size = 18),
                axis.text = element_text(size = 14, colour = "black"),
                axis.title.y = element_text(size = 18),
                legend.key.width = unit(1.5,"cm"),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                legend.position = "right",
                plot.title = element_text(size = 20, colour = "black"),
                panel.background = element_rect(fill="grey90")) +
          scale_x_continuous(limits = c(-1, 1), 
                             breaks = seq(-1, 1, 0.2)) +
  labs(x = "Correlation estimates", y = "Density") +
  ggtitle("Correlation sampling distributions") +
  guides(colour = guide_legend(override.aes = list(size=3), # make thicker legend lines
        title="Sample size")) # change legend title
p
```

![](corr_pval_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# save figure
ggsave(filename='./figures/figure_rpval_ori_04.png',width=9,height=5) 
```

#### Mean and median correlations

``` r
load("./data/rpval_04.RData")
round(apply(res.r, 2, mean), digits = 3)
```

    ##  [1] 0.381 0.391 0.394 0.395 0.396 0.397 0.398 0.397 0.398 0.397 0.399
    ## [12] 0.400 0.400

``` r
round(apply(res.r, 2, median), digits = 3)
```

    ##  [1] 0.423 0.409 0.407 0.404 0.403 0.403 0.402 0.402 0.402 0.401 0.402
    ## [12] 0.401 0.400

``` r
round(apply(res.r, 2, Mode), digits = 3)
```

    ##  [1] 0.490 0.348 0.609 0.508 0.403 0.254 0.408 0.269 0.379 0.591 0.374
    ## [12] 0.368 0.462

### Correlation estimates p&lt;=0.05

``` r
# get data
load("./data/rpval.kde_04.RData")

# make data frame
fm <- array(0, dim = c(Np, Nn+1)) # make full matrix
fm[,1] <- pts
fm[,2:(Nn+1)] <- kde.r.cond
colnames(fm) <- c("x",nseq)
df <- as_tibble(fm)
df <- tidyr::gather(df, SS, Density,2:(Nn+1))
df[[2]] <- as.character(df[[2]])
df[[2]] <- factor(df[[2]], levels=unique(df[[2]]))

# make plot
p <- ggplot(df, aes(x, Density)) + theme_classic() +
          geom_line(aes(colour = SS), size = 1)  + 
          scale_color_viridis(discrete = TRUE) + 
          theme(axis.title.x = element_text(size = 18),
                axis.text = element_text(size = 14, colour = "black"),
                axis.title.y = element_text(size = 18),
                legend.key.width = unit(1.5,"cm"),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                legend.position = "right",
                plot.title = element_text(size = 20, colour = "black"),
                panel.background = element_rect(fill="grey90")) +
          scale_x_continuous(limits = c(-1, 1), 
                             breaks = seq(-1, 1, 0.2)) +
  labs(x = "Correlation estimates", y = "Density") +
  ggtitle("Correlation sampling distributions") +
  guides(colour = guide_legend(override.aes = list(size=3), # make thicker legend lines
        title="Sample size")) # change legend title
p
```

![](corr_pval_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
# save figure
ggsave(filename='./figures/figure_rpval_cond_04.png',width=9,height=5) 
```

#### Mean and median correlations

``` r
load("./data/rpval_04.RData")

mean.res <- vector(mode = "numeric", length = Nn) 
median.res <- vector(mode = "numeric", length = Nn) 
mode.res <- vector(mode = "numeric", length = Nn) 

for(iter.n in 1:Nn){

mean.res[iter.n] <- round(mean(res.r[res.p[,iter.n]<=0.05,iter.n]), digits=3)
median.res[iter.n] <- round(median(res.r[res.p[,iter.n]<=0.05,iter.n]), digits=3)
mode.res[iter.n] <- round(Mode(res.r[res.p[,iter.n]<=0.05,iter.n]), digits=3)
}

mean.res
```

    ##  [1] 0.723 0.567 0.496 0.457 0.434 0.420 0.412 0.406 0.404 0.401 0.399
    ## [12] 0.400 0.400

``` r
median.res
```

    ##  [1] 0.721 0.553 0.485 0.449 0.429 0.417 0.411 0.406 0.404 0.402 0.402
    ## [12] 0.401 0.400

``` r
mode.res
```

    ##  [1] 0.644 0.457 0.609 0.508 0.403 0.442 0.408 0.269 0.379 0.591 0.374
    ## [12] 0.368 0.462
