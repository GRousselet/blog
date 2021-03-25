Simulate correlations using g & h distributions
================
Guillaume A. Rousselet
2019-07-24

# Dependencies

``` r
library(tibble)
library(ggplot2)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

``` r
source("./functions/theme_gar.txt")
source("./functions/ghmult.txt")
```

The code to generate *g\&h* distributions is from Rand Wilcox’s
`Rallfun-v35.txt`, which is available in full
[here](http://dornsife.usc.edu/labs/rwilcox/software/). A subset of
functions was copied into the `ghmult.txt` file for convenience. The
function `ghmul` is used to generate bivariate normal distributions,
then each marginal distribution is transformed to a *g\&h* distribution.
After transformation, the correlation is not the same as the desired
one, and the modification varies as a function of *g*, *h*, and *rho*.
So here we use the function `rngh.sub` to determine the actual value of
rho needed to obtained the desired value after *g\&h* transformation. It
uses *n*=one million observations to check on the actual value of rho.

# Examples: rho = 0.5

## Get samples

``` r
rho.pop <- 0.5
gseq <- seq(0, 1, 0.1)
hseq <- seq(0, 0.4, 0.1)
gn <- length(gseq)
hn <- length(hseq)
n <- 200
m1 <- array(0, dim = c(n, 11, 5))
m2 <- array(0, dim = c(n, 11, 5))
for(G in 1:length(gseq)){
  for(H in 1:length(hseq)){
    g <- gseq[G]
    h <- hseq[H]
    rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
    cmat <- matrix(c(1,rho,rho,1),2,2)    
    out <- ghmul(n,g,h,p=2,cmat=cmat)
    m1[, G, H] <- out[,1]
    m2[, G, H] <- out[,2]
  }
}
```

    ## Loading required package: MASS

## Make figure

``` r
df <- tibble(x = as.vector(m1),
             y = as.vector(m2),
             g = factor(rep(rep(gseq, each = n), hn)),
             h = factor(rep(hseq, each = n*gn)))

ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(g), cols = vars(h)) +
  coord_cartesian(xlim = c(-5, 20), ylim = c(-5, 20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())
```

![](g_h_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# ggsave("./figures/gh_rho05_examples.pdf", width = 40, height = 30, units = "cm")
```

# Examples: rho = 0

## Get samples

``` r
rho.pop <- 0
gseq <- seq(0, 1, 0.1)
hseq <- seq(0, 0.4, 0.1)
gn <- length(gseq)
hn <- length(hseq)
n <- 200
m1 <- array(0, dim = c(n, 11, 5))
m2 <- array(0, dim = c(n, 11, 5))
for(G in 1:length(gseq)){
  for(H in 1:length(hseq)){
    g <- gseq[G]
    h <- hseq[H]
    rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
    cmat <- matrix(c(1,rho,rho,1),2,2)    
    out <- ghmul(n,g,h,p=2,cmat=cmat)
    m1[, G, H] <- out[,1]
    m2[, G, H] <- out[,2]
  }
}
```

## Make figure

``` r
df <- tibble(x = as.vector(m1),
             y = as.vector(m2),
             g = factor(rep(rep(gseq, each = n), hn)),
             h = factor(rep(hseq, each = n*gn)))

ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(g), cols = vars(h)) +
  coord_cartesian(xlim = c(-5, 20), ylim = c(-5, 20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())
```

![](g_h_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# ggsave("./figures/gh_rho0_examples.pdf", width = 40, height = 30, units = "cm")
```

# Examples: g = 0, h = 0

## Get samples

``` r
g <- 0 
h <- 0
rhoseq <- seq(0, 0.95, 0.05)
rn <- length(rhoseq)
n <- 200
m1 <- array(0, dim = c(n, rn))
m2 <- array(0, dim = c(n, rn))
for(R in 1:rn){
  rho <- rngh.sub(50,g,h,rhoseq[R])$rho.adjusted
  cmat <- matrix(c(1,rho,rho,1),2,2)    
  out <- ghmul(n,g=g,h=h,p=2,cmat=cmat)
  m1[, R] <- out[,1]
  m2[, R] <- out[,2]
}
```

## Make figure

``` r
df <- tibble(x = as.vector(m1),
             y = as.vector(m2),
             rho = factor(rep(rhoseq, each = n)))

ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(rho)) +
  coord_cartesian(xlim = c(-5, 20), ylim = c(-5, 20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())
```

![](g_h_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# ggsave("./figures/gh_g0h0_examples.pdf", width = 40, height = 30, units = "cm")
```

# Examples: g = 1, h = 0

## Get samples

``` r
g <- 1
h <- 0
rhoseq <- seq(0, 0.95, 0.05)
rn <- length(rhoseq)
n <- 200
m1 <- array(0, dim = c(n, rn))
m2 <- array(0, dim = c(n, rn))
for(R in 1:rn){
  rho <- rngh.sub(50,g,h,rhoseq[R])$rho.adjusted
  cmat <- matrix(c(1,rho,rho,1),2,2)
  out <- ghmul(n,g=g,h=h,p=2,cmat=cmat)
  m1[, R] <- out[,1]
  m2[, R] <- out[,2]
}
```

## Make figure

``` r
df <- tibble(x = as.vector(m1),
             y = as.vector(m2),
             rho = factor(rep(rhoseq, each = n)))

ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(rho)) +
  coord_cartesian(xlim = c(-5, 20), ylim = c(-5, 20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())
```

![](g_h_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# ggsave("./figures/gh_g1h0_examples.pdf", width = 40, height = 30, units = "cm")
```

# Examples: g = 0, h = 0.2

## Get samples

``` r
g <- 0 
h <- 0.2
rhoseq <- seq(0, 0.95, 0.05)
rn <- length(rhoseq)
n <- 200
m1 <- array(0, dim = c(n, rn))
m2 <- array(0, dim = c(n, rn))
for(R in 1:rn){
  rho <- rngh.sub(50,g,h,rhoseq[R])$rho.adjusted
  cmat <- matrix(c(1,rho,rho,1),2,2)
  out <- ghmul(n,g=g,h=h,p=2,cmat=cmat)
  m1[, R] <- out[,1]
  m2[, R] <- out[,2]
}
```

## Make figure

``` r
df <- tibble(x = as.vector(m1),
             y = as.vector(m2),
             rho = factor(rep(rhoseq, each = n)))

ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(rho)) +
  coord_cartesian(xlim = c(-5, 20), ylim = c(-5, 20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())
```

![](g_h_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# ggsave("./figures/gh_g0h02_examples.pdf", width = 40, height = 30, units = "cm")
```

# Examples: g = 1, h = 0.2

## Get samples

``` r
g <- 1
h <- 0.2
rhoseq <- seq(0, 0.95, 0.05)
rn <- length(rhoseq)
n <- 200
m1 <- array(0, dim = c(n, rn))
m2 <- array(0, dim = c(n, rn))
for(R in 1:rn){
  rho <- rngh.sub(50,g,h,rhoseq[R])$rho.adjusted
  cmat <- matrix(c(1,rho,rho,1),2,2)
  out <- ghmul(n,g=g,h=h,p=2,cmat=cmat)
  m1[, R] <- out[,1]
  m2[, R] <- out[,2]
}
```

## Make figure

``` r
df <- tibble(x = as.vector(m1),
             y = as.vector(m2),
             rho = factor(rep(rhoseq, each = n)))

ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(rho)) +
  coord_cartesian(xlim = c(-5, 20), ylim = c(-5, 20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())
```

![](g_h_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# ggsave("./figures/gh_g1h02_examples.pdf", width = 40, height = 30, units = "cm")
```

# Summary figure: g/h = 0/0.2

n = 5000 g = 0 / 0.2 h = 0 / 0.2

## rho = 0

### Get samples

``` r
set.seed(21)

n <- 5000
rho.pop <- 0

g <- 0
h <- 0
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x1 <- out[,1]
y1 <- out[,2]

g <- 0.2
h <- 0
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x2 <- out[,1]
y2 <- out[,2]

g <- 0
h <- 0.2
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x3 <- out[,1]
y3 <- out[,2]

g <- 0.2
h <- 0.2
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x4 <- out[,1]
y4 <- out[,2]

df <- tibble(x = c(x1, x2, x3, x4),
             y = c(y1, y2, y3, y4),
             g = factor(rep(c("g=0","g=0.2","g=0","g=0.2"), each = n)),
             h = factor(rep(c("h=0","h=0","h=0.2","h=0.2"), each = n)))
```

### Make figure

``` r
pA <- ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.1) +
  # geom_smooth(method = "lm") +
  facet_grid(rows = vars(h), cols = vars(g)) +
  coord_cartesian(xlim = c(-10, 25), ylim = c(-10, 25)) +
  theme(axis.title = element_blank())
  # theme(panel.grid.minor = element_blank())
pA
```

![](g_h_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# ggsave("./figures/gh_rho0.png", width = 20, height = 20, units = "cm")
```

## rho = 0.5

### Get samples

``` r
set.seed(21)

n <- 5000
rho.pop <- 0.5

g <- 0
h <- 0
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x1 <- out[,1]
y1 <- out[,2]

g <- 0.2
h <- 0
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x2 <- out[,1]
y2 <- out[,2]

g <- 0
h <- 0.2
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x3 <- out[,1]
y3 <- out[,2]

g <- 0.2
h <- 0.2
rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
cmat <- matrix(c(1,rho,rho,1),2,2)
out <- ghmul(n,g,h,p=2,cmat=cmat)
x4 <- out[,1]
y4 <- out[,2]

df <- tibble(x = c(x1, x2, x3, x4),
             y = c(y1, y2, y3, y4),
             g = factor(rep(c("g=0","g=0.2","g=0","g=0.2"), each = n)),
             h = factor(rep(c("h=0","h=0","h=0.2","h=0.2"), each = n)))
```

### Make figure

``` r
pB <- ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.1) +
  # geom_smooth(method = "lm") +
  facet_grid(rows = vars(h), cols = vars(g)) +
  coord_cartesian(xlim = c(-10, 25), ylim = c(-10, 25)) +
  theme(axis.title = element_blank())
  # theme(panel.grid.minor = element_blank())
pB
```

![](g_h_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# ggsave("./figures/gh_rho05.png", width = 20, height = 20, units = "cm")
```

## Combine panels

``` r
cowplot::plot_grid(pA, pB,
                   labels = c("A", "B"))
ggsave("./figures/gh_summary.pdf", width = 40, height = 20, units = "cm")
```

# Summary figure: vary g, h = 0

n = 5000 g = 0.1 -\> 0.9 h = 0

## rho = 0

### Get samples

``` r
set.seed(21)
n <- 5000
rho.pop <- 0
h <- 0
gseq <- seq(0, 1, 0.1)
ng <- length(gseq)
res <- array(NA, dim = c(n, 2, ng))
for(G in 1:ng){
  g <- gseq[G]
  rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
  cmat <- matrix(c(1,rho,rho,1),2,2)
  res[,,G] <- ghmul(n,g=g,h=h,p=2,cmat=cmat)
}

df <- tibble(x = as.vector(res[,1,]),
             y = as.vector(res[,2,]),
             g = factor(rep(paste0("g=",gseq), each = n)))
```

### Make figure

``` r
pA <- ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.1) +
  # geom_smooth(method = "lm") +
  facet_wrap(vars(g), ncol = 3) +
  coord_cartesian(xlim = c(-5, 25), ylim = c(-5, 25)) +
  theme(axis.title = element_blank())
  # theme(panel.grid.minor = element_blank())
pA
```

![](g_h_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# ggsave("./figures/gvar_h0_rho0.png", width = 20, height = 20, units = "cm")
```

## rho = 0.5

### Get samples

``` r
set.seed(21)
n <- 5000
rho.pop <- 0.5
h <- 0
gseq <- seq(0, 1, 0.1)
ng <- length(gseq)
res <- array(NA, dim = c(n, 2, ng))
for(G in 1:ng){
  g <- gseq[G]
  rho <- rngh.sub(50,g,h,rho.pop)$rho.adjusted
  cmat <- matrix(c(1,rho,rho,1),2,2)
  res[,,G] <- ghmul(n,g=gseq[G],h=h,p=2,cmat=cmat)
}

df <- tibble(x = as.vector(res[,1,]),
             y = as.vector(res[,2,]),
             g = factor(rep(paste0("g=",gseq), each = n)))
```

### Make figure

``` r
pB <- ggplot(df, aes(x=x, y=y)) + theme_gar +
  geom_point(size = 0.5, alpha = 0.1) +
  facet_wrap(vars(g), ncol = 3) +
  coord_cartesian(xlim = c(-5, 25), ylim = c(-5, 25)) +
  theme(axis.title = element_blank())
pB
```

![](g_h_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
# ggsave("./figures/gvar_h0_rho05.png", width = 20, height = 20, units = "cm")
```

## Combine panels

``` r
cowplot::plot_grid(pA, pB,
                   labels = c("A", "B"))
ggsave("./figures/gh_summary.pdf", width = 40, height = 20, units = "cm")
```
