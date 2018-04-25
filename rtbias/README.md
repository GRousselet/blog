R code for the series of blog posts on median bias:
# Reaction times and other skewed distributions: problems with the mean and the median

In this series of 4 posts, I replicate, expand and discuss the results from
> Miller, J. (1988) A warning about median reaction time. J Exp Psychol Hum Percept Perform, 14, 539-543.

## Part 1 = replicate Miller’s simulations + apply bootstrap bias correction
[Blog post](https://garstats.wordpress.com/2018/02/02/rtbias1/)

[R code](docs/miller1988.md)

## Part 2 = expand Miller’s simulations to group comparison

[Blog post](https://garstats.wordpress.com/2018/02/08/rtbias2/)

[R code](docs/bias_diff.md)

## Part 3 = problems with the mean

[Blog post](https://garstats.wordpress.com/2018/03/30/rtbias3/)

[R code](docs/samp_dist.md)

## Part 4 = application to a large dataset

[Blog post](https://garstats.wordpress.com/2018/04/25/rtbias4/)

R code:

- [Description of dataset](docs/flp_illustrate_dataset.md)

- [Simulation](docs/flp_simulation.md)

## Extra post on Cohen's *d*

[Blog post](https://garstats.wordpress.com/2018/04/04/dbias/)

[R code](docs/miller1988_cohend.md)

## Dependencies

Packages:
`ggplot2`
`tibble`
`tidyr`
`cowplot`
`retimes`
`viridis`
`knitr`
`HDInterval`

`Rallfun-v34.txt` contains Rand Wilcox’s functions and is also available [here](http://dornsife.usc.edu/labs/rwilcox/software/).

`akerd.txt` contains a subset of Rand Wilcox's functions to compute adaptive kernel density estimates.

`HDIofMCMC.R` contains a function to compute highest density intervals.
The code is from John K. Kruschke’s [`Doing bayesian data analysis book`](https://sites.google.com/site/doingbayesiandataanalysis/).

## Data

The `.RData` files contain simulation results, to speed things up if you want to interact with the code without running slow chunks.

Part 4 uses data from the French lexicon project. Data and code to turn them into a data frame are available [here](https://github.com/GRousselet/blog/tree/master/10000). Or directly on the project's [webpage](https://sites.google.com/site/frenchlexicon/results).
