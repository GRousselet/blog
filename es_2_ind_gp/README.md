Matlab code for my blog post on robust measures of effect size for two independent groups:
<https://garstats.wordpress.com/2016/05/02/robust-effect-sizes-for-2-independent-groups/>

`es_2_ind_gp.m` generates all the figures presented in the post, and a few extra ones.
`varres.mat` & `sysres.mat` contain the results of two simulations described in the m file. 
The data can be generated from the m file, but using the mat files will save a few minutes.

## Dependencies
To run the code, you will need functions from these two toolboxes:

<https://github.com/robince/gcmi>

<https://github.com/GRousselet/matlab_stats>

## Main functions

Here is a table of the main functions, from the `matlab_stats` toolbox, relevant to this post, and their R equivalent, 
available on Rand Wilcox’s [website](http://dornsife.usc.edu/labs/rwilcox/software/) unless indicated otherwise. 
Mutual information is calculated using functions from the `gcmi` toolbox, which also contains Python code.

|Matlab function|What it computes|R function|
|-----|-----|-----| 
|`cohend`|Cohen’s *d*|available in packages listed below|
|`cid`|Cliff’s *delta* + probabilities + confidence interval|`cid`|
|`cliffdelta`|for convenience, only outputs Cliff’s *delta*|na|
|`ksstat`|Kolmogorov-Smirnov test statistic, no p value is returned|`ks` implements regular & weighted versions of the KS test|
|`ksstat_fig`|illustration of the Kolmogorov-Smirnov test statistic|na|
|`qhat`|Wilcox & Muska’s *Q* non-parametric measure of effect size|`qhat`|
|`mapd`|median of all pairwise differences|`wmwloc`|
|`l2dci`|bootstrap confidence interval for mapd|`l2dci`|
|`akerd`|adaptive kernel density estimate|`akerd`|

## Other R packages

In addition to Wilcox’s code, these R packages compute Cohen’s *d*, Cliff’s *delta* and other measures of effect sizes:
[effsize](https://cran.r-project.org/web/packages/effsize/) 
[compute.es](https://cran.r-project.org/web/packages/compute.es/)
[orddom](https://cran.r-project.org/web/packages/orddom/)
