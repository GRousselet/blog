Matlab code for my blog post on the Harrell-Davis estimator:
<https://garstats.wordpress.com/2016/06/09/the-harrell-davis-quantile-estimator/>

`hd_demo.m` generates all the figures presented in the post.
The file `simres.mat` contains the results of the simulations described in the blog.
It is used by `hd_demo.m` to reproduce the main results without having to run the simulation.

## Dependencies
To run the code, you will need functions from these toolboxes:

<https://github.com/GRousselet/matlab_stats>
<https://github.com/GRousselet/matlab_visualisation>

And the archive `data.zip` from this data package:

<https://figshare.com/articles/Face_noise_ERP_onsets_from_194_recording_sessions/1588513/1>

The relevant files `age.txt` and `onset.txt` are included in the blog folder for convenience. 
So if you have the Matlab toolboxes in your path, and you cd to the blog folder, you should be good to go. 

## Abbreviations
pb = percentile bootstrap

ci = confidence interval

se = standard error

hd = Harrell-Davis quantile estimator

## Matlab and R functions
Table of the main functions, from the `matlab_stats` toolbox, relevant to this post.
R functions are available on Rand Wilcox’s [website](http://dornsife.usc.edu/labs/rwilcox/software/).

|Matlab function|What it computes|R function
|-----|-----|-----|
|`hd`|Harrell-Davis estimator|`hd`|
|`hdci`|1 quantile and its ci using a pb estimation of the se of hd|`hdci`|
|`hdpbci`|1 quantile and its ci using a pb of hd|`qcipb`|
|`deciles`|deciles, using hd|`deciles`|
|`decilesci`|deciles and their ci using a pb estimation of the se of hd|not available|
|`decilespbci`|deciles and their ci using a pb of hd|not available
|`pbci`|pb ci of any estimator, including hd|`onesampb`|
|`bootse`|bootstrap estimate of the se of any estimator, including hd|`bootse`|

See my previous [post](https://garstats.wordpress.com/2016/05/27/the-percentile-bootstrap/) on the percentile bootstrap for a list of other Matlab and R functions to compute confidence intervals.

You can see the Harrell-Davis estimator in action here for instance:
http://onlinelibrary.wiley.com/doi/10.1111/ejn.13100/full