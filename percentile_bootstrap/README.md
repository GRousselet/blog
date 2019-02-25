Matlab code for my blog post on the percentile bootstrap:
<https://garstats.wordpress.com/2016/05/27/the-percentile-bootstrap/>

`pb_demo.m` generates all the figures presented in the post.

## Dependencies
To run the code, you will need functions from these toolboxes:

<https://github.com/GRousselet/matlab_stats>

<http://www.mathworks.com/matlabcentral/fileexchange/54243-univarscatter>

However, for the purpose of this post, I have modified the main function of the `univarscatter` toolbox
and renamed it `UnivarScatter_nofig`. So I provide the original toolbox + modified function in
the folder for this post.

## Abbreviations
pb = percentile bootstrap
ci = confidence interval
se = standard error
hd = Harrell-Davis quantile estimator
M = M-estimator of location based on Huber’s *Psy*
tm = trimmed mean

## Matlab functions
Table of the main functions, from the `matlab_stats` toolbox, relevant to this post.

|Matlab function|What it computes|
|-----|-----|
|`pbci`|pb ci of any estimator|
|`bootse`|bootstrap estimate of the se of any estimator|
|`pb2ig`|pb ci of the difference between two estimators of independent observations|
|`pb2dg`|pb ci of the difference between two estimators of paired observations|

## R functions
Table of the main functions, relevant to this post, and available on Rand Wilcox’s [website](http://dornsife.usc.edu/labs/rwilcox/software/).

|R function|What it computes|
|-----|-----|
|`onesampb`|pb ci of any estimator|
|`mestci`|pb ci of M|
|`momci`|pb ci of modified one-step M-estimator of location|
|`bootse`|bootstrap estimate of the se of any estimator|
|`hdseb`|bootstrap estimate of the se of hd|
|`mestseb`|bootstrap estimate of the se of M|
|`pb2gen`|pb ci of the difference between two estimators of independent observations|
|`m2ci`|ci for the difference between two M-estimators|
|`comvar2`|ci for the difference between two variances|
|`bootdpci` with option dif=F|pb ci of the difference between two estimators of paired observations|
|`bootdpci` with option dif=T|pb ci of the estimator of the paired differences|
|`trimpb2`|ci of the difference between 2 tm of independent observations|
|`trimci`|ci of the tm|
