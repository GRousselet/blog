R code and notebook for my blog post on bootstrap bias correction:
<https://garstats.wordpress.com/2018/01/24/10000/>

`10000_experiments.Rmd` generates all the figures presented in the post, and the notebook `10000_experiments.pdf`.

`french_lexicon_project_rt_data.RData` contains the raw data from the French lexicon project. They were organised into a data frame using `getflprtdata.Rmd`.

`akerd.txt` contains functions to generate kernel density estimates.
The code is from Rand Wilcox’s much larger `Rallfun-v34.txt`, available [here](http://dornsife.usc.edu/labs/rwilcox/software/).

`HDIofMCMC.R` contains a function to compute highest density intervals.
The code is from John K. Kruschke’s Doing bayesian data analysis [book](https://sites.google.com/site/doingbayesiandataanalysis/).

The various `.RData` files contain simulation results, to speed things up if you want to interact with the code without running slow chunks.