R code for my blog post on difference distributions:
<https://garstats.wordpress.com/2016/07/19/typical-differences/>

`typical_difference_demo.R` generates all the figures presented in the post.

All the dependencies are listed at the beginning of the file: 
in addition to packages `ggplot2`, `cowplot` & `tidyr`, a few house-made toolboxes are provided as text files.

The main R functions are available on Rand Wilcox’s [website](http://dornsife.usc.edu/labs/rwilcox/software/) - I’ve included the official current version, as of writing, as `Rallfun-v30.txt`. Most functions were covered in previous posts, except for these functions on the asymmetry of difference distributions:
- *independent groups*
`cbmhd`
`qwmwhd`
- *dependent groups*
`Dqdif`
`difQpci`
`difQplot`

Datasets are either generated in the demo file, or loaded from an `.rds` file - ensure all the files in this folder are in your working directory.