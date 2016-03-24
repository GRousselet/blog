# onset_priors
# Guillaume Rousselet <Guillaume.Rousselet@glasgow.ac.uk>
# https://garstats.wordpress.com
# University of Glasgow
# version: 2016-03-17

# The code below is mostly adapted from:
# BEST1G.R & BEST1Gexample.R - version of Dec 02, 2015
# BESTgamma.R & BESTgamma_example - version of 4/28/15
# John K. Kruschke
# johnkruschke@gmail.com
# http://www.indiana.edu/~kruschke/BEST/

# OPTIONAL: Clear R's memory and graphics:
rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off() # This closes all of R's graphics windows.

# set working directory to the onset_priors folder
setwd("/onset_priors/")
# declare functions:
# these files should be in the onset_priors folder
# the first two files are from Kruschke's BEST toolbox, but with small cosmetic modifications
# http://www.indiana.edu/~kruschke/BEST/
# https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/tree/master/2e
# the third file is Rand Wilcox's R toolbox
# http://dornsife.usc.edu/labs/rwilcox/software/
source("BEST1G.txt")
source("DBDA2E-utilities.txt")
source("Rallfun-v29.txt")

# ------------------------------------------------------------------------
# First, create fake data to check that the procedure works
# ------------------------------------------------------------------------
set.seed(2) # set the seed to get same data - but MCMC results will differ
targetsd = 5 # our fake data will have sd of 5
targetmean = 100 # and mean of 100
nobs = 20 # number of observations
y = rt(nobs,df=30) # we sample from a t distribution
y = (y-mean(y))/sqrt(mean((y-mean(y))^2))*targetsd + targetmean

# data description
akerd(y) # adaptive kernel density estimate
summary(y)
sd(y)

# Run Bayesian analysis using BEST defaults:
#   here I've modified the original BEST1Gmcmc (in BEST1G.txt) so it returns codaSamples instead of mcmcChain
#   this is useful to run diagnostic tools on the chains
codaSamples = BEST1Gmcmc( y , numSavedSteps=10000, thinSteps=1, showMCMC=FALSE )
mcmcChain = as.matrix(codaSamples)
# Display the results:
BEST1Gplot( y , mcmcChain , compValm=100 )
# save the figure:
saveGraph( file="test" , type="eps" )

# call diagnostic plotting function
diagMCMC( codaObject=codaSamples, parName = "mu", saveName="test")
# the 3 chains appear to be well mixed and to be sampling from a high-probability region.
# however there is some autocorrelation, and the ESS is <7000.
# to see if this is a concern, we could run the MCMC analysis again, and compare the 95% HDIs.
# we could also thin the chains, or increase the burn in period, or increase the number of steps,
# or any combination of these actions.

# With mild thinning:
codaSamples = BEST1Gmcmc( y , numSavedSteps=10000, thinSteps=5, showMCMC=FALSE )
diagMCMC( codaObject=codaSamples, parName = "mu", saveName="test_with_thinning")
# the results show improved HDI consistency across the 3 chains, from already very consistent results.
# in this case, mild thinning is sufficient to get an ESS around the expected 10,000.
# However, in this case, if our main interest is the mode and the 95% HDI of the mean,
# thinning makes no noticeable difference. The required precision of the estimation depends
# on the application of course.

# ***************************************************
# add outliers to check robustness
openGraph(width=7,height=7)
y <- c(y,150,200)
akerd(y) # adaptive kernel density estimate
summary(y)
sd(y)

# Run Bayesian analysis again:
codaSamples = BEST1Gmcmc( y , numSavedSteps=10000, thinSteps=1, showMCMC=FALSE )
mcmcChain = as.matrix(codaSamples)
# Display the results:
BEST1Gplot( y , mcmcChain , compValm=100 )
# save the figure:
saveGraph( file="test_with_outliers" , type="eps" )
# MCMC with a t distribution is robust:
# the estimation of the mean is barely different from the original,
# and within random fluctuations expected from running the same MCMC several times.
# ***************************************************

# Run Bayesian analysis using uniform prior on mean and gamma prior on sd,
# instead of BEST's default settings:
y <- y[1:20]
source("BEST1Gunif.txt")
# we expect onsets to be contained in the broad interval 0-200 ms:
unif_low <- 0
unif_high <- 200
codaSamples = BEST1Gunifmcmc( y , unif_low , unif_high )
mcmcChain = as.matrix(codaSamples)
# Display the results:
BEST1Gplot( y , mcmcChain , compValm=100 , ROPEm=c(-5,5) , pairsPlot=FALSE )
saveGraph( file="test_with_unif_prior" , type="eps" )
# We get very similar results, as expected because the original BEST1Gmcmc function
# uses very broad priors informed by the data, so that the algorithm does not
# start too far off from the high-probability region.

# Run Bayesian analysis with strong priors on mean:
source("BEST1Gpriors.txt")
# Here we expect onsets to be arbitrarily around 80 ms, with little dispersion around that value.
# I'm using 100,000 steps because the steps are strongly correlated in this example.
# You can try different combinations of steps and thinning to see how that affects autocorrelation.
codaSamples = BEST1Gpriorsmcmc( y , muPriorMean=80 , muPriorSD=3 ,
              sigmaPriorMode=sd(y) , sigmaPriorSD=sd(y)*5 ,
              nuPriorMean=30 , nuPriorSD=30,
              numSavedSteps=100000, thinSteps=1, showMCMC=FALSE )
# in that case the posterior predictive samples for the mean look
# very different from those obtained using broad priors:
openGraph(width=7,height=7)
plotPost( codaSamples[,"mu"] , main="mean", cenTend="mode",
                     compVal=100, ROPE=c(95,105), credMass=0.95, HDItextPlace=0.7,
                     xlab=bquote(mu) , xlim=c(80,120))
saveGraph( file="test_with_strong_priors" , type="eps" )
diagMCMC( codaObject=codaSamples, parName = "mu")
# in that case the 95% HDI does not include 100, the sample mean.
# the HDI does overlap with a proposed (95,105) ROPE (region of practical equivalence),
# so we should reserve judgement.
# more about using ROPEs here:
# http://doingbayesiandataanalysis.blogspot.co.uk/2013/08/how-much-of-bayesian-posterior.html

# ------------------------------------------------------------------------
# Now that we've gained some confidence in the tools,
# let's do Bayes on our real dataset of ERP onsets:
# ------------------------------------------------------------------------
x <- read.table("onset.txt")
ses1 <- x[,3]
x <- read.table("onset2.txt")
ses2 <- x[,3]
# ses1 contains onsets from 120 participants
# ses2 contains onsets from 74 participants who also provided ses1 onsets
# although this is a paired design, and the results can be used to investigate
# test-retest reliability (see our EJN 2015 paper), here we treat the two sessions
# as independent

# kernel density estimation for the two sessions
# openGraph(width=7,height=7)
g2plot(ses1,ses2,op=4,xlab="Onsets in ms",ylab="Density")
legend('topright',legend=c('session 1', 'session 2'),
       col=c('black', 'black'),
       lty=c(1,2))
summary(ses1)
sd(ses1)

# ------------------------------------------------------------------------
# Run Bayesian analysis on session 1 using BEST defaults:
codaSamples = BEST1Gmcmc( ses1 , numSavedSteps=10000 )

# call diagnostic plotting function
diagMCMC( codaObject=codaSamples, parName = "mu", saveName="ses1")
# the diagnostic tools suggests the chain are well mixed and seem to
# be sampling from a high-probability region.
# However the ESS is ~6000, which calls for some thinning or more steps.
# So let's try again with 20,000 steps:
codaSamples = BEST1Gmcmc( ses1 , numSavedSteps=20000, thinSteps=1 )
diagMCMC( codaObject=codaSamples, parName = "mu", saveName="ses1_with_more_steps")
# ESS is now ~10000, and all the diagnostic indicators look fine.
# It seems we've got a good sampling of the posterior distribution.

# Display the results, with a 100 ms mean onset for comparison:
mcmcChain = as.matrix(codaSamples)
BEST1Gplot( ses1, mcmcChain, compValm=100, pairsPlot=FALSE  )
saveGraph( file="ses1" , type="eps" )
# here I get a mode of 92.6, 95% HDI [88.4, 96.7]
# If we had the hypothesis of a mean onset of 100 ms, the results suggest
# a credible mean lower than 100 ms.

# plot only mean posterior samples + ROPE:
openGraph(width=7,height=7)
plotPost( codaSamples[,"mu"] , main="mean", cenTend="mode",
          compVal=100, ROPE=c(95,105), credMass=0.95, HDItextPlace=0.7,
          xlab=bquote(mu) , xlim=c(80,120))
saveGraph( file="ses1_mean_posterior" , type="eps" )
# a 95-105 ms ROPE overlaps with the 95% HDI, so we conclude that 100 ms
# is still a credible value.

# ------------------------------------------------------------------------
# What happens is we use more informed priors?
# For instance, the literature suggests that face onsets are around 50-150 ms.
source("BEST1Gunif.txt")
unif_low <- 50
unif_high <- 150
codaSamples1 = BEST1Gunifmcmc( ses1 , unif_low , unif_high , numSavedSteps=20000, thinSteps=1 )
mcmcChain = as.matrix(codaSamples1)
diagMCMC( codaObject=codaSamples1, parName = "mu", saveName="ses1_unifprior")

# Display the results:
BEST1Gplot( ses1 , mcmcChain , compValm=100 , pairsPlot=TRUE )
saveGraph( file="ses1_unifprior" , type="eps" )

openGraph(width=7,height=7)
plotPost( codaSamples1[,"mu"] , main="mean", cenTend="mode",
          compVal=100, ROPE=c(95,105), credMass=0.95, HDItextPlace=0.7,
          xlab=bquote(mu) , xlim=c(80,120))
saveGraph( file="ses1_unifprior_mean_posterior" , type="eps" )
# as expected with uniform priors, the results changed very little at all
# compared to using BEST default options.

# save chains for re-analysis:
save( codaSamples1, file="ses1_codasamples.Rdata" )

# ------------------------------------------------------------------------
# With our new posterior estimates, we can now study results from session 2
# using more informed priors. We plugin estimates of the posterior samples
# from session 1 as priors for session 2. This makes a lot of sense given that
# session 1 has a large sample size (for this type of research) and is therefore
# our best description of the population we're trying to estimate.
source("BEST1Gpriors.txt")

mupm = mean( mcmcChain[,"mu"] )
mupsd = sd( mcmcChain[,"mu"] )
mcmcDensity = density(mcmcChain[,"sigma"])
sigmapm = mcmcDensity$x[which.max(mcmcDensity$y)] # mode
sigmapsd = sd( mcmcChain[,"sigma"] )
nupm = mean( mcmcChain[,"nu"] )
nupsd = sd( mcmcChain[,"nu"] )

codaSamples2 = BEST1Gpriorsmcmc( ses2 , muPriorMean=mupm , muPriorSD=mupsd ,
                              sigmaPriorMode=sigmapm , sigmaPriorSD=sigmapsd ,
                              nuPriorMean=nupm , nuPriorSD=nupsd ,
                              numSavedSteps=20000 , thinSteps=1)
mcmcChain = as.matrix(codaSamples2)
diagMCMC( codaObject=codaSamples2, parName = "mu", saveName="ses2_ses1prior")
save( codaSamples2, file="ses2_codasamples.Rdata" )

# Display the results:
BEST1Gplot( ses2 , mcmcChain , compValm=100 , pairsPlot=TRUE )
saveGraph( file="ses2_ses1prior" , type="eps" )

openGraph(width=7,height=7)
plotPost( codaSamples2[,"mu"] , main="mean", cenTend="mode",
          compVal=100, ROPE=c(95,105), credMass=0.95, HDItextPlace=0.7,
          xlab=bquote(mu) , xlim=c(80,120))
saveGraph( file="ses2_ses1prior_mean_posterior" , type="eps" )
# here we've gained information:
# the ROPE is now completely outside the 95% HDI.
# ROPE = [95, 105], 95% HDI = [87.3, 93.4], mode of the mean = 90.3
# So 100 ms is not a credible onset value anymore.

# What do we get if we pretend that session 1 does not exist?
codaSamples = BEST1Gmcmc( ses2 , numSavedSteps=20000, thinSteps=1 )
diagMCMC( codaObject=codaSamples, parName = "mu", saveName="ses2_broadpriors")
mcmcChain = as.matrix(codaSamples)
BEST1Gplot( ses2, mcmcChain, compValm=100, pairsPlot=FALSE  )
saveGraph( file="ses2_broadpriors" , type="eps" )
# ses2 with ses1 priors: 95% HDI = [87.3, 93.4], mode of the mean = 90.3
# ses2 with braod priors: 95% HDI = [83.2, 92.2], mode of the mean = 88.1
# because the results are very similar across the two sessions,
# using broad priors instead of session 1 priors give very similar results.
# However, ses2 estimates using broad priors are shifted to the left, because
# the session 2 data distribution suggests a slightly lower mean compared to session 1.

# Comparison to a new point estimate
# Let say someone finds an onset at 50 ms, probably using group statistics.
# Before writing a paper declaring "face processing takes 50 ms", is it credible given
# our estimation from our 2 large samples?
testpe <- 50
rside <- 10 # 10 ms ROPE boundaries -> 20 ms ROPE
load("ses2_codasamples.Rdata" )

plotPost( codaSamples2[,"mu"] , main="mean", cenTend="mode",
          compVal=testpe, ROPE=c(testpe-rside,testpe+rside), credMass=0.95, HDItextPlace=0.7,
          xlab=bquote(mu) , xlim=c(0,100))
# Well, as you've guessed, the answer is no, it is not statistically credible!
# It is also not physiologically credible!
