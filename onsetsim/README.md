Can we use cluster-based permutation tests to estimate MEG/EEG onsets?

R and Matlab code for my [blog post][1] on onset estimation:

## R code | Notebook | Content || ----- | ----- ||`onsetsim_eeg.Rmd`| Simulate one participant using EEG-like noise with the approach from Yeung et al. (2004).||`onsetsim_eeg_group.Rmd`| Simulate 20 participants with random onsets.||`onsetsim_1overf.Rmd`| Simulate one participant with 1/f noise.||`examples.Rmd`| Illustrate noise, signal and statistical methods.|## Matlab code| Notebook | Content || ----- | ----- ||`onsetsim.m`| Illustrate 1/f noise & methods and perform simulations||`onsetsim_eegnoise.m`| Illustrate EEG-like noise (Yeung et al., 2004) & methods and perform simulations|## DependenciesExtra code dependencies are in the `code` folder for R and in the `functions` folder for Matlab.## Simulation resultsAll simulation results are in the `data` folder.

[1]:	https://garstats.wordpress.com/2023/11/01/onsetsim/