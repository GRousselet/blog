Data and Matlab code for my blog post on ERP differences
<https://garstats.wordpress.com/2016/04/02/simple-steps-for-more-informative-erp-figures/>

The file `erp_differences.m` generates the main figures, providing you have installed the toolbox `limo eeg`, available on github.
<https://github.com/LIMO-EEG-Toolbox/limo_eeg>

`erp_differences` uses two files, both containing data from one electrode only. I used electrode B8 from the Biosemi 128-channel layout. It is a right posterior electrode which tends to record the largest face responses across participants. No bad trial rejection was done. 

- `erps.mat` contains the variable `erps`, a 20 participants x 2 conditions cell array. Each cell contains a matrix with dimensions 1 electrode x 451 time points x 192 trials. The 451 time points correspond to a 500 Hz sampling from -300 ms to 600 ms post-stimulus onset.
- `erpm.mat` contains the variable `erpm`, which is a matrix of 451 time points x 20 participants x 2 conditions of mean averages across trials.

The file `fake_erp_differences.m` demonstrates a few examples of fake ERPs and their differences. The point is that it is sometimes difficult to form a mental image of the difference from the original conditions alone.
