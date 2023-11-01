# Adapted from Matlab code:
# https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator

# function signal = noise (frames, epochs, srate)
# 
# % function signal = noise (frames, epochs, srate)
# %
# % Function generates noise with the power spectrum of human EEG
# % Inputs:
#   %  frames - number of signal frames per each trial
# %  epochs - number of simulated trials
# %  srate - sampling rate of simulated signal
# % Output:
#   %  signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
# % Implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002
# 
# load meanpower  
# sumsig = 50;	%number of sinusoids from which each simulated signal is composed of
# 
# signal = zeros (1, epochs * frames);
# for trial = 1:epochs
# freq=0;
# range = [(trial-1)*frames+1:trial*frames];
# for i = 1:sumsig
# freq = freq + (4*rand(1));
# freqamp = meanpower(min (ceil(freq), 125)) / meanpower(1);
# phase = rand(1)*2*pi;
# signal (range) = signal (range) + sin ([1:frames]/srate*2*pi*freq + phase) * freqamp;
# end
# end

# Generate noise with the power spectrum of human EEG
# Inputs:
#   frames - number of signal frames per each trial
#   epochs - number of simulated trials
#   srate - sampling rate of simulated signal
# Output:
#   signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
# Implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002
# R version + control variance: GAR, University of Glasgow, Oct 2023
eeg_noise <- function(frames = 51, srate = 100, outvar = 1, meanpower){
  sumsig <- 50 # number of sinusoids from which each simulated signal is composed of
  out <- vector(mode= "numeric", length = frames)
    freq <- 0
    for(i in 1:sumsig){
      freq <- freq + (4*runif(1))
      freqamp <- meanpower[min(ceiling(freq), 125)] / meanpower[1]
      phase <- runif(1)*2*pi
      out <- out + sin((1:frames)/srate*2*pi*freq + phase)*freqamp
    }
    out <- out - mean(out);
    out <- out / sd(out);
    out <- out * sqrt(outvar);
    return(out)
}
