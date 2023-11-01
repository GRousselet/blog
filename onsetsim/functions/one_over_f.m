% Matlab version of the R function one_over_f
% https://rdrr.io/cran/primer/src/R/one_over_f.R

% #' Function to generate 1/f noise.
% #'
% #' Generates 1/f noise with a specified power or amplitude.
% #' @param gamma spectral power, numeric, where 0 generates a white noise time series, 2 generates reddened noise. Defaults to 1 (pink).
% #' @param N length of the time series.
% #' @keywords noise 1/f color spectra
% #' @export
% #' @author Hank Stevens
% #' @references
% #' J. M. Halley. Ecology, evolution and 1/f-noise. Trends in Ecology & Evolution, 11:33-37, 1996.
% #' O. L. Petchey, A. Gonzalez, and H. B. Wilson. Effects on population persistence: the interaction between environmental noise colour, intra-specific competition and space. Proceedings of the Royal Society of London Series B, 264:1841-1847, 1997.
% #' J. E. Cohen, C. M. Newman, A. E. Cohen, O. L. Petchey, and A. Gonzalez. Spectral mimicry: a method of synthesizing matching time series with different Fourier spectra. Circuits, Systems and Signal Processing, 18:431-442, 1999.
% #' @seealso
% #' [spec_mimic()] to rearrange one vector, X, to mimic the spectrum of another vector, Y; [plot_f()] to plot the time series and the spectrogram of the series.
% #' @keywords 1/f color noise spectra
% #' @examples
% #'
% #' set.seed(1)
% #' time.series <- one_over_f(gamma=2, N=50)
% #' plot(1:50, time.series, type='l', main="Reddened noise")
% #' time.series <- one_over_f(gamma=0, N=50)
% #' plot(1:50, time.series, type='l', main="White noise")
% #' one_over_f()
% #'
% #' @export one_over_f
% one_over_f <- function(gamma=1, N=200){
%   ## Generate 1/f noise with power = gamma
%   ## after Petchey SAS code, etc.
%   N2 <- N/2
%   sine.waves <- matrix(NA, nrow=N, ncol=N2)
%   steps=2*pi*(1:N)/N
%   phase <- stats::runif(N2, 0, 2*pi)
%   for(i in 1:N2) {
%     freq <- i
%     weight <- 1/(freq^gamma)
%     y <- weight*sin(freq*steps+phase[i])
%     sine.waves[,i] <- y
%
%   }
%   out <- rowSums(sine.waves)
%   # force sd = 1 and mean = 0
%   out <- out - mean(out)
%   out <- out / sd(out)
%   return(out)
% }

function out = one_over_f(gamma, N, outvar)

if nargin < 3
   outvar = 1; % output variance 
end

% ensure that the M is even
if rem(N,2)
    M = N+1;
else
    M = N;
end
N2 = M/2;
sinewaves = zeros(N, N2);
steps = 2*pi*(1:N)/N;
phase = unifrnd(0,2*pi, N2, 1);
for F = 1:N2
    freq = F;
    weight = 1/(freq^gamma);
    sinewaves(:,F) = weight*sin(freq*steps+phase(F));
end
out = sum(sinewaves, 2);
% force sd = 1 and mean = 0
out = out(1:N);
out = out - mean(out);
out = out / std(out);
out = out .* sqrt(outvar);

end
