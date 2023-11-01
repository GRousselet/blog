function onset = find_onset(sigmask,Xf,rmzero)

% FIND_ONSET gets the latency of the first value larger than zero in a vector.
%
% FORMAT    onset = find_onset(sigmask,Xf,rmzero)
%
% INPUTS:
%           sigmask = a vector of 0s and 1s (output from limo_ecluster_test for instance).
%           Xf      = a vector of time points, matching sigmask (typically in ms, from EEG.times for instance).
%           rmzero  = 0/1 option to discard an onset if it belongs to a cluster that starts in the baseline; default = 1.
%
% OUTPUT:
%           onset   = latency of the first cluster in the units of Xf, typically in ms.
%
% v1 Guillaume Rousselet, University of Glasgow, May 2022
% ----------------------------------------------------------------------------

onset = [];
if nargin < 3
    rmzero = 1;
end

if rmzero == 0 % do not remove baseline onsets
    onset = find(sigmask > 0);
    onset = Xf(onset(1));
end

if rmzero == 1 % ignore baseline onsets
    [L,NUM] = bwlabeln(sigmask); % find clusters
    maskedsig = sigmask;
    if NUM~=0
        tmp = zeros(1,NUM);
        for C = 1:NUM % check if cluster includes baseline
            sigtimes = Xf(L == C);
            if min(sigtimes) <= 0
                maskedsig(L == C) = 0;
            end
        end
    end
    onset = find(maskedsig > 0);
    onset = Xf(onset(1));
end
