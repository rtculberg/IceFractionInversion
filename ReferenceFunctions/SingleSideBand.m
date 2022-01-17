%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 12/7/2020
%
% This function models the transmitted signal from a chirped pulse radar
% system.
%
% Inputs:
%   slope - frequency slope of the chirp, calculated as bandwidth (in Hz)
%   divided by pulse length (in seconds)
%   length - pulse length in seconds
%   rate - sample rate in Hz
%   center_freq - center frequency in Hz
%
% Outputs:
%   data - modeled chirp (1xM complex double vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = SingleSideBand(slope, length, rate, center_freq)
    npts = length*rate;
    t = linspace(0,length,npts);
    phase =(pi*slope*(t.^2))+(2*pi*center_freq*t);
    data = exp(1i.*phase);
end
