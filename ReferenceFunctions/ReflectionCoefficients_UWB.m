%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This function calculates the effective reflectivity of a stack of layers
% of differing refractive index, taking into account the full bandwidth of
% the radar system.
% 
% Input Variables:
% fc - radar center frequency in Hz
% BW - radar bandwidth in Hz
% n - 1 x N vector of refractive indices for each layer
% n_depth - 1 x N depth axis in meters for the refractive index profile
% ref - frequency spectrum of the transmitted chirp (1 x M vector)
% freq_win - frequency window used in matched filtering (1 x M vector)
% freq_axis - frequency domain axis for the frequency window (1 x M vector)
% 
% Output Variables:
% reflectivity - effective reflectivity (complex scalar double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reflectivity = ReflectionCoefficients_UWB(fc, BW, n, n_depth, ref, freq_win, freq_axis)
    
    c = 299792458;   % speed of light in a vacuum m/s
    % Set up the spacing of frequency samples for calculating the
    % reflectivity
    f = (fc - BW/2):10e6:(fc + BW/2);

    core_diff = diff(n_depth);
    R = zeros(length(f), 1);

    % For every frequency, use the transfer matrix method to calculate the
    % effective reflectivity of the layer stack
    for m = 1:length(f)
        lambda = c/f(m); 
        tmp = TransferMatrix(lambda, 0, n, core_diff(1:end-1));
        R(m) = tmp(1);    
    end
    
    % Generate the radar system impulse response
    ipr = fftshift(ifft(conj(ref).*ref.*freq_win));
    freq_reflec = interp1(f - (fc - BW/2), R, freq_axis);
    
    % Set the effective reflectivity at frequencies outside the radar
    % bandwidth
    freq_reflec(isnan(freq_reflec)) = mean(abs(R));
    freq_reflec = [freq_reflec mean(abs(R))*ones(1, length(ref) - length(freq_reflec))];
    
    % Weight the transmitted signal by the refletivity in the frequency
    % domain, matched filter, and rescale by max value of IPR to get true
    % reflectivity
    reflec_spectrum = ref.*freq_reflec;
    pc = ifft(fftshift(conj(ref).*reflec_spectrum.*freq_win));
    reflectivity = max(pc)./max(ipr);
      
end