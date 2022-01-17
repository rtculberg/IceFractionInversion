%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This function simulates a random density profiles based on an ARMA model
% trained on high resolution density measurement from Greenland firn cores.
% 
% Input Variables:
% M - melt feature percentage (MFP) - fraction of range resolution occupied 
% by refrozen ice, value from 0 to 1
% mean_density - mean density of the firn in g/cm^3
% coefficients - ARMA model coefficients
% variability - variability of the firn density in g/cm^3
% depth - depth axis in meters
% 
% Output Variables:
% density - 1 x length(depth) vector of densities in g/cm^3
% M - melt feature percentage
% N - number of layers
% T - 1 x N vector of ice layer thickness
% G - 1 x N vector of gap layer thickness
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [density, M, N, T, G] = DensityEmulator(M, mean_density, coefficients, variability, depth)
    
    spacing = mean(diff(depth));
    
    % Simulate a random density profiles using the ARMA model, mean
    % density, and variability
    density_variance = simulate(coefficients, length(depth));
    density_variance = density_variance';
    density = (density_variance.*variability) + mean_density;
    
    D = 1 - M;  % dry firn layer percentage 
    
    M_remainder = M*depth(end);
    D_remainder = D*depth(end);
    T = [];
    G = [];
    N = 0;
    stack_pointer = 1;
    while M_remainder > 0.01 && D_remainder > 0.01
        % Draw a random ice layer and gap thickness that exceed the minimum
        % vertical sampling rate of the density profile
        T1 = (M_remainder-0.01)*rand + 0.01;
        G1 = (D_remainder-0.01)*rand + 0.01;
        N = N + 1;
        % Calculate the remaining vertical distance in which we can add ice
        % layers
        M_remainder = M_remainder - T1;
        D_remainder = D_remainder - G1;
        T = [T T1];
        G = [G G1];
        % Convert depth in meters to samples in the density profile
        ice_samples = round(T1/spacing);
        gap_samples = round(G1/spacing);
        % Generate a random ice layer density profile centered on a density
        % of 0.873 g/cm^3 based on observations in Machguth (2016)
        rho_ice = 0.005*randn(1,ice_samples+1) + 0.873; 
        % Make sure not samples exceed the density of solid ice
        ind = find(rho_ice > 0.917);
        rho_ice(ind) = 0.917;
        % Add ice layer to density profile and move stack pointer to below
        % the gap
        density(stack_pointer:stack_pointer+ice_samples) = rho_ice;
        stack_pointer = stack_pointer + ice_samples + gap_samples;
    end
    
    % If we have not yet met the desired MFP after layer tiling, add a
    % single ice layer of the thickness needed to meet the MFP
    T = [T M_remainder];
    ice_samples = round(T(end)/spacing);
    if ice_samples ~= 0
        rho_ice = 0.005*randn(1,ice_samples+1) + 0.873;
        ind = find(rho_ice > 0.917);
        rho_ice(ind) = 0.917;
        density(stack_pointer:stack_pointer+ice_samples) = rho_ice;
        N = N + 1;
    end
    
    % Split the density profile at a random point and swap halves to
    % further randomize the layer thickness and positioning
    if ~isempty(T)
        split_sample = randi(length(depth));
        density = [density(split_sample:end) density(1:split_sample-1)];
    end
   
end