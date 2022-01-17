%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This function runs a Markov Chain Monte Carlo inversion with
% Metropolis-Hastings sampling to invert radar reflectivity in terms of
% firn density, firn density variability, melt feature percentage, and
% fractional firn area in the radar footprint. 
% 
% Input Variables:
% data - mean radar reflectivity in dB (value to be inverted) (scalar double)
% sigma - standard deviation of radar reflectivity distribution in dB
% (uncertainty in value to be inverted) (scalar double)
% Niter - number of iterations for MCMC run (scalar double)
% filename - name of output file to save results (string)
% 
% Output Variables:
% x_keep - 4 x Niter matrix of accepted model parameters
% L_keep - 1 x Niter vector of the likelihood of the accepted model
% parameters at that step
% count - number of accepted proposals (count/Niter = acceptance rate)
% reflectivity - 1 x Niter vector of modeled reflectivity from accepted
% parameters (can be compared to data distributions to sanity check the
% model runs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function RunFacetMCMC(data, sigma, Niter, filename)

addpath('RefData');
addpath('DemoData');
addpath('ReferenceFunctions');
addpath('SensitivityTests');

depth = 0:0.005:0.43;                % bounded by range resolution of system             
load('DensityCoefficients.mat');     % ARMA coefficients for density model
coefficients.Constant = 3.0051e-5;   % ARMA constant
coefficients.Variance = 0.0106;      % ARMA variance

% Starting model parameters for MCMC
% x0(1) = mean firn density in g/cm^3
% x0(2) = firn density variability in g/cm^3
% x0(3) = melt feature percentage (MFP) - fraction of range resolution
% occupied by refrozen ice, value from 0 to 1
% x0(4) = fractional firn area - fraction of horizontal footprint covered
% by relict firn patches, value from 0 to 1
x0 = [0.6 0.02 0.5 0.5]; 
% Range of acceptable values for model parameters
xbnds = [[0.45 0.01 0 0]; [0.75 0.04 1 1]];  
% Step size for random walk
xstep = [0.01 0.001 0.03 0.03];

% Run MCMC inversion
% Edit this function call to run the sensitivity tests (mcmc_gauss_fresnel,
% mcmc_gauss_phaseerror, mcmc_gauss_swap_firn)
[x_keep, L_keep, count, reflectivity] = mcmc_gauss(data,x0,xstep,xbnds,sigma,Niter,coefficients,depth,filename);

% Save final inversion results in the results directory
save(strcat('./DerivedData/', filename), 'data', 'sigma', 'Niter', 'x_keep', 'L_keep', 'count', 'reflectivity');

end