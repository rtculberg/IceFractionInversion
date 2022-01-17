%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% Subroutine for Markov Chain Monte Carlo inversion with Metropolis-Hastings 
% sampling with a Gaussian proposal distribution and Gaussian likelihood
% function
% 
% Input Variables:
% data - mean radar reflectivity in dB (value to be inverted) (scalar double)
% x0 - vector of initial model parameters
% xstep - vector of step sizes for model parameters
% xbnds - matrix of upper and lower bounds for model parameters
% sigma - standard deviation of radar reflectivity distribution in dB
% (uncertainty in value to be inverted) (scalar double)
% Niter - number of iterations for MCMC run (scalar double)
% coefficients - ARMA model coefficients for density emulation
% depth - vector of depths within the range bin for density simulation
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

function [x_keep, L_keep, count, reflectivity, accept] = mcmc_gauss_swapfirn(data,x0,xstep,xbnds,sigma,Niter,coefficients,depth,filename)

    h = 435;         % average flight altitude
    n_ice = 1.78;    % index of refraction in solid ice
    fc = 300e6;      % radar center frequency
    lambda = 3e8/(fc*n_ice);                % radar wavelength in ice
    A_tot = 0.5*pi*h*lambda;    % illuminated area
    A_min = pi*(2*lambda).^2;               % minimum patch area

    Fs = 3e9;           % Sampling frequency (Hz)
    Fc = 0;             % baseband frequency
    fc = 300e6;         % radar center frequency
    BW = 300e6;         % radar bandwidth
    tau = 1e-6;         % radar pulse length
    slope = BW/tau;     % Chirp slope

    % Generate transmit chirp
    chirp = SingleSideBand(slope, tau, Fs, Fc);
    time_window = tukeywin(length(chirp),0.2);
    chirp = chirp.*time_window';

    % Set up frequency window
    freq_axis = linspace(0, Fs/2, ceil(length(chirp)/2));
    [~, ind] = min(abs(freq_axis - BW));
    freq_win = hann(ind)';
    freq_win = [freq_win zeros(1,length(chirp) - ind)];

    % Frequency spectrum of transmitted chirp
    ref = fft(chirp, length(chirp));

    % Initial storage variable for MCMC run
    count = 0;
    L_keep = zeros(Niter,1);
    x_keep = zeros(Niter+1,length(x0));
    reflectivity = zeros(1,Niter);
    
    % Run radar model on initial parameters as starting point for MCMC
    % random walk
    % Choose the number and area of firn patches
    [A, ~] = DeterministicTilePatches(x0(4), A_tot, A_min);
    R = zeros(size(A));
    for m = 1:length(A)
        % Generate a random density profile for each firn patch
        [density, ~, ~, ~, ~] = DensityEmulator(x0(3), x0(1), coefficients, x0(2), depth);
        % Convert density to refractive index with Kovacs relationship
        n = 1 + 0.845*density;
        % Calculate effective reflectivity for each stack of density layers
        R(m) = ReflectionCoefficients_UWB(fc, BW, n, depth, ref, freq_win, freq_axis);
    end
    % Calculate background area
    A = [A 1 - sum(A)];
    % Generate a density profile for the background firn
    [density, ~, ~, ~, ~] = DensityEmulator(0, x0(1), coefficients, x0(2), depth);
    % Convert density to refractive index with Kovacs relationship
    n = 1 + 0.845*density;
    % Calculate effective reflectivity for each stack of density layers
    R(m+1) = ReflectionCoefficients_UWB(fc, BW, n, depth, ref, freq_win, freq_axis);
    % Calculate the effective reflectivity of radar footprint as the
    % coherent, area weighted sum of all the reflectors within the
    % footprint
    reflectivity(1) = 10*log10(abs(sum(A.*R)).^2);
    % Save the model parameters
    x_keep(1,:) = x0;
    
    % Run the desired number of iternations in the random walk
    for k = 1:Niter
        
        % Save an intermediate results file every 10,000 iterations
        if mod(k,10000) == 0
            save(strcat('./DerivedData/',filename), 'x_keep', 'L_keep', 'count', 'reflectivity');
        end
        
        % Print the percent completion every 1000 iterations
        if mod(k,1000) == 0
            fprintf('%f percent\n', 100*(k/Niter));
        end
        
        % Propose the next set of model parameters drawing from a
        % multivariate gaussian distribution whose mean is set by the last
        % set of accepted model parameters and whose standard deviation is
        % set by the model parameter step size
        xprop = mvnrnd(x_keep(k,:),diag(xstep.^2),1);
        flag = 0;
        % Check that the proposed model parameters are within the
        % predefined bounds
        for m = 1:length(xprop)
            if xprop(m) <= xbnds(1,m) || xprop(m) >= xbnds(2,m)
                flag = 1;
            end
        end
        
        if flag == 0  % if proposed model is within the bounds,
           
            % Get current reflectivity
            dcurr = reflectivity(k);
            
            % Run radar model with newly proposed model parameters
            [A, ~] = DeterministicTilePatches(xprop(4), A_tot, A_min);
            R = zeros(size(A));
            for m = 1:length(A)
                [density, ~, ~, ~, ~] = DensityEmulator(xprop(3), xprop(1), coefficients, xprop(2), depth);
                n = 1 + 0.845*density;
                R(m) = ReflectionCoefficients_UWB(fc, BW, n, depth, ref, freq_win, freq_axis);
            end
            A = [A 1 - sum(A)];
            [density, ~, ~, ~, ~] = DensityEmulator(0, xprop(1), coefficients, xprop(2), depth);
            n = 1 + 0.845*density;
            R(m+1) = ReflectionCoefficients_UWB(fc, BW, n, depth, ref, freq_win, freq_axis);
            dprop = 10*log10(abs(sum(A.*R)).^2);
            
            % Calculate the likelihood ratio of the newly proposed model
            p_d_xprop = -0.5*sum(((data - dprop).^2)/(sigma^2));
            p_d_x = -0.5*sum(((data - dcurr).^2)/(sigma^2));
            ratio = p_d_xprop - p_d_x;
            % If new model is more likely than previous model, accept
            if ratio > 0
                count = count + 1;
                L_keep(k,:) = p_d_xprop;
                x_keep(k+1,:) = xprop;
                reflectivity(k+1) = dprop;
            else
                % If new model is less likely, accept randomly based on
                % whether ratio exceeds a random number between 0 and 1
                u = log(rand(1,1));
                if ratio > u
                    count = count + 1;
                    L_keep(k,:) = p_d_xprop;
                    x_keep(k+1,:) = xprop;
                    reflectivity(k+1) = dprop;
                else
                    x_keep(k+1,:) = x_keep(k,:);
                    L_keep(k,:) = p_d_x;
                    reflectivity(k+1) = reflectivity(k);
                end
            end
        else  
            % if proposed model is outside the bounds, keep the previous
            % model
            x_keep(k+1,:) = x_keep(k,:);
            reflectivity(k+1) = reflectivity(k);
            dcurr = reflectivity(k);
            L_keep(k,:) = -0.5*sum(((data - dcurr).^2)/(sigma^2));
        end
    end
    accept = count/Niter;
    fprintf("Acceptance Ratio: %f\n", accept);
end