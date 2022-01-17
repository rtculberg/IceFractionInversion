%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This interactive script takes the results files from an MCMC inversion
% run (output of FacetModelMCMC.m) and calculates the fractional ice
% content. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
addpath('DemoData');

% Load the results of teh MCMC inversion for each material type
a = load('MCMC_Firn_Final.mat');
b = load('MCMC_Ice_Final.mat');
c = load('MCMC_Ridge_Final.mat');

% Physical and System Constants
h = 435;          % Average aircraft altitude in m
c_light = 299792458;    % Speed of light in a vacuum in m/s
fc = 300e6;       % Radar center frequency in Hz
lambda = fc/c_light;    % Radar wavelength
n_ice = 1.78;     % Refractive index of ice
del_x = 2.5;      % Radar azimuth resolution in m
del_r = (1.53*c_light)/(2*fc*n_ice);  % Radar range resolution in m
r_f = sqrt((lambda*h)/2);       % Radar fresnel zone diameter

% Fractional ice content in firn
% x_keep(:,3) = MFP - percentage of each firn patch filled with ice
% x_keep(:,4) = fractional firn area - percentage of footprint filled with
% firn patches
% del_r*del_x*2*r_f = volume of footprint (range res x az res x cross-track
% res)
v_f = (1-a.x_keep(:,3)).*(a.x_keep(:,4))*del_r*del_x*2*r_f;
v_tot = del_r*del_x*2*r_f;
percent_ice_firn = (v_tot-v_f)./v_tot;

% Fractional ice content in refrozen ice shell
v_f = (1-b.x_keep(:,3)).*(b.x_keep(:,4))*del_r*del_x*2*r_f;
percent_ice_ice = (v_tot-v_f)./v_tot;

% Fractional ice content in sub-ridge material
v_f = (1-c.x_keep(:,3)).*(c.x_keep(:,4))*del_r*del_x*2*r_f;
percent_ice_ridge = (v_tot-v_f)./v_tot;

figure;
histogram(percent_ice_ridge, 200, 'Normalization', 'pdf', 'FaceColor',(1/256)*[118 173 48]);
hold on;
histogram(percent_ice_ice, 200, 'Normalization', 'pdf', 'FaceColor',(1/256)*[238 178 32]);
hold on;
histogram(percent_ice_firn, 200, 'Normalization', 'pdf', 'FaceColor',(1/256)*[0 114 190]);
xlabel('Percent Ice by Volume');
ylabel('PDF');
set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
legend('Ice below ridge', 'Ice next to ridge', 'Firn');

fprintf('Firn Mean: %f\n', mean(percent_ice_firn));
fprintf('Firn Median: %f\n', median(percent_ice_firn));
fprintf('Ice Mean: %f\n', mean(percent_ice_ice));
fprintf('Ice Median: %f\n', median(percent_ice_ice));
fprintf('Ridge Mean: %f\n', mean(percent_ice_ridge));
fprintf('Ridge Median: %f\n', median(percent_ice_ridge));


%% Fractional Ice Content Calculation Method for Swap Firn Version

% This calculates the same results as above, but for the case where the
% footprint background is assumed to be firn rather than ice. 

% Load resuls from firn background inversion sensitivity test
a = load('MCMC_Firn_SwapFirn.mat');
b = load('MCMC_Ice_SwapFirn.mat');
c = load('MCMC_Ridge_SwapFirn.mat');

% Physical and System Constants
h = 435;                % Average aircraft altitude in m
c_light = 299792458;    % Speed of light in a vacuum in m/s
fc = 300e6;             % Radar center frequency in Hz
lambda = fc/c_light;    % Radar wavelength
n_ice = 1.78;           % Refractive index of ice
del_x = 2.5;            % Radar azimuth resolution in m
del_r = (1.53*c_light)/(2*fc*n_ice);  % Radar range resolution in m
r_f = sqrt((lambda*h)/2);             % Radar fresnel zone diameter

% Fractional ice content in firn
% x_keep(:,3) = MFP - percentage of each firn patch filled with ice
% x_keep(:,4) = fractional firn area - percentage of footprint filled with
% firn patches
% del_r*del_x*2*r_f = volume of footprint (range res x az res x cross-track
% res)
v_f = ((1-a.x_keep(:,3)).*a.x_keep(:,4) + (1-a.x_keep(:,4)))*del_r*pi*r_f^2;
v_tot = del_r*pi*r_f^2;
percent_ice1 = 1-(v_f)./v_tot;

% Fractional ice content in refrozen ice shell
v_f = ((1-b.x_keep(:,3)).*b.x_keep(:,4) + (1-b.x_keep(:,4)))*del_r*pi*r_f^2;
percent_ice2 = 1-(v_f)./v_tot;

% Fractional ice content in sub-ridge material
v_f = ((1-c.x_keep(:,3)).*c.x_keep(:,4) + (1-c.x_keep(:,4)))*del_r*pi*r_f^2;
percent_ice3 = 1-(v_f)./v_tot;

figure;
histogram(percent_ice3, 200, 'Normalization', 'pdf', 'FaceColor',(1/256)*[126 47 142]);
hold on;
histogram(percent_ice2, 200, 'Normalization', 'pdf', 'FaceColor',(1/256)*[238 178 32]);
hold on;
histogram(percent_ice1, 200, 'Normalization', 'pdf', 'FaceColor',(1/256)*[218 83 25]);
xlabel('Percent Ice by Volume');
ylabel('PDF');
set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
legend('Ice below ridge', 'Ice next to ridge', 'Firn');

fprintf('Firn Mean: %f\n', mean(percent_ice1));
fprintf('Firn Median: %f\n', median(percent_ice1));
fprintf('Ice Mean: %f\n', mean(percent_ice2));
fprintf('Ice Median: %f\n', median(percent_ice2));
fprintf('Ridge Mean: %f\n', mean(percent_ice3));
fprintf('Ridge Median: %f\n', median(percent_ice3));

toc