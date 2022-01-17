%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radar Horizontal Phase Coherence Calculations
%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This interactive scripts calculates the horizontal phase coherence of
% pulse compressed radar data at varying length scales and finds the mean
% coherence of different material types (ice, firn, and noise/clutter). 
% 
% Input Data:
% pc - pulse compressed radargram with complex voltage (must not be power
% detected!)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Data and Set Variables

clear;

addpath('RefData');

% Load the pulse compressed complex radargram to be analyzed
% Update the file path to point to the location where the radargram of
% interest is stored
I = ncread('.\DempData\PC_20160510_01_037.nc', 'Data_I');
Q = ncread('.\DempData\PC_20160510_01_037.nc', 'Data_Q');
pc.Data = I + 1i*Q;
pc.Surface = ncread('.\DempData\PC_20160510_01_037.nc', 'Surface');
pc.Time = ncread('.\DempData\PC_20160510_01_037.nc', 'Time');
pc.Latitude = ncread('.\DempData\PC_20160510_01_037.nc', 'Latitude');
pc.Longitude = ncread('.\DempData\PC_20160510_01_037.nc', 'Longitude');

% Constant variables
c = 299792458;           % speed of light in vacuum
fc = 300e6;              % radar center frequency in Hz

%% Radargram Motion Compensation and Flatenning

lambda = c/fc;   % radar wavelength

phase_corrected = zeros(size(pc.Data));
for k = 1:length(pc.Surface)
    [~, surf_ind] = min(abs(pc.Time - pc.Surface(k)));
    seg = pc.Data(surf_ind:end,k);
    % Find surface elevation of current trace relative to first trace
    % (accounts for changes in flight altitude and surface topography)
    del_clearance = 0.5*c*(pc.Surface(k) - pc.Surface(1));
    % Apply phase correction to trace based on change in flight clearance
    del_phi = ((4*pi)/lambda)*del_clearance;
    phase_corrected(1:length(seg),k) = seg.*exp(1i*del_phi);
end

% Plot motion compensated and flattened radargram
figure;
imagesc(10*log10(abs(phase_corrected).^2));

%% Pick the Boundaries for Materials of Interest

% % These sections allow you to manually pick the regions of the radargram
% % representing materials of interest - in this case, ice, firn, and noise
% % To run this script with pre-defined boundaries provided in the demo data,
% % skip to line 155
% 
% % Pick upper ice boundary
% figure;
% imagesc(10*log10(abs(phase_corrected).^2));
% ylim([0 150]);
% caxis([-110 -40]);
% [x1,y1] = getpts;
% 
% % Interpolate to trace spacing
% x = 1:1:size(phase_corrected,2);
% x1 = [x(1) x1' x(end)];
% y1 = [y1(1) y1' y1(end)];
% y_upper = interp1(x1, y1, x);
% 
% % Pick lower ice boundary
% figure;
% imagesc(10*log10(abs(phase_corrected).^2));
% ylim([0 150]);
% caxis([-110 -40]);
% [x2,y2] = getpts;
% 
% % Interpolate to trace spacing
% x2 = [x(1) x2' x(end)];
% y2 = [y2(1) y2' y2(end)];
% y_lower = interp1(x2, y2, x);
% 
% % Save ice boundaries
% ice.upper = y_upper;
% ice.lower = y_lower;
% 
% % Pick upper firn boundaries
% imagesc(10*log10(abs(phase_corrected).^2));
% ylim([0 150]);
% caxis([-110 -40]);
% [x1,y1] = getpts;
% 
% % Interpolate to trace spacing
% x = [1:1:4767 9880:1:size(phase_corrected,2)];
% x1 = [x(1) x1' x(end)];
% y1 = [y1(1) y1' y1(end)];
% y_upper = interp1(x1, y1, x);
% 
% % Pick lower firn boundaries
% figure;
% imagesc(10*log10(abs(phase_corrected).^2));
% ylim([0 150]);
% caxis([-110 -40]);
% [x2,y2] = getpts;
% 
% % Interpolate to trace spacing
% x2 = [x(1) x2' x(end)];
% y2 = [y2(1) y2' y2(end)];
% y_lower = interp1(x2, y2, x);
% 
% % Save firn boundaries
% firn.upper = y_upper;
% firn.lower = y_lower;
% 
% % Pick upper noise boundary
% figure;
% imagesc(10*log10(abs(phase_corrected).^2));
% ylim([0 500]);
% caxis([-110 -40]);
% [x1,y1] = getpts;
% 
% % Interpolate to trace spacing
% x = 1:1:size(phase_corrected,2);
% x1 = [x(1) x1' x(end)];
% y1 = [y1(1) y1' y1(end)];
% y_upper = interp1(x1, y1, x);
% 
% % Pick lower noise boundary
% figure;
% imagesc(10*log10(abs(phase_corrected).^2));
% ylim([0 500]);
% caxis([-110 -40]);
% [x2,y2] = getpts;
% 
% % Interpolate to trace spacing
% x2 = [x(1) x2' x(end)];
% y2 = [y2(1) y2' y2(end)];
% y_lower = interp1(x2, y2, x);
% 
% % Save noise boundaries
% noise.upper = y_upper;
% noise.lower = y_lower;

%% Total Coherence Metric for Image

lambda = c/fc;                    % radar wavelength
rf = 2*sqrt((lambda*500)/2);      % radar fresnel zone diameter

% Along-track spacing
[at, ~, ~, ~] = geodetic_to_along_track(pc.Latitude, pc.Longitude, []);
spacing = mean(diff(at));

% Length scale over which to calculate phase coherence - edit this number
% to change
length_scale = rf;
stacks = round(length_scale./spacing);  % convert meters to trace count

% Calculate phase coherence for every horizontal block at length scale of
% interest
coherence = zeros(size(phase_corrected));
for m = 1:size(phase_corrected,2)-stacks
    seg = phase_corrected(:,m:m+stacks);
    coherence(:,m) = (abs(sum(seg,2)))./sum(abs(seg),2);
end

%% Plots Image Coherence 

[at, ~, ~, ~] = geodetic_to_along_track(pc.Latitude, pc.Longitude, []);
depth = 0.5*(1:1:size(phase_corrected,1))*mean(diff(pc.Time))*(c/1.78);

figure;
imagesc(at/1000, depth, coherence);
cmocean('thermal');
h = colorbar;
ylim([0 90]);
title('Coherence - Fresnel Zone');
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
set(get(h, 'label'), 'string', 'Coherence');
xlabel('Distance Along-Track (km)');
ylabel('Approximate Depth (m)');

%% Histograms of Coherence

% Load predefined boundaries for ice, firn, and noise in teh demo image
load('IceCoherenceBounds.mat');
load('FirnCoherenceBounds.mat');
load('NoiseCoherenceBounds.mat');

% Extract image coherence within the boundaries of each material type
ice_coh = [];
for k = 1:length(ice.upper)
    ice_coh = [ice_coh coherence(round(ice.upper(k)):round(ice.lower(k)),k)'];
end

firn_coh = [];
for k = 1:length(firn.upper)
    firn_coh = [firn_coh coherence(round(firn.upper(k)):round(firn.lower(k)),k)'];
end

noise_coh = [];
for k = 1:length(noise.upper)
    noise_coh = [noise_coh coherence(round(noise.upper(k)):round(noise.lower(k)),k)'];
end

% Plot histograms of image coherence within each material type
figure;
histogram(ice_coh);
hold on;
histogram(firn_coh);
hold on;
histogram(noise_coh);

%% Descriptive Statistics

% Print the descriptive statistics for the image coherence of each material
fprintf('Ice Median: %f\n', median(ice_coh));
fprintf('Ice Interquartile Range: %f\n', prctile(ice_coh,75)-prctile(ice_coh,25));
fprintf('Firn Median: %f\n', median(firn_coh));
fprintf('Firn Interquartile Range: %f\n', prctile(firn_coh,75)-prctile(firn_coh,25));
fprintf('Noise Median: %f\n', median(noise_coh));
fprintf('Noise Interquartile Range: %f\n', prctile(noise_coh,75)-prctile(noise_coh,25));

