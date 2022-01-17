%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This script calculates the reflectivity distributions for various
% material types based on predefined portions of a radargram. The outputs
% of this interactive script can be used to run the MCMC inversion from
% RunFacetMCMC.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Data

addpath('RefData');
addpath('DemoData');

clear;

% Data directory and loading - update these paths to wherever the radargram
% under analysis is stored
I = ncread('D:\Data\Greenland\FirnStructure\IceBlobs\DoubleRidgeData\SAR_20160510_01_037.nc', 'Data_I');
Q = ncread('D:\Data\Greenland\FirnStructure\IceBlobs\DoubleRidgeData\SAR_20160510_01_037.nc', 'Data_Q');
new.Data = I + 1i*Q;
new.Time = ncread('D:\Data\Greenland\FirnStructure\IceBlobs\DoubleRidgeData\SAR_20160510_01_037.nc', 'Time');
new.Surface = ncread('D:\Data\Greenland\FirnStructure\IceBlobs\DoubleRidgeData\SAR_20160510_01_037.nc', 'Surface');
new.Latitude = ncread('D:\Data\Greenland\FirnStructure\IceBlobs\DoubleRidgeData\SAR_20160510_01_037.nc', 'Latitude');
new.Longitude = ncread('D:\Data\Greenland\FirnStructure\IceBlobs\DoubleRidgeData\SAR_20160510_01_037.nc', 'Longitude');

% Physical constants
c = 299792458;   % speed of light in a vacuum 
n_ice = 1.78;     % refractive index of ice 

%% Correct Data for Geometric Spreading

% Convert radargram to dB and trim to depth of interest (y pixel limits
% should be updated if a different radargram is being used)
radar = 10*log10(new.Data(950:1300,:));
time = new.Time(950:1300);

% Apply range squared geometric corrections following Equation 1 in Peters,
% et al "Analysis techniques for coherent airborne radar sounding: 
% Application to West Antarctic ice streams" (2005)
h = zeros(1,size(new.Data,2));
for k = 1:size(new.Data,2)
    [~, surf_ind] = min(abs(time - new.Surface(k)));
    for m = 1:size(radar,1)
        if m <= surf_ind
            h(k) = time(m)*0.5*c;
            radar(m,k) = radar(m,k) + 20*log10(2*h(k));
        else
            h(k) = time(surf_ind)*0.5*c;
            z = (time(m) - time(surf_ind))*0.5*(c/n_ice);
            radar(m,k) = radar(m,k) + 20*log10(2*(h(k) + z/n_ice));
        end
    end
end

%% Apply the Absolute Calibration Constant

radar = radar - 16.2632;   

figure;
imagesc(radar);
colorbar;
cmocean('matter');

%% Correct Data for Attenuation

load('MARAttenuation.mat');

% Alternate attenuation corrections that assume a constant temperature
% profile of either -14 C or 0 C
% load('MaxAttenuation.mat');
% load('MinAttenuation.mat');

for k = 1:size(new.Data,2)
    [~, surf_ind] = min(abs(time - new.Surface(k)));
    for m = 1:size(radar,1)
        if m > surf_ind
            z = (time(m) - time(surf_ind))*0.5*(c/n_ice);
            [~, ind] = min(abs(d - z));
            radar(m,k) = radar(m,k) + atten(ind);
        end
    end
end

figure;
imagesc(radar);
colorbar;
cmocean('matter');

%% Load and Interpolate Boundaries for Sub-Ridge Material

load('PlugIceBoundaries.mat');

start = max(min(top_x), min(bot_x));
stop =  min(max(top_x), max(bot_x));
x_plug = ceil(start):1:floor(stop);

uy_plug = ceil(interp1(top_x, top_y, x_plug));
by_plug = floor(interp1(bot_x, bot_y, x_plug));

plug_ice = [];
for k = 1:length(x_plug)
    plug_ice = [plug_ice radar(uy_plug(k):by_plug(k),x_plug(k))'];
end

%% Load and Interpolate Boundaries for Firn Material

load('FirnBoundaries.mat');

start = max(min(top_x), min(bot_x));
stop =  min(max(top_x), max(bot_x));
x_firn = ceil(start):1:floor(stop);

uy_firn = ceil(interp1(top_x, top_y, x_firn));
by_firn = floor(interp1(bot_x, bot_y, x_firn));

firn = [];
for k = 1:length(x_firn)
    firn = [firn radar(uy_firn(k):by_firn(k),x_firn(k))'];
end

load('FirnBoundaries2.mat');

start = max(min(top_x), min(bot_x));
stop =  min(max(top_x), max(bot_x));
x_firn2 = ceil(start):1:floor(stop);

uy_firn2 = ceil(interp1(top_x, top_y, x_firn2));
by_firn2 = floor(interp1(bot_x, bot_y, x_firn2));

for k = 1:length(x_firn2)
    firn = [firn radar(uy_firn2(k):by_firn2(k),x_firn2(k))'];
end

%% Load and Interpolate Boundaries for Refrozen Ice Shell

load('IceBoundaries.mat');

start = max(min(top_x), min(bot_x));
stop =  min(max(top_x), max(bot_x));
x_ice = ceil(start):1:floor(stop);

uy_ice = ceil(interp1(top_x, top_y, x_ice));
by_ice = floor(interp1(bot_x, bot_y, x_ice));

ice = [];
for k = 1:length(x_ice)
    ice = [ice radar(uy_ice(k):by_ice(k),x_ice(k))'];
end

load('IceBoundaries2.mat');

start = max(min(top_x), min(bot_x));
stop =  min(max(top_x), max(bot_x));
x_ice2 = ceil(start):1:floor(stop);

uy_ice2 = ceil(interp1(top_x, top_y, x_ice2));
by_ice2 = floor(interp1(bot_x, bot_y, x_ice2));

for k = 1:length(x_ice2)
    ice = [ice radar(uy_ice2(k):by_ice2(k),x_ice2(k))'];
end

%% Plot the Selected Regions Each Type of Material on Radargram

[at, ~, ~, ~] = geodetic_to_along_track(new.Latitude, new.Longitude, []);

figure;
imagesc(radar);
hold on;
patch([x_ice fliplr(x_ice)], [by_ice fliplr(uy_ice)], 'r', 'FaceColor', 'none', 'EdgeColor', (1/256)*[238 178 32], 'LineWidth', 1.5);
hold on;
patch([x_plug fliplr(x_plug)], [by_plug fliplr(uy_plug)], 'r', 'FaceColor', 'none', 'EdgeColor', (1/256)*[126 47 142], 'LineWidth', 1.5);
hold on;
patch([x_firn fliplr(x_firn)], [by_firn fliplr(uy_firn)], 'r', 'FaceColor', 'none', 'EdgeColor', (1/256)*[218 83 25], 'LineWidth', 1.5);
hold on;
patch([x_ice2 fliplr(x_ice2)], [by_ice2 fliplr(uy_ice2)], 'r', 'FaceColor', 'none', 'EdgeColor', (1/256)*[238 178 32], 'LineWidth', 1.5);
hold on;
patch([x_firn2 fliplr(x_firn2)], [by_firn2 fliplr(uy_firn2)], 'r', 'FaceColor', 'none', 'EdgeColor', (1/256)*[218 83 25], 'LineWidth', 1.5);
colormap gray;
caxis([-55 -10]);
h = colorbar;
xlabel('Along-track Sample');
ylabel('Fast-time Sample');
set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
set(get(h, 'label'), 'string', 'Calibrated Power (dB)');

%% Plot the Log Power Histograms for Each Region

% Extract peak powers within each region (to avoid skewing statistics with
% samples taken on the rising or falling edge of a reflection peak)
[firn_peaks, firn_locs] = findpeaks(firn);
[ice_peaks, ice_locs]  = findpeaks(ice);
[plug_ice_peaks, plug_ice_locs] = findpeaks(plug_ice);

figure;
histogram(firn_peaks, 35, 'Normalization', 'pdf', 'FaceColor',(1/256)*[218 83 25]);
hold on;
histogram(ice_peaks, 35, 'Normalization', 'pdf', 'FaceColor',(1/256)*[238 178 32]);
hold on;
histogram(plug_ice_peaks, 35, 'Normalization', 'pdf', 'FaceColor',(1/256)*[126 47 142]);
legend('Firn', 'Refrozen Ice Shell', 'Sub-Ridge Material');
xlabel('Calibrated Power (dB)');
ylabel('PDF');
set(gca, 'FontSize', 15, 'FontWeight', 'bold');

%% Descriptive Statistics

fprintf('Firn Mean: %f\n', mean(firn_peaks));
fprintf('Firn STD: %f\n', std(firn_peaks));
fprintf('Firn Kurtosis: %f\n', kurtosis(firn_peaks));
fprintf('Firn Skewness: %f\n', skewness(firn_peaks));

fprintf('Ice Mean: %f\n', mean(ice_peaks));
fprintf('Ice STD: %f\n', std(ice_peaks));
fprintf('Ice Kurtosis: %f\n', kurtosis(ice_peaks));
fprintf('Ice Skewness: %f\n', skewness(ice_peaks));

fprintf('Ridge Mean: %f\n', mean(plug_ice_peaks));
fprintf('Ridge STD: %f\n', std(plug_ice_peaks));
fprintf('Ridge Kurtosis: %f\n', kurtosis(plug_ice_peaks));
fprintf('Ridge Skewness: %f\n', skewness(plug_ice_peaks));


%% Bootstrap Confidence Intervals on Mean

% This portion of the script can be used to calculate the 95% confidence
% interval on the mean reflectivities listed above. This was used to
% confirm that the different attenuation corrections (line 71) do not lead 
% to statistically different distributions of reflectivity for each
% material type. 

num_trials = 10000;
firn_avg = zeros(1, num_trials);
ice_avg = zeros(1, num_trials);
plug_avg = zeros(1, num_trials);
for k = 1:num_trials
    for m = 1:length(firn_peaks)
        pick = randi(length(firn_peaks));
        firn_avg(k) = firn_avg(k) + firn_peaks(pick);
    end
    firn_avg(k) = firn_avg(k)/length(firn_peaks);
    
    for m = 1:length(ice_peaks)
        pick = randi(length(ice_peaks));
        ice_avg(k) = ice_avg(k) + ice_peaks(pick);
    end
    ice_avg(k) = ice_avg(k)/length(ice_peaks);
    
    for m = 1:length(plug_ice_peaks)
        pick = randi(length(plug_ice_peaks));
        plug_avg(k) = plug_avg(k) + plug_ice_peaks(pick);
    end
    plug_avg(k) = plug_avg(k)/length(plug_ice_peaks);
end

figure;
histogram(firn_avg, 'Normalization', 'pdf', 'FaceColor',(1/256)*[218 83 25]);
hold on;
histogram(ice_avg, 'Normalization', 'pdf', 'FaceColor',(1/256)*[238 178 32]);
hold on;
histogram(plug_avg, 'Normalization', 'pdf', 'FaceColor',(1/256)*[126 47 142]);
legend('Firn', 'Blob Ice', 'Plug Ice');
xlabel('Power (dB)');
ylabel('Probability');
title('Bootstrap Confidence Intervals for Mean');
set(gca, 'FontSize', 15, 'FontWeight', 'bold');

fprintf('Firn 2.5: %f\n', prctile(firn_avg, 2.5));
fprintf('Firn 97.5: %f\n', prctile(firn_avg, 97.5));
fprintf('Ice 2.5: %f\n', prctile(ice_avg, 2.5));
fprintf('Ice 97.5: %f\n', prctile(ice_avg, 97.5));
fprintf('Ridge 2.5: %f\n', prctile(plug_avg, 2.5));
fprintf('Ridge 97.5: %f\n', prctile(plug_avg, 97.5));
