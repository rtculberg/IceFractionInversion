%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This interactive script calculates attenuation profiles based on firn
% core conductivity measurements from the North Greenland Traverse and
% modeled temperatures from MARv3.5.2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load conductivity data

clear;

addpath('RefData');

% Load conductivity measurements from North Greenland Traverse Firn Cores
% (original data available at pangea.de)
core1 = load('B26_conductivity.mat');
core2 = load('B27_conductivity.mat');
core3 = load('B28_conductivity.mat');
core4 = load('B29_conductivity.mat');

%% Interpolate mean conductivity profile

% Set up a common depth axis with 5mm sampling and 50 m depth cutoff
d = 0:0.005:50;

% Extract unique depth measurements from each core
[~, ind1, ~] = unique(core1.conductivity(:,1));
[~, ind2, ~] = unique(core2.conductivity(:,1));
[~, ind3, ~] = unique(core3.conductivity(:,1));
[~, ind4, ~] = unique(core4.conductivity(:,1));

% Interpolate conductivity to common depth axis
c1_interp = interp1(core1.conductivity(ind1,1), core1.conductivity(ind1,2), d);
c2_interp = interp1(core2.conductivity(ind2,1), core2.conductivity(ind2,2), d);
c3_interp = interp1(core3.conductivity(ind3,1), core3.conductivity(ind3,2), d);
c4_interp = interp1(core4.conductivity(ind4,1), core4.conductivity(ind4,2), d);

% Calculate mean conductivity across all for cores and convert units to
% microSiemens/meter
mean_cond = 1e-6*(0.25*(c1_interp+c2_interp+c3_interp+c4_interp));


%% Check Temperature Profile from MAR

% Latitude and longitude of double ridge
lat = 74.56234505;
lon = -54.05583724;

% Convert to x,y coordinates in EPSG 3413
[x,y] = ll2psn(lat, lon);

%% Load Variables from MAR Data

% Load data
load('./RefData/MARTemperatureProfile.mat');

% Interpolate to common depth axis
temp_interp = interp1(depth, temp_profile, d);

% Load the 1980-2011 mean annual temperature data and extract the mean
% annual temperature at the double ridge location
averages = load('MAR_20km_MeanVals_1980-2011.mat');
[mx, my] = ll2psn(averages.mlat, averages.mlon);

opt = abs(mx - x).^2 + abs(my - y).^2;
[value, index] = min(opt(:));
[row, col] = ind2sub(size(opt), index);

mean_temp = averages.mean_temp(row,col);

%% Temperature Correct the Conductivity

% Source of the Dome C temperature correction from Wolff (1999) which 
% references Wolff (1995) for specific procedures
% Wolff (1995) Long-term changes in the acid and salt concentrations of the 
% Greenland Ice Core Project ice core from electrical stratigraphy
% Wolff (1999) Comparison of Holocene Electrial records from Dome C and
% Vostok, Antarctica

E_p = 0.5;             % Wolff, 1995
E_i = 0.22;            % Wolff, 1995
R1 = 8.62173303e-5;
T0 = -15 + 273.15;     % Wolff, 1995
T1 = temp_interp + 273.15;                            % MAR temperature profile correction
T2 = (mean_temp + 273.15)*(ones(size(temp_interp)));  % Mean annual temperature correction
T3 = 273.15*(ones(size(temp_interp)));                % 0C temperature correction
sigma_pure = 6.1e-6;

% Arrhenius correction from MacGregor, 2007

% MAR temperature profile correction
sigma_p1 = sigma_pure./exp((E_p./R1).*((1./T1) - (1./T0)));                  % background pure ice component
sigma_i1 = (mean_cond - sigma_pure)./exp((E_i./R1).*((1./T1) - (1./T0)));      % impurity components
rescaled1 = sigma_p1 + sigma_i1;

% Mean annual temperature correction
sigma_p2 = sigma_pure./exp((E_p./R1).*((1./T2) - (1./T0)));                  % background pure ice component
sigma_i2 = (mean_cond - sigma_pure)./exp((E_i./R1).*((1./T2) - (1./T0)));      % impurity components
rescaled2 = sigma_p2 + sigma_i2;

% 0C temperature correction
sigma_p3 = sigma_pure./exp((E_p./R1).*((1./T3) - (1./T0)));                  % background pure ice component
sigma_i3 = (mean_cond - sigma_pure)./exp((E_i./R1).*((1./T3) - (1./T0)));      % impurity components
rescaled3 = sigma_p3 + sigma_i3;


%% Calculate Attenuation

e0 = 8.854e-12;          % Permittivity of a vacuum
e_real = 3.15;           % Real component of permittivity
mu = 1.26e-6;            % Magnetic permeability of a vacuum

% Fill in data gap near surface with first data point from subsurface
rescaled1(1:49) = rescaled1(50);
rescaled2(1:49) = rescaled2(50);
rescaled3(1:49) = rescaled3(50);

% Calculate attenuation from temperature corrected conductivity following
% Equation 1 in Zirizzotti, et al "Electromagnetic ice absorption rate at 
% Dome C, Antarctica" (2014)
A1 = 8.686.*sqrt(mu/(e0*e_real)).*(rescaled1./2);
A2 = 8.686.*sqrt(mu/(e0*e_real)).*(rescaled2./2);
A3 = 8.686.*sqrt(mu/(e0*e_real)).*(rescaled3./2);
del_x = horzcat(0.005, diff(d));

% Integrate over profile to find total attenuation as a function of depth

atten1 = zeros(size(A1));
atten1(1) = 2*A1(1)*del_x(1);
for p = 2:length(A1)
    atten1(p) = atten1(p-1) + 2*A1(p)*del_x(p);
end

atten2 = zeros(size(A2));
atten2(1) = 2*A2(1)*del_x(1);
for p = 2:length(A2)
    atten2(p) = atten2(p-1) + 2*A2(p)*del_x(p);
end

atten3 = zeros(size(A3));
atten3(1) = 2*A3(1)*del_x(1);
for p = 2:length(A3)
    atten3(p) = atten3(p-1) + 2*A3(p)*del_x(p);
end

%% Plot different total attenuation plots

figure;
plot(d, atten1, 'LineWidth', 2);
hold on;
plot(d, atten2, 'LineWidth', 2);
hold on;
plot(d, atten3, 'LineWidth', 2);
xlabel('Depth (m)');
ylabel('Total Attenuation (dB)');
legend('MAR Temperature', 'Mean Annual Temperature', 'Temperate Ice');
set(gca, 'FontSize', 15, 'FontWeight', 'bold');

