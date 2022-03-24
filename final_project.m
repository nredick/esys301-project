%% note: work in Gt and convert to ppm for graphing (CO2)

% Constants 
S_o = 1367; % solar constant; W/m^2
sigma = 5.67*10^(-8); % stefan-boltzmann const; Js^-1m^-2K^-4
rho_w = 1000; % water density; kg/m^3
rho_a = 1.225; % air density; kg/m^3
Cp_w = 4184; % specific heat of water; J/kgK
Cp_a = 1003.5; % specific heat of air; J/kgK
H_ml = 100; % height of mixing layer (active) of ocean; m
H_a = 10000; % height of the atmosphere; m
H_ocn = 1000; % deep ocean depth; m
k = 5.55*10^(-5); % piston velocity; m/s
S_ocn = 34; % salinity of the ocean surface; ppm
P_atm = 10^5; % atmospheric pressure; Pa
CO2_s = 10.5; % dissolved CO2 in the surface ocean; umol/kg
A_earth = 4*pi*6371e03^2; %Earth surface area; m^2

% Time steps and stuff
dt = 10^(-5); % timestep in seconds *todo CHANGE THIS TO WHAT IT ACTUALLY IS
Ntot = 0; % number of timesteps needed *todo CHANGE THIS TO WHAT IT ACTUALLY IS
 
% Conversion factors
s2y = 31536000; % seconds to years
% todo: get consts for the following equation
GT2g = 1e15; %Gt to gram conversion; g/Gt
mol2umol = 1e06; %mol to umol conversion; umol/mol
mu_CO2 = 44; %CO2 molecular weight; g/mol

GT2umolperkg = GT2g*mol2umol / (mu_CO2 * rho_w * Aearth * Hocn);  % umol/kg to GT conversion umol/kg/GT
Gt2ppm = 0.1292; % Gt to ppm conversion; Gt/ppmv

% Parameters we will change for experiments 
F_CO2 = 0; % Amount of CO2 released by eruption *CHANGE TO REAL VALUE
F_a = 0; % Amount of aerosols released by eruption *CHANGE TO REAL VALUE
% todo: something for how explosive it is which will change the time constant
 
% Initialize & preallocate arrays
T_e = nan(1, Ntot+1); % temperature of earth; deg K
T_a = nan(1, Ntot+1); % temperature of atmosphere; deg K
CO2_atm = nan(1, Ntot+1); % conc. of CO2 in atmosphere; Gt 
M_a = nan(1, Ntot+1); % mass of aerosols in atmosphere; kg
 
time = nan(1, Ntot+1); % time; s
e_a = nan(1, Ntot+1); % atmospheric emissivity (longwave)
a = nan(1, Ntot+1); % solar constant; W/m^2
 
% Define initial conditions
T_e(1) = 280;
T_a(1) = 230;
CO2_atm(1) = 280;
M_a(1) = 0;

time(1) = 0;
e_a(1) = 0.8;
 
% Time step
 
% Equations
e_a(t) = 0.8*(1 + 0.054*log(CO2_atm(t)/280));
 
% Plot
