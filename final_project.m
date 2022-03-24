% note: work in Gt and convert to ppm for graphing (CO2)

% Constants 
S_o = 1367; % solar constant; W/m^2
sigma = 5.67*10^(-8); % stefan-boltzmann const; Js^-1m^-2K^-4
rho_w = 1000; % water density; kg/m^3
rho_a = 1.225; % air density; kg/m^3
Cp_w = 4184; % specific heat of water; J/kgK
Cp_a = 1003.5; % specific heat of air; J/kgK
H_ml = 100; % height of mixing layer (active) of ocean; m
H_a = 10000; % height of the atmosphere; m
k = 5.55*10^(-5); % piston velocity; m/s
S_ocn = 34; % salinity of the ocean surface; ppm
P_atm = 10^5; % atmospheric pressure; Pa
CO2_s = 10.5; % dissolved CO2 in the surface ocean; umol/kg
 
% Time steps and stuff
dt = 10^(-5); % timestep in seconds *todo CHANGE THIS TO WHAT IT ACTUALLY IS
Ntot = 0; % number of timesteps needed *todo CHANGE THIS TO WHAT IT ACTUALLY IS
 
% Conversion factors
s2y = 31536000; % seconds to years
% todo: get consts for the following equation
GT2umolperkg = GT2g*mol2umol / (mu_CO2 * rhow * Aearth * Hocn);  % umol/kg to GT conversion umol/kg/GT
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
CO2_atm(1) = 280; % ppm 
M_a(1) = 0;
 
time(1) = 0;
e_a(1) = 0.8;
 
% Time step
dt = 0.24 / s2y; % time step (atmospheric temp); secs 
n = 100; % number of time steps 

%% Run model 

for t = 1 : n 
    % earth temp
    T_e(t+1) = dt * ( ((S_o/4)*(1-a(t))+e_a(t)*sigma*(T_a(t)^4)- ...
        sigma*(T_e(t)^4) ) / (rho_w*Cp_w*H_ml)) + T_e(t);

    % emissivity
    e_a(t) = e_a(1) * (1+);

    % time
    time(t+1) = time(t) + dt;
end
 
%% Plot
