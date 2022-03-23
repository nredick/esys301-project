% Constants 
a = 0.3; % albedo; no units 
sigma = 5.67*10^(-8); % stefan-boltzmann const; W/m^2
rho_w = 1000; % water density; kg/m^3
rho_a = 1.225; % air density; kg/m^3
Cp_w = 4184; % specific heat of water; J/kgK
Cp_a = 1003.5; % specific heat of air; J/kgK
H_ml = 100; % height of mixing layer (active) of ocean; m
H_a = 10000; % height of the atmosphere; m
k = 5.55*10^(-5); % m/s
S_ocn = 34; % salinity of the ocean; ppm
P_atm = 10^5; % atmospheric pressure; Pa
CO2_s = 10.5; % dissolved CO2 in the surface ocean; umol/kg
 
% Timesteps and stuff
dt = 10^(-5); % timestep in seconds *CHANGE THIS TO WHAT IT ACTUALLY IS
Ntot = 0; % number of timesteps needed *CHANGE THIS TO WHAT IT ACTUALLY IS
 
% Conversion factors
s2y = 31536000; % seconds to years
GT2umolperkg = GT2g*mol2umol / (mu_CO2 * rhow * Aearth * Hocn);  % umol/kg to GT conversion umol/kg/GT
     
% something for umol/kg to ppm
 
% Parameters we will change
F_CO2 = 0; % Amount of CO2 released by eruption *CHANGE TO REAL VALUE
F_a = 0; % Amount of aerosols released by eruption *CHANGE TO REAL VALUE
% something for how explosive it is which will change the time constant
 
% Initialize variables, preallocate
T_e = nan(1, Ntot+1); % K
T_a = nan(1, Ntot+1); % K
CO2_atm = nan(1, Ntot+1); % ppm?
M_a = nan(1, Ntot+1); % kg
 
time = nan(1, Ntot+1); % s
e_LW = nan(1, Ntot+1);
S_o = nan(1, Ntot+1); % W/m^2
 
% Define initial conditions
T_e(1) = 280;
T_a(1) = 230;
CO2_atm(1) = 280;
M_a(1) = 0;
 
time(1) = 0;
e_LW(1) = 0.8;
S_o = 1367; 
 
% Time step
 
% Equations:
e_LW(t) = 0.8*(1 + 0.054*log(CO2_atm(t)/280));
S_o = 1367 - 1.25*10^(-10)*M_a(t);
 
% Plot
