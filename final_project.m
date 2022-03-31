%% note: work in Gt and convert to ppm for graphing (CO2)
clear all, close all, clc
% Constants 
S_o = 1367; % solar constant; W/m^2
sigma = 5.67e-8; % stefan-boltzmann const; Js^-1m^-2K^-4
rho_w = 1000; % water density; kg/m^3
rho_a = 1.225; % air density; kg/m^3
Cp_w = 4184; % specific heat of water; J/kgK
Cp_a = 1003.5; % specific heat of air; J/kgK
H_ml = 100; % height of mixing layer (active) of ocean; m
H_a = 10000; % height of the atmosphere; m
H_ocn = 1000; % deep ocean depth; m
k = 5.55e-5; % piston velocity; m/s
S_ocn = 34; % salinity of the ocean surface; ppt??
P_atm = 10^5; % atmospheric pressure; Pa
CO2_s = 242.7; % dissolved CO2 in the surface ocean; Gt
A_earth = 4*pi*6371e03^2; % Earth surface area; m^2
c_emissivity = 0.054; % relating CO2 concentration to atmospheric emissivity
c_aerosol = 3.65e-13; % relating aerosol concentration to earth albedo; kg^-1
 
% Conversion factors
s2y = 1/31536000; %seconds to year conversion; year/second
% todo: get consts for the following equation
GT2g = 1e15; %Gt to gram conversion; g/Gt
mol2umol = 1e06; %mol to umol conversion; umol/mol
mu_CO2 = 44; %CO2 molecular weight; g/mol
GT2umolperkg = GT2g*mol2umol / (mu_CO2 * rho_w * A_earth * H_ocn);  % GT to umol/kg conversion umol/kg/GT
Gt2ppm = 1/0.1291; % Gt to ppm conversion; ppmv/Gt
atm2Pa = 101325; %Pa/atm

% Time step
dt = (0.24/100)/s2y; % time step (atmospheric temp); secs 
n = 10000;%round((100*s2y)/dt); % number of time steps needed to it runs for 100 years. we can change this
%made it three just for testing

% Model parameters we will change for experiments (aka eruptions) 
F_CO2 = 0; % Amount of CO2 released by eruption *CHANGE TO REAL VALUE
F_aero = 0; % Amount of aerosols released by eruption *CHANGE TO REAL VALUE
tau_aero = 1; % residence time for aerosols in the atmosphere *CHANGE TO REAL VALUE
% todo: something for how explosive it is which will change the time constant
 
% Initialize & preallocate arrays
T_e = nan(1, n+1); % temperature of earth; deg K
T_a = nan(1, n+1); % temperature of atmosphere; deg K
CO2_atm = nan(1, n+1); % conc. of CO2 in atmosphere; Gt 
M_a = nan(1, n+1); % mass of aerosols in atmosphere; kg
 
time = nan(1, n+1); % time; s
e_a = nan(1, n+1); % atmospheric emissivity (longwave)
a = nan(1, n+1); % solar constant; W/m^2
 
% Define initial conditions
T_e(1) = 280;
T_a(1) = 230;
CO2_atm(1) = 280 / Gt2ppm; % Gt 
M_a(1) = 0;

time(1) = 0;
e_a(1) = 0.8; % reference emissivity
a(1) = 0.3; % reference albedo 

%% Run model 

for t = 1 : n 
    % co2 pre boom
    k0 = k0calc(T_e(t), S_ocn);
    P_CO2 = (CO2_atm(t)*Gt2ppm)*10^(-6)*atm2Pa;
    CO2_atm(t+1) = CO2_atm(t) + (dt*((-k/(H_ml))*((k0*P_CO2)-CO2_s)));
    
    % emissivity
    e_a(t) = e_a(1) * (1 + c_emissivity * log(CO2_atm(t) / CO2_atm(1)));

    % albedo
    a(t) = a(1) + (c_aerosol * M_a(t));

    % earth temp
    T_e(t+1) = (dt * (((S_o/4)*(1-a(t))+e_a(t)*sigma*(T_a(t)^4) - ...
        sigma*(T_e(t)^4)) / (rho_w*Cp_w*H_ml))) + T_e(t);

    % atmospheric temp
    T_a(t+1) = (dt * ((e_a(t)*sigma*(T_e(t)^4)-2*e_a(t)*sigma*(T_a(t)^4)) / ...
        (rho_a*Cp_a*H_a))) + T_a(t);

    % todo: calc when co2 reaches equilibrium then add boom 

    % aerosols 
    M_a(t+1) = (dt * (M_a(t)/tau_aero)) + M_a(t);

    % time
    time(t+1) = time(t) + dt;
end
 
%% Plot
