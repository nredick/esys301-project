% Program OcnCarbonClimateModel
%
% Simple bio-geochemmical Box Model.
% This model includes, one atmospheric layer, one surface ocean layer, and
% one deep ocean layer. We are using the conservation of energy (Ta, To1, 
% To2)and conservation of Carbon (DIC, Dissolved Inorganic Carbon) and 
% Nutrient (Ph) to calculate the temporal evolution of the atmospheric and 
% ocean temperature, nutrient and DIC concentrations, and CO2 atmospheric
% concentration. The Alkalinity is kept constant and no nutrient from 
% rivers enters the model which implies that the PO4 reaches steady-states
% and stops evolving. For this reason, I calculate the stead-state for PO4
% analytically and I use those values as initial conditions.
%
% References:
% Sarmiento and Gruber, Book.
% Bjerum diagram: http://biocycle.atmos.colostate.edu/shiny/carbonate/
%
% Simple Atmospheric CO2 model:
% -----------------------------
% Calculates atmospheric [CO2] given anthopogenic emission.
% d[CO2]/dt = FCO2_in - FCO2_atm-ocn
% Coupled model: algebraic growth of CO2 eission
%
% A few basic facts regarding CO2(yr) and FCO2_in(yr):
% [CO2](1850)   = 280 ppmv
% [CO2](2015)   = 400 ppmv
% FCO2_in(1850) =  1.83 GT/year
% FCO2_in(2015) = 36.00 GT/year
%
% Simple Ocean Model - Energy Balance:
% ------------------------------------
% We only consider a surface mixed-layer for the energy balance.  
%
% Simple Surface Ocean Carbon Cycle Model:
% ----------------------------------------
%
% [DIC] = [H2CO3] + [HCO3-] + [CO3]-2
% [DIC] = Carbonic Acid + Bicarbonate + Carbonate
%
% d[DIC]s/dt = k / Hml ([CO2]sat - [CO2]ocn) + nu/Hml ([DIC]d - [DIC]s)
%                                            - Phi_P (C:P)
% where
% [DIC]s = Dissolved Inorganic Carbon in the surface ocean [~2000 umol/kg]
% [DIC]d = Dissolved Inorganic Carbon in the deep ocean [~2500 umol/kg]
% k: piston velocity from k0calc [~ 20cm/s or 5.5 e-05 m/s]
% H_ml : Mixedlayer depth [100 m]
% [CO2]sat: Ocn concentration at saturation [~10 umol]/kg]
% [CO2]ocn: Ocn concentration - from the Bjerrum diagram [~10 umol/kg]
% nu: volume flux between the deep ocean and the surface
% Phi_P: Biological pump [umol/kg]
% (C:P): Redfield ratio [1:106]
%
% Simple Deep Ocean Carbon Cycle Model:
% -------------------------------------
%
% d[DIC]s/dt = nu/Hocn ([DIC]d - [DIC]s) - Phi_P (C:P) * Hml/Hocn
% 
% where
% Hocn: Deep Ocean depth [m]
% Note: Phi_P is multiplied by Hml/Hocn since it is equal to that of the
% surface layer and this insures conservation of mole of DIC.
%
% Equilibrium Solution: Useful to set some model parameters
%
% DIC:
% -----
% At equilibium during pre-industrial (no net CO2 flux between the
% atmophere and the ocean), the balance is between the surface to deep
% ocean DIC fluxes and the Biological Pump. I.e. 
% nu/Hml ([DICd] - [DICs]) - Phi_P C:P = 0
% For a 300 umol/kg difference between DIC surface and deep (as per
% Sarmiento Gruber above), the Biological Pump should be equal to:
% Phi_P = nu/Hml ([PO4d]-[PO4s]) / C:P = 3 m/yr / 100m * (300 umol/kg)/106
%       = 0.085 umol/kg/yr
% PO4:
% ----
% At equilibrium during the pro-industrial (no net CO2 flux between the
% atmosphere and the ocean), the balance is between the surface to deep PO4
% fluxes and the Biological Pump. I.e.:
% nu/Hml ([PO4d] - [PO4s]) = Phi_P = Vmax [PO4s] / ([PO4s] + Kp) = 0.085
% umol/kg/yr
% I am assuming here that I know Phi_P, Vmax and [PO4d] and calculate the
% [PO4s] that gives me the correct Phi_P. Note that I set the
% Phi_P_equibrium in order to get the correct DICs and DICd (see just
% above). This is much simpler than what I have done in class where I have
% fixed [PO4d] and then calculated [PO4s] from the steady solution of the
% PO4 equation. But this gave me no garantee that my [DICs] and [DICd]
% would be as per observations. Since the steady value of DICs and DICd are
% very sensitive to the exact choice of PO4d or PO4s, it is better to work
% this way.
% Solving for [PO4s], we get: [PO4s] = Phi_P Kp / (Vmax - Phi_P)
% [PO4s] = 0.085 umol/kg/yr * 0.4 umol/kg / (12 umol/kg/yr - 0.085 umol/kg/yr)
%        = 0.0029 umol/kg
% From the equation above, we can then calculate the [PO4d] at steady
% state:
% nu/Hml ([PO4d] - [PO4s]) = Phi_P ==> [PO4d] = Phi_P Hml/nu + [PO4s]
% [PO4d] = 0.085umol/kg/yr * 100m / 3m/yr + 0.0029 umol/kg
%        = 2.8362 umol/kg
%
% Derivation of Equations:
%
% Atmospheric Carbon Budget
%
% Units of each flux term is mol/yr and eta_co2 is mol
%
% We like however the following units for each of the term
%
% [CO2] : [umol CO2/mol air] --> typical magnitude: 280 -- 400
% F'_co2 : [GT co2/yr]
% 
% d(eta_co2)/dt = F_co2 - Fao
% 1st term:  d(M_air*kg2g/mu_air [CO2]/mol2umol)/dt = M_air*kg2g/(mu_air*mol2umol) d([CO2])/dt = ...
%                                                   = P_atm*A/g * kg2g/(mu_air*mol2umol)d([CO2])/dt = ...
% 2nd term:  F_co2 = F'_co2 * GT2g / mu_co2
% 3rd term:  Fao = k/H_ml * ([co2]sat - [co2]ocn) * rhow*A*H_ml/mol2umol
% All terms: d([CO2])/dt = F'_co2 * GT2ppmv - k/H_ml * ([co2]sat - [co2]ocn) * umolperkgperyr2ppmv
% REWRITE THOSE USING VARIABLES DEFINED HERE
% where GT2ppmv = GT2g * g * (mu_air*mol2umol) / (mu_co2 * Patm * A * kg2g)
%               =  1e15 * 10 * (29 * 1e06) / (44 * 1e05 * 4 * pi * 6371e03^2 * 1000)
%               = 0.1292 GT of CO2 / ppmv
% where umolperkg2ppmv = (rhow * H_ml * g * mu_air) /  (Patm * kg2g)
%                           = (1e03 * 100 * 10 * 29) / (1e05 * 1e03)
%                           = 29e12 / 1e014 = 0.29
%
% Sink at the bottom of the ocean must be equal to the input of CO2 in the
% pre-industrial era for an equilibrium to exist.
% 1e15 g/GT * 1e06 umol/mol / (nu_co2 * rhow * A * Hd)


clear all, close all

%
% Run parameters
%

tstart          = 0 ;           % Initial time [year AD]
tend            = 4000 ;        % End time [year AD]
dt              = 14 ;          % delta t [days] 
pctCO2_emiss    = 1.8 ;         % percentage increase in CO2 emission [%/year]

%
% Define constants
%

Rearth      = 6371e03 ;         % Earth radius [m]
Aearth      = 4*pi*Rearth^2 ;   % Earth surface area [m2]

g           = 9.81 ;            % gravitational acceleration [m/s2]

Ha          = 10e03 ;           % Atmosphere scale height [m]
Hml         = 100 ;             % Surface Ocean depth [m]
Hocn        = 1000 ;            % Deep Ocean depth [m]

rhoa        = 1 ;               % air density [kg/m3]
rhow        = 1030 ;            % water density [kg/m3]

Cpa         = 1000 ;            % Specific Heat of air [J/kg/K]
Cpw         = 4000 ;            % Specific Heat of water [J/kg/K]

Patm        = 1e05 ;            % Mean surface Atmospheric pressure [Pa]

So          = 1360 ;            % Solar constant [W/m2]
A           = 0.3 ;             % Planetary albedo []
sigma       = 5.67e-08 ;        % Stefan Boltzmann constant [W/K^4/m2]
epsa0       = 0.8 ;             % reference atmospheric emissivity
a           = 0.1 ;             % atmospheric absorptivity to SW []

mu_CO2      = 44 ;              % CO2 moelcular weight [g/mol]
mu_O2       = 32 ;              % O2 molecular weight [g/mol]
mu_N2       = 28 ;              % N2 molecular weight [g/mol]
O2frac      = 0.21 ;            % Mole fraction of O2 in air []
N2frac      = 0.79 ;            % Mole fraction of N2 in air []
mu_air      = O2frac*mu_O2 ...  % Air molecular weight [g/mol]
             +N2frac*mu_N2 ;     
         
k           = 20 ;              % piston velocity [cm/hr]
nu          = 3 ;               % Volume exchange [m/yr --> m/s) 

% Vmac: 1 umol/kg of [PO4] is consumed in 1 month. 
% Phi_P = Vmax [PO4] / (Kp + [PO4]_mean) 

Vmax        = 3.8052e-7 ;       % Maximum uptake rate [umol/kg/sec] ** changed units to umol/kg
Kp          = 0.4 ;             % [umol/kg] [0.3-0.5 mmol/m3, Eric Galbraith, pers comm] ** changed units without changing actual #
CPratio     = 106 ;             % Redfield Ratio C:N:P = 106:16:1

Socn        = 34 ;              % Global mean surface ocean salinity [ppm]

%
% Conversion factors
%

day2sec     = 24*60*60 ;        % day to sec conversion [sec/day]
sec2yr      = 1/day2sec/365 ;   % sec to year conversoin [year/sec]
yr2sec      = 60*60*24*365 ;    % yr to sec conversion [sec/yr]
hr2sec      = 60*60 ;           % sec to hour conversion [sec/hr]
sec2mo      = 1/day2sec/30 ;    % second to month conversion [mo/sec]

kg2g        = 1000 ;            % kg to gram conversion [g/kg]
g2Pg        = 1e-15 ;           % g to Petagram conversion [Pg/g]
GT2g        = 1e15 ;            % GT to gram conversion [g/GT]

oneppmv     = 1*1e-06;          % One ppmv Atmospheric concentration [ppmv]

cm2m        = 1e-02 ;           % cm to meter conversion [m/cm]

umol2mmol   = 1e-03 ;           % micromol to milimol [mmol/umol]
GT2g        = 1e15 ;            % GT to gram conversion [g/GT]
mol2umol    = 1e06 ;            % mol to umol conversion [umol/mol]

T0          = 273 ;             % 0 Celsius in Kelvin [K]


% ppmv to GT of Carbon conversion - We dont use it but useful
% [CO2] = n_co2 (mol) / n_air (mol) = (M_co2/mu_co2) / (M_air/mu_air)
% Patm  = M_air g / Aearth ==> M_air = Patm Aearth / g
% [CO2] = M_co2 / (Patm Aearth / g) * mu_air/mu_co2 
%       = M_co2 * mu_air/mu_co2 * g / (Patm Aearth)
%       = M_co2 GT * mu_air/mu_co2 [] * g / (Patm Aearth) [kg] * 1e12 kg/GT * 1e06 []/ppmv
% [CO2] = 0.1261 GT/ppmv

ppmv2mol     = Patm*Aearth/g * kg2g/(mu_air*mol2umol) ; % ppmv to mole conversion [mol/ppmv]
GT2mol       = GT2g / mu_CO2 ;                          % GT to mole conversion [mol/GT]
mol2kgmol    = rhow*Aearth*Hml/mol2umol ;               % mol/kg to umol/kg conversion [mole/(umol/kg)]
GT2umolperkg = GT2g*mol2umol / ...                      % umol/kg to GT conversion umol/kg/GT
               (mu_CO2 * rhow * Aearth * Hocn) ;
%
% Convert Model Parameters to the correct units
%

dt          = dt * day2sec ;       % time step [day --> sec]
k           = k*cm2m/hr2sec;       % piston velocity [cm/hr --> m/s]
nu          = nu*sec2yr ;         % Volume exchange [m/yr --> m/s) 

%
% Initialize variables
%
                                  % Set vector length
N        = length(tstart:dt*sec2yr:tend) ; 
time     = nan(1,N) ;             % time [year]
Tw       = nan(1,N) ;             % ocean temperature [K] - from spin up
Ta       = nan(1,N) ;             % atmosphere temperature [K] - from spin up
Ew       = nan(1,N) ;             % Energy of water [W/m2]
Ea       = nan(1,N) ;             % Energy of air [W/m2]
epsa     = nan(1,N) ;             % atmospheric emissivity 

CO2      = nan(1,N) ;             % [CO2] concentration in 1850
FCO2_in  = nan(1,N) ;             % Pre-industrial CO2 emission [GT/year] --> [GT/sec]

DICs     = nan(1,N) ;             % Surface ocean DIC
DICd     = nan(1,N) ;             % Deep ocean DIC
Alks     = nan(1,N) ;             % Surface ocean alkalinity
PO4s     = nan(1,N) ;             % Surface ocean phosphate
PO4d     = nan(1,N) ;             % Deep ocean phosphate
Phi_P    = nan(1,N) ;             % Biological Pump

%
% Define Initial Conditions
%

time(1)     = tstart ;            % Initial time [year]
Tw(1)       = 287.8846 ;          % Initial ocean temperature [K] - Equilibrium from spin up
Ta(1)       = 249.9117 ;          % Initial atmosphere temperature [K] - Equlibrium from spin up
Ew(1)       = rhow*Cpw*Hml*Tw(1) ;% Initial Energy of water [W/m2]
Ea(1)       = rhoa*Cpa*Ha*Ta(1) ; % Initial Energy of air [W/m2]
epsa(1)     = epsa0 ;             % Initial atmospheric emissivity 

CO2(1)      = 280 ;               % Initial [CO2] concentration in 1850
FCO2_in(1)  = 1.83 / yr2sec ;     % Pre-industrial CO2 emission [1.83 GT/year] --> [GT/sec]

DICs(1)     = 1500 ;              % DIC surface [umol/kg] - [DICs]ss = 2000 umol/kg S&G, page 204 fig 8.1.2
DICd(1)     = 2000 ;              % DIC deep [umol/kg]    - [DICd]ss = 2300 umol/kg       "   "
Alks(1)     = 1920 ;              % Alk surface [umol/kg] - Alk = 2000 umol/kg            "   "
% http://sam.ucsd.edu/sio210/gifimages/A16_PHSPHT.gif
PO4s(1)     = 0.0029 ;            % Phosphate surface [umol/kg] - see section below for #estimate
PO4d(1)     = 2.8362 ;            % Phosphate deep [umol/kg]  == we dont use this val
Phi_P(1)    = Vmax * PO4s(1) ...
              / (Kp + PO4s(1)) ;  % Biological pump

%
% Cacluate the solution at future time
%

t = 0 ;                     % Time counter

for tt = tstart : dt*sec2yr : tend
    t = t + 1 ;  

% Energy Ocean 

% 0.043 gives epsa = 0.83 for double [CO2] (=560 ppmv) 
    epsa(t+1)   = epsa0 + 0.043*log(CO2(t)/CO2(1)) ; 

    Ew(t+1)     = Ew(t) + dt * (So*(1-A)*(1-a)/4 + epsa(t)*sigma*Ta(t)^4 ...
                                - sigma*Tw(t)^4) ; 
    Tw(t+1)     = Ew(t+1)/(rhow*Cpw*Hml) ;
%    Tw(t+1)     = Tw(t) ;  % Keep Tw constant to debug carbon component of the model
    
% Energy Atmosphere 
    Ea(t+1)     = Ea(t) + dt * (So*(1-A)*a/4 + epsa(t)*sigma*Tw(t)^4 ...
                          - 2*epsa(t)*sigma*Ta(t)^4) ; 
    Ta(t+1)     = Ea(t+1)/(rhoa*Cpa*Ha) ;
%    Ta(t+1)     = Ta(t) ; % Keep Ta constant to debug carbon component of the model

% Carbon Atmosphere
    pCO2atm     = CO2(t) * 1e-06 * 1                      % Partial pressure of CO2 [1 atm]
    CO2sat(t)   = k0calc(Tw(t), Socn) * pCO2atm;           % k0: output:[umol/kg/atm]; input{[K, psu]; 
    k000 = k0calc(Tw(t), Socn)
    [CO2ocn(t),pH,HCO3,CO3] = co2calc(Tw(t)-T0,Socn,DICs(t),Alks(t)) ; % [umol/kg,[],umol/kg,umol/kg] = [C, psu, umol/kg, umol/kg]
    
    FCO2_ao(t)  = k/Hml * (CO2sat(t) - CO2ocn(t)) ;        % CO2 flux in surface ocean [umol/kg/sec]

    CO2(t+1)    = CO2(t) + dt * (FCO2_in(t) * GT2mol/ppmv2mol - ...
                                 FCO2_ao(t) * mol2kgmol/ppmv2mol) ; 
%    CO2(t+1)    = CO2(t) ;   % KEEP constant for debugging purposes
%    CO2(t+1)    = CO2(t) + dt * (FCO2_in(t) - ...
%                                 FCO2_in(1)*sec2yr) * 0.6 * GT2ppmv ;   % 40% absorption by ocn - 1.83 GT/yr

% Ocean Phosphate - assume steady state

    Phi_P(t+1)  = Phi_P(t) ;                            % Biological pump
%    Phi_P(t+1)  = Vmax * PO4s(t) / (Kp + PO4s(t)) ;    % Biological pump
%    Upwelling   = nu * (PO4d(t) - PO4s(t)) ;           % Upwelling of nutrient from depth
    
    PO4s(t+1)   = PO4s(t) ;                             % Keep PO4 in the surface/deep ocean constant   
    
    PO4d(t+1)   = PO4d(t) ; 
    
% Surface Ocean DIC 

    FDIC_do_so           = nu * (DICd(t)   - DICs(t)) ;             % DIC flux from deep to surface ocean [mol/m3/sec]
    
    DICs(t+1)            = DICs(t) + dt * (FCO2_ao(t) ...
                                         + FDIC_do_so / Hml ...     % Idivide by Hml and by Hocn in the DICd
                                         - Phi_P(t) * CPratio); 
    
% Deep Ocean DIC

    DICd(t+1) = DICd(t) + dt * (-FDIC_do_so / Hocn + ...
                                 Phi_P(t) * CPratio * Hml/Hocn - ...
                                 FCO2_in(1) * GT2umolperkg );     %
%    DICd(t+1) = DICd(t) ;              % for debugging purposes

                                     
% Alkalinity Equation - we are not solving for the alkalinity - important
% for longer time scale
 
    Alks(t+1) = Alks(t) ;           % Keep alkalinity constant [mol/m3]


%    Ideas about how to code it. See notes from Eric.
%    CO3 = Alk(t) - DIC(t) ;                 % this is an approximation, but is it contradicting Bjerrum???
%    CaCO3_sinking_diss = 0.2 * Phi_P(t) ; 
%    CO3_compensation   = CO3_target - CO3 / tau ;   
%    Alk(t+1) = Alk(t) + dt * (CaCO3_sinking_diss + CO3_compensation) ;
%    

%
% Calculate next time step CO2 emission
%
    if time(t) <= 1850
        FCO2_in(t+1)    = FCO2_in(t) ; 
    elseif (time(t) > 1850 & time(t) < 2050)
        % this transform a Q%/yr increase into a %/dt increase
        FCO2_in(t+1) = FCO2_in(t) * (1+pctCO2_emiss/100*sec2yr*dt) ; 
    elseif time(t) > 2030
        FCO2_in(t+1)    = FCO2_in(1) ;
    end

    time(t+1)   = time(t) + dt*sec2yr ; 

end

%%
% Plot Results
%

% Water and Atmosphere Temperature temporal evolution

figure(1), clf
plot(time, Tw-T0, 'r-')
%hold on
%plot(time, Ta-T0, 'b-') ;
xlabel('time (years)')
ylabel('Temperature (C)')
title('Earth Surface and Atmosphere Temperature - Ocean-Atm Model')
grid on 

% DIC temporal evolution

figure(2), clf
plot(time, DICs, 'r-')
hold on
plot(time, DICd, 'b-')
xlabel('time (years)')
ylabel('DICs [umol/kg]')
title('Dissolved Inorganic Carbon - Surface/Deep Ocean - Temporal Evolution')
legend('DICs','DICd')
grid on

% Surface and Deep Ocean CO2 concentration temporal evolution

figure(3), clf
plot(time(2:end), CO2sat, 'r-')
hold on
plot(time(2:end), CO2ocn, 'b-')
xlabel('time (years)')
ylabel('[CO2] [umol/kg]')
title('Surface Ocean Carbin Fluxes - Temporal Evolution')
legend('[CO2]sat', '[CO2]_ocn')
hold on
grid on

% Atmopsheric [CO2] concentration temporal evolution

figure(4), clf
plot(time, CO2, 'b-')
hold on
xlabel('time (years)')
ylabel('[CO2] (ppmv)')
title('CO2 Concentration Temporal Evolution')
grid on

% Atmospheric emissivity temporal evolution

figure(5), clf
plot(time, epsa, 'b-')
hold on
xlabel('time (years)')
ylabel('emissivity []')
title('Atmospheric Emissivity Temporal Evolution')
hold on
grid on


% CO2 emission temporal evolution: 
% pre-industrial: 1.83 GT of CO2 / year
% Today: 36 GT of CO2 / year

figure(6), clf
plot(time, FCO2_in/sec2yr, 'r-')
hold on
%plot(time, FCO2_out/sec2yr, 'b-')
xlabel('time (years)')
ylabel('CO2 Emission [GT/year]')
title('CO2 Emission Temporal Evolution [GT or Pg]')
grid on

% Atmosphere-Ocean [CO2] flux temporal evolution

figure(7), clf
plot(time(1:end-1), FCO2_ao/sec2yr, 'r-')
hold on
xlabel('time (years)')
ylabel('CO2 Flux surface [mol/m3/yr]')
title('Surface Ocean Carbin Fluxes - Temporal Evolution')
hold on
grid on




