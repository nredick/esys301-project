function k0_Gt_Pa = k0calc(t,s)
%
% Calculate solubility coefficient of H2CO3 K0 as a function of
% T (K) and S (psu). 
%
% Input: 
% T: ocean temperature [K]
% S: ocean salinity [psu]
%
% Output:
% ##K0 = Solubility coefficient [mol/kg/atm]
% K0 = Solubility coefficient [umol/kg/atm] - Changed to micromol for
% consistency with co2calc and with Sarmiento
%
% Reference: Weiss 1974.
%
% NOTE: I checked the weiss reference and the output is definitely in
% mol/kg/m3. For instance, Weiss gives a k0 = 7e-02 mol/liter/atm = 70 
% mol/m3/atm or 70,000 umol/kg/atm at T=273K and S=30psu, and the routine gives a k= of ~7e-02
% mol/kg/atm or ~70 mol/m3/atm.
%
% BT: T used to be input in C. Now it is in K.
% To this end, I commented out the first line of code. 
% Old Commment: Convert from C to K and calculate some terms
% This is not true: t is not used anywherem, only tk is used.

%tk = 273.15 + t;
tk = t ;
tk100 = tk/100.0;
tk1002=tk100*tk100;


k0_umol_kg_atm = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) + ...
		s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002));
    
% BTU: Change units from mol/kg/atm --> micromol/kg/atm

k0_umol_kg_atm = k0_umol_kg_atm * 1e06 ; 

%Convert units so it works with our calculations
GT2g = 1e15; %Gt to gram conversion; g/Gt
mol2umol = 1e06; %mol to umol conversion; umol/mol
mu_CO2 = 44; %CO2 molecular weight; g/mol
rho_w = 1000; % water density; kg/m^3
A_earth = 4*pi*6371e03^2; % Earth surface area; m^2
H_ocn = 1000; % deep ocean depth; m

GT2umolperkg = GT2g*mol2umol / (mu_CO2 * rho_w * A_earth * H_ocn);  % umol/kg to GT conversion umol/kg/GT
atm2Pa = 101325; %Pa/atm

k0_umol_kg_Pa = (k0_umol_kg_atm/atm2Pa);
k0_Gt_Pa = k0_umol_kg_Pa/GT2umolperkg;
