function [co2,pH,hco3,co3]=co2calc(t,s,dic,alk)
%
%
% This routine computes various constants as functions of T and S 
% required by the carbon chemistry solver. 
%
%
% BTU - ORIGINAL UNITS
% Input: 
% T   : Temperature [C]
% S   : Salinity [psu] 
% DIC : Dissolved Inorganic Carbon [mol/kg]
% Alk : Alkalinity [mol/kg]

% Output:
% CO2 : CO2 [mol/kg]
% pH  : []
% HCO3: Bicarbonate [mol/kg]
% CO3 : Carbonate [mol/kg]
%
% Example:
% 
% [aco2,pH,hco3,co3]=co2calc(10.9, 34.78, 2.102/1000, 2.000/1000) -- old
% routine before conversion
% [aco2,pH,hco3,co3]=co2calc(10.9, 34.78, 2102, 2000) -- new
% routine after conversion
%
%
%  aco2 = 129.98 umol/kg; pH = 7.1612; hco3 = 1954.2 umol/kg; co3 =  17.86 umol/kg;
%
% This gives the same results as the online calculator at:
% http://biocycle.atmos.colostate.edu/shiny/carbonate/
% Units for this calculator: T in C, S in psu, DIC and Alk in micromol/kg
% See also ScreenShot in ESYS301/Matlab-Model in my dropbox
% COMMENT ON THE http://biocycle.atmos.colostate.edu/shiny/carbonate/ WED
% SITE!!
% 
% Sea Water Properties:
% Tw (C); Sw (per mil); Alkalinity (uequiv/kg); DIC (umol/kg); pH;
% UNKNOWN; aco2 (umol/kg); hco3 (umol/kg); hco3 (umol/kg); co3 (umol/kg);
% UNKNOWN; UNKNOWN; UNKNOWN
% 
%
%-----------------------------------------------------------------
% Note: mol/kg are actually what the body of this routine uses 
% for calculations.  
%---------------------------------------------------------------------
%
% UNITS AFTER CONVERSION:
%
% BTU
% Input: 
% T   : Temperature [C]
% S   : Salinity [psu] 
% DIC : Dissolved Inorganic Carbon [umol/kg]
% Alk : Alkalinity [umol/kg]

% Output:
% CO2 : CO2 [umol/kg]
% pH  : []
% HCO3: Bicarbonate [umol/kg]
% CO3 : Carbonate [umol/kg]
%
%



% BTU: Here I convert the input from umol/kg to mol/kg - required by this
% routine.

dic = dic * 1e-06 ;
alk = alk * 1e-06 ;

% Derive simple terms used more than once

tk = 273.15 + t;
tk100 = tk/100.0;
tk1002=tk100*tk100;
invtk=1.0/tk;
dlogtk=log(tk);
s2=s*s;
sqrts=sqrt(s);
s15=s^1.5;
scl=s/1.80655;

%-------------------------------------------------------------
% K0 from Weiss 1974
%

 k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) + ...
 		s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002));

%-------------------------------------------------------------
%
% k1 = [H][HCO3]/[H2CO3]
% k2 = [H][CO3]/[HCO3]
%
% Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
%

%k1=10^(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk - ...
%		0.0118 * s + 0.000116*s2));

%   Mehrbach et al (1973) refit by Lueker et al. (2000).

 k1 = 10^(-1*(3633.86*invtk - 61.2172 + 9.6777.*dlogtk - 0.011555.*s + 0.0001152.*s2));

%
%-------------------------------------------------------------

%k2=10^(-1*(1394.7*invtk + 4.777 - ...
%		0.0184*s + 0.000118*s2));

%   Mehrbach et al. (1973) refit by Lueker et al. (2000).

 k2 = 10^(-1*(471.78*invtk + 25.9290 - 3.16967.*dlogtk - 0.01781.*s + 0.0001122.*s2));


%
% kb = [H][BO2]/[HBO2]
%
% Millero p.669 (1995) using data from Dickson (1990)
%

kb=exp((-8966.90 - 2890.53*sqrts - 77.942*s +  ...
		1.728*s15 - 0.0996*s2)*invtk + ...
		(148.0248 + 137.1942*sqrts + 1.62142*s) + ...
		(-24.4344 - 25.085*sqrts - 0.2474*s) * ...
		dlogtk + 0.053105*sqrts*tk);
%------------------------------------------------------------------------
% kw = [H][OH]
%
% Millero p.670 (1995) using composite data
%

kw = exp(-13847.26*invtk + 148.96502 - 23.6521 * dlogtk + ...
		(118.67*invtk - 5.977 + 1.0495 * dlogtk) * ...
		sqrts - 0.01615 * s);
%------------------------------------------------------------------------
% Calculate concentrations for borate, sulfate, and fluoride
%
% Uppstrom (1974)

bor = 1.*(416.*(s/35.))* 1.e-6;   % (mol/kg), DOE94
%bor = 0.000232 * scl/10.811;

 p5  = -1.;
 p4  = -alk-kb-k1;
 p3  = dic*k1-alk*(kb+k1)+kb*bor+kw-kb*k1-k1*k2;
 tmp = dic*(kb*k1+2.*k1*k2)-alk*(kb*k1+k1*k2)+kb*bor*k1;
 p2  = tmp+(kw*kb+kw*k1-kb*k1*k2);
 tmp = 2.*dic*kb*k1*k2-alk*kb*k1*k2+kb*bor*k1*k2;
 p1  = tmp+(+kw*kb*k1+kw*k1*k2);
 p0  = kw*kb*k1*k2;
 p   = [p5 p4 p3 p2 p1 p0];
 hroots   = roots(p);
 hplus=max(hroots);
 pco2=dic/k0*(hplus*hplus)/(hplus*hplus+k1*hplus+k1*k2);
 co2=pco2*k0;


pH=-log10(hplus);
hco3=co2*k1/hplus;
co3=hco3*k2/hplus;
%dic=co2+hco3+co3;
%   disp(['K0 = ' num2str(k0)]); 
%   disp(['K1 = ' num2str(k1)]);
%   disp(['K2 = ' num2str(k2)]); 
%   disp(['Kw = ' num2str(kw)]); 
%   disp(['Kb = ' num2str(kb)]);


% BTU: Convert in units that makes more sense :

co2     = co2 * 1e06 ;
hco3    = hco3 * 1e06 ;
co3     = co3 * 1e06 ;

 return
 end





