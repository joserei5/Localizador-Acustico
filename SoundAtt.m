function dbm = SoundAtt(Psin, Tin, hrin, f)
% Function for calculate the attenuation in the air because only the atmospheric
% conditions in dB/m
%
%   dbm = SoundAtt(Psin, Tin, hrin, f)
%
% Source:
%   Nathan Burnside 10/5/04
%   AerospaceComputing Inc.
%   nburnside@mail.arc.nasa.gov
%
% INPUT
%       Psin - atmospheric pressure in atm;
%       Tin  - temperature in ºC;
%       hrin - hunidity of the air in %;
%       f    - central frequency.
%
%  Copyright 2008 LocUS Project.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 1.0$    $Date: 2008/03/17 20:39$
%
T = Tin + 273.15; % temp input in K
To1 = 273.15; % triple point in K
To = 293.15; % ref temp in K

Ps = Psin; % static pressure in atm
Pso = 1; % reference static pressure

F = f./Ps; % frequency per atm


% calculate saturation pressure
Psat = 10^(10.79586*(1-(To1/T))-5.02808*log10(T/To1)+1.50474e-4*(1-10^(-8.29692*((T/To1)-1)))-4.2873e-4*(1-10^(-4.76955*((To1/T)-1)))-2.2195983);


h = hrin*Psat/Ps; % calculate the absolute humidity 

% Scaled relaxation frequency for Nitrogen
FrN = (To/T)^(1/2)*(9+280*h*exp(-4.17*((To/T)^(1/3)-1)));

% scaled relaxation frequency for Oxygen
FrO = (24+4.04e4*h*(.02+h)/(.391+h));

% attenuation coefficient in nepers/m
alpha = Ps.*F.^2.*(1.84e-11*(T/To)^(1/2) + (T/To)^(-5/2)*(1.275e-2*exp(-2239.1/T)./(FrO+F.^2/FrO) + 1.068e-1*exp(-3352/T)./(FrN+F.^2/FrN)));

dbm = 10*log10(exp(2*alpha));