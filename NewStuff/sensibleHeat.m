function H = sensibleHeat(wT,P,T,Pv)
% function H = sensibleHeat(wT,P,T,Pv)
% Converts a sensible heat flux from kinematic units (K m/s) to energy units (W/m^2).
% INPUTS:
% wT: kinematic heat flux in K m/s. Can be scalar, vector or matrix.
% P: average air pressure, mbar.
% T: average temperature, K.
% Pv: average water vapor partial pressure, mbar.
%
% OUTPUT:
% H: sensible heat flux in W/m^2.
%
% 20140612 GMW
% 20170907 GMW  Commented out additional term from water vapor. This is only necessary when using
%               virtual temperature...I think.

Md  = 0.028964; % mass of dry air, kg/mol
Mv = 0.018016; % mass of water, kg/mol
Cpd = 1004;     %specific heat of dry air, J/kg*K
Cpv = 1952; %specific heat of water vapor, J/kg*K (note: no clear consensus on this #)

R   = 8.314e-2; %gas constant, mbar*m^3/mol*K

% Pd = P - Pv; %partial pressure of dry air
% rhod = Pd.*Md./R./T; %dry air density, kg/m^3
% rhov = Pv.*Mv./R./T; %vapor density
% 
% H = wT.*(Cpd.*rhod + Cpv.*rhov);

rhod = P.*Md./R./T;
H = wT.*Cpd.*rhod;

