function L = obukhov(ustar,Tv0,H0)
% function L = obukhov(ustar,Tv0,H0)
% Calculates Obukhov Length.
% INPUTS:
% ustar: friction velocity in m s^-1
% Tv0: surface virtual potential temperature in K
% H0: kinematic heat flux at surface, K m/s
%
% OUTPUT is monin-obukhov length in m.
%
% 20071204 GMW

k = 0.41; %von Karman constant
g = 9.81; %m s^-2

L = -(ustar.^3)./(k.*g.*H0./Tv0);