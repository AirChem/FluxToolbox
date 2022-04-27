function E = latentHeat(wH2O,T)
% function E = latentHeat(wH2O,T)
% Converts water vapor fluxes from density units (molec/m^2/s) to Energy units (W/m^2).
%
% INPUTS:
% wH2O: water flux in molec/m^2/s. Can be scalar, vector or matrix.
% P: average air pressure, mbar.
% T: average temperature, K.
%
% OUTPUT:
% E: latent heat flux in W/m^2.
%
% 20140612 GMW

% heat of vaporization, J/g
% taken from Wiki article "Latent Heat"
% see also Table 2.1. R. R. Rogers & M. K. Yau (1989). A Short Course in Cloud Physics (3rd ed.). Pergamon Press. p. 16. ISBN 0-7506-3215-1
TC = T-273.15; %celcius
Lv = 2500.8 - 2.36*TC + 0.0016*TC.^2 - 6e-5*TC.^3;

mv = 18.02; %g/mol
A = 6.022e23; %molec/mol

E = wH2O.*mv./A.*Lv;

