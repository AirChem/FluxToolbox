function [C_H2O,C_T] = WPL(P,T,Nv,Nc,wT,wH2O)
% function [C_H2O,C_T] = WPL(P,T,Nv,Nc,wT,wH2O)
% Calculates the correction terms for the effect of density fluctuations on eddy covariance fluxes.
% Reference: Webb, Pearman and Leuning, "Correction of flux measurements for density effects due to
%            heat and water vapour transfer," Quart. J. R. Met. Soc. (1980), 106, pp.85-100.
% 
% Notes:
% 1) The use of one or both corrections depends on the nature of the flux measurement.
%    - If the sample is dried and brought to constant T, no correction is needed.
%    - If the sample is not dried but is at constant T, only the H2O correction is needed.
%    - If the sample is at ambient conditions, both corrections are needed.
% 2) Corrections may require post-calculation unit conversion (e.g. from number density to mixing ratio).
% 3) The original WPL equations are given in mass density, but this function uses number densities (cleaner math).
% 4) If needed, use ConvertHumidity.m to get H2O in correct units for input (molec/cm^3).
%
% INPUTS:
% P:    pressure, mbar
% T:    temperature, K
% Nv:   water vapor concentration, molec/cm^3
% Nc:   scalar concentration, molec/cm^3
% wT:   kinematic heat flux, K m/s
% wH2O: uncorrected water flux, molec/m^2/s
%
% OUTPUTS:
% C_H2O: correction term for water vapor effects, molec/m^2/s
% C_T:   correction term for temperature effects, molec/m^2/s
%
% 20140617 GMW

N = NumberDensity(P,T); %concentration of moist air, molec/cm^3
Na = N - Nv; %dry air

%Eq. 24 from WPL80
C_H2O = Nv./Na.*wH2O;
C_T = (1 + Nv./Na).*Nc*1e6./T.*wT;


