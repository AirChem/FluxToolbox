function D = D_H2O(T,P)
% Parameterization for diffusivity of water vapor in air.
% Found on the internet, with a reference to Bolz and Tuve (1976).
% Diffusivity for other molecules may be estimated by scaling with ratio of square root of masses:
% Dm = Dh2o * sqrt(18./Mm)
%
% INPUTS:
% T: Temperature, K
% P: pressure, mbar
%
% OUTPUTS:
% D: diffusivity of H2O, m^2/s
%
% 20140321 GMW

a = -2.775e-6;
b=4.479e-8;
c=1.656e-10;

D = a + b*T + c*T.^2;

D = D.*1013./P;