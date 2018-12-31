function Rb = Rb_JH95(ustar,P,lai,leafWidth,D)
% function Rb = Rb_JH95(ustar,P,lai,leafWidth,D)
% calculates laminar sublayer resistance for a forest canopy following the parameterization of 
% Jensen and Hummelshoj (1995,1997).
% INPUTS:
% ustar: friction velocity, m/s
% P: surface pressure, mbar
% lai: leaf area index, m^2/m^2
% leafWidth: characteristic leaf width, m
% D: moleular diffusion coefficient, m^2/s
%
% OUTPUT, Rb, is resistance in s/m.
%
% 20140321 GMW

nu = 1.46e-5.*1013./P; %air viscosity
c = 100; %tunable constant

Rb = (nu./D).*(c.*leafWidth.*ustar./lai.^2./nu).^(1/3)./ustar;

