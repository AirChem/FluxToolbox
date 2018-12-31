function Ra = Ra_SP(ustar,Z,d,z0,L)
% function Ra = Ra_SP(ustar,Z,d,z0,L)
% calculates aerodynamic resistance to surface transfer in the boundary layer.
% Equations from Seinfeld and Pandis, Section 19.3.
% INPUTS:
% ustar: friction velocity, cm/s
% Z: height, m
% d: displacement height, m (typically ~70% of height of roughness elements)
% z0: surface roughness length, m
% L: obukhov length, m
%
% OUTPUT is resistance in s/m.
%
% 2014321 GMW

%get started
k = 0.41; %von Karman constant
r0 = 1./(k.*ustar); %neutral resistance
zeta = (Z-d)./L; %stability parameter
zeta0 = z0./L;

%neutral
Ra = r0.*log((Z - d)./(z0));

%stable
j = zeta>0;
Ra(j) = Ra(j) + r0(j).*4.7.*(zeta(j) - zeta0(j));

%unstable
j = zeta<0;
n = (1 - 15*zeta(j)).^0.25;
n0 = (1 - 15*zeta0(j)).^0.25;
Ra(j) = Ra(j) + r0(j).*(log((n0.^2 + 1).*(n0 + 1).^2./(n.^2 + 1)./(n + 1).^2) + 2.*(atan(n) - atan(n0)));


% older method
%stability corrections, from Arya (2002)
% zeta = (Z-d)./L;
% phiH = nan(size(zeta));
% phiM = nan(size(zeta));
% i = find(zeta<0);%unstable
% x = (1 - 15*zeta(i)).^0.25;
% phiH(i) = 2*log((1 + x.^2)./2);
% phiM(i) = log(((1 + x.^2)./2).*((1 + x)./2).^2) - 2*atan(x) + pi/2;
% i = find(zeta>=0);%stable
% phiH(i) = -5*zeta(i);
% phiM(i) = -5*zeta(i);
% 
% Ra = ubar./ustar.^2 - (phiH - phiM)./(k.*ustar);

