function E = FluxErrors_flight(z,zbl,L,T,I,tau)
% Based on formulae found in Lenschow et al., J. Atmos. Ocean. Technol. (1994).
% See also Karl et al., J. Atmos. Sci. (2013)
%
% INPUTS:
% z: measurement height, m above ground level
% zbl: boundary layer height, m
% L: length of flight leg, m
% T: time-width of flight leg, s
% I: effective sampling interval (time between points), s
% tau: integral time scale for w'c', s. Only needed for disjunct errors.
%
% OUTPUTS:
% All outputs are in a structure, E.
% SE: systematic error, %
% RE: random error, %
% SE_dj: systematic error from disjunct sampling, %
% RE_dj: random error from disjunct sampling, %
% 
%
% 20140512 GMW

% sampling errors
SE = 2.2.*zbl.*sqrt((z./zbl))./L; %Eq. 64
RE = 1.75.*(z./zbl).^0.25.*sqrt(zbl./L); %Eq. 66
E.SE = SE*100;
E.RE = RE*100;

%disjunct time series errors
if nargin>3
    SE_dj = I./T.*(coth(I./tau./2) - (I./T).*(1 - exp(-T./tau))./(2.*sinh(I./tau./2).^2)); %Eq. 55
    RE_dj = I./tau./2.*coth(I./tau./2) - 1; %Karl, Eq. 5
    E.SE_dj = SE_dj.*100;
    E.RE_dj = RE_dj.*100;
end


% continuous formulae; can be used if you have all integral timescales
% tau_w: integral timescale of w, seconds
% tau_s: integral timescale of scalar s, seconds
% tau_ws: integral timescale of w-s correlation, seconds
% tau_f: integral timescale of w*s timeseries, seconds
% r_ws: correlation coefficient of w and s (at optimum lag)

% SE = 2.*tau_ws./T; %eq. 27
% SEalt = 2.*sqrt(tau_w.*tau_s)./r_ws./T; %eq. 29
% 
% RE = sqrt(w.*tau_f./tau_w).*sqrt((1+r_ws.^2)./r_ws.^2); %eq. 48
% REalt = 2.*sqrt(min([tau_w,tau_s])./T)./r_ws; %eq. 49


