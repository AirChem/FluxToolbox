function og = ogive(f,s,absFlag)
% function og = ogive(f,s,absFlag)
% calculates an ogive (cumulative integral) for a spectrum.
% currently starts integration at low-frequency end.
%
% INPUTS:
% f: frequency, Hz. Should be column vector.
% s: spectra power density (variance per unit frequency).
% absFlag: optional flag to integrate absolute value of spectral power density.
%
% OUTPUTS:
% og: ogive, scaled from 0 to 1. Same size as s.
%
% 20071128 GMW
% 20090711 GMW  removed for-loop.
% 20131008 GMW  modified to use absolute value of cospectral power.
% 20161130 GMW  modified to accept matrix input for s, added input dimension checking
%               converted from summation to trapezoidal integration.
% 20161214 GMW  Removed absolute value usage.
% 20170131 GMW  added absFlag input, b/c I can't seem to decide what I want.

% massage inputs
N = length(f);
if size(f,2)>1, f=f(:); end %turn to columns if necessary

if size(s,1)==N
    srot = 0;
elseif size(s,2)==N
    s=s';
    srot = 1;
else
    error('Ogive: spectra input must have one dimension equal to size of freqeuncy input');
end
s(isnan(s)) = 0;

% ensure integration from low to high freqency
if f(2)<f(1)
    f = flipud(f);
    s = flipud(s);
    fflip = 1;
else
    fflip = 0;
end

% do integral
if nargin==3 && absFlag
    og = cumtrapz(f,abs(s)); %absolute value
else
    og = cumtrapz(f,s); %or not
end
og = og./repmat(og(end,:),N,1); %normalize to total area (variance)

% adjust outputs if needed
if fflip, og = flipud(og); end
if srot, og = og'; end


