function [cohr,phase] = coherenceNphase(psdx,psdy,co,qu)
% function [cohr,phase] = coherenceNphase(psdx,psdy,co,qu)
% Calculates coherence and phase angle for a set of co-spectra.
% For reasons not well understood by the author, this MUST be done with bin-averaged spectra.
% INPUTS:
% psdx: power spectral density of x
% psdy: power spectral density of y
% co: cospectrum of x and y
% qu: quadrature spectrum of x and y
%
% OUTPUTS:
% cohr: estimate of spectral correlation between x and y
% phase: phase angle between x and y, in degrees
%
% 20140321 GMW

cross = co + qu*1i; %cross spectrum

cohr   = abs(cross).^2./(psdx.*psdy); 
phase = atan2(co,qu)*180/pi;
