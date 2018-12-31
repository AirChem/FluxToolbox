function W = WaveletSpectralCorrection(W,tau)
% function W = WaveletSpectralCorrection(W,tau)
% Applies a high-frequency spectral correction of the form (1 + (2*pi*tau*f)^2) to a wavelet power
% spectrum and calculates new fluxes, cospectra and ogives.
% INPUTS:
% W: wavelet structure as output by WaveletFlux.m
% tau: instrument response time, seconds. Typically derived from inspection of cospectral shape
% against a another faster-response instrument.
%
% OUTPUT, is the W input structure with a few more fields:
% flux_hfc: fluxes after high-frequency correction
% co_hfc: global cospectrum afer high-frequency correction
% og_hfc: global ogive after high-frequency correction
%
% 20170109 GMW

Tr = (1 + (2.*pi.*tau.*W.freq).^2); %transfer function for high-frequency loss
power_hfc = W.power.*repmat(Tr,1,length(W.flux));
co_hfc = nanmean(power_hfc,2); % time average [Eqn(22)], but in units of variance and bias-corrected
% df = -ndiff(W.freq);
% W.og_hfc = ogive(W.freq,W.co_hfc./df);
W.quality.HFcorr = sum(co_hfc)./sum(W.co); %multiplication factor

% W.flux_hfc = sum(power_hfc)';   % scale average [Eqn(24)] (divide by scale done above)
% W.flux_hfc = W.flux.*W.quality.HFcorr;


