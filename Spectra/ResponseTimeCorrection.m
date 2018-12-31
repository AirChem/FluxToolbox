function RTcorr = ResponseTimeCorrection(freq,co,tau)
% function RTcorr = ResponseTimeCorrection(freq,co,tau)
% Applies a high-frequency spectral correction of the form (1 + (2*pi*tau*f)^2) to a cospectrum.
%
% INPUTS:
% freq: frequency, Hz
% co: co-spectrum of w and scalar. same length as freq.
% tau: instrument response time, seconds. Typically derived from inspection of cospectral shape
% against a another faster-response instrument.
%
% OUTPUT RTcorr is a correction factor to be applied to the flux.
%
% 20170109 GMW

Tr = (1 + (2.*pi.*tau.*freq).^2); %transfer function for high-frequency loss
coTr = co.*Tr;
RTcorr = sum(coTr)./sum(co); %multiplication factor



