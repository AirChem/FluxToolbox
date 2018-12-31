function [flux_err,lagcovs] = WaveletFluxError(x_wave,w_wave,scale,freq,param,errMinFreq,errMaxLag)
% function [flux_err,lagcovs] = WaveletFluxError(x_wave,w_wave,scale,freq,param,errMinFreq,errMaxLag)
% Calculates statistical errors in wavelet fluxes based on an extension of the "variance of
% covariance" method discussed in Finklestein and Sims (2001).
%
% INPUTS:
% x_wave, w_wave: raw wavelet coefficients from WaveletFlux function. Can be time-averaged.
% scale: wavelet scale vector
% freq: wavelet frequency vector
% param: structure of wavelet parameters. Must include fields: dj, dt, Cdelta.
% errMinFreq: OPTIONAL lowest frequency to use for error calculations. Default = 0.
% errMaxLag: OPTIONAL maxmimum lag for lag-covariances. Default = floor(N/2).
%
% OUTPUT:
% flux_err: time series of flux errors in flux units.
% lagcovs: 5-column matrix giving lags and covariances for ww, xx, wx, xw.
%
% 20170801 GMW
% 20170913 GMW Changed last equation to not divide by N.

% check inputs
N = size(w_wave,2);
if nargin<6, errMinFreq = 0; end
if nargin<7, errMaxLag = floor(N/2); end

pnames = {'dj','dt','Cdelta'};
for i=1:length(pnames)
    assert(isfield(param,pnames{i}),'Input "param" must include field "%s".',pnames{i});
end

% variance scaling
scale_big  = repmat(scale,1,N);
w_wave = w_wave.*sqrt(param.dj.*param.dt./param.Cdelta./scale_big);
x_wave = x_wave.*sqrt(param.dj.*param.dt./param.Cdelta./scale_big);

% calculate covariances at all lags
k = freq>errMinFreq;
wxcov = real(w_wave(k,:)'*x_wave(k,:)); %not sure why we only need to take the real part, but this agrees with brute-force method
xxcov = real(x_wave(k,:)'*x_wave(k,:));
wwcov = real(w_wave(k,:)'*w_wave(k,:));

% extract subset of lags and rearrange
lags = -errMaxLag:errMaxLag;
xwcov = spdiags(wxcov,lags); %exploiting symmetry of xw and wx
wxcov = spdiags(wxcov',lags);
xxcov = spdiags(xxcov',lags);
wwcov = spdiags(wwcov',lags);

% bias correction
bc = N./repmat(N-abs(lags),N,1); %correction factor
xwcov = xwcov.*bc;
wxcov = wxcov.*bc;
xxcov = xxcov.*bc;
wwcov = wwcov.*bc;

% compute error based on variance of covariance
flux_err = sqrt(abs(nansum(wwcov.*xxcov + wxcov.*xwcov,2))./1); %divide by 1 instead of N b/c considering each point

% additional outputs: time-averaged lag covariances
if nargout>1
    lagcovs = [lags', nanmean(wwcov,1)', nanmean(xxcov,1)', nanmean(wxcov,1)', nanmean(xwcov,1)'];
end

