function [sigxy,lags] = lagCovFFT(x,y,maxlag,plotme,Fcut,Fsample)
% function [sigxy,lags] = lagCov(x,y,maxlag,plotme,Fcut,Fsample)
% Calculates a lagged covariance between two vectors.
% The function lags y relative to x from -maxlag:maxlag.
% Also makes a lag-correlation plot.
%
% INPUTS:
% x,y: vectors to be compared. Should be column vectors.
% maxlag: maximum number of points to lag y relative to x. Default is length(x)/2.
% plotme: flag for making lag-correlation plot. default is 0 (no).
% Fcut: optional minimum frequency, Hz. Frequencies lower than or equal to Fcut will be removed before
%   calculating autocovariance.
% Fsample: sample frequency for x and y, Hz. Required if Fcut is specified.
%
% OUTPUTS:
% sigma_xy: covariance of x and y at each lag.
% lags: # of lag points for each r. Positve values correspond to y being shifted forward.
%
% Adapted from lagCov.
%
% 20170106 GMW
% 20170727 GMW  Fixed error on line 57. df = Fsample./N should read Fsample./nfft.

%%%%%CHECK INPUTS%%%%%
% assert(~any(isnan(x)),'input x cannot contain NaN.')
% assert(~any(isnan(y)),'input y cannot cantain NaN.')
assert(length(x)==length(y),'inputs x and y must be same length.')
N = length(y);

% denan
% note, nearest neighbor interp seems to best recover true covariance for both intermittent and large gaps
t = (1:N)'; %fake time
n = isnan(x);
if any(n), x = interp1(t(~n),x(~n),t,'nearest','extrap'); end
n = isnan(y);
if any(n), y = interp1(t(~n),y(~n),t,'nearest','extrap'); end

x = x - mean(x);
y = y - mean(y);

%%%%%DEFAULTS%%%%%
if nargin<3 || isempty(maxlag), maxlag = fix(N/2); end
if nargin<4, plotme = 0; end
if nargin<5, Fcut=[]; end

%%%%%CALCULATE LAGS%%%%%
lags = (-maxlag:maxlag)';

%get cross spectrum
nfft = 2^nextpow2(2*N-1); %number of points
xfft = fft(x,nfft);
yfft = fft(y,nfft);
xyCr = xfft.*conj(yfft); %cross-spectrum

% optional detrending
if ~isempty(Fcut) && exist('Fsample','var')
    df = Fsample./nfft; %frequency spacing for fft
    icut = floor(Fcut./df); %index for highest frequency to cut
    icut = [2:icut+1 nfft-icut+1:nfft]; %all frequencies to cut (start at 2 b/c 1 is direct current)
    xyCr(icut)=0; %snip
end

% inverse fft to get lag covariance
sigxy = ifft(xyCr);
sigxy = [sigxy(end-maxlag+1:end); sigxy(1:maxlag+1)]; %rearrange
sigxy = sigxy./(N-abs(lags)); %un-bias

%%%%%PLOT%%%%%
if plotme
    figure;
    plot(lags,sigxy,'-')
    xlabel('y lag (# of points)')
    ylabel('Cov(x,y)')
    grid on
end

