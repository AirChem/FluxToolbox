% function [sigma_xy,lags,dt] = lagCov(x,y,maxlag,plotme,t)
% Calculates a lagged covariance between two vectors.
% The function lags y relative to x from -maxlag:maxlag.
% Also makes a lag-correlation plot.
% INPUTS:
% x,y: vectors to be compared. Should be column vectors.
% maxlag: maximum number of points to lag y relative to x.
% plotme: flag for making lag-correlation plot. default is 1 (yes).
% t: time stamps for x and y, used to determine time-step of each lag (optional).
%    t MUST be evenly-spaced for this to work properly.
% OUTPUTS:
% sigma_xy: covariance of x and y at each lag.
% lags: # of lag points for each r. Positve values correspond to y being shifted forward.
% dt: corresponding delta-t of each lag.
%
% Adapted from lagcorr.
%
% 20130904 GMW

function [sigma,lags,dt] = lagCov(x,y,maxlag,plotme,t)

Ly = length(y);

%%%%%DEFAULTS%%%%%
if nargin<4
    plotme = 1;
end

if nargin<5
    t = nan(Ly,1);
end

%%%%%CALCULATE LAGS%%%%%
lags = (-maxlag:maxlag)';
sigma = nan(length(lags),1);
dt = nan(length(lags),1);
xp = x - nanmean(x);
yp = y - nanmean(y);
for i=1:length(lags)
    l = lags(i);
    if l>0
        yl = [nan(l,1); yp(1:Ly-l)];
        dt(i) = t(l+1) - t(1);
    elseif l<0
        l=-l;
        yl = [yp(l+1:end); nan(l,1)];
        dt(i) = t(end-l) - t(end);
    else
        yl = yp;
        dt(i) = 0;
    end
    
    sigma(i) = nanmean(xp.*yl);
end

if plotme
    [~,imax] = max(abs(sigma));
    figure;
    if nargin==5
        plot(dt,sigma,'-')
        xlabel('\Delta t_y')
    else
        plot(lags,sigma,'-')
        xlabel('y lag (# of points)')
    end
    ylabel('Cov(x,y)')
    grid on
end

