function [tlag,T,Talt] = integralTimeScale(t,x,y,maxlag,plotme)
% function T = integralTimeScale(t,x,y,maxlag,,plotme)
% Calculates and plots an "integral time scale", which is a measure of the time over which a
% variable is correlated with itself (or with another variable).
% Technially this is defined in the limit of integration out to infinite lag times,
% so users should use their best judgement when interpreting integrals of actual data.
%
% INPUTS:
% t: time vector, usually in seconds
% x: first data vector
% y: second data vector
% maxlag: number of points to lag by. Default is half of data length.
% plotme: flag for plotting results. default is 1 (yes).
%
% OUTPUTS:
% tlag: lag times for each value of T
% T: integral timescale over all lag times
% Talt: different way of calculating T (from autocorrelations of x and y)
% 
% 20140210 GMW

% default values
if nargin<5
    plotme=1;
end

if nargin<4
    maxlag = floor(length(t)./2);
end

%calculate lagged correlation coefficients
[rlag,~,tlag] = lagcorr(x,y,maxlag,0,t);

%fold lags and average
l = length(tlag);
mid = l./2 + 0.5;
thalf = tlag(mid:end);
rhalf = rlag(mid:end);
rhalf(2:end) = (rhalf(2:end) + flipud(rlag(1:mid-1)))./2; %average + and - lags

%calculate cumulative integral
T = cumtrapz(thalf,rhalf);

%alternative method: T<sqrt(Tx*Ty)/rxy
xlag = lagcorr(x,x,maxlag,0);
ylag = lagcorr(y,y,maxlag,0);

xhalf = xlag(mid:end);
yhalf = ylag(mid:end);

Tx = cumtrapz(thalf,xhalf);
Ty = cumtrapz(thalf,yhalf);

Talt = real(sqrt(Tx.*Ty)./rhalf(1));

%plot
if plotme
    figure
    subplot(121)
    plot(tlag,rlag,'.-')
    xlabel('\Delta t_y')
    ylabel('r_x_y')
    grid on
    
    subplot(122)
    plot(thalf,T,'b-',thalf,Talt,'r--')
    xlabel('\Delta t_y')
    ylabel('\int R(t)dt')
    legend('x-y','xx-yy')
    grid on
end

