function [x_dt, xtrend] = detrend_flux(x,t,trendType,frameSize,plotMe)
% function x_dt = detrend(x,t,trendType,frameSize,plotMe)
% Removes a trend from a piece of data. Meant to be used with ECFluxToolBox.
%
% INPUTS:
% x: variable to be detrended. 1-D vector.
% t: time vector for x.
% trendType: string specifier for detrending. options are "mean", "linear" (default) and "smooth."
% frameSize: number of points to use for smoothing. Default is 100.
% plotMe: flag for generating plots. Default is 0 (no).
%
% OUTPUTS:
% x_dt is x with the trend subtracted.
% xtrend is the trend that was subtracted from x.
%
% 20130915 GMW
% 20140206 GMW  Added xtrend output.
% 20190410 GMW  Renamed from detrend.m to avoid conflict with MATLAB canned routine.

%defaults
if nargin<3, trendType = 'linear'; end
if nargin<4, frameSize = 100; end
if nargin<5, plotMe = 0; end


%calculate trends
xmean = mean(x,'omitnan').*ones(size(x)); %mean
fit = nanpolyfit(t,x,1);
xlin = polyval(fit,t); %linear

if strcmp(trendType,'smooth') || plotMe %skip this if not needed (expensive)
    xsmooth = smooth(x,frameSize,'mean'); %smooth
end

% detrend
switch trendType
    case 'mean',   xtrend = xmean;
    case 'linear', xtrend = xlin;
    case 'smooth', xtrend = xsmooth;
    case 'none',   xtrend = 0;
    otherwise
        error(['detrend: trendType "' trendType '" not recognized. '...
            'Valid values are "mean", "linear" and "smooth".'])
end
x_dt = x - xtrend;

%make plots if desired
if plotMe
    figure
    hold all
    plot(t,x,'c.-',...
        t,xmean,'k-',...
        t,xlin,'r-',...
        t,xsmooth,'b-')
    xlabel('Time')
    ylabel('x')
    legend('Data','Mean','Linear',['Smooth ' num2str(frameSize) 'p'])
    box on
    
    figure
    plot(t,x_dt,'b-')
    xlabel('Time')
    ylabel('x detrended')
end


