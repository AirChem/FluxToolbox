function [ErrSub,Fsub,ErrRun,Frun] = stationarityTest(t,wp,xp,nsub,runflag,plotme)
% function [error,Fsub,Frun] = stationarityTest(t,wp,xp,nsub,runflag,plotme)
% Calculates eddy covariance fluxes from subsets of an EC flux leg.
% There are two ways to go about doing this:
% 1) divide time series into non-overlapping subsets and calculate covariance for each
% 2) calculate a "running covariance" over whole flux leg.
% This is often used as an indicator of flux quality.
% For more info, see Foken and Wichura, Agricultural Forest Meteorology (1996).
%
% INPUTS:
% t:        time stamp, s
% wp:       vertical wind speed, rotated so that mean is 0
% xp:       de-trended and lagged scalar variable
% nlegs:    number of subset to divide data set into. Typically 3 - 5.
% runflag:  flag to also do the running covariance. Default is 0 (no).
% plotme:   flag for plotting results. Default is 0 (no).
%
% OUTPUTS:
% ErrSub: percentage difference between mean of subset fluxes and actual flux.
% Fsub:   Fluxes calculated from subsets of full time series.
% ErrRun: percentage difference between mean of running covariance and actual flux.
% Frun:   Fluxes calculated from running covariance. Same length as input time series.
%
% 20140310 GMW

if nargin<6, plotme = 0;  end
if nargin<5, runflag = 0; end

L = length(wp); %length of whole time series
l = floor(L./nsub); %window size

wp = wp - mean(wp,'omitnan'); %remove means
xp = xp - mean(xp,'omitnan');

wpxp = wp.*xp; %individual covariances
F = mean(wpxp,'omitnan'); %total flux

% option 1: discrete binned covariance
wp_short = wp(1:l*nsub); %remove remainder points off the end
wp_mat = reshape(wp_short,l,nsub);
wp_mean = mean(wp_mat,1,'omitnan');

xp_short = xp(1:l*nsub);
xp_mat = reshape(xp_short,l,nsub);
xp_mean = mean(xp_mat,1,'omitnan');

Fsub = mean(wp_mat.*xp_mat,1,'omitnan') - wp_mean.*xp_mean; %sub-fluxes
tsub = mean(reshape(t(1:l*nsub),l,nsub),1,'omitnan'); %center of each window

ErrSub = 100*(1 - mean(Fsub,'omitnan')./F);

% option 2: running covariance
if runflag
    Frun = smooth(wpxp,l,'mean') - smooth(wp,l,'mean').*smooth(xp,l,'mean'); %see FW96, Eq. 13
    ErrRun = 100*(1 - mean(Frun,'omitnan')./F);
else
    Frun = nan*t;
    ErrRun = 0;
end

% plot
if plotme
    figure
    hold on
    h1 = plot(t,F*ones(L,1),'k-','LineWidth',2);
    h2 = plot(t,0.7*F*ones(L,1),'k--','LineWidth',2);
    h3 = plot(t,1.3*F*ones(L,1),'k--','LineWidth',2);
    h4 = plot(t,Frun,'b-','LineWidth',3);
    h5 = plot(tsub,Fsub,'rd','LineWidth',4,'MarkerSize',18);
    xlabel('Time (s)')
    ylabel('Flux')
    box on
    xlim([min(t) max(t)])
    
    legend([h1 h2 h5 h4],...
        'Flux',...
        '+/- 30%',...
        ['F_{sub} (' num2str(ErrSub,'%2.0f') '%)'],...
        ['F_{run} (' num2str(ErrRun,'%2.0f') '%)'])
    
    dt = median(diff(t));
    lt = num2str(l*dt);
    title(['Staionarity Test, window = ' lt 's'])
end


