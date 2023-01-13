function [x_ds,spikePct,spikeFlag] = despike(x,q,fill)
% function [x_ds,spikePct,spikeFlag] = despike(x,q,fill)
% Removes spikes from a data series.
% Spike detection follows the "mean absolute deviation" method outlined in
% Mauder et al., Agricultural and Forest Meteorology (2013).
%
% INPUTS:
% x:    time series of variable to despike. Must be a vector.
% q:    optional cutoff threshold. Default value is 7 following Mauder (2013).
% fill: optional specification for how to fill in spikes.
%       Can include "nan" (default) or "interp" for linear interpolation.
%
% OUTPUTS:
% x_ds:     despiked x.
% spikePct: percentage of points ID'd as spikes.
% spikeFlag: logical flag for ID'd spikes.
%
% 20140312 GMW
% 20170106 GMW  Added spikeFlag output.

if nargin<3, fill = 'nan'; end
if nargin<2, q = 7; end

x_ds = x;

%ID spikes
m = nanmedian(x);
MAD = nanmedian(abs(x - m)); %median absolute deviation
spike = abs(x - m) >= q.*MAD./0.6745;

spikePct = sum(spike)./length(x).*100;
spikeFlag = spike;

if spikePct==0, return; end %if no spikes, no need to fill

%fill 'em in
switch fill
    case 'nan'
        x_ds(spike) = nan;
    case 'interp'
        i = find(spike);
        xint = interp1(x,i,'linear');
        x_ds(spike) = xint;
    otherwise
        warning(['despike: fill option ' fill ' not recognized. Filling with NaNs.'])
        x_ds(spike) = nan;
end


