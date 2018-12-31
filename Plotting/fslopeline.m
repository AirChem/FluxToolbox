function [fslope,h] = fslopeline(f,y1,slope,plotme)
% function [fslope,h] = fslopeline(f,yscale,,slope,plotme)
% generates a line with a user-defined slope in log space.
% Typically plotted alongside turbulence spectra.
%
% INPUTS:
% f: frequencies (x values) for line
% y1: y value for first point of line
% slope: value of slope, e.g. -5/3
% plotme: flag for plotting in current figure. 1=yes, 0=no (default).
%
% OUTPUT:
% fslope: line with specified slope and scaled to y1. Same size as f.
% h: handle for plot
%
% 20131010 GWW

if nargin<4, plotme=0; end

fslope = f.^(slope);
fslope = fslope./fslope(1).*y1;

if plotme
    hold on
    h = plot(f,fslope,'k--');
end