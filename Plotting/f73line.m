function f73 = f73line(f,y1,plotme)
% function f73 = f73line(f,yscale,plotme)
% generates a line with a slope of -7/3 in log space.
% Typically plotted alongside turbulence spectra.
%
% INPUTS:
% f: frequencies (x values) for line
% y1: y value for first point of line
% plotme: flag for plotting in current figure. 1=yes, 0=no (default).
%
% OUTPUT:
% f73: line with -7/3 slope and scaled to y1. Same size as f.
%
% 20131010 GWW

if nargin<3, plotme=0; end

f73 = f.^(-7./3);
f73 = f73./f73(1).*y1;

if plotme
    hold on
    plot(f,f73,'k--')
end