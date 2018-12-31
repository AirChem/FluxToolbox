function f53 = f53line(f,y1,plotme)
% function f53 = f53line(f,yscale,plotme)
% generates a line with a slope of -5/3 in log space.
% Typically plotted alongside turbulence spectra.
%
% INPUTS:
% f: frequencies (x values) for line
% y1: y value for first point of line
% plotme: flag for plotting in current figure. 1=yes, 0=no (default).
%
% OUTPUT:
% f53: line with -5/3 slope and scaled to y1. Same size as f.
%
% 20131010 GWW

if nargin<3, plotme=0; end

f53 = f.^(-5./3);
f53 = f53./f53(1).*y1;

if plotme
    hold on
    plot(f,f53,'k--')
end