function [wdir,wspd]=UV2Wind(U,V,offset)
% function [wdir,wspd]=UV2Wind(U,V,offset)
% Converts cartesian wind vectors from a 3-D wind probe into polar wind direction
% and speed.
%
% INPUTS:
% U:       wind speed along x principal axis (if already compass-corrected, should be + in E->W)
% V:       wind speed along y principle axis (if already compass-corrected, should be + in N->S)
% offset:  Compass orientation of positive wind probe x-axis, relative to TRUE NORTH (not magnetic north).
%          Not needed if data is compassed corrected (e.g. flight data).
%
% OUTPUTS:
% wdir: wind direction in degrees East of true North
% wspd: horizontal wind speed (same units as U and V)
%
% 20140206 GMW

if nargin<3
    offset=0;
end

[rwdir,wspd] = cart2pol(V,U);

wdir = rwdir*180/pi + 180 + offset; %convert to degrees and 

wdir(wdir>=360) = wdir(wdir>=360) - 360;

