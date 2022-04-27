function [x_2d_rot,y_2d_rot,f_2d_rot,flag_err] = calc_footprint_FFP_rotated(x_2d,y_2d,f_2d,wind_dir)
                                  
% function [x_2d_rot,y_2d_rot,f_2d_rot,flag_err] = calc_footprint_FFP_rotated(x_2d,y_2d,f_2d,wind_dir)
% Rotate the FFP footprint into mean wind direction of selected time stamp (e.g. 30-min mean)
% 
% See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
% The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
% Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
% contact: n.kljun@swansea.ac.uk
%
%
% Input for FFP rotated
%    x_2d,y_2d,f_2d = FFP output of calc_footprint_FFP.m or calc_footprint_FFP_umean.m
%    wind_dir       = wind direction in degrees (of 360)      
%
% Output of FFP rotated
%    x_2d_rot  = rotated x-grid of 2-dimensional footprint [m]
%    y_2d_rot  = rotated y-grid of 2-dimensional footprint [m]
%    f_2d_rot  = rotated footprint function values of 2-dimensional footprint [m-2]
%    flag_err  = 1 in case of error, 0 otherwise
%
% created: 15 April 2015 natascha kljun
% version: 1.01
% last change: 18/11/2015 natascha kljun
%
% Copyright (C) 2015, Natascha Kljun

%--------------------------------------------------------------------
% Check input variables
%--------------------------------------------------------------------

flag_err   = 0;
ind_return = 0;

if nargin ~= 4
    display('wrong number of input arguments')
    ind_return = 1;
elseif size(x_2d)~=size(y_2d) | size(x_2d)~=size(f_2d)
    display('problem with footprint input')
    ind_return = 1;
elseif wind_dir>360
    display('problem with wind direction')
    ind_return = 1;
elseif wind_dir<0
    display('problem with wind direction')
    ind_return = 1;
end

if ind_return
    x_2d_rot = NaN;
    y_2d_rot = NaN;
    f_2d_rot = NaN;
    flag_err = 1;
    return
end %if

%--------------------------------------------------------------------------
% rotate 3d footprint
%--------------------------------------------------------------------------
    
wind_dir_rad = wind_dir .* pi ./180;

len = size(x_2d,2);
x_2d_rot = NaN .*ones(size(x_2d));
y_2d_rot = NaN .*ones(size(y_2d));
f_2d_rot = NaN .*ones(size(f_2d));

for i = 1:len
    x_sel = x_2d(:,i);
    y_sel = y_2d(:,i);
    dist  = sqrt(x_sel.^2 + y_sel.^2);
    angle = atan2(y_sel,x_sel);

    x_2d_rot(:,i) = dist .* sin(wind_dir_rad - angle);
    y_2d_rot(:,i) = dist .* cos(wind_dir_rad - angle);
    f_2d_rot(:,i) = f_2d(:,i);
end %for
   
end

