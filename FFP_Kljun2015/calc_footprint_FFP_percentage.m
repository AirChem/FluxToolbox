function [FFP,flag_err] = calc_footprint_FFP_percentage(x_2d,y_2d,f_2d,r)
                                  
% function [FFP,flag_err] = calc_footprint_FFP_percentage(x_2d,y_2d,f_2d,r)
% Derive source area of R% of the flux footprint
% 
% See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
% The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
% Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
% contact: n.kljun@swansea.ac.uk
%
%
% Input for FFP percentage
%    x_2d,y_2d,f_2d = output of calc_footprint_FFP.m or calc_footprint_FFP_umean.m 
%    r              = percentage of footprint, i.e. a value between 10 and 90. 
%                     Can be either a single value (e.g., "80") or an array of 
%                     increasing percentage values (e.g., "[10:10:80]") 
%
% Output of FFP percentage
%    FFP            = structure array with footprint contour lines
%    FFP.r          = percentage of footprint as in input
%    FFP.f          = footprint value at r
%    FFP.x          = x-array for contour line of r
%    FFP.y          = y-array for contour line of r
%                     For array of percentage values, structure entries can be accessed 
%                     as FFP(1).r, FFP(1).x, etc.
%    flag_err       = 1 in case of error, 0 otherwise
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
elseif max(r(:))>90
    display('R must be <= 90')
    ind_return = 1;
elseif min(r(:))<10
    display('R must be >= 10')
    ind_return = 1;
end

%--------------------------------------------------------------------
% Create output structure array
%--------------------------------------------------------------------

FFP = struct('r',[],'f',[],'x',[],'y',[]);

if ind_return
    FFP.r    = NaN;
    FFP.f    = NaN;
    FFP.x    = NaN;
    FFP.y   = NaN;
    flag_err = 1;
    return
end %if

%--------------------------------------------------------------------
% Derive footprint ellipsoid incorporating R% of the flux
%--------------------------------------------------------------------

df = 0.01; 

dx = abs(x_2d(11,10)-x_2d(10,10));
dy = dx;

% Calculate integral of f_2d starting at peak value until R% are reached
fmax = max(f_2d(:));
fr   = fmax;
for i = 1:length(r)
    rloop = r(i)/100;
    if rloop >= 0.8
        dfloop = df./10;
    else
        dfloop = df;
    end
    
    
    fr_sum     = 0;
    fr_sum_old = 0;
    while fr_sum<=rloop
        fr_sum_old      = fr_sum;
        fr_old          = fr;
        fr              = fr-fmax.*dfloop;
        fdummy          = f_2d;
        ind_xyr         = f_2d < fr;
        fdummy(ind_xyr) = NaN;
        fr_sum          = nansum(nansum(fdummy)).*dx.*dy;
    end %while
      
    % Select nearest f_int
    diff_1 = abs(rloop-fr_sum_old);
    diff_2 = abs(rloop-fr_sum);
    if diff_2 < diff_1
        contour_level = fr;
    else
        contour_level = fr_old;
    end %if
      
    contour_r = contourc(x_2d(:,1),y_2d(1,:),f_2d',[contour_level contour_level]);

    % Decrease number of digits and sort/unique
    contour_r = round(contour_r(:,2:end),1);
      
    if max(contour_r(2,:))>max(y_2d(:))
        flag_err  = 1;
        contour_r = NaN;
    end
    
    % Fill output structure
    FFP(i).r = r(i);
    FFP(i).f = contour_level;
    FFP(i).x = contour_r(1,:);
    FFP(i).y = contour_r(2,:);

end %for i


end

