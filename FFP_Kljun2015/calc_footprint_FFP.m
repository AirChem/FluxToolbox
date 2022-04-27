function [x_ci_max,x_ci,f_ci,x_2d,y_2d,f_2d,flag_err]=calc_footprint_FFP(zm,z0,umean,h,ol,sigmav,ustar)

% [x_ci_max,x_ci,f_ci,x_2d,y_2d,f_2d,flag_err]=calc_footprint_FFP(zm,z0,umean,h,ol,sigmav,ustar)
% Derive a flux footprint estimate based on the simple parameterisation FFP
% 
% See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
% The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
% Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
% contact: n.kljun@swansea.ac.uk
%
%
% FFP Input
%    zm        = Measurement height above displacement height (i.e. z-d) [m]
%    z0        = Roughness length [m] - enter [NaN] if not known 
%    umean     = Mean wind speed at zm [ms-1] - enter [NaN] if not known 
%                Either z0 or umean need to be entered. If both are given,
%                umean is selected to calculate the footprint
%    h         = Boundary layer height [m]
%    ol        = Obukhov length [m]
%    sigmav    = standard deviation of lateral velocity fluctuations [ms-1]
%    ustar     = friction velocity [ms-1]
%
% FFP output
%    x_ci_max  = x location of footprint peak (distance from measurement) [m]
%    x_ci      = x values of crosswind integrated footprint [m]
%    f_ci      = footprint function values of crosswind integrated footprint [m-1] 
%    x_2d      = x-grid of 2-dimensional footprint [m]
%    y_2d      = y-grid of 2-dimensional footprint [m]
%    f_2d      = footprint function values of 2-dimensional footprint [m-2]
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

if nargin ~= 7
    display('wrong number of input arguments')
    ind_return = 1;
elseif zm <= 0
    display('zm needs to be larger than 0')
    ind_return = 1;
elseif h < 10
    display('h needs to be larger than 10 m')
    ind_return = 1;
elseif (zm>h) && (ol>0)
    display('zm needs to be smaller than h')
    ind_return = 1;
elseif (zm>0.8*h) && (ol<0)
    display('for convective stratification, zm needs to be smaller than the entrainment height')
    ind_return = 1;
elseif zm/ol<=-15.5
    display('zm/L needs to be equal or larger than -15.5')
    ind_return = 1;
end %if
if ~isnan(z0)
    if z0 < 0
        display('z0 needs to be larger than 0')
        ind_return = 1;
    elseif zm<12.5*z0
        %changed to lowest limit of roughness sublayer definition
        display('zm needs to be above roughness sublayer')
        ind_return = 1;
    end
elseif ~isnan(umean)
    if umean < 0
        display('umean needs to be larger than 0')
        ind_return = 1;
    end
elseif isnan(z0) && isnan(umean)
    display('enter either z0 or umean')
    ind_return = 1;

end %if

if ind_return
    x_ci_max = NaN;
    x_ci     = NaN;
    f_ci     = NaN;
    x_2d     = NaN;
    y_2d     = NaN;
    f_2d     = NaN;
    flag_err = 1;
    return
end %if

%--------------------------------------------------------------------
% Initialize variables
% Selection of nx has impact on accuracy and on output file size, 
% decrease for speed, increase for accuracy
% (nx=3000 ideal but slow, nx=600 fast but may not resolve details)
%--------------------------------------------------------------------

xstar_end = 30;
nx        = 800;

a = 1.4524;
b = -1.9914;
c = 1.4622;
d = 0.1359;

ac = 2.17; 
bc = 1.66;
cc = 20.0;

%limit for neutral scaling
ol_n = 5000;

%von Karman
k = 0.4;

%--------------------------------------------------------------------
% Create scaled X* for crosswind integrated footprint
%--------------------------------------------------------------------

xstar_ci_param = linspace(d,xstar_end,nx+2);
xstar_ci_param = xstar_ci_param(2:end);
         
%--------------------------------------------------------------------
% Calculate crosswind integrated scaled F* 
%--------------------------------------------------------------------

fstar_ci_param = a.*(xstar_ci_param-d).^b .* exp(-c./(xstar_ci_param-d));
     
ind_notnan     = ~isnan(fstar_ci_param);
fstar_ci_param = fstar_ci_param(ind_notnan);
xstar_ci_param = xstar_ci_param(ind_notnan);

%--------------------------------------------------------------------
% Calculate scaled sig_y*
%--------------------------------------------------------------------

sigystar_param = ac.*sqrt(bc.*(xstar_ci_param).^2 ./ (1+cc.*(xstar_ci_param)));

%--------------------------------------------------------------------
% Calculate real scale x and f_ci
%--------------------------------------------------------------------

if ~isnan(umean)
    x = xstar_ci_param.*zm ./ (1-(zm./h)) .* (umean./ustar.*k);

    if (umean/ustar)>0
        x_ci = x;
        f_ci = fstar_ci_param./zm .* (1-(zm./h)) ./ (umean./ustar.*k);
    else
        x_ci_max = NaN;
        x_ci     = NaN;
        f_ci     = NaN;
        x_2d     = NaN;
        y_2d     = NaN;
        f_2d     = NaN;
        flag_err = 1;
    end

else
    if ol <=0 || ol >=ol_n
        xx  = (1 - 19.0.*zm./ol).^0.25;
        psi_f = log((1+xx.^2)./2) + 2.*log((1+xx)./2) - 2.*atan(xx) + pi./2;
    elseif ol > 0 && ol < ol_n
        psi_f = -5.3.*zm./ol;
    end
    
    x = xstar_ci_param.*zm ./ (1-(zm./h)) .* (log(zm./z0)-psi_f);
    if (log(zm./z0)-psi_f)>0
        x_ci = x;
        f_ci = fstar_ci_param./zm .* (1-(zm./h)) ./ (log(zm./z0)-psi_f);
    else
        x_ci_max = NaN;
        x_ci     = NaN;
        f_ci     = NaN;
        x_2d     = NaN;
        y_2d     = NaN;
        f_2d     = NaN;
        flag_err = 1;
    end
end %if

if size(x_ci) == size(x)
%--------------------------------------------------------------------
% Calculate maximum location of influence (peak location)
%--------------------------------------------------------------------

   xstarmax = -c./b+d;
   if ~isnan(umean)
       x_ci_max = xstarmax.*zm ./ (1-(zm./h)) .* (umean./ustar.*k);
   else
       x_ci_max = xstarmax.*zm ./ (1-(zm./h)) .* (log(zm./z0)-psi_f);
   end %if

%--------------------------------------------------------------------
% Calculate real scale sigy
%--------------------------------------------------------------------

   if abs(ol) >ol_n
       ol = -1000000;
   end
   if ol <=0 %convective
       scale_const = 1E-5.*abs(zm./ol).^(-1)+0.8;
   elseif ol > 0  %stable
       scale_const = 1E-5.*abs(zm./ol).^(-1)+0.55;
   end %if
   if scale_const>1
       scale_const  = 1.0;
   end %if
   
   sigy         = sigystar_param./scale_const .*zm .*sigmav./ustar;
   sigy(sigy<0) = NaN;

%--------------------------------------------------------------------
% Calculate real scale f(x,y)
%--------------------------------------------------------------------

   dx    = x_ci(3)-x_ci(2);
   y_pos = 0:dx:(length(x_ci)/2)*dx*1.5;

   f_pos = NaN.*ones(length(f_ci),length(y_pos));
   for i=1:length(f_ci)
       f_pos(i,:) = f_ci(i) * 1./(sqrt(2.*pi).*sigy(i)) .* exp(-y_pos.^2./(2.*sigy(i).^2));
   end

%--------------------------------------------------------------------
% Complete footprint for negative y (symmetrical)
%--------------------------------------------------------------------

   y_neg = -fliplr(y_pos); 
   f_neg = fliplr(f_pos);
   y     = [y_neg(1:end-1) y_pos];
   f     = [f_neg(:,1:end-1) f_pos];

%--------------------------------------------------------------------
% Matrixes for output
%--------------------------------------------------------------------

   x_2d = repmat(x',1,length(y));
   y_2d = repmat(y,length(x),1);
   f_2d = f;


end

end
