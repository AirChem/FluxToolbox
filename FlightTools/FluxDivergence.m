function [div,div_err,fits,usedpts] = FluxDivergence(flux,RE,SE,lat,lon,alt,leg,divLegs)
% function [div,div_err,fits,usedpts] = FluxDivergence(flux,RE,SE,lat,lon,alt,leg,divLegs)
% Calculates divergence correction factors for airborne fluxes.
% INPUTS:
% flux: time series of fluxes
% RE: random error in fluxes (in flux units)
% SE: systematic error in fluxes (in fractional units)
% lat: latitude
% lon: longitude
% alt: altitude, m above ground
% leg: index for flight legs
% divLegs: 2-column cell array specifying legs to use when calculating divergence corrections.
%          First column identifies which legs to use for slope calculation.
%          Second column identifes which legs the correction applies to.
%          divLegs can contain multiple rows, one for each set of calculations.
%
% OUTPUTS:
% div: divergence correction factors, same size as input flux
% div_err: uncertainty in divergence correction factors (in same units as div)
% fits: 5-column matrix of fit parameters for flux vs alt, [m b r sm sb]
% usedpts: flag indicating which flux points were used for divergence calculation.
%           this can be a subset of total due to the search for the "box of overlap."
%
% 20170526 GMW
% 20170712 GMW  Added usedpts output.
% 20170801 GMW  Corrected error averaging to use RMS.
% 20170814 GMW  Switched overlap determination to using "lloverlap" function with circles
% 20171026 GMW  Added separate input of RE and SE and adjusted error averaging accordingly
%               Included nan exclusion within iavg

% initialize some stuff
Npts = length(flux);
Ndiv = size(divLegs,1);
fits = nan(Ndiv,5);
M = nan(Npts,1); B = M; SM = M; SB = M; %fit coefficients
allLegs = unique(leg);
allLegs(isnan(allLegs))=[];
usedpts = false(Npts,1);

% get divergence
for i = 1:Ndiv
    
    % break out indices
    fitLegs = divLegs{i,1}; %legs to use for slope calculation
    corLegs = divLegs{i,2}; %legs to apply correction to
    icorr = ismember(leg,corLegs); %index for corrected legs
    
    % define region of overlap among legs
    radius = 1; %km
    iavg = lloverlap(lat,lon,leg,fitLegs,'circle',radius);
    iavg = iavg & ~isnan(flux);
    usedpts = (usedpts | iavg); %accumulate

    % average legs within box
    alt_avg = BinAvg(leg(iavg),alt(iavg),allLegs); alt_avg = alt_avg(fitLegs);
    flux_avg = BinAvg(leg(iavg),flux(iavg),allLegs); flux_avg = flux_avg(fitLegs);
    
    %average errors
    % REavg = sqrt(sum(RE^2))/N
    [RE_avg,~,Navg] = BinAvg(leg(iavg),RE(iavg).^2,fitLegs);
    RE_avg = sqrt(RE_avg./Navg);
    SE_avg = abs(BinAvg(leg(iavg),SE(iavg),fitLegs).*flux_avg); %now in flux units
    flux_err_avg = sqrt(RE_avg.^2 + SE_avg.^2);
    
    % error-weighted fit flux vs alt
    [m b r sm sb] = lsqfityz(alt_avg,flux_avg,flux_err_avg); 
    fits(i,:) = [m b r sm sb]; %save fit coefficients
    
    % distribute fit coefficients along legs
    M(icorr) = m; B(icorr) = b; SM(icorr) = sm; SB(icorr) = sb;
end

% calculate corrections
div = B./(M.*alt + B); %multiplicative correction factor

% error in correction factor
dcdm = -B.*alt./(M.*alt + B).^2; %partial derivative of F(0)/F(z) wrt m
dcdb = M.*alt./(M.*alt + B).^2; %wrt b
div_err = sqrt((dcdm.*SM).^2 + (dcdb.*SB).^2);


