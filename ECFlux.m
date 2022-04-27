function F = ECFlux(data,options)
% function F = ECFlux(data,options)
% Runs a series of calculations to derive a flux using the eddy covariance method.
% The order of operations is as follows:
% 1) rotate 3-D winds into natural wind coordinate (such that mean cross and vertical wind speeds are 0)
% 2) despike and detrend scalar variable x
% 3) lag scalar and vertical wind vectors to optimize covariance
% 4) calculate covariance and exchange velocity
% 5) plot quadrant data
% 6) plot frequency spectra
% 7) Quality analysis
%
% INPUTS
% data: structure containing the necessary data. Note that all vectors must be the same length.
%   t: time vector, typically in seconds. Must have constant spacing.
%   x: scalar vector (e.g. concentration).
%      This can also be a string for momentum fluxes: 'u' (mean wind) or 'v' (cross wind).
%   u: horizontal wind speed
%   v: cross-wind speed
%   w: vertical wind speed
%   optional fields for flight errors:
%   z: sampling altitude, m agl
%   zi: boundary layer depth, m agl
%   Lleg: flight leg length, m
% options: structure containing parameters that specify operations and output.
%   name:            name of flux chunk for labeling figures                Default: []
%   despikeFlag:     flag for despiking x data following Mauder (2013)      Default: 0
%   xTrendType:      detrending method: 'mean', 'linear', 'smooth', 'none'. Default: 'mean'
%   frameSize:       frame size for trend if 'smooth' is selected.          Default: 100
%   nRot:            number of wind rotations.                              Default: 2
%   plotX:           flag for plotting x data.                              Default: 0
%   plotW:           flag for plotting wind data.                           Default: 0
%   plotLag:         flag for lag-covariance plot.                          Default: 0
%   plotQuad:        flag for quadrant plot.                                Default: 0
%   plotSpectra:     flag for frequency-spectrum plots.                     Default: 0
%   lagWindow:       Range of lags to look over for max covariance.         Default: [-100 100];
%   xLag:            lag to apply to x (if not determined automatically).   Default: []
%   nStat:           number of subsets for stationarity test.               Default: 5
%   errMinFreq       lowest frequency to use for error calculations.        Default: 0
%   errMaxLag        maximum lag to use in error calculation.               Default: floor(N/2)
%   tauResponse      response time for high-frequency correction.           Default: 0
%   plotWave         flag for wavelet plots.                                Default: 0
%   waveParam        wavelet parameters. See WaveletFlux.m for details.     Default: []
%   gapFill          wavelet gap-filling method: 'none','stitch','interp'.  Default: 'none'
%   scrubGapBound    flag to scrub wavelet fluxes around data gaps.         Default: 0
%   waveError        flag to calculate wavelet flux errors (slow).          Default: 0
%
% OUTPUTS
% F: a structure containing the following fields:
%   data:   same as data input
%   options: same as options input (but xLag replaced with calculated lag if it was empty)
%   flux:   eddy covariance flux, <w'x'>
%   vex:    exchange velocity, <w'x'>/mean(x))
%   x_dt:   x after detrending
%   x_dtl:  detrended x, lagged to optimize covariance with w
%   xtrend: trend subtracted from x.
%   u_r:    rotated horizonal wind speed
%   v_r:    rotated cross wind speed
%   w_r:    rotated vertical wind speed
%   angles: structure of wind rotation angles:
%       eta:    rotation angle around the z1 axis
%       theta:  rotation angle around the y1 axis
%       beta:   rotation angle around x2 axis to make mean(w2'v2') = 0.
%   cov_wx: covariance of w' and x' at various lags
%   cov_wxdt: covariance of w' and x' at various lags with errMinFreq removed.
%   lags:   # of x lag points relative to w. Corresponds to cov_wx.
%   spectra:    structure containing spectra and cospectra:
%       f:      frequency, Hz
%       psdx:   x power spectrum
%       psdxn:  variance-normalized x power spectrum
%       ogx:    x ogive 
%       psdw:   w power spectrum
%       psdwn:  variance-normalized w power spectrum
%       ogw:    w ogive
%       co:     w-x cospectrum
%       con:    covariance-normalized w-x cospectrum
%       ogwx:   w-x ogive
%       qu:     quadrature spectrum (imaginary cospectrum)
%       coherence: coherence for cospectrum
%       phase:  phase angle for cospectrum
%   quality: structure containing various quality criteria
%       spikePct: percentage of x-values identified as spikes
%       spikeFlag: logical flag marking location of identified spikes in input x.
%       statTest: results of stationarity test (% deviation from mean flux for data subsets)
%       xvarS2N: signal/noise ratio for variance in x.
%       xstd_noise: standard deviation of white noise in input x.
%       REnoise: flux random error due to instrument noise (Lenschow 2000). 
%                 Given in absolute (flux) units.
%       REfs01: total flux random error (Finklestein & Sims 2001).
%                 Given in absolute (flux) units.
%       REturb: flux random error due to flight sampling of stochastic turbulence (Lenschow 1994).  
%                  Given in fractional units.
%       SEturb: flux systematic error due to flight under-sampling (Lenschow 1994).
%                  Given in fractional units.
%       SErt_ft: systematic error fraction due to response time, based on FFT cospectrum.
%       SErt_wv: systematic error fraction due to response time, based on wavelet cospectrum.
%   stats: statistics for various variables, including
%           x, w, x_dtl ,w_r, and the product w_r*x_dtl (the sum of which give the flux).
%           Also includes correlation coefficient of x_dtl and w_r.
%   wave: structure of wavelet outputs. See WaveletFlux.m for details.
%
% 20130920 GMW
% 20131008 GMW  Added option to input data.x as a string for momentum fluxes.
% 20140324 GMW  Added despiking option and stats output. Also started adding quality stuff.
% 20140812 GMW  Added white-noise flux detection limit calculation.
% 20160920 GMW  Added wavelet transforms
% 20170106 GMW  Replaced lagCov with lagCovFFT
%               Replaced "nLags" input with "lagWindow"
%               Removed white-noise DL code (it's bogus)
%               Added calculations for random error (instrument, flight, total) and systematic error (flight)
%               Added tauResponse option and ResponseTimeCorrection call.
% 20170621 GMW  Modified all errors to output as a fraction of total flux.
%               Changed to detrended lag covariance for peak determination.
% 20170801 GMW  Propagated wavelet error flags
%               Renamed "REtotal" to "REfs01"
% 20170817 GMW  Changed REfs01 and REnoise from fractional to flux units
% 20220423 GMW  Wavelet call: added removal of nans from lag-shifting. Not sure why this is broken
%               now and not previously.
    

%% GET STARTED
% assign default options
defaults = {...
   'name'               '';...
   'despikeFlag'        0;...
   'xTrendType'         'mean';...
   'frameSize'          100;...
   'nRot'               2;...
   'plotX'              0;...
   'plotW'              0;...
   'plotLag'            0;...
   'plotQuad'           0;...
   'plotSpectra'        0;...
   'lagWindow'          [-100 100];...
   'xLag'               [];...
   'nStat'              5;...
   'errMinFreq'         0;...
   'errMaxLag'          floor(length(data.t)/2);...
   'tauResponse'        0;...
   'gapFill'            'none';...
   'plotWave'           0;...
   'waveParam'          [];...
   'scrubGapBound'      0;...
   'waveError'          0;...
   };

Onames = fieldnames(options);
check = ismember(defaults(:,1),Onames);
imiss = find(~check);
for i=1:length(imiss)
    options.(defaults{imiss(i),1}) = defaults{imiss(i),2};
end

% break out structures
struct2var(data)
struct2var(options)
sampleFreq = 1./mean(diff(t)); %data frequency, Hz

%% BASIC EDDY COVARIANCE
% rotate winds
[u_r,v_r,w_r,angles] = natWindRot(u,v,w,t,nRot,plotW);

% assign x in the case of momentum fluxes
if ischar(x)
    switch x
        case 'u', x = u_r;
        case 'v', x = v_r;
        otherwise
            error(['data.x input "' x '" not recognized. ' ...
                'Valid string inputs are "u" and "v".'])
    end
    data.x = x; %overwrite string input
end

% despike
if despikeFlag
    [x,quality.spikePct,quality.spikeFlag] = despike(x);
end

% detrend x
[x_dt,xtrend] = detrend_flux(x,t,xTrendType,frameSize,plotX);

% calculate flux and lag
[cov_wx,lags] = lagCovFFT(w_r,x_dt,[],plotLag);
cov_wxdt      = lagCovFFT(w_r,x_dt,[],0,errMinFreq,sampleFreq); %detrended
if ~isempty(xLag)
    flux = cov_wx(lags==xLag);
else
    ilook = lags>=lagWindow(1) & lags<=lagWindow(2);
    cov_wx_mask = cov_wxdt;
    cov_wx_mask(~ilook)=nan; %mask points to ignore
    [~,imax] = max(abs(cov_wx_mask)); %max covariance method
%     [~,imax] = max(abs(gradient(gradient(cov_wx_mask)))); %max second derivative method
    flux = cov_wx(imax);
    xLag = lags(imax);
    options.xLag = xLag; %replace empty value
    disp(['Best lag is ' num2str(xLag) ' points.'])
end
vex = flux./mean(x,'omitnan'); %exchange velocity
x_dtl = lagVar(x_dt,xLag); %lag x

if plotLag
    hold on
    plot(xLag,flux,'rp') %mark lag peak
end

% more plots
if plotQuad
    quadPlot(w_r,x_dtl);
end

% spectra
nbin = 50; %number of log bins
spectra = allSpectra(w_r,x_dtl,sampleFreq,nbin,plotSpectra);

%% WAVELETS
good = ~isnan(x_dtl); %nans from lag-shifting
wdata.t = data.t(good);
wdata.x = x_dtl(good);
wdata.w = w_r(good);
wdata.name = options.name;
if isfield(data,'speed'); wdata.speed = data.speed; end
wave = WaveletFlux(wdata,gapFill,plotWave,waveParam,scrubGapBound,waveError,errMinFreq);

%% QUALITY AND CORRECTIONS
% stationarity
quality.statTest = stationarityTest(data.t,w_r,x_dtl,nStat);

% get autocovariances (for error calculations below)
[xvar,Elags] = lagCovFFT(x_dtl,x_dtl,[],0,errMinFreq,sampleFreq);
wvar = lagCovFFT(w_r,w_r,[],0,errMinFreq,sampleFreq);
wxcv = lagCovFFT(w_r,x_dtl,[],0,errMinFreq,sampleFreq);
xwcv = lagCovFFT(x_dtl,w_r,[],0,errMinFreq,sampleFreq);

% instrument white noise, Lenschow (2000)
N = sum(~isnan(w_r + x_dtl)); % # of good points
lags2fit = (1:5)';
xvar2fit = xvar(ismember(Elags,lags2fit));
xvar_sig = polyval(polyfit(lags2fit,xvar2fit,1),0); %extrapolate trend
xvar_noise = xvar(Elags==0) - xvar_sig; %noise variance
if xvar_noise<0, xvar_noise = 0; end %might be negative if relatively no noise
quality.xvarS2N = xvar_sig./xvar_noise; %Signal/Noise Ratio for variance
quality.xstd_noise = sqrt(xvar_noise); %white noise standard deviation
quality.REnoise = sqrt(xvar_noise.*std(w_r,'omitnan').^2./N);
% quality.REnoise = quality.REnoise./abs(wave.flux_avg); %fraction

% Total Random Error, Finkelstein and Sims (2001) and Mauder (2013)
i = abs(Elags)<=errMaxLag;
quality.REfs01 = sqrt(abs(sum(xvar(i).*wvar(i) + wxcv(i).*xwcv(i)))./N);
% quality.REfs01 = quality.REfs01./abs(wave.flux_avg); %fraction

% flight relative and systematic error, Lenshow (1994)
if exist('z','var') && exist('zi','var') && exist('Lleg','var')
    quality.SEturb = 2.2.*zi.*sqrt((z./zi))./Lleg; % Systematic Error, Eq. 64
    quality.REturb = 1.75.*(z./zi).^0.25.*sqrt(zi./Lleg); %Random Error, Eq. 66
end

% spectral corrections
RTcorr_ft = ResponseTimeCorrection(spectra.f,spectra.co,tauResponse);
RTcorr_wv = ResponseTimeCorrection(wave.freq,wave.co,tauResponse);

quality.SErt_ft = RTcorr_ft - 1; %fraction
quality.SErt_wv = RTcorr_wv - 1;

% statistics
stats.info  = 'median mean std min max N';
stats.x     = allStats(x);
stats.w     = allStats(w);
stats.x_dtl = allStats(x_dtl);
stats.w_r   = allStats(w_r);
stats.flux  = allStats((x_dtl-mean(x_dtl,'omitnan')).*(w_r-mean(w_r,'omitnan')));
stats.r_wx  = flux./(stats.w_r(3).*stats.x_dtl(3)); %correlation coefficient

%% OUTPUT
F.data      = data;
F.options   = options;
F.flux      = flux;
F.vex       = vex;
F.x_dt      = x_dt;
F.u_r       = u_r;
F.v_r       = v_r;
F.w_r       = w_r;
F.angles    = angles;
F.cov_wx    = cov_wx;
F.cov_wxdt  = cov_wxdt;
F.lags      = lags;
F.x_dtl     = x_dtl;
F.spectra   = spectra;
F.xtrend    = xtrend;
F.quality   = quality;
F.stats     = stats;
F.wave      = wave;

F = orderfields(F);


