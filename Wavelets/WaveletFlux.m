function W = WaveletFlux(data,gapFill,plotMe,param,scrubGapBound,getError,errMinFreq,errMaxLag)
% function W = WaveletFlux(data,gapFill,plotMe,param,scrubGapBound,getError,errMinFreq,errMaxLag)
% Calculates the wavelet co-spectrum and related properties for two time series.
% Code is based on scripts from P. Misztal and T. Karl, which is based on WAVELET.m.
% For gory details, see: Torrence and Compo, "A Practical Guide to Wavelet Analysis," BAMS (1998).
% 
% NOTE that input time series, x and w, must be monotonized, despiked, detrended, deNaNed, rotated,
% etc.. There are some built-in methods to deal with gaps.
%
% INPUTS:
% data: structure containing the following:
%       x: First time series (e.g. temperature or species mixing ratio)
%       w: Second time series (e.g. vertical wind speed in m/s)
%       t: time vector, s
%       speed: platform speed, typically in m/s. used to calculate distance.
%       name: optional string name to append to plots and warnings.
% gapFill: specifier for how to deal with gaps in time series.
%          'none': do nothing. In this case, if there are gaps, you will get an error.
%          'stitch': remove gaps before transform. better for rare large random gaps.
%          'linear': linear interpolate. better for short, repetative gaps.
%          'nearest': nearest neighbor interpolation. Better for gaps near ends of dataset.
%          'covfill': predict missing data using local covariance with other time series.
%                     See covfill.m for more info on this method. better for everything.
% plotMe: OPTIONAL flag for generating plots. Default = 0.
% param: OPTIONAL structure containing optional inputs for WAVELET.m.
%        These inlude: padflag, dj, s0, j1, mother, Cdelta. For more info, see WAVELET.m documentation.
% scrubGapBound: OPTIONAL flag to remove wavelet-derived parameters in regions around data gaps, +/-
%                   the gap width. Default = 0 (no).
% getError: OPTIONAL flag for calculating variance-based uncertainty in flux time series.
%           Calculation is based on Finklestein and Sims (2001) and is VERY SLOW for fast data. Default = 0 (no).
% errMinFreq: lowest frequency to use for error calculations. Default = 0.
% errMaxLag: maximum lag to include in error calculation. Default = floor(N/2).
%
% OUTPUT is a structure, W, containing the following fields:
%
% W.param       wavelet parameters
% W.data        same as input structure
% W.time        time vector (same as input t)
% W.dist        distance vector (based on dt and speed)
% W.period      wavelet fourier period, s
% W.scale       wavelet scale
% W.freq        wavelet fourier freqency, Hz
% W.coi         cone of influence at each time
% W.w_wave      raw wavelet transform of w
% W.x_wave      raw wavelet transform of x
% W.w_power     2-D cospectral power of w, variance units
% W.x_power     as above but for x
% W.w_psd       time-averaged power spectral density of w
% W.x_psd       time averaged power spectral density of x
% W.w_var       time series of w variance
% W.x_var       time series of x variance
% 
% W.power       2-D co-spectral power of w and x
% W.power_sig   95% signficance levels for 2-D cospectrum
% W.co          time-averaged wavelet w-x cospectrum
% W.co_fft      fourier transform w-x cospectrum
% W.co_sig      95% confidence level for co
% W.flux        scale-averaged flux, units of w*x
% W.flux_avg    mean flux over whole leg
% W.flux_sig    95% confidence level for fluxes
% W.flux_err    1-sigma error in flux, same units as flux
% W.og          ogive of cospectrum
% W.og_fft      ogive of fft cospectrum
% W.flux_nocoi  flux without coi included
% W.co_nocoi    cospectrum without coi included
% W.og_nocoi    ogive for above
% W.qcoi     fraction of flux within coi
% W.quality     structure of quality markers for transform:
%               xvar,wvar: % deviation of reconstructed variance (Eq. 14)
%               xtime,wtime: % RMS deviation of reconstructed time series (Eq. 11)
%               flux: % deviation of wavelet-average and EC flux.
%               gapfraction: fraction of missing x-w data pairs
%
% 20140515 GMW
% 20140625 GMW  added functionality to deal with gaps in time series.
%               Also changed inputs.
% 20160921 GMW  added ogive output.
% 20160929 GMW  Modified for vector speed input.
% 20161128 GMW  Multiple modifications:
%               - Tweaked plotting to show period instead of distance, and moved to separate fn.
%               - Changed cospectral power to variance units and corrected for scale bias (Liu 2007)
%               - removed the "lag1" input. Now using w*x time series auto-correlation.
% 20161213 GMW  Altered times series and variance reconstruction warnings
%               Added x and w power spectra
%               Added _nocoi outputs
%               Added quality output
% 20170501 GMW  Added w_var and x_var outputs.
% 20170610 GMW  Added covfill and nearest gap-filling methods
%               Added scrubGapBound option
%               Added gapfraction quality output
% 20170801 GMW  Added error estimation


%% INPUT CHECKING AND DEFAULTS
if nargin<2, gapFill = 'none'; end
if nargin<3, plotMe = 0;       end
if nargin<4 || isempty(param), param = struct;   end
if nargin<5, scrubGapBound = 0; end
if nargin<6, getError = 0; end
if nargin<7, errMinFreq = 0; end
if nargin<8, errMaxLag = floor(length(data.t)/2); end

struct2var(data)
if ~isfield(data,'name'), name=''; end

if length(t)~=length(x) || length(t)~=length(w)
    error([name ' WaveletFlux: Inputs x and w are different lengths.'])
end

if ~isfield(data,'speed')
    speed = ones(size(t));
elseif ~any(length(speed)==[1 length(t)])
     error([name ' WaveletFlux: Input "speed" must be scalar or same length as t.'])
end

% fill gaps if needed
badx = isnan(x) | isinf(x);
badw = isnan(w) | isinf(w);
bad = badx | badw; %either missing
badxw = badx & badw; %both missing
quality.gapfraction = sum(bad)./length(x);
switch gapFill
    case 'none'
        if any(badx), error([name ' WaveletFlux: Input x contains NaNs or Inf.']); end
        if any(badw), error([name ' WaveletFlux: Input w contains NaNs or Inf.']); end
    case 'stitch'
        x(bad)=[]; w(bad)=[]; t(bad)=[];
    case 'linear'
        if any(badx), x = interp1(t(~badx),x(~badx),t,'linear','extrap'); end
        if any(badw), w = interp1(t(~badw),w(~badw),t,'linear','extrap'); end
    case 'nearest'
        if any(badx), x = interp1(t(~badx),x(~badx),t,'nearest','extrap'); end
        if any(badw), w = interp1(t(~badw),w(~badw),t,'nearest','extrap'); end
    case 'covfill'
        if any(badx), x = covfill(x,w); end
        if any(badw), w = covfill(w,x); end
        
        %might've missed some if both missing. clean up with interp
        if any(isnan(x)), x = interp1(t(~isnan(x)),x(~isnan(x)),t,'linear','extrap'); end
        if any(isnan(w)), w = interp1(t(~isnan(w)),w(~isnan(w)),t,'linear','extrap'); end    
        
    otherwise
        error([name ' WaveletFlux: gapFill input "' gapFill '" not recognized'])
end

% Ensure no gaps before proceeding
if any(isnan(x+w))
    error([name ' WaveletFlux: gaps not filled, cannot proceed.'])
end

time = t;
N = length(time);
dt = nanmedian(diff(time));
param.dt = dt; %save for later
dist = cumsum(dt.*speed); %cumulative distance

%% DEFAULT WAVELET PARAMETERS
defaults = {...
   'padflag'                    1;... % pad the time series with zeroes (recommended)
   'dj'                      0.25;... % this will do 1/dj sub-octaves per octave
   's0'                      2*dt;... % smallest scale.
   'j1'                        -1;... % # of scales - 1. Default is J1 = LOG2(N DT/S0)/DJ.
   'mother'              'Morlet';... % mother wavelet
   };

Pnames = fieldnames(param);
check = ismember(defaults(:,1),Pnames);
imiss = find(~check);
for i=1:length(imiss)
    param.(defaults{imiss(i),1}) = defaults{imiss(i),2};
end
struct2var(param) %break out parameters

if j1 == -1 %have to do this separately b/c it is a function of other params
    j1 = ceil(log2(N*dt/s0)/dj); %same as WAVELET.m default
    param.j1 = j1;
end

% get reconstruction factor (TC98, Table 2)
switch upper(mother)
    case 'MORLET'
        Cdelta = 0.776;
    case 'DOG'
        Cdelta = 3.541; %for m = 2 only ( Mexican hat)
    otherwise
        error('WaveletFlux:InvalidInput','No Cdelta specified for mother wavelet type "%s".',mother)
end
param.Cdelta = Cdelta; %store for posterity

%% GET SOME STATISTICS
w_var_avg = std(w).^2;
x_var_avg = std(x).^2;
wx = (w-mean(w)).*(x-mean(x));
wx_cov = mean(wx); %EC-equivalent flux

%% DO TRANSFORMS
[w_wave,period,scale,coi] = WAVELET(w,dt,padflag,dj,s0,j1,mother);
x_wave                    = WAVELET(x,dt,padflag,dj,s0,j1,mother);

scale = scale';
period = period';
freq = 1./period; %fourier frequency
df = -ndiff(freq);
scale_big  = repmat(scale,1,N);
period_big = repmat(period,1,N);
coi_big    = repmat(coi,length(scale),1);

%% APPLY TRANSFER FUNCTION TO W (experimental)
% load W858toRadTransferFunction.mat %fTr, Tr
% Tr = interp1(fTr,Tr,freq);
% Tr(isnan(Tr)) = 1; %lower frequencies
% w_wave = w_wave.*repmat(sqrt(Tr),1,N);

%% POWER SPECTRA
w_power = dj.*dt./Cdelta.*abs(w_wave).^2 ./scale_big; %variance units (Eq. 14) and bias-rectified
x_power = dj.*dt./Cdelta.*abs(x_wave).^2 ./scale_big;

w_psd = mean(w_power,2); %power spectra
x_psd = mean(x_power,2);

w_var = sum(w_power)'; % variance time series
x_var = sum(x_power)'; 

%% CROSS POWER SPECTRUM, TIME SERIES AND GLOBAL COSPECTRUM
power = real(w_wave).*real(x_wave) + imag(w_wave).*imag(x_wave); %cospectral power
power = dj.*dt./Cdelta.*power./scale_big; %bias rectification (divide by scale, see Liu (2007)) and variance scaling
flux = sum(power)';   % scale average [Eqn(24)] (divide by scale done above)
co = mean(power,2); % time average [Eqn(22)], but in units of variance and bias-corrected
og = ogive(freq,co./df); %ogive, /df to get spectral density

%% SIGNIFICANCE TESTING
% NOTE: using the lag-1 for wx may not be strictly correct (see e.g. Eqn (31)), but seems like a fair guess.
% I'm not smart enough to figure out how to get significance levels for the cross-spectrum.
lag1_wx = corrcoef(wx(1:end-1),wx(2:end)); lag1_wx = lag1_wx(1,2);

power_sig = wave_signif(wx_cov,dt,scale,0,lag1_wx,-1,-1,mother).*dj.*dt./Cdelta./scale'; %significance level
flux_sig  = wave_signif(wx_cov,dt,scale,2,lag1_wx,-1,[min(scale) max(scale)],mother)';
co_sig    = wave_signif(wx_cov,dt,scale,1,lag1_wx,-1,N - scale',mother)'.*dj.*dt./Cdelta./scale; % the -scale corrects for padding at edges (which determines coi)
% co_sig    = co_signif(wx_cov,dt,scale,lag1_wx,-1,mother,param);

%% COI INFLUENCES
icoi = period_big>coi_big; %periods greater than coi subject to edge effects

flux_nocoi = sum(power.*~icoi)';   % scale average [Eqn(24)] (divide by scale done above)
p = power; p(icoi)=nan;
co_nocoi = nanmean(p,2); % time average [Eqn(22)], but in units of variance and bias-corrected
og_nocoi = ogive(freq,co_nocoi./df); %ogive, /df to get spectral density

% og_big = ogive(freq,power.*period_big,1); %ogive time series (using absolute value)
% qcoi = max(og_big.*icoi)';

% qcoi = 1 - flux_nocoi./flux; %fraction of flux within the coi

og_abs = ogive(freq,co./df,1); %now with absolute value
qcoi = interp1(period,og_abs,coi)'; %fraction of global-average power within coi at each point in time

% og_nocoi_abs = ogive(freq,co_nocoi./df,1); %now with absolute value
% qcoi = interp1(period,og_nocoi_abs,coi)'; %fraction of global-average power within coi at each point in time

%% FOURIER TRANSFORM COSPECTRUM
[freq_fft co_fft] = cospectra(1./dt,x,w);
co_fft = flipud(BinAvg(log2(freq_fft),co_fft,flipud(log2(freq)))).*df; %*df gives units of variance
og_fft = ogive(freq,co_fft./df); %ogive,/df to get spectral density

%% RECONSTRUCT TIME SERIES AND VARIANCE
psi0 = pi.^-0.25; %for MORLET wavelet (Table 2)

%NOTE: WAVELET.m subtracts mean before doing transform
w_recon = dj.*sqrt(dt)./Cdelta./psi0.*sum(real(w_wave)./sqrt(scale_big))' + mean(w); %Eq. 11
% w_var_recon = dj.*dt./Cdelta./N.*sum(sum(abs(w_wave).^2./scale_big)); % Eq. 14
w_var_recon = mean(w_var);

x_recon = dj.*sqrt(dt)./Cdelta./psi0.*sum(real(x_wave)./sqrt(scale_big))' + mean(x); %Eq. 11
% x_var_recon = dj.*dt./Cdelta./N.*sum(sum(abs(x_wave).^2./scale_big)); % Eq. 14
x_var_recon = mean(x_var);

%% ISSUE WARNINGS
threshold_time = 30; %percent time series deviation
threshold_var  = 30; %percent variance deviation
threshold_flux = 30; %percent flux deviation

quality.wvar  = abs(1 - w_var_recon./w_var_avg)*100;
quality.xvar  = abs(1 - x_var_recon./x_var_avg)*100;
quality.wtime = sqrt(mean((w_recon - w).^2))./sqrt(mean(w.^2))*100; %have to rms w because its near 0
quality.xtime = sqrt(mean((x_recon - x).^2))./sqrt(mean(x.^2))*100;
quality.flux  = abs(1 - mean(flux)./wx_cov)*100;

% if quality.wvar > threshold_var, disp([name ' WaveletFlux: original and reconstructed w variance differ by ' num2str(quality.wvar,'%3.0f') '%.']); end
% if quality.xvar > threshold_var, disp([name ' WaveletFlux: original and reconstructed x variance differ by ' num2str(quality.xvar,'%3.0f') '%.']); end
% if quality.wtime > threshold_time, disp([name ' WaveletFlux: original and reconstructed w time series differ by ' num2str(quality.wtime,'%3.0f') '%.']); end
% if quality.xtime > threshold_time, disp([name ' WaveletFlux: original and reconstructed x time series differ by ' num2str(quality.xtime,'%3.0f') '%.']); end
if quality.flux > threshold_flux, disp([name ' WaveletFlux: EC and Wavelet flux differ by ' num2str(quality.flux,'%3.0f') '%.']); end
    
%% RE-INSERT GAPS
switch gapFill
    case 'stitch'
        N = length(data.t); %original length
        Ns = length(scale);
        good = ~bad;
        
        time = data.t;
        w = data.w;
        x = data.x;
        
        temp = nan(N,1);    temp(good) = coi;           coi = temp;
        temp = nan(N,1);    temp(good) = qcoi;       qcoi = temp;
        temp = nan(Ns,N);   temp(:,good) = w_wave;      w_wave = temp;
        temp = nan(Ns,N);   temp(:,good) = x_wave;      x_wave = temp;
        temp = nan(Ns,N);   temp(:,good) = w_power;     w_power = temp;
        temp = nan(Ns,N);   temp(:,good) = x_power;     x_power = temp;
        temp = nan(Ns,N);   temp(:,good) = power;       power = temp;
        temp = nan(N,1);    temp(good) = w_var;         w_var = temp;
        temp = nan(N,1);    temp(good) = x_var;         x_var = temp;
        temp = nan(N,1);    temp(good) = flux;          flux = temp;
        temp = nan(N,1);    temp(good) = flux_nocoi;    flux_nocoi = temp;
        temp = nan(N,1);    temp(good) = flux_sig;      flux_sig = temp;
        
    case {'linear','covfill','nearest'}
        w_power(:,badw) = nan;
        x_power(:,badx) = nan;
        power(:,bad) = nan;
        w_var(badw) = nan;
        x_var(badw) = nan;
        flux(bad) = nan;
        flux_nocoi(bad) = nan;
        flux_sig(bad) = nan;
        x_wave(:,badx) = nan;
        w_wave(:,badx) = nan;
end

%% ADDITIONAL REMOVAL AROUND GAPS
if scrubGapBound
    
    %w
    g = chunker(find(badw)); %start and stop indices for gap chunks
    n = diff(g,1,2)+1; %gap width
    r = g + [-n n]; %region around gaps
    r(r<1)=1; r(r>N)=N; %boundaries
    for i=1:length(n)
        badw(r(i,1):r(i,2)) = 1; %add to bad
    end
    
    %x
    g = chunker(find(badx)); %indices for gap chunks
    n = diff(g,1,2)+1; %gap width
    r = g + [-n n]; %region around gaps
    r(r<1)=1; r(r>N)=N; %boundaries
    for i=1:length(n)
        badx(r(i,1):r(i,2)) = 1; %add to bad
    end
    
    bad = badx | badw;
    
    %scrubadub
    w_wave(:,badw) = nan;
    x_wave(:,badx) = nan;
    w_power(:,badw) = nan;
    x_power(:,badx) = nan;
    power(:,bad) = nan;
    w_var(badw) = nan;
    x_var(badw) = nan;
    flux(bad) = nan;
    flux_nocoi(bad) = nan;
    flux_sig(bad) = nan;
end

%% FLUX ERROR
% This code is based on the work of Finklestein and Sims (2001), extended to Wavelet land.
if getError
    flux_err = WaveletFluxError(x_wave,w_wave,scale,freq,param,errMinFreq,errMaxLag);
else
    flux_err = nan(size(flux));
end

%% OUTPUTS

W = struct;

W.data      = data;
W.time      = time;
W.dist      = dist;
W.period    = period;
W.scale     = scale;
W.freq      = freq;
W.coi       = coi;
W.qcoi   = qcoi;
W.w_power   = w_power;
W.x_power   = x_power;
W.w_psd     = w_psd;
W.x_psd     = x_psd;
W.w_var     = w_var;
W.x_var     = x_var;

W.power     = power;
W.power_sig = power_sig;
W.co        = co;
W.co_fft    = co_fft;
W.co_sig    = co_sig;
W.flux      = flux;
W.flux_sig  = flux_sig;
W.flux_avg  = nanmean(flux);
W.flux_err  = flux_err;
W.param     = param;
W.og        = og;
W.og_fft    = og_fft;
W.flux_nocoi= flux_nocoi;
W.flux_avg_nocoi  = nanmean(flux_nocoi);
W.co_nocoi  = co_nocoi;
W.og_nocoi  = og_nocoi;
W.quality   = quality;

W.x_wave = x_wave;
W.w_wave = w_wave;

% W = orderfields(W);

%% PLOTS
if plotMe
    WaveSummaryPlot(W,name);
end


