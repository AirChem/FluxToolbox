function W = WaveletFlux(data,gapFill,plotMe,param)
% function W = WaveletFlux(data,gapFill,plotMe,param)
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
%       speed: aircraft speed, m/s. Optional. If not specified, 1 m/s assumed.
% gapFill: specifier for how to deal with gaps in time series.
%          'none': do nothing. In this case, if there are gaps, you will get an error.
%          'stitch': remove gaps before transform. better for rare large random gaps.
%          'interp': interpolate. better for short, repetative gaps.
% plotMe: OPTIONAL flag for generating plots. Default = 0.
% param: OPTIONAL structure containing optional inputs for WAVELET.m.
%        These inlude: pad, dj, s0, j1, mother. For more info, see WAVELET.m documentation.
%        Can also specify the following:
%        lag1: lag-1 autocorrelation for noise, used in the wave_signif function.
%              Default is 0.72 (red noise).
%        Cdelta: empirical time series coefficient. Default is 0.776 for the MORLET wavelet
%
% OUTPUT is a structure, W, containing the following fields:
%
% W.param       wavelet parameters
% W.x           input x
% W.w           input w
% W.time        time vector (based on dt and # of points)
% W.dist        distance vector (based on time and speed)
% W.period      wavelet fourier period, s
% W.scale       wavelet scale, s
% W.freq        wavelet freqency, Hz
% W.coi         cone of influence at each time/dist
% W.w_wave      raw wavelet transform of w
% W.x_wave      raw wavelet transform of x
% 
% W.power       2-D co-spectral power of w and x
% W.power_sig   95% signficance levels for 2-D cospectrum
% W.co          time-averaged wavelet w-x cospectrum
% W.co_fft      fourier transform w-x cospectrum
% W.co_sig      95% confidence level for co
% W.flux        scale-averaged flux, units of w*x
% W.flux_sig    95% confidence level for fluxes
% W.og          ogive of cospectrum
%
% 20140515 GMW
% 20140625 GMW  added functionality to deal with gaps in time series.
%               Also changed inputs.
% 20160921 GMW  added ogive output.
% 20160929 GMW  Modified for vector speed input.

%% INPUT CHECKING AND DEFAULTS
if nargin<2, gapFill = 'none'; end
if nargin<3, plotMe = 0;       end
if nargin<4 || isempty(param), param = struct;   end

struct2var(data)
if length(t)~=length(x) || length(t)~=length(w)
    error('WaveletFlux: Inputs x and w are different lengths.')
end

if ~exist('speed','var')
    speed = ones(size(t));
end

% fill gaps if needed
switch gapFill
    case 'none'
        if any(isnan(x) | isinf(x)), error('WaveletFlux: Input x contains NaNs.'); end
        if any(isnan(w) | isinf(w)), error('WaveletFlux: Input w contains NaNs.'); end
    case 'stitch'
        bad = isnan(x+w);
        x(bad)=[]; w(bad)=[]; t(bad)=[];
    case 'interp'
        badx = isnan(x);
        if sum(badx), x = interp1(t(~badx),x(~badx),t,'linear','extrap'); end
        badw = isnan(w);
        if sum(badw), w = interp1(t(~badw),x(~badw),t,'linear','extrap'); end
    otherwise
        error(['WaveletFlux: gapFill input "' gapFill '" not recognized'])
end

time = t;
N = length(time);
dt = nanmedian(diff(time));
mean_speed = nanmean(speed);
if length(speed)>1
dist = [0; cumsum(diff(time).*speed(1:end-1))];
else
    dist = [0; cumsum(diff(time).*mean_speed)];
end

%% DEFAULT WAVELET PARAMETERS
defaults = {...
   'pad'                        1;... % pad the time series with zeroes (recommended)
   'dj'                      0.25;... % this will do x sub-octaves per octave
   's0'                      2*dt;... % smallest scale.
   'j1'                        -1;... % # of scales - 1. Default is J1 = LOG2(N DT/S0)/DJ.
   'mother'              'Morlet';... % mother wavelet (they are asexual)
   'lag1'                    0.72;... % lag-1 autocorrelation for red noise background
   'Cdelta'                 0.776;... % empirical scaling factor. Default is 0.776 for Morlet.
   };

Pnames = fieldnames(param);
check = ismember(defaults(:,1),Pnames);
imiss = find(~check);
for i=1:length(imiss)
    param.(defaults{imiss(i),1}) = defaults{imiss(i),2};
end
struct2var(param) %break out parameters

if j1 == -1 %have to do this separately b/c it is a function of other params
    j1 = round(log2(N*dt/s0)/dj); %same as WAVELET.m default
    param.j1 = j1;
end

%% GET SOME STATISTICS
w_std = nanstd(w);
w_var = w_std.^2;
x_std = nanstd(x);
x_var = x_std.^2;
wx_cov = nanmean(w.*x); %EC-equivalent flux

%% DO TRANSFORMS
[w_wave,period,scale,coi] = WAVELET(w,dt,pad,dj,s0,j1,mother);
x_wave                    = WAVELET(x,dt,pad,dj,s0,j1,mother);
period = period'; %make column vector and flip to make trapz integration give the right sign
scale = scale';
freq = 1./period; %fourier frequency

%% COSPECTRUM
power = real(w_wave).*real(x_wave) + imag(w_wave).*imag(x_wave); %cospectral power
power_sig = wave_signif(wx_cov,dt,scale,0,lag1,-1,-1,mother); %significance level

%% SCALE AVERAGE
scale_big = scale*ones(1,N); %expand size to x_wave and w_wave
flux = dj.*dt./Cdelta.*sum(power./scale_big)';   % [Eqn(24)]
flux_sig = wave_signif(wx_cov,dt,scale,2,lag1,-1,[min(scale) max(scale)],mother);

% removal of portion below COI (not tested)
% coi_big = repmat(coi,length(scale),1);
% icoi = scale_big<=coi_big;
% f0 = power./scale_big;
% f0(~icoi) = nan;
% flux_nocoi = dj.*dt./Cdelta.*nansum(f0)';   % [Eqn(24)]

%% TIME AVERAGE
co = sum(power,2)./N;
co_norm = nanmean(flux)./-trapz(freq,co); %scaling factor for flux units
co = co.*co_norm; %normalize for flux units
og = flipud(ogive(flipud(freq),flipud(co))); %ogive
dof = N - scale';  % the -scale corrects for padding at edges
co_sig = wave_signif(wx_cov,dt,scale,1,lag1,-1,dof,mother);
co_sig = co_sig'.*co_norm;

%% FOURIER TRANSFORM COSPECTRUM
[freq_fft co_fft] = cospectra(1./dt,x,w);
co_fft = BinAvg(log10(freq_fft),co_fft,flipud(log10(freq)));
og_fft = ogive(flipud(freq),co_fft); %ogive
co_fft = flipud(co_fft);
og_fft = flipud(og_fft);

%% RECONSTRUCT TIME SERIES AND VARIANCE
psi0 = pi.^-0.25; %for MORLET wavelet (Table 2)

w_recon = dj.*sqrt(dt)./Cdelta./psi0.*sum(real(w_wave)./sqrt(scale_big))'; %Eq. 11
w_var_recon = dj.*dt./Cdelta./N.*sum(sum(abs(w_wave).^2./scale_big)); % Eq. 14
w_mse = nanmean((w_recon - w).^2); %mean square error

x_recon = dj.*sqrt(dt)./Cdelta./psi0.*sum(real(x_wave)./sqrt(scale_big))'; %Eq. 11
x_var_recon = dj.*dt./Cdelta./N.*sum(sum(abs(x_wave).^2./scale_big)); % Eq. 14
x_mse = nanmean((x_recon - x).^2);

%% ISSUE WARNINGS
threshold_mse  = 1; %percent
threshold_var  = 10; %percent
threshold_flux = 10; %percent

check_wvar = abs(1 - w_var_recon./w_var)*100;
check_xvar = abs(1 - x_var_recon./x_var)*100;
check_wmse = w_mse./w_var*100;
check_xmse = x_mse./x_var*100;
check_flux = abs(1 - nanmean(flux)./wx_cov)*100;

% if check_wvar > threshold_var, disp(['WaveletFlux: original and reconstructed w variance differ by ' num2str(check_wvar,'%3.0f') '%.']); end
% if check_xvar > threshold_var, disp(['WaveletFlux: original and reconstructed x variance differ by ' num2str(check_xvar,'%3.0f') '%.']); end
% if check_wmse > threshold_mse, disp(['WaveletFlux: original and reconstructed w time series differ by ' num2str(check_wmse,'%3.0f') '%.']); end
% if check_xmse > threshold_mse, disp(['WaveletFlux: original and reconstructed x time series differ by ' num2str(check_xmse,'%3.0f') '%.']); end
if check_flux > threshold_flux, disp(['WaveletFlux: EC and Wavelet flux differ by ' num2str(check_flux,'%3.0f') '%.']); end

%% UN-STITCH, IF NEEDED
if strcmp(gapFill,'stitch')
    N = length(data.t); %original length
    Ns = length(scale);
    good = ~bad;
    
    time = data.t;
    w = data.w;
    x = data.x;
    
    temp = nan(N,1);
    temp(good) = coi;
    coi = temp;
    
    temp = nan(Ns,N);
    temp(:,good) = w_wave;
    w_wave = temp;
    
    temp = nan(Ns,N);
    temp(:,good) = x_wave;
    x_wave = temp;
    
    temp = nan(Ns,N);
    temp(:,good) = power;
    power = temp;
    
    temp = nan(N,1);
    temp(good) = flux;
    flux = temp;
    
    temp = nan(N,1);
    temp(good) = flux_sig;
    flux_sig = temp;
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
W.w_wave    = w_wave;
W.x_wave    = x_wave;

W.power     = power;
W.power_sig = power_sig;
W.co        = co;
W.co_fft    = co_fft;
W.co_sig    = co_sig;
W.flux      = flux;
W.flux_sig  = flux_sig;
W.param     = param;
W.og        = og;
W.og_fft    = og_fft;

W = orderfields(W);

%% PLOTS
if plotMe
    
    h1pos = [0.07 0.70 0.85 0.29];
    h2pos = [0.07 0.40 0.85 0.29];
    h3pos = [0.07 0.10 0.85 0.29];
    
    % Time series
    w_norm = (w - nanmean(w))./w_std;
    x_norm = (x - nanmean(x))./x_std;

    figure
    h1 = subplot('position',h1pos);
    plot(dist,w_norm,'k',dist,x_norm,'c')
    xlim([min(dist) max(dist)])
    ylim([min([w_norm;x_norm]) max([w_norm;x_norm])])
    ylabel('(x - <x>)/\sigma_x')
    legend('w','x','Location','northwest')
    legend boxoff
    set(gca,'XTickLabel',[])
    
    % Contour plot wavelet power spectrum
    h2 = subplot('position',h2pos);
    pnorm = max(max(abs(power)));
    power_norm = power./pnorm; %scale from -1 to 1
    power_norm = power_norm - 1; %must adjust for significance contours to plot properly
    surface(dist,log10(period),power_norm);
    shading interp
    ylabel('Period (s)')
    xlim([min(dist) max(dist)])
    ylim(log10([min(period) max(period)]))
    set(gca,'YDir','reverse')
    Yticks = 10.^(floor(log10(min(period))):ceil(log10(max(period))))';
    set(gca,'YTick',log10(Yticks(:)),'YTickLabel',Yticks)
    set(gca,'XTickLabel',[])
%     c = interp1([1;20;44;64],[1 1 1;0 0 1;1 0 0;0 0 0],(1:64)'); %custom colormap, w-b-r-k
    c = interp1([1;32;64],[0 0 1; 1 1 1; 1 0 0],(1:64)'); %b-w-r
    colormap(c)
    
    hold on
    power2sig = power./(power_sig'*(ones(1,N)));  % where ratio > 1, power is significant
    power2sig(power2sig<=0) = 0;
    contour(dist,log10(period),power2sig,[1 1],'k','LineWidth',2)  % 95% significance contour
    caxis([-2 0])
    box on
    
    plot3(dist,log10(coi),ones(1,N),'--','LineWidth',5,'Color',[0 0.6 0]) % cone-of-influence
    
    text(0.01,0.85,'+ power','Color','r')
    text(0.01,0.75,'- power','Color','b')
    
    %right ticks
    space = period.*mean_speed/1000; %spatial scale, km
    YtickRight = 10.^(floor(log10(min(space))):ceil(log10(max(space))))'; %powers of 10
    LinkAxisData('y',log10(YtickRight*1000/mean_speed),YtickRight,'Length (km)');
    
    % scale-average time series
    h3 = subplot('position',h3pos);
    plot(dist,flux,'b');
    xlim([min(dist) max(dist)])
    xlabel('Distance (km)')
    ylabel('Flux')
    hold on
    plot(xlim,wx_cov + [0 0],'r-')
    text(0.01,0.8,['<wave>: ' num2str(nanmean(flux),'%3.2g')],'Color','b')
    text(0.01,0.66,['<w''x''>: ' num2str(wx_cov,'%3.2g')],'Color','r')
    text(0.01,0.52,['Ratio: ' num2str(nanmean(flux)./wx_cov,'%2.2f')])
    
    linkaxes([h1 h2 h3],'x')
    
    % cospectrum
    figure
    semilogx(freq,freq.*co,'b-','LineWidth',3)
    hold on
    semilogx(freq,freq.*co_fft,'r-','LineWidth',3)
    xlabel('Frequency (Hz)')
    ylabel('Co(wx)*f')
    legend('Wave','FFT')
    legend boxoff
    xlo = 10.^floor(log10(min(freq)));
    xhi = 10.^ceil(log10(max(freq)));
    xlim([xlo xhi])
    set(gca,'XTick',10.^(log10(xlo):log10(xhi)))
    plot([xlo xhi],[0 0],'k:') %zero line
    
    topTickLoc = get(gca,'XTick');
    topTickLabel = mean_speed./topTickLoc/1000; %km
%     LinkTopAxisData(topTickLoc,topTickLabel,'Length (km)');
    LinkAxisData('x',topTickLoc,topTickLabel,'Length (km)');
    
    set(gca,'Position',[0.15 0.15 0.76 0.72])
end


