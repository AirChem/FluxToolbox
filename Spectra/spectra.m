function [f psdx psdxnorm] = spectra(Fs,x)
%function [f psdx psdxnorm] = spectra(Fs,x)
%calculates spectrum and power spectrum for a data time series.
%INPUTS:
%x: data vector as a column
%Fs: Sampling Frequency in Hz
%OUTPUTS:
%f: frequency vector for spectra
%psdx: power spectral density of x
%psdxnorm: psdx normalized to total variance
%Lifted from Ian's spec_all.m file.
%071128 GMW

x = x - mean(x,'omitnan');

%% CALCULATE FREQEUNCY
N=length(x);
df=Fs/N; %spacing
Ny = floor(N/2); %Nyquist length
f =(0:df:Fs)';
f = f(2:Ny);
 
%% X POWER SPECTRA
xx   = fft(x,N);            % spectrum
pxx  = xx.*conj(xx)/N^2;    % power spectrum
psdx = 2*pxx(2:Ny)./df;     % power spectral density 
psdx = psdx(:);

%% COMPARE VARIANCES
varx = (std(x))^2;
var2x = sum(psdx)*df;
reldif=abs(varx-var2x)/varx;
if reldif > .01
   warning(['Relative difference between calculated and psd-estimated '...
       'variance is:  ', num2str(reldif*100), ' %'])
end

%normalize by variance
psdxnorm = psdx./var2x;

