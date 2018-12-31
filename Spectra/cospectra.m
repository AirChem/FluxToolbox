function [f co co_norm qu qu_norm] = cospectra(Fs,x,y)
% function [f co co_norm qu qu_norm coherence phase] = cospectra(Fs,x,y)
% calculates cospectrum for two data time series.
% INPUTS:
% x:  first data vector as a column
% y:  second data vector
% Fs: Sampling Frequency in Hz
%
% OUTPUTS:
% f: frequency vector for spectra
% co: x-y cospectral density
% co_norm: co normalized to total covariance absolute value
% qu: quadrature spectral density (complex conjugate of co)
% qu_norm: qu normalized by its integral absolute value
%
% 20140319 GMW

x = x - nanmean(x);
y = y - nanmean(y);

%% CALCULATE FREQEUNCY
N=length(x);
df=Fs/N; %spacing
Ny = floor(N/2); %Nyquist length
f =(0:df:Fs)';
f = f(2:Ny);

%%%%%CALCULATE FOURIER TRANSFORMS AND SPECTRA%%%%%
xx = fft(x,N);              % spectrum
pxx  = xx.*conj(xx)/N^2;    % power spectrum
psdx = 2*pxx(2:Ny)./df;     % power spectral density 
psdx = psdx(:);

yy = fft(y,N);
pyy  = yy.*conj(yy)/N^2;    % power spectrum
psdy = 2*pyy(2:Ny)./df;     % power spectral density 
psdy = psdy(:);

cross = xx(2:Ny).*conj(yy(2:Ny))/(N^2*df);
cross = cross(:);

co = 2*real(cross);
qu = 2*imag(cross);

%%%%%COMPARE COVARIANCES%%%%%
covxy=nanmean(x.*y);
cov2xy=sum(co)*df;
reldif=abs(covxy-cov2xy)/covxy;
if reldif > .01
    warning(['Relative difference between calculated and psd-estimated '...
        'covariance is:  ', num2str(reldif*100), ' %'])
end

%normalize by variance
co_norm = co./abs(cov2xy);
qu_norm = qu./abs((sum(qu)*df));

