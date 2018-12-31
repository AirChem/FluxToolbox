%CO_SIGNIF  Significance testing for the 2D Wavelet transform
%
%   [SIGNIF,FFT_THEOR] = ...
%      wave_signif(Y1,Y2,DT,SCALE,SIGLVL,DOF,MOTHER,PARAM)
%
% INPUTS:
%
%    Y = the time series, or, the VARIANCE of the time series.
%        (If this is a single number, it is assumed to be the variance...)
%    DT = amount of time between each Y value, i.e. the sampling time.
%    SCALE = the vector of scale indices, from previous call to WAVELET.
%
%
% OUTPUTS:
%
%    SIGNIF = significance levels as a function of SCALE
%    FFT_THEOR = output theoretical red-noise spectrum as fn of PERIOD
%
%
% OPTIONAL INPUTS:
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    LAG1 = LAG 1 Autocorrelation for Y1, used for SIGNIF levels. Default is 0.0
%    LAG2 = LAG 1 Autocorrelation for Y2
%
%    SIGLVL = significance level to use. Default is 0.95
%
%    DOF = degrees-of-freedom for signif test.
%         DOF = 2 (or 1 for MOTHER='DOG')
%
%
%%%% ADAPTED FROM TORRENCE AND COMPO WAVE_SIGNIF.M by RAH 20170728 --------


function [signif] = co_signif(Y1,Y2,dt,scale1,siglvl,dof,mother,param)

if (nargin < 8), param = -1; end
if (nargin < 7), mother = -1; end
if (nargin < 6), dof = -1; end
if (nargin < 5), siglvl = -1; end
if (nargin < 4)
	error('Must input a vector Y, sampling time DT, and SCALE vector')
end

n1 = length(Y1);
n2 = length(Y2);
J1 = length(scale1)-1;

scale(1:J1+1) = scale1;
s0 = min(scale);
dj = log(scale(2)/scale(1))/log(2.);

if (n1 == 1)
	var1 = Y1;
elseif (n2 ==1)
    var2 = Y2;
else
	var1 = std(Y1)^2;
    var2 = std(Y2)^2;
end

if (siglvl == -1), siglvl = 0.95; end
if (mother == -1), mother = 'MORLET'; end

mother = upper(mother);

% get the appropriate parameters [see Table(2)]
if (strcmp(mother,'MORLET'))  %----------------------------------  Morlet
	if (param == -1), param = 6.; end
	k0 = param;
	fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
	empir = [2.,-1,-1,-1]; Zv = 3.999;
	if (k0 == 6), empir(2:4)=[0.776,2.32,0.60]; end
elseif (strcmp(mother,'PAUL'))  %--------------------------------  Paul
	if (param == -1), param = 4.; end
	m = param;
	fourier_factor = 4*pi/(2*m+1);
	empir = [2.,-1,-1,-1]; Zv = 3.999;
	if (m == 4), empir(2:4)=[1.132,1.17,1.5]; end
elseif (strcmp(mother,'DOG'))  %---------------------------------  DOG
	if (param == -1), param = 2.; end
	m = param;
	fourier_factor = 2*pi*sqrt(2./(2*m+1));
	empir = [1.,-1,-1,-1]; Zv = 2.182;
	if (m == 2), empir(2:4) = [3.541,1.43,1.4]; end
	if (m == 6), empir(2:4) = [1.966,1.37,0.97]; end
else
	error('Mother must be one of MORLET,PAUL,DOG')
end

period = scale.*fourier_factor;
dofmin = empir(1);     % Degrees of freedom with no smoothing
Cdelta = empir(2);     % reconstruction factor
gamma_fac = empir(3);  % time-decorrelation factor
dj0 = empir(4);        % scale-decorrelation factor

freq = dt ./ period;   % normalized frequency

lag1 = corrcoef(Y1(1:end-1),Y1(2:end)); lag1 = lag1(1,2);
lag2 = corrcoef(Y2(1:end-1),Y2(2:end)); lag2 = lag2(1,2);

fft_theor1 = (1-lag1^2) ./ (1-2*lag1*cos(freq*2*pi)+lag1^2); % [Eqn(16)]
fft_theor2 = (1-lag2^2) ./ (1-2*lag2*cos(freq*2*pi)+lag2^2);

fft_theor  = var1*var2*sqrt(fft_theor1.*fft_theor2);
signif = fft_theor;

if (dof == -1), dof = dofmin; end

dof = dofmin;
signif = (Zv./dof)*fft_theor;  % [Eqn(31)]

return

