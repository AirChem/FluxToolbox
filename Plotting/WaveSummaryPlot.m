function h = WaveSummaryPlot(W,name,Tunit)
% function h = WaveSumarymPlots(W,name,Tunit)
% Generates a single-panel visualization of important flux parameters, including:
% Time series of c' and w', wavelet cospectra and fluxes
% lag covariance
% leg-averaged cospectra and ogives
% INPUTS:
% W: structure containing variables as output from WaveletFlux.
% name: optional name to give to figure.
% Tunit: optional specifier for independent coordinate ('time' or 'dist')
%
% OUTPUT h is the set of handles for the figure and all subplots.
%
% 20160929 GMW
% 20161128 GMW  Modified to accept wavelet structure instead of flux structure.
% 20170108 GMW  Added Tunit input
%               Added power spectra

if nargin<2, name=''; end
if nargin<3, Tunit = 'time'; end

N = length(W.time);

switch Tunit
    case 'dist'
        T = W.dist;
        Tlab = 'Distance';
        F = W.freq./mean(W.data.speed);
        Flab = 'Wave#';
        COI = 1./(W.coi.*mean(W.data.speed)); %from period to wvnm
    otherwise %default is to use time
        T = W.time;
        Tlab = 'Time';
        F = W.freq;
        Flab = 'Freq';
        COI = 1./W.coi; %from period to freq
end

% Time series
w_plot = (W.data.w - mean(W.data.w,'omitnan'))./std(W.data.w,'omitnan');
x_plot = (W.data.x - mean(W.data.x,'omitnan'))./std(W.data.x,'omitnan');

% flux variables
[cov_wx,lags] = lagCovFFT(W.data.w,W.data.x,[]); %note, x and w inputs already lagged
ECflux = cov_wx(lags==0);

% position info
fpos = [0.01 0.05 0.96 0.86]; %figure
h1pos = [0.05 0.70 0.62 0.29]; %inputs
h2pos = [0.05 0.40 0.62 0.29]; %wavelet power
h3pos = [0.05 0.10 0.62 0.29]; %flux time series

h4pos = [0.74 0.10 0.24 0.18]; %lag
h5pos = [0.74 0.70 0.24 0.29]; %power spectra
h6pos = [0.74 0.40 0.24 0.29]; %cospectrum

%start figure
figure('name',name)
set(gcf,'units','normalized','position',fpos)

% time series
h1 = subplot(3,3,[1 2]);
set(gca,'units','normalized','position',h1pos);
plot(T,x_plot,'k',T,w_plot,'m')
xlim([min(T) max(T)])
ylim([min([w_plot;x_plot]) max([w_plot;x_plot])])
ylabel('(x - <x>)/\sigma_x')
text(0.01,0.98,'w''','color','m')
text(0.05,0.98,'x''','color','k')
set(gca,'XTickLabel',[])

% Contour plot wavelet power spectrum
pnorm = max(max(abs(W.power)));
power_norm = W.power./pnorm; %scale from -1 to 1
power_norm = power_norm - 1; %must adjust for significance contours to plot properly
y = log10(F);

h2 = subplot(3,3,[4 5]);
set(gca,'units','normalized','position',h2pos);
surface(T,y,power_norm);
shading interp
ylabel(Flab)
xlim([min(T) max(T)])
ylim([min(y) max(y)])
Yticks = floor(min(y)):ceil(max(y));
set(gca,'YTick',Yticks(:),'YTickLabel',10.^Yticks)
set(gca,'XTickLabel',[])
c = interp1([1;32;64],[0 0 1; 1 1 1; 1 0 0],(1:64)'); %b-w-r
colormap(c)
caxis([-2 0])
box on

%contours
hold on
power2sig = W.power./(W.power_sig'*(ones(1,N)));  % where ratio > 1, power is significant
power2sig(power2sig<=0) = 0;
contour(T,y,power2sig,[1 1],'k','LineWidth',2)  % 95% significance contour
plot3(T,log10(COI),ones(1,N),'--','LineWidth',4,'Color',[0.5 0 0.5]) % cone-of-influence

%legend
text(0.40,0.14,'+ power','Color','r')
text(0.53,0.14,'- power','Color','b')

% scale-average time series
h3 = subplot(3,3,[7 8]);
set(gca,'units','normalized','position',h3pos);
plot(T,W.flux,'b','LineWidth',3);
xlim([min(T) max(T)])
ylim([min(W.flux) max(W.flux)])
xlabel(Tlab)
ylabel('Flux')
hold on
plot(xlim,ECflux + [0 0],'c--','LineWidth',3)
text(0.40,0.95,['<wave>: ' num2str(mean(W.flux,'omitnan'),'%3.2g')],'Color','b')
text(0.40,0.86,['<w''x''>: ' num2str(ECflux,'%3.2g')],'Color','c')
text(0.40,0.75,['Ratio: ' num2str(mean(W.flux,'omitnan')./ECflux,'%2.2f')])

i=W.qcoi>0.5;
plot(T(i),W.flux(i),'kx') %mark fluxes with majority influence w/in coi
text(0.85,0.95,'x: COI>50%')

linkaxes([h1 h2 h3],'x')

% lag covariance
h4 = axes('units','normalized','position',h4pos);
plot(lags,cov_wx,'-')
hold on
plot([0 0],[min(cov_wx) max(cov_wx)],'k:')
xlabel('Lag points')
ylabel('<w''x''>')
axis tight

% power spectra
h5 = axes('units','normalized','position',h5pos);
loglog(F,W.w_psd./sum(W.w_psd),'m-','LineWidth',3)
hold on
loglog(F,W.x_psd./sum(W.x_psd),'k-','LineWidth',3)
ylabel('Pw(x'' or w'')*f')
xlo = 10.^floor(log10(min(F)));
xhi = 10.^ceil(log10(max(F)));
axis tight
xlim([xlo xhi])
set(gca,'XTick',10.^(log10(xlo):log10(xhi)),'xticklabel',[])
y1 = max([W.w_psd./sum(W.w_psd);W.x_psd./sum(W.x_psd)]); %starting position for f-slope line
f23 = fslopeline(xlim,y1,-2/3,0);
plot(xlim,f23*2.5,':','color',[0.6 0.6 0.6])
plot(xlim,f23*5,':','color',[0.6 0.6 0.6])
plot(xlim,f23*10,':','color',[0.6 0.6 0.6])
text(0.65,0.17,'f^{-2/3}','color',[0.6 0.6 0.6])
text(0.45,0.15,'w''','color','m')
text(0.55,0.15,'x''','color','k')

% cospectrum
h6 = axes('units','normalized','position',h6pos);
semilogx(F,W.co./max(abs(W.co)),'b-','LineWidth',3)
hold on
semilogx(F,W.co_fft./max(abs(W.co)),'c-','LineWidth',3)
ylabel('Co(w''x'')*f')
xlo = 10.^floor(log10(min(F)));
xhi = 10.^ceil(log10(max(F)));
axis tight
xlim([xlo xhi])
plot([xlo xhi],[0 0],'k:') %zero line

% ogive
semilogx(F,W.og,'b--','LineWidth',2)
semilogx(F,W.og_fft,'c--','LineWidth',2)
xlabel(Flab)
set(gca,'XTick',10.^(log10(xlo):log10(xhi)))
text(0.45,0.3,'wave','color','b')
text(0.45,0.15,'FFT','color','c')

linkaxes([h5 h6],'x')

%output
if nargout
    h = [gcf h1 h2 h3 h4 h5 h6]; %handles
end

