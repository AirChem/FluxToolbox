function DivergencePlot(alt,lat,lon,leg,time,F,Ferr,div,derr,fit,name)
% function DivergencePlot(alt,lat,lon,leg,time,F,Ferr,div,derr,fit,name)
% generates some plots for evaluating flux divergence.
% note, map currently uses USEast shapefile.
% INPUTS:
% alt: altitude, m asl
% lat: latitude, deg
% lon: longitude, deg
% leg: leg index
% time: time, s
% F: flux time series
% Ferr: flux error time series
% div: optional divergence correction factors
% derr: error in div
% name: optional name of species
% fit: optional n x 5 matrix of fit coefficients, with columns [m b r sm sb]
%
% 20170511 GMW
% 20170801 GMW Corrected error averaging to RMS.

% defaults
if ~exist('name','var'), name = ''; end
if ~exist('fit','var'),  fit = []; end
if ~exist('div','var'),  div = ones(size(alt)); end
if ~exist('derr','var'), derr = zeros(size(alt)); end

% average data over each leg
legB = unique(leg);
legB(isnan(legB))=[];
x = [alt lat lon time F];
xavg = BinAvg(leg,x,legB);
altB = xavg(:,1); latB = xavg(:,2); lonB = xavg(:,3);
timeB = xavg(:,4); FB = xavg(:,5);
N = length(legB);
c = colormap(jet(N));

FerrB = sqrt(BinAvg(leg,Ferr.^2,legB));

% initialize figure
figure('name',[name ' Divergence'],...
    'position',[0.1 0.1 0.8 0.8],...
    'units','normalized',...
    'defaultaxesunits','normalized',...
    'defaulttextunits','normalized');
s1 = subplot(221); hold on; box on; grid on %vert profile
s2 = subplot(223); hold on; box on; grid on % map
s3 = subplot(222); hold on; box on; grid on % time series of fluxes
s4 = subplot(224); hold on; box on; grid on % divergence corrections

%initialize map
mlat = IBWread('USEast_outline_lat.ibw'); mlat=mlat.y;
mlon = IBWread('USEast_outline_lon.ibw'); mlon=mlon.y;
plot(s2,mlon,mlat,'-','color',[0.6 0.6 0.6])
plot(s2,lon,lat,'k-')
set(s2,'xlim',[min(lon) max(lon)],'ylim',[min(lat) max(lat)])

for i=1:N
    j = leg==i;
    
    % alt profile
    subplot(s1);
    h = ploterr(FB(i),altB(i),FerrB(i),[],'o');
    set(h,'color',c(i,:))
    text(FB(i),altB(i)*1.01,['L' num2str(i)],'color',c(i,:),'Units','data')
    
    % map
    subplot(s2)
    plot(s2,lon(j),lat(j),'color',c(i,:))
    text(nanmean(lon(j)),nanmean(lat(j)),['L' num2str(i)],'color',c(i,:),'Units','data')
    
    % time series
    subplot(s3)
    plot(s3,time(j),F(j),'color',c(i,:))
    plot(s3,time(j),F(j).*div(j),'--','color',c(i,:))
    text(nanmean(time(j)),max(F),['L' num2str(i)],'color',c(i,:),'Units','data')
    
    % correction factors
    subplot(s4)
    plot(time(j),div(j),'-','color',c(i,:))
    plot(time(j),div(j)+derr(j),':','color',c(i,:))
    plot(time(j),div(j)-derr(j),':','color',c(i,:))
    
end
xlabel(s1,'Flux')
ylabel(s1,'Alt (m)')
xlabel(s2,'Lon')
ylabel(s2,'Lat')
xlabel(s3,'Time')
ylabel(s3,'Flux')
xlabel(s4,'Time')
ylabel(s4,'F(0)/F(z)')

linkaxes([s3 s4],'x')

subplot(s3)
text(0.6,0.3,'- raw flux')
text(0.6,0.15,'-- corrected')

% add fits
if ~isempty(fit)
    subplot(s1)
    z = [0 max(altB)];
    for i=1:size(fit,1)
        plot(fit(i,1).*z + fit(i,2),z,'k-','linewidth',3)
    end
end


