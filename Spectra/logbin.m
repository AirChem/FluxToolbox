function [newx, newy, bin] = logbin(xdata, ydata, numint)
%function [newx, newy, bin] = logbin(xdata, ydata, numint)
% this function bins vectors so that they are plotted with even spacing on a log-scale...
% USE: [x,y] = logbin(xdata, ydata, numint);
% Optional output bin is matrix showing which points go into which bins.
% DKF 040126
% Stolen and modified, 071210 GMW.


logx = log(xdata);
logx(isinf(logx)) = nan; %deal with 0 in xdata
logmin = min(logx);
logmax = max(logx);
interval = abs(logmax - logmin)/numint;

intvec = exp(logmin + interval.*(0:numint)');
newx = nan(numint,1);
newy = nan(numint,1);
bin = zeros(length(xdata),numint);
for j = 1:numint
    a = find(xdata>intvec(j) & xdata<(intvec(j) + intvec(j+1)));
    newx(j) = mean(xdata(a),'omitnan');
    newy(j,:) = mean(ydata(a,:),1,'omitnan');
    
    bin(a,j) = 1;
end

% alternate code
% logxb = (logmin:interval:logmax)' + interval/2; %bin centers
% logxb = logxb(1:numint);
% newx = exp(logxb);
% newy = BinAvg(logx,ydata,logxb);
