function Xd = damp(X,t,alpha,tau)
% function Xd = damp(X,t,alpha,tau)
% Applies a correction to scalar measurements to simulate slow signal response time (i.e. due to
% stickiness on inlets).
% Taken from Nyugen et al., PNAS (2014).
%
% INPUTS:
% X: raw scalar (e.g. mixing ratio).
% t: time vector, seconds. Code assumes this is evenly spaced.
% alpha: pre-exponential factor for max-normalized pulse-response.
% tau: decay constant for max-normalized pulse-response.
%
% OUTPUT, Xd, is damped scalar values.
%
% 20140807 GMW

L = length(t);
nt = ceil(2*tau/nanmedian(diff(t))); %number of points in window

%pad front end
tp = [nan(nt-1,1); t(:)];
Xp = [nan(nt-1,1); X(:)];
Lp = length(Xp);
 
%build column index for chunks of data
i = repmat((1:Lp)',[1 nt]);
i = spdiags(i,0:-1:-Lp+nt); %grab diagonals
tbig = tp(i);
Xbig = Xp(i);
t0 = repmat(max(tbig),nt,1); %time for last point of each chunk

%calculate damped time series
e = exp(-(t0 - tbig)./tau);
dt = gradient(tbig);
dt(isnan(Xbig))=nan;
bot = nansum(dt.*e)';
top = nansum(dt.*e.*Xbig)';

Xd = X.*(1-alpha) + alpha.*top./bot;


% BRUTE FORCE METHOD
% L = length(t);
% Xd = X;
% for i=1:L
%     tnow = t(i);
%     Xnow = X(i);
%     j = t>=tnow-2*tau & t<=tnow & ~isnan(X);
%     if sum(j)<2, continue; end
%     
%     e = exp(-(tnow - t)./tau);
%     top = trapz(t(j),X(j).*e(j));
%     bot = trapz(t(j),e(j));
%     Xd(i) = Xnow.*(1-alpha) + alpha.*top./bot;
% end
