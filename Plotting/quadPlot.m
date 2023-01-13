function [p,f] = quadPlot(wprime,xprime)
% function [p,f] = quadPlot(wprime,xprime)
% Generates a quadrant plot for evaluating scalar covariance.
%
% INPUTS:
% wprime: instantaneous deviation of vertical wind speed from mean.
% xprime: instantaneous deviation of scalar from mean (i.e. detrended).
%
% OUTPUTS:
% p: percentage of points in each quadrant
% f: flux calculated in each quadrant
%
% 20130920 GMW

%index quadrants
q1=wprime>0 & xprime>0;
q2=wprime<0 & xprime>0;
q3=wprime<0 & xprime<0;
q4=wprime>0 & xprime<0;

%calculate some statistics
n = [sum(q1) sum(q2) sum(q3) sum(q4)];
p = n./sum(n)*100; %percentage of points in each quadrant
f = [mean(wprime(q1).*xprime(q1),'omitnan'),... %average fluxe from each quadrant
    mean(wprime(q2).*xprime(q2),'omitnan'),...
    mean(wprime(q3).*xprime(q3),'omitnan'),...
    mean(wprime(q4).*xprime(q4),'omitnan')]; 

%plot
wl = max(abs(wprime))*1.5;
xl = max(abs(xprime))*1.5;

figure
hold on
box on
plot(wprime,xprime,'.')
plot([-wl wl],[0 0],'k-')
plot([0 0],[-xl xl],'k-')
xlabel('w''')
ylabel('x''')
xlim([-wl wl])
ylim([-xl xl])

fform = '%6.3g';
text(0.99,0.99,['p = ' num2str(p(1)) '%'],'horizontalAlignment','right','verticalAlignment','top','fontSize',14)
text(0.99,0.9,['w''x'' = ' num2str(f(1),fform)],'horizontalAlignment','right','verticalAlignment','top','fontSize',14)

text(0.01,0.99,['p = ' num2str(p(2)) '%'],'horizontalAlignment','left','verticalAlignment','top','fontSize',14)
text(0.01,0.9,['w''x'' = ' num2str(f(2),fform)],'horizontalAlignment','left','verticalAlignment','top','fontSize',14)

text(0.01,0.1,['p = ' num2str(p(3)) '%'],'horizontalAlignment','left','verticalAlignment','bottom','fontSize',14)
text(0.01,0.01,['w''x'' = ' num2str(f(3),fform)],'horizontalAlignment','left','verticalAlignment','bottom','fontSize',14)

text(0.99,0.1,['p = ' num2str(p(4)) '%'],'horizontalAlignment','right','verticalAlignment','bottom','fontSize',14)
text(0.99,0.01,['w''x'' = ' num2str(f(4),fform)],'horizontalAlignment','right','verticalAlignment','bottom','fontSize',14)

