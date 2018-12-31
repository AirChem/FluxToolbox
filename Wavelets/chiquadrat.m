function xchi=chiquadrat(xi,layers)
%idealized lai profile cumulative vai
x=[1:20]; n=7;
yc=((x.^(n/2-1)).*exp(-x./2)./(gamma(n/2)*2^(n/2)));
xc=[1:20]*layers/20;
xchi=interp1(xc,yc,xi);%./sum(interp1(xc,yc,[1:layers]));
%xchi=xyci(xi);