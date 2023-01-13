function p = nanpolyfit(x,y,n)
%function p = nanpolyfit(x,y,n)
%Just like polyfit, but removes the nans first.
%090212 GMW

i = find(isnan(x+y));
x(i)=[];
y(i)=[];

p = polyfit(x,y,n);