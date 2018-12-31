function xlagged = lagVar(x,lag)
% function xlagged = lagVar(x,lag)
% Lags a variable forward or backward by a specified number of points.
% 
% INPUTS:
% x: variable to lag. must be a column vector.
% lag: number of points to lag. Positive is forward, negative is backward.
%
% OUTPUTS:
% xlagged: lagged variable x, padded by nans.
%
%20130915 GMW

x = x(:); %convert to column if needed

n = nan(abs(lag),1); %filler
if lag<0
    xlagged = [x(-lag+1:end); n]; %shift backward
elseif lag>0
    xlagged = [n; x(1:end-lag)]; %shift forward
else
    xlagged = x;
end


