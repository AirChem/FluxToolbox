function s = allStats(x)
% function s = allStats(x)
% Calculates statistics for a variable or variables.
% 
% INPUTS:
% x: data vector or matrix. If x is a matrix, each column is assumed to be a variable. 
%
% OUTPUTS:
% s is an N x 6 matrix (where N is # of variables) containing the following:
% median mean std min max #pts
%
% 20140324 GMW

[nrow,ncol] = size(x);

if isempty(x)==1
    s = nan(nrow,6);
    n = nan(nrow,1);
else
    s(:,1)  = nanmedian(x);
    s(:,2)  = nanmean(x);
    s(:,3)  = nanstd(x);
    s(:,4)  = min(x);
    s(:,5)  = max(x);
    s(:,6)  = sum(~isnan(x),1);
end


