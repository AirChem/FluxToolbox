function xfill = covfill(x,w)
% function xfill = covfill(x,w)
% Uses the covariance of one variable with another to fill in gaps in the former.
% Specifically, the algorithm fills in gaps in x such that the w-x covariance in the gap is the same
% as that in the region outside the gap. Only works in regions where this is a gap in x but NOT w.
% Originally designed to fill time series before performing wavelet transforms.
% Steps:
% 1) identify gaps in x
% 2) for each gap, define a window bounded by +/- the gap size and calculate w-x covariance for good
% data pairs in that window.
% 3) calculate the expected variance in x in the gap, given the covariance outside of the gap and variance of w within the gap.
% 4) scale w within the gap to create predicted values of x with the "correct" variance.
%
% INPUTS:
% x: vector containing NaNs (gaps).
% w: reference vector, presumed to co-vary with x.
% OUPUT, xfill, is x with the gaps filled by scaling w.
% Note xfill may still contain NaN if there were any data pairs in x and w that were both NaN.
%
% 20170608 GMW

c = chunker(find(isnan(x))); %indices for gap chunks
n = diff(c,1,2)+1; % # of points outside of gap to consider. For now assume same as gap width
xfill = x;
for g = 1:size(c,1)
    
    %define indices
    gap = c(g,1):c(g,2);
    reg = gap(1)-n(g) : gap(end)+n(g); % index for region to consider
    reg(reg<=0 | reg>length(x)) = []; %out of bounds
    reg(ismember(reg,gap))      = []; %exclude gap
    reg(isnan(x(reg) + w(reg))) = []; %exclude points with nan w or x
    
    %interpolate if gap is small or if too many adjacent gaps
    %the limits chosen here are somewhat arbitrary
    minGapSize = 3;
    maxMissData = n(g); %must have at least 50% good data
    if n(g)<minGapSize || length(reg)<maxMissData
        t = [gap(1) gap(end)] + [-1 1]; %index for nearest non-nans
        if     t(1)<=0,        xfill(gap) = x(t(2));
        elseif t(2)>length(x), xfill(gap) = x(t(1));
        else                   xfill(gap) = interp1(t,x(t),gap);
        end
        continue
    end
    
    % covariance for region around gap
    xr = x(reg); xr = xr-nanmean(xr);
    wr = w(reg); wr = wr-nanmean(wr);
    covr = nanmean(xr.*wr);
    
    % expected std of x, using cov = r*wstd*xstd
    r = 1; %assume perfect correlation b/c we are scaling w below
    wstd = nanstd(w(gap));
    xstd = covr./r./wstd;
    
    % scale wind data to give a predicted value of x
    xfill(gap) = nanmean(x(reg)) + (w(gap)-nanmean(w(gap)))./wstd.*xstd;
    
end


