% white noise detection limit
% DELETED FROM ECFLUX.M, 20170106 GMW
nrep = 10; %number of times to repeat (b/c using random time series)
std_noise = nan(nrep,1);
for i=1:nrep
    x_noise = nanmean(x_dtl) + nanstd(x_dtl).*randn(size(x_dtl)); %white noise time series
    x_noise(isnan(x_dtl)) = nan;
    cov_noise = lagCov(w_r,x_noise,nLags,0);
    std_noise(i) = nanstd(cov_noise);
end
quality.fluxDL = S2Nlimit.*nanmean(std_noise);