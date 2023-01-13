function S = compileStats(varargin)

nF = length(varargin);
n = nan(nF,6);
S.x = n; S.w = n; S.x_dtl = n; S.w_r = n; S.flux = n;
S.r_wx = n(:,1);

for i=1:nF
    F = varargin{i};
    
    S.x(i,:)     = F.stats.x;
    S.w(i,:)     = F.stats.w;
    S.x_dtl(i,:) = F.stats.x_dtl;
    S.w_r(i,:)   = F.stats.w_r;
    S.flux(i,:)  = F.stats.flux;
    S.r_wx(i)    = F.stats.r_wx;
end

S.info = F.stats.info;


