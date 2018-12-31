function compareLags(varargin)
% function compareLags(varargin)
% Creates a plot to compare lag-covariances for multiple flux calculations.
% INPUTS are strucutures containing the fields "lags" and "cov_wx."
% Typically, these are from the output of the ECFlux function.
%
% 20140217 GMW

nF = length(varargin);
L = cell(nF,1); %legend strings

figure
hold all
for i=1:nF
    F = varargin{i};
    L{i} = inputname(i); %for legend
    plot(F.lags,F.cov_wx)
end
xlabel('lag points')
ylabel('<w''x''>')
legend(L)
grid on
box on


