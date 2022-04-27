function X = undamp(Xd,t,alpha,tau,tol)
% function Xd = damp(X,t,alpha,tau,tol)
% Applies a correction to scalar measurements to "undampen" a slow signal response time (i.e. due to
% stickiness on inlets).
% Calculations are done iteratively:
% 1) take a guess at the undampened variable, X
% 2) call the "damp" function to dampen it
% 3) compare to input Xd and adjust guess
% 4) repeat until convergence.
%
% INPUTS:
% Xd: raw dampened scalar (e.g. mixing ratio).
% t: time vector, seconds. Code assumes this is evenly spaced.
% alpha: pre-exponential factor for max-normalized pulse-response.
% tau: decay constant for max-normalized pulse-response.
% tol: OPTIONAL tolerance for convergence, representing fractional difference between guessed and
%            input Xd. Default is 0.01.
%
% OUTPUT, X, is undamped scalar values.
%
% 20140807 GMW

if nargin<5
    tol = 0.01; %defualt convergence criteria
end

L = length(t);
X = Xd; %first guess
Xdg = Xd/2; %dummy
while any(abs(1 - Xdg./Xd)>tol)
    Xdg = damp(X,t,alpha,tau); %damped version of guess X
    X = X + (Xd - Xdg); %correct to new guess
end


