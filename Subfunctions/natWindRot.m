function [u_r,v_r,w_r,angles] = natWindRot(u,v,w,t,nRot,plotMe)
% function [u_r,v_r,w_r,angles] = natWindRot(u,v,w,t,thirdRotation,plotMe)
% Rotates wind vectors such that vavg=wavg=0 for a chunk of data.
% Routine has been taken from "Handbook of Micrometeorology" by Lee, Massman and Law.
% Designed for use with ECFluxToolBox.
%
% INPUTS:
% u,v,w: x,y and z wind vectors in the instrument coordinate
% t: time vector for winds.
% nRot: # of rotations to do:
%       0: none
%       1: wavg = 0 only
%       2: wavg = 0 and vavg = 0 (default)
%       3: wavg = 0 and vavg = 0 and <w'v'> = 0
% plotMe: flag for generating time series plot. default is 0 (no).
%
% OUTPUTS:
% u_r,v_r,w_r: wind vectors rotated into the natural wind coordinate.
%              u_r: horizontal wind speed in mean wind direction.
%              v_r: cross-wind speed
%              w_r: vertical wind speed
% angles: structure containing rotation angles (in degrees):
%   eta: rotation angle around the z1 axis
%   theta: rotation angle around the y1 axis
%   beta: rotation angle around x2 axis to make cross-wind momentum flux (w2'v2')bar = 0
%
% 20130915 GMW
% 20140729 GMW  changed "thirdRotation" input to "nRot" and added option to only do first rotation.
% 20141118 GMW  added short-circuit if nRot=0.
% 20220316 GMW  removed nanmean calls and added de-naning at top of function.

%defaults
if nargin<6
    plotMe = 0;
end

if nargin<5
    nRot = 2;
end

if nRot==0
    u_r = u;
    v_r = v;
    w_r = w;
    angles.eta=0;
    angles.theta=0;
    angles.beta=0;
    return
end

% denan
i = isnan(u+v+w);
u(i)=[]; v(i)=[]; w(i)=[]; t(i)=[];

%Define some quantities
u1=u; v1=v; w1=w;
u1m=mean(u1); v1m=mean(v1); w1m=mean(w1);
u1m2=u1m.^2; v1m2=v1m.^2; w1m2=w1m.^2;
rms_uv=sqrt(u1m2+v1m2); rms_uvw=sqrt(u1m2+v1m2+w1m2);

%cosines and sines
CE = u1m./rms_uv;
SE = v1m./rms_uv;
CT = rms_uv./rms_uvw;
ST = w1m./rms_uvw;

%first two rotations
u2 = u1*CT*CE + v1*CT*SE + w1*ST;
v2 = v1*CE - u1*SE;
w2 = w1*CT - u1*ST*CE - v1*ST*SE;

%third rotation (optional)
v2p = v2-mean(v2); w2p = w2-mean(w2);
B = 0.5*atan(2*mean(v2p.*w2p)./(mean(v2p.^2)-mean(w2p.^2)));
CB = cos(B); SB = sin(B);
v3 = v2*CB + w2*SB;
w3 = w2*CB - v2*SB;

% output wind vectors
switch nRot
    case 1
        u_r = u1;
        v_r = v1;
        w_r = w2;
    case 2
        u_r = u2;
        v_r = v2;
        w_r = w2;
    case 3
        u_r = u2;
        v_r = v3;
        w_r = w3;
    otherwise
        warning('natWindRot: nRot input must be 1,2 or 3. Using default value of 2.')
        u_r = u2;
        v_r = v2;
        w_r = w2;
end

%rotation angles
angles.eta   = asin(SE)*180/pi;
angles.theta = asin(ST)*180/pi;
angles.beta  = B*180/pi;

%plot if desired
if plotMe
    o = ones(size(t));
    
    % U
    figure
    plot(t,u,'g-',...
        t,u_r,'c-',...
        t,mean(u)*o,'k--',...
        t,mean(u_r)*o,'b--')
    xlabel('Time')
    ylabel('U wind speed')
    legend('raw','rot','raw mean','rot mean')
    box on
    
    % V
    figure
    plot(t,v,'g-',...
        t,v_r,'c-',...
        t,mean(v)*o,'k--',...
        t,mean(v_r)*o,'b--')
    xlabel('Time')
    ylabel('V wind speed')
    legend('raw','rot','raw mean','rot mean')
    box on
    
    % W
    figure
    plot(t,w,'g-',...
        t,w_r,'c-',...
        t,mean(w)*o,'k--',...
        t,mean(w_r)*o,'b--')
    xlabel('Time')
    ylabel('W wind speed')
    legend('raw','rot','raw mean','rot mean')
    box on
    
end

