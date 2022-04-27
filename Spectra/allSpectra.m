function S = allSpectra(w,x,F,nbin,plotSpectra)
% function S = allSpectra(w,x,F,nbin,plotSpectra)
% Calls subroutines to calculate spectra and cospectra (i.e. fourier transforms) of data used in
% eddy covariance calculations.
% This code borrows heavily from the script spec_All.m, provided by Delphine Farmer and Ian Faloona.
% IMPORTANT NOTE: input data w and x must be constantly spaced in time!
%                 NaNs will be filled with 0s.
%
% INPUTS:
% w: vertical wind speed
% x: scalar vector (e.g. concentration, mixing ratio, etc).
% F: frequency of observations, Hz.
% nbin: number of bins for log-averaging spectra. Default is 100.
% plotSpectra: flag for plotting spectra. 0=no, 1=yes.
%
% OUTPUTS:
% S: a structure containing all of the computed spectra.
%       f:      frequency, Hz
%       psdx:   x power spectrum
%       psdxn:  variance-normalized x power spectrum
%       ogx:    x ogive 
%       psdw:   w power spectrum
%       psdwn:  variance-normalized w power spectrum
%       ogw:    w ogive
%       co:     w-x cospectrum
%       con:    covariance-normalized w-x cospectrum
%       ogwx:   w-x ogive
%       qu:     quadrature spectrum (imaginary cospectrum)
%       qun:    qu normalized by area
%       bin:    a sub-structure containing log-binned versions of spectra and cospectra.
%               also contains the following:
%               cohr:   coherence for cospectrum
%               phase:  phase angle for cospectrum
%
% 20130921 GMW
% 20131007 GMW  Added bin outputs
% 20140310 GMW  Added gap-filling
% 20140319 GMW  Added quadrature, coherence and phase angle cospectra.

%% INITIALIZE VARIABLES

%defaults
if nargin<4, nbin        = 100; end
if nargin<5, plotSpectra = 0;   end

%subtract means
x = x - mean(x,'omitnan');
w = w - mean(w,'omitnan');

%zero fill
n = isnan(w+x);
w(n) = 0;
x(n) = 0;

%% X POWER SPECTRA
[f psdx psdxn]  = spectra(F,x);            % power spectrum
ogx             = ogive(f,psdx);           % integrated spectrum
[fb, psdxb]     = logbin(f, psdx, nbin);   % bin-average
[~, psdxnb]     = logbin(f, psdxn, nbin);

% plot if desired
if plotSpectra
    [m,i] = max(psdxn);
    f53 = f.^(-5/3)./f(i).^(-5/3).*m;

    figure
    loglog(f,psdxn,'c-')
    hold on
    loglog(fb,psdxnb,'b-','lineWidth',3)
    loglog(f,f53,'k--')
    xlabel('Frequency (Hz)')
    ylabel('x Spectral Power')
    legend('Data','Avg','-5/3')
    
    figure
    semilogx(f,psdxn.*f,'c-')
    hold on
    semilogx(fb,psdxnb.*fb,'b-','lineWidth',3)
    xlabel('Frequency (Hz)')
    ylabel('SP_x*F')
    
    figure
    semilogx(f,ogx,'b-')
    xlabel('Frequency (Hz)')
    ylabel('x Ogive')
end

%% power spectra w
[~, psdw, psdwn] = spectra(F,w);
ogw              = ogive(f,psdw);
[~, psdwb]       = logbin(f, psdw, nbin);
[~, psdwnb]      = logbin(f, psdwn, nbin);

if plotSpectra
    [m,i] = max(psdwn);
    f53 = f.^(-5/3)./f(i).^(-5/3).*m;

    figure
    loglog(f,psdwn,'c-')
    hold on
    loglog(fb,psdwnb,'b-','lineWidth',3)
    loglog(f,f53,'k--')
    xlabel('Frequency (Hz)')
    ylabel('w Spectral Power')
    legend('Data','Avg','-5/3')
    
    figure
    semilogx(f,psdwn.*f,'c-')
    hold on
    semilogx(fb,psdwnb.*fb,'b-','lineWidth',3)
    xlabel('Frequency (Hz)')
    ylabel('SP_w*F')
    
    figure
    semilogx(f,ogw,'b-')
    xlabel('Frequency (Hz)')
    ylabel('w Ogive')
end

%% cospectra
[~,co,con,qu,qun] = cospectra(F,x,w);
ogwx = ogive(f,co);

[~, cob] = logbin(f, co, nbin);
[~, conb] = logbin(f, con, nbin);
[~, qub] = logbin(f, qu, nbin);
[~, qunb] = logbin(f, qun, nbin);

%Plot if desired
if plotSpectra
    [m,i] = max(abs(conb));
    f53 = abs(f.^(-5/3)./f(i).^(-5/3).*m);
    
    %separate positive and negative
    pos = con>0;
    neg = con<=0;
    
    posb = conb>0;
    negb = conb<=0;

    figure
    h1 = loglog(f(pos),con(pos),'c-');
    hold on
    h2 = loglog(f(neg),-con(neg),'cx');
    h3 = loglog(fb(posb),conb(posb),'bo','lineWidth',2);
    h4 = loglog(fb(negb),-conb(negb),'rd','lineWidth',2);
    h5 = loglog(f,f53,'k--');
    xlabel('Frequency (Hz)')
    ylabel('w-x Spectral Power')
    legend([h1 h3 h4 h5],'Data','+ Avg','- Avg','-5/3')
    
    figure
    semilogx(f,con.*f,'c-')
    hold on
    semilogx(fb,conb.*fb,'b-','lineWidth',3)
    xlabel('Frequency (Hz)')
    ylabel('SP_w_-_x*F')
    
    figure
    semilogx(f,ogwx,'b-')
    xlabel('Frequency (Hz)')
    ylabel('w-x Ogive')
end

%% Accumulate structure
S.f = f;

S.psdx = psdx;
S.psdxn = psdxn;
S.ogx = ogx;

S.psdw = psdw;
S.psdwn = psdwn;
S.ogw = ogw;

S.co = co;
S.con = con;
S.ogwx = ogwx;

S.qu = qu;
S.qun = qun;

S.bin.f = fb;
S.bin.psdx = psdxb;
S.bin.psdxn = psdxnb;
S.bin.psdw = psdwb;
S.bin.psdwn = psdwnb;
S.bin.co = cob;
S.bin.con = conb;
S.bin.qu = qub;
S.bin.qun = qunb;

[S.bin.cohr,S.bin.phase] = coherenceNphase(psdxb,psdwb,cob,qub);

