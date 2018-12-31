function h = compareSpectra(varargin)
% function h = compareSpectra(varargin)
% Allows comparison of spectra and cospectra from multiple flux calculations.
% INPUTS are strucutures containing from the output of the ECFlux function.
% One of the inputs can also be an index vector with values of 1-8 specifying which spectra to plot:
%   1: x power, binned and power-normalized
%   2: x power, f-weighted, binned and power-normalized
%   3: x ogive
%   4: cospectrum, binned and power-normalized
%   5: cospectrum, f-weighted, binned and power-normalized
%   6: cospectral ogive
%   7: coherence
%   8: phase angle
%   9: wavelet cospectrum
%   10: wavelet ogive
% OUTPUT h give handles for axes.
%
% 20140217 GMW
% 20140407 GMW  Added option to specify subset of plots to make using an index.
% 20160921 GMW  Added axes handle output.
% 20161207 GMW  Added wavelet stuff

%% sub-select plots if desired
varCheck = cellfun(@isstruct,varargin);
if all(varCheck)
    whichPlots = 1:10;
else
    whichPlots = varargin{~varCheck};
    varargin = varargin(varCheck);
end

xpw     = any(whichPlots==1);
xpwf    = any(whichPlots==2);
xog     = any(whichPlots==3);
co      = any(whichPlots==4);
cof     = any(whichPlots==5);
coog    = any(whichPlots==6);
coh     = any(whichPlots==7);
phs     = any(whichPlots==8);
wvco    = any(whichPlots==9);
wvog    = any(whichPlots==10);
    
nF = length(varargin);
L = cell(nF,1); %legend strings

h = []; %axes handles

%% x Power Spectra
if xpw
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        loglog(F.spectra.bin.f,F.spectra.bin.psdxn)
    end
    xlabel('Frequency (Hz)')
    ylabel('psd(x)')
    legend(L)
    set(gca,'XScale','log','YScale','log')
    box on
    f53line([0.2 2],10,1);
    h(end+1) = gca;
end

%% x Power Spectra, f-weighted
if xpwf
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.spectra.bin.f,F.spectra.bin.psdxn.*F.spectra.bin.f)
    end
    xlabel('Frequency (Hz)')
    ylabel('psd(x)*f')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

%% x Ogives
if xog
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.spectra.f,F.spectra.ogx)
    end
    xlabel('Frequency (Hz)')
    ylabel('x Ogive')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

%% w-x CoSpectra
if co
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        loglog(F.spectra.bin.f,abs(F.spectra.bin.con))
    end
    
    %deal with negative points
    l = get(gca,'children'); %line handles
    l = flipud(l);
    for i=1:nF
        F = varargin{i};
        neg = F.spectra.bin.con<=0;
        c = get(l(i),'Color'); %color of last line
        loglog(F.spectra.bin.f(neg),-F.spectra.bin.con(neg),'x','Color',c) %mark negative values
    end
    
    xlabel('Frequency (Hz)')
    ylabel('Co(w-x)')
    legend(L)
    set(gca,'XScale','log','YScale','log')
    box on
    f53line([0.2 2],10,1);
    h(end+1) = gca;
end

%% w-x CoSpectra, f-weighted
if cof
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.spectra.bin.f,F.spectra.bin.con.*F.spectra.bin.f)
    end
    xlabel('Frequency (Hz)')
    ylabel('Co(w-x)*f')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

%% w-x Ogives
if coog
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.spectra.f,F.spectra.ogwx)
    end
    xlabel('Frequency (Hz)')
    ylabel('w-x Ogive')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

%% Coherence
if coh
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.spectra.bin.f,F.spectra.bin.cohr)
    end
    xlabel('Frequency (Hz)')
    ylabel('w-x Coherence')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

%% Phase angle
if phs
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.spectra.bin.f,F.spectra.bin.phase)
    end
    xlabel('Frequency (Hz)')
    ylabel('w-x Phase Angle (degrees)')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

%% Wavelet cospectra
if wvco
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.wave.freq,F.wave.co)
    end
    xlabel('Frequency (Hz)')
    ylabel('Wave Co(w-x)')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

%% wavelet Ogives
if wvog
    figure
    hold all
    for i=1:nF
        F = varargin{i};
        L{i} = inputname(i); %for legend
        semilogx(F.wave.freq,F.wave.og)
    end
    xlabel('Frequency (Hz)')
    ylabel('Wave w-x Ogive')
    legend(L)
    set(gca,'XScale','log')
    box on
    h(end+1) = gca;
end

