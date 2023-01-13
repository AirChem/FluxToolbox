% exampleFluxCalculation.m
% An example of how to use the flux toolbox to calculate an eddy covariance flux.
% The data used in this example was taken from the MMS instrument onboard the DC-8 aircraft during SEAC4RS.
% The selected data is a 20Hz subset from a level flight leg (~400 m above ground level) over a forest.
% For this case, "data.x" is the MMS-measured air temperature.
% Flux and associated diagnostics are computed for both temperature (i.e. sensible heat) and
% horizontal wind (i.e. momentum).
% For more information, see the documentation for ECFlux.m
%
% 20130921 GMW

%% STEP 1: GET DATA
% Here we load a structure, "data", that contains all of the variables needed for the flux
% calculation: time (t), 3-D winds (u, v, w) and scalar values (x).
% The time vector must have constant spacing, and all vectors must be the same length.

load exampleData.mat

%% STEP 2: DEFINE OPTIONS
% The options stucture defines various parameters for how the calculation is done and what plots are
% generated. For more info, see Readme.txt or comments in ECFlux.m.

options.despikeFlag     = 0;        %flag for despiking x
options.xTrendType      = 'smooth'; %mean, linear or smooth
options.frameSize       = 400;     %frame size for smoothing, # of points
options.thirdRotation   = 0;        %flag for third wind rotation
options.plotX           = 1;        %flag for plotting x data
options.plotW           = 1;        %flag for plotting wind data
options.plotLag         = 1;        %flag for lag-covariance plot
options.plotQuad        = 1;        %flag for quadrant plot
options.plotSpectra     = 1;        %flag for frequency-spectrum plots
options.nLags           = 1000;     %number of lag points for lag-covariance
options.xLag            = [];       %lag to apply to x (if not determined automatically)
options.nStat           = 3;        %number of sub-sets for stationarity test
options.plotWave        = 1;        % plot wavelet fluxes

%% STEP 3: CALCULATE FLUXES
% Here, "data" and "options" are used as inputs to the ECFlux function, which performs all necessary
% calculations and generates plots.
% Fluxes are calculated for both sensible heat (i.e. temperature) and momentum (i.e. horizontal wind). 

TorMOM = 'T'; %pick one
switch TorMOM
    case 'T' % heat flux
        F_heat = ECFlux(data,options); % heat flux
        
    case 'MOM' % momentum flux
        data.x = 'u'; %overwrite for new calculation
        options.xTrendType = 'mean';
        F_mom = ECFlux(data,options); % momentum flux
end


