
Flux Tool Box
v1.0
Last Updated 20130921
Glenn Wolfe
zlonewolfe@gmail.com


A Few Words of Caution
    Eddy covariance (EC) is a statistical method for assessing surface-atmosphere interactions in the turbulent boundary layer.
    At first glance it seems very straight-forward: acquire high-frequency observations of vertical wind and concentrations, 
    apply a few corrections and voila! Users should keep in mind, however, that EC is subject to a number of assumptions regarding site location, 
    measurement setup and instrument performance. Results should thus be interpreted with due caution. This tool box is meant to provide 
    a starting point for those interested in developing their own EC measurements. Serious users should expend some effort on understanding 
    the physical underpinnings of EC. In particular, the "Handbook of Micrometeorology" by Lee, Massman and Law is a great introduction.

    As always, feedback is welcome and can be sent to Glenn at the above e-mail address.


Summary
    This toolbox is a set of functions and scripts designed to aid in calculation of turbulent vertical fluxes using the eddy 
    covariance (EC) method. All scripts are written for use in MATLAB, with algorithms borrowed or derived from various sources. 
    It is expected that users are familiar with the theory behind EC. If not, we suggest you start with Wikipedia and go from there.

    Before applying these tools, data must be appropriately reduced. Often, 3-D wind data and concentrations are collected on 
    different systems and/or with different timebases. To run EC calculations, data must be averaged or binned to the slowest 
    timebase. It is up to the user to decide the most appropriate method of data preparation.

    The primary function, ECflux.m, runs all of the necessary calculations and has a number of user-specified options. 
    This function is designed to work with a single "chunk" of data, i.e. it will calculate one EC flux for a given time series. 
    Specific operations include:
    1) detrend scalar variable x
    2) rotate 3-D winds into natural wind coordinate (such that mean cross and vertical wind speeds are 0) 
    3) lag scalar and vertical wind vectors to optimize covariance
    4) calculate covariance and exchange velocity
    5) plot quadrant data
    6) plot frequency spectra
    More info can be found in the documentation for this function.

    An example (exampleFluxCalculation.m) is included to demonstrate how these functions are utilized.
    Please read through this example and documentation on other functions before requesting assistance.
    Happy data crunching!


Included Functions and Scripts
    NOTE: some of these functions utilize routines (e.g. nanmean) from the MATLAB statistics toolbox, 
          which is not included with the base MATLAB package. Please contact Glenn if you require 
          a work-around for this issue.

    ECFlux           Primary function. Performs all calculations and makes plots.
    detrend          Calculates and removes a trend from a data vector.
    natWindRot       Applies a 2 or 3-step rotation to 3-D wind data, such that mean cross and vertical winds are 0. 
    lagCov           Calculates covariance between two variables at multiple time lags.
    lagVar           Applies a set time lag to a vector.
    quadPlot         Visualize statistics detailing contributions of + and - fluctuations to total flux.
    allSpectra       Calculates and plots spectra and cospectra for a piece of EC data.
    spectra          Calculates power spectrum (fourier transform) for a vector.
    ogive            Calculates a normalized cumulative integral of a power or cospectrum.
    cospectra        Calculates a co-spectrum (product) for two co-varying variables (typically x and w).
    struct2var       Takes all fields in a structure and makes them variables in the workspace.
    smooth           Applies a running mean/median filter to a piece of data.
    logbin           Bins and averages data in base-10 log space.
    f53line          Calculates and plots a line with a -5/3 slope in log-log space.


Improvements to Come
    Quality analysis
    Uncertainty estimation
    Spectral corrections
    Density (WPL) corrections

