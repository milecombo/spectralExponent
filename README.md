# spectralExponent
this code allows to compute the spectral exponent of the resting EEG, based on the Power Spectral Density (PSD), over a given scaling region.

The spectral exponent describes the decay of the PSD. it is computed as the slope of an OLS line, fit on log-freq vs log-PSD, excluding oscillatory peaks (and their base).

% USAGE EXAMPLE: 
% % first compute the PSD
% epLen= 2* sRate; epShift= 1*srate;numFFT=[];
%  [myPSD,frex]= pwelch( myEEGch  , epLen, epShift,numFFT, sRate); 
%  frBand=[1 40];
%  frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
% XX= frex(frBins);
% YY= myPSD(frBins);
% robRegMeth= 'ols';
% doPlot= 1;
% thisCol= [0 0 1];
%  [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
%%%%%%%%% INPUT
% XX ---> are the frequency bins of interest (hz)
% YY ---> are the PSD relative to the bins of interest (mV^2 /Hz)
% robRegMeth --->is the method for robust regression fit (default: ols)
%                type help robustfit for more details
% doPlot  --->    0 (default) means no plot ; 
%                 1 is plot XX and YY on log-log scale;  
%                 2 is plot log(XX) and log(YY) on linear scale; 
% thisCol --->    if plotting, is the rgb triplet
%%%%%%%%% OUTPUT
% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit
