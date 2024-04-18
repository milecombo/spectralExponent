# Spectral Exponent
### a.k.a. spectral slope, power-law slope, PSD slope, slope of the aperiodic component, 1/f exponent


this code allows to compute the spectral exponent of the resting EEG, based on the Power Spectral Density (PSD), over a given scaling region.

The spectral exponent describes the decay of the PSD. it is computed as the slope of an OLS line, fit on log-freq vs log-PSD, excluding oscillatory peaks (and their base).

![figure 1 fit steps](https://user-images.githubusercontent.com/6671316/49792598-72717e00-fd33-11e8-993c-0430313dee43.png)

(Note: while the Y axis label reads as mV^2/Hz, it should really read µ^2/Hz. Thus, the proper version should be microVolt (and not milliVolt) squared over Hz.

## USAGE EXAMPLE:
````matlab
% you should have in your workspace
% sRate: sampling Rate, e.g. 1450 for Nexstim
% myEEGch: a vector of datapoints.  

% here, we generate 5 minutes of dummy data, with a 1/f^2 decay, alpha and beta oscillations
sRate= 1450;
dataPoints= 1: sRate*5*60 +2*sRate;
myEEGch=  randn(1,dataPoints(end)); 
myEEGch= smooth( myEEGch, sRate*2)'; 
for ff= [ 9:0.05:11  18:0.05:24] 
    myEEGch= myEEGch+ sin( 2*pi*ff/sRate.*dataPoints)/20/(ff^2); 
end 
myEEGch= myEEGch(sRate+1:end-sRate);


% first compute the PSD
 epLen= 3* sRate; epShift= 1*sRate;numFFT=[];
 [myPSD,frex]= pwelch( myEEGch  , epLen, epShift,numFFT, sRate); 
 
 frBand=[1 40];
 frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
 XX= frex(frBins);
 YY= myPSD(frBins);
 robRegMeth= 'ols'; % method to perform linear regression. see >> help robustfit
 doPlot= 1; figure;
 thisCol= [0 0 1];
 [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
 spectralExponent= intSlo(2);
 
 % repeat for every electrode to obtain the average spectral exponent across the scalp

%% fitPowerLaw3steps %%

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
% intSlo(2) slope of 2nd (final) powerLaw Fit. i.e. the spectral exponent
% stat is a structure with informations on the 2nd fit, (see 2nd argument of robustfit)

% Pows is a structure with
%           pred: 1*length(frBins) doubles, predicted PSD by the power-law fit
%            obs: 1*length(frBins) doubles, observed PSD    = YY
%            res: 1*length(frBins) doubles, residuals (pred-obs) PSD
%           frex: 1*length(frBins) doubles, observed frex   = XX
%       meanPred: 1 double, mean of the predicted PSD
%        meanObs: 1 double, mean of the observed PSD
%        meanRes: 1 double, mean of the residuals of the PSD
%     meanResPos: 1 double, mean of the residuals of the PSD (exclude negative values)
%       cenFrObs: 1 double, central Frequency of the observed PSD
%       cenFrRes: 1 double, central Frequency of the residual PSD

% Deviants is a structure with 
%       res2: 1*length(frBins) doubles, 
%        res: 1*length(frBins) doubles, 
%        rej:  length(frBins)*1 logicals, 1 if bins are rejected, 0 if kept   for fitting the 2nd power-law
%         fr: min and max freqs considered
%       thre: threshold used for bins adjacent to the peaks
%     clusts: structure with information on the contiguous freq. bins 

% stat0 is a structure with informations on the 1st fit, (see 2nd argument of robustfit)

% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit
````

## HOW DOES IT WORK?
#### • log10 transform of both PSD and Frequency, and resample evenly in the log-log graph
First transform PSD (YY) and frequencies (XX) each by log10, 
resample the values using logarithmically spaced bins (upsample 4x to obtain PSD bins evenly sampled in the loglog space)

#### • fit a first-guess/naive line 
A 1st order ordinary least square, or "best fitting straight line", is fit under the log-log axes

#### • Discard frequency bins with large peaks and adjacent bins above the first line
Find all the residual >0 (i.e. bins emerging from the naive power-law fit). Only large peaks (whose residuals are larger than 1*mad(residuals)) count as proper peaks.
Find clusters of bins with residuals >0, and consider only those clusters of frequency bins where there is a peak, for subsequent rejection.
i.e. reject large enough peaks (large residuals well above the naive power-law trend, where residuals > mad(residuals) ) 
along with the base of the peaks (adjacent bins till above the naive power-law trend, where residuals >0)

#### • fit a second and final line, on bins sufficiently resembling a power-law 
Perform the fit only on the remaining residuals (those closer to a power-law)

#### • the slope of this second line is the spectral exponent
indexing the decay of the PSD over frequencies.
statistics relative to the fit are also assessed


## SCIENTIFIC REFERENCES
Scientific References where we employed the present code to estimate the spectral exponent. 
#### • General Anesthesia
the first article of the series, the EEG spectral slope reflects the presence/absence of consciousness (according to delayed reports)
across anesthetic agents all tailored to reach complete behavioural unresponsiveness,
discriminating cases of unconsciousness (during propofol and xenon anesthesia) 
against baseline wakefulness and cases of sensorimotor disconnection from the environment (during ketamine)
```` 
Colombo, M. A., Napolitani, M., Boly, M., Gosseries, O., Casarotto, S., Rosanova, M., ... & Sarasso, S. (2019).
The spectral exponent of the resting EEG indexes the presence of consciousness
during unresponsiveness induced by propofol, xenon, and ketamine.
NeuroImage, 189, 631-644.
````
[link to Neuroimage article](https://doi.org/10.1016/j.neuroimage.2019.01.024).



#### • Severe Brain Injury
The EEG spectral exponent (dubbed spectral gradient) reflects the presence/absence of consciousness
in severe brain-injury with Disorders of Cosnciousness, in conscious patients with brain injury,
and even in cases of sensorimotor disconnection from the environment: 
those patients with unresponsive wakefulness syndrome (vegetative state) that nonetheless show capacity for consciousness,
as shown by complex reactivity to direct cortical perturbation,
(i.e. high values of Perturbational Complexity Index (PCI), assessed from TMS-EEG)

```` 
Colombo, M. A., Comanducci, A., Casarotto, S., Derchi, C. C., Annen, J., Viganò, A., ... & Rosanova, M. (2023).
Beyond alpha power: EEG spatial and spectral gradients robustly stratify disorders of consciousness.
Cerebral cortex, 33(11), 7193-7210.
````
[link to Cerebral cortex article](https://doi.org/10.1093/cercor/bhad031).


#### • Human development
Wake and light-sleep EEG aperiodic activity diverge from deep sleep activity over age from infants to adolescents.
EEG Spectral exponent flattens from infants to adolescent during wakefulness and light sleep, 
while the spread between wake and sleep values become progressively larger;
further, the hot-spot where the PSD is the slowest migrates following a postero-anterior gradient over age.

````
Favaro, J., Colombo, M. A., Mikulan, E., Sartori, S., Nosadini, M., Pelizza, M. F., ... & Toldo, I. (2023). 
The maturation of aperiodic EEG activity across development reveals a progressive differentiation of wakefulness from sleep. 
NeuroImage, 277, 120264.
````
[link to Neuroimage article](10.1016/j.neuroimage.2023.120264)


#### • Cortical Stroke
EEG Spectral exponent in patients with stroke (in a regimen of values compatible with consciousness) 
discriminate the lesioned from the contralateral hemisphere,
and correlate with clinical indexes of stroke recovery

````
Lanzone, J., Colombo, M. A., Sarasso, S., Zappasodi, F., Rosanova, M., Massimini, M., ... & Assenza, G. (2022).
EEG spectral exponent as a synthetic index for the longitudinal assessment of stroke recovery. 
Clinical Neurophysiology, 137, 92-101.
````
[link to Clinical Neurophysiology article](https://doi.org/10.1016/j.clinph.2022.02.022).

## OTHER RESOURCES
A good place explaining the basics; a brief description of the distinction between aperiodic and periodic activity

[link to blog](https://sapienlabs.org/lab-talk/eeg-and-depth-of-anesthesia/)

A youtube video where I present the findings of the EEG spectral exponent during various anesthetic protocols (Colombo et al., 2019)

[link to video](https://www.youtube.com/watch?v=d7csRxa1nCI)

