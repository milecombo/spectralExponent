# Spectral Exponent
### a.k.a. spectral slope, power-law slope, PSD slope, slope of the aperiodic component, 1/f exponent

## UPDATE >>> PYTHON TRANSLATION NOW AVAILABLE<<<

this code allows to compute the spectral exponent of the resting EEG, based on the Power Spectral Density (PSD), over a given scaling region.

The spectral exponent describes the decay of the PSD. it is computed as the slope of an OLS line, fit on log-freq vs log-PSD, excluding oscillatory peaks (and their base).

![figure 1 fit steps](https://user-images.githubusercontent.com/6671316/49792598-72717e00-fd33-11e8-993c-0430313dee43.png)

(Note: while the Y axis label reads as mV^2/Hz, it should really read µ^2/Hz. Thus, the proper version should be microVolt (and not milliVolt) squared over Hz.

## USE EXAMPLE (MATLAB)
 you should have in your workspace:
 
     sRate: sampling Rate, e.g. 1450 for Nexstim amplifier
    signal: one EEG channel, a vector of datapoints.  
     % use one eeg channel after preprocessing and filtering (but be aware of how filters impact the psd)
     % if high-pass filter (butter order 5) is set at 0.5 Hz, then begin fitting at 1 Hz
     % if low-pass filter (butter order 5) is set at 60 Hz, then stop fitting at approximately 45 Hz. 
     % Better even, do NOT use any low-pass

````matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%% FAKE SIGNAL GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate signal data (a sum of sinusoids with random phase) with a given PSD exponent (-1,5)
sRate= 1000;
freqs= linspace( .5, sRate/2, 10000);
time= (1 : 60*sRate )./sRate; %1: sRate*5*60 +2*sRate;
exp= -1.5;
signal= zeros(size(time));
for ii = 1: length(freqs)
     pwr = freqs(ii)^exp;
     signal= signal + sqrt(pwr) *sin(2*pi*freqs(ii)*time + rand(1)*2*pi ) ; % + randn(size(time)).*sqrt(pwr)
end

%%% add to the signal alpha and beta oscillations , leading to two bell-shaped peaks in the PSD
pks= [8 13; 16 24];
for fb= 1 : size(pks,1)
    ampFac= hamming(200);
    ffs= linspace(pks(fb,1) , pks(fb,2) ,200);
    for ii= 1 : 200;  
        amp = 60/( ffs(ii)^ 1.75)* ampFac(ii); % 60 1.75  % 90 2
        signal= signal+ sin( 2*pi* ffs(ii)*time + rand(1)*2*pi )*  amp;
    end
end
% figure, plot(signal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
````
now that you have generated a fake signal (or imported a real signal), you can compute the PSD
and estimate the spectral exponent in a given frequency range (1-40 Hz)

````matlab
% first compute the PSD
 epLen= 3* sRate; epShift= 1*sRate; numFFT=[];
 [myPSD,frex]= pwelch( signal  , epLen, epShift,numFFT, sRate); 
 
 frBand=[1 40];
 frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
 XX= frex(frBins);
 YY= myPSD(frBins);
 robRegMeth= 'ols'; % method to perform linear regression. see >> help robustfit
 doPlot= 1; figure;
 thisCol= [0 0 1];
 [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
 spectralExponent= intSlo(2);
 
 % repeat for every scalp electrode, and compute the average spectral exponent across the scalp

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

## USE EXAMPLE (PYTHON)
First, generate an EEG-looking signal
````python
# # HOW TO USE THIS CODE ON A GENERATED EEG SIGNAL

from scipy.signal import welch
def dsearchn(x, value):
    dist_from_value = np.abs(x-value)
    return np.where(dist_from_value == dist_from_value.min())[0][0]

################################### SIGNAL GENERATION ###################################

# generate noise data (a sum of sinusoids with random phase) with a given PSD exponent (-1.5)
fs = 1000 #sampling rate in Hz
time = np.arange(60, step=1/fs)
freqs = np.linspace(0.5,fs/2,10000)
exponent=-1.5

signal = np.zeros(time.shape)
for f in range(len(freqs)):
    power = freqs[f]**(exponent)
    signal += np.sqrt(power)*np.sin(2*np.pi*freqs[f]*time + np.random.uniform()*2*np.pi)

# add to the noise alpha and beta oscillations, leading to two bell-shaped peaks in the PSD
pks= np.array([[8, 13], [16, 24]])
for fb in range(len(pks)):
    ampFac= np.hamming(200)
    ffs= np.linspace(pks[fb,0] , pks[fb,1], 200)
    for ii in range(len(ffs)):
        amp = 60/(ffs[ii]**1.75)*ampFac[ii] #90 **2
        signal += amp * np.sin(2*np.pi*ffs[ii]*time + np.random.uniform()*2*np.pi)
        
########## Compute the PSD and estimate the exponent in the range 1 to 40 Hz ########## 
````
now that you have generated a fake signal (or imported a real signal), you can compute the PSD and estimate the spectral exponent in a given frequency range (1-40 Hz)

````python
#PSD
frex, psd_1chan = welch(signal, fs=fs, nperseg=fs*3, noverlap=fs*2, detrend='linear')
frband = [1, 40]
frBins = np.arange(dsearchn(frex, frband[0]), dsearchn(frex, frband[1])+1)
XX = frex[frBins]
# spectral exponent
YY = psd_1chan[frBins]
plt.figure(dpi=300, figsize=(4,4))
slope, inter, stats, vectors = compute_SpectralExponent(XX, YY, do_plot=True)

````



![PSD of fake signal](Figure_1_PSD_exponent fit0_and_fit1.png)
### PSD and estimates of the spectral exponent (preliminary vs excluding peaks) of a generated EEG-looking signal
This generated EEG signal has a theoretical PSD spectral exponent of -1.5. 
As you can see, the first preliminary/naive fit (in the legend> slope fit0, grey dotted line) of the whole PSD is biased by the peaks in the alpha and beta frequency bands.
The second and final fit (in the legend> slope 1-40 Hz, blue dotted line) , by excluding most of those peaks, well adheres to the aperiodic PSD background,
 and its slope is much closer to the theoretical value (not exact, each simulation gives slightly different values).


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

