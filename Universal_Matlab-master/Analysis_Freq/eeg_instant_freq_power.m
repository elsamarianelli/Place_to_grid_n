function [ freqEEG50hz, ampEEG50hz, phaseEEG, ampEEG, freqEEG, filtEEG] = ...
    eeg_instant_freq_power( eeg, filterRange, varargin)
%BASIC_EEG_ANALYSIS Filters EEG & extracts power & instant frequency
% Takes EEG sequence and uses analytic function to calculate instantaneous frequency and
% power within a certain frequency range (i.e. filtered)

% --- ARGUMENTS --- two or three
% eeg           eeg sequence in full. Should already be converted to volts 
%               for meaningful results
%
% filterRange   range in which to filter eeg - lower then upper value in Hz
%               [1x2] e.g. [6,12]
%
% eegFreq       Specify frequency of eeg signal in Hz- if not supplied 
%               defaults to 250Hz (i.e. freq of standard LFP) hence if using
%               EGF file need to specify. 
%
%
% --- RETURNS ---
% freqEEG50hz   Frequency of eeg at each position point (NB. downsampled
%               from 250Hz to 50Hz to match .pos
%
% ampEEG50hz    Amplitude of eeg at each position point
%
% phaseEEG      Phase of signal at each eeg sample point (normally 250hz) Note this
%               monotonically increasing. Phase is not reset to 0 at 2pi but keeps
%               increasing. Peak of EEG is defined as 0 phase.
%
% ampEEG        Amplitude of eeg at each eeg sample point 
%
% freqEEG       Instantaneous frequency of eeg at each sample point
%
% filtEEG       Eeg timeseries filtered  according to input variables


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impose an filter on low amplitude values - set the phase of these values
%to nan as they are unreliable and worse than that often get assigned to a
%phase around 4.2 rads (odd). So values below the value ampFilt set to nan
ampFilt         =1E-7; 



%%% PARSE VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==2
    
%If samp freq not specified assuming standard EEG i.e. 250Hz
eegSampFreq=250; %EEG freq in Hz

elseif nargin==3
    eegSampFreq=varargin{1};
else
    error('eeg_instant_freq_power takes 2 or 3 arguments.');
end


%%% FILTER AND PROCESS EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%EEG as returned from read_eeg_file is a two column matrix, first column
%being the voltage data and the second being timestamps. This second column
%is removed by eeg_convert_to_volts.m but double check here that it is not
%present. First check eeg is column form, then just take first column.
if size(eeg,2)>size(eeg,1)
    eeg         =eeg';
end
eeg         =eeg(:,1);

%Detrend - remove best linear fit to take out DC shift & effectivly mean
%normalise
eeg             =detrend(eeg);

%Need to filter the EEG (bandpass on theta range) and find analytic function -
%this allows us to calculate power, freq, amp etc.
window = blackman(eegSampFreq+1); %Define window of interest - 1 second winowd here
nqLim=eegSampFreq/2;
EEGFilter = fir1(eegSampFreq,[filterRange(1)/nqLim filterRange(2)/nqLim],window); % Denominator - 125 is half the sampling rate - Nyquist frequency.
filtEEG = filtfilt(EEGFilter, 1, eeg); % The numerators are then the real actual frequencies you want to band pass.
analyticEEG = hilbert(filtEEG); %Hilbert transform itself - returns a complex var.
clear eeg

%Process analytic function to get phase, frequency and power
%(actually using true amp instead of power - power would just be that
% squared)
phaseEEG=angle(analyticEEG); %Phase (radians) is the angle of the complex variable at each time point
phaseEEG=unwrap(phaseEEG); %Adds 2pi if necessary to render plot of phase smooth
ampEEG=abs(analyticEEG); %Modulus of analytic function i.e. instantaneous amplitude

%Impose a filter based on very low amp points (actual value defined
%in ampFilt at top of this function). These very low amp points have
%spurious phase and freq - worse than that don't tent to take a random
%phase but typically take a value around 4.2rads. Can be seen as a large
%peak in hist(phaseEEG)
percLowAmp              =sum(ampEEG<ampFilt)/length(ampEEG)*100;
if percLowAmp >10 %Warning if many pts low
    warning(['>10% (' num2str(percLowAmp,2) '%) of EEG  points have v.low amplitude - set to nan.']);
end
clear percLowAmp
phaseEEG(ampEEG<ampFilt)=nan;
ampEEG(ampEEG<ampFilt)  =nan;


%Calculate instantaneous freq - note this vector is length of EEG-1 use
%interp2 to add back the extra point (interp2 needs to work on a mat not a
%vector). NB previously just added a nan to the end which resulted in a
%small but avoidable slippage.
freqEEG             =diff(phaseEEG)*(eegSampFreq/(2*pi));
freqEEG             =...
    interp2([freqEEG, ones(size(freqEEG))], ...
    ones(size(freqEEG,1)+1,1), ...
    [0.5:1:size(freqEEG,1)+0.5]');


%Sometimes Hilbert assigns regressing phase which end up as negative
%frequency, these are very likly to be spurious. Set negative freqs to nan,
%most will be lost in the final nanmean below which reduces sample rate to
%that of the pos file (50Hz). Also set the matching phase values to nan
%since these are likely to be spurious too [might consider setting values
%on either side to nan too but that could be too conservative]
phaseEEG(freqEEG<0)     =nan;
ampEEG(freqEEG<0)       =nan;
freqEEG(freqEEG<0)      =nan;



%To compare against position (smapled at 50hz) need to down sample - so
%take mean accross 5 bins at a time for 250hz eeg
downSampConst=eegSampFreq/50;
freqEEG50hz=nanmean(reshape(freqEEG, downSampConst, []))';
ampEEG50hz=nanmean(reshape(ampEEG, downSampConst, []))';


