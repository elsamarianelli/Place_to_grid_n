function [ power, freqs, maxPower, maxFreq, s2n ] = ...
    eeg_power_spectra( eeg, pos2use, smthKern, thetaRange, upperFrq, s2nWdth, varargin )
%EEG_POWER_SPECTRA2 FFT of LFP for specified positional points.
%
% --- ARGS ---
% eeg       [nEEG pts, 1] eeg sequence in full. Should already be converted 
%           to volts for meaningful results
%
% pos2use   [nPts2use,1] or [empty] index of pos samples to analyse. If 
%           an empty matrics is passed in all points will be used to make
%           the power spectra.
%           (note pos is sampled at 50Hz eeg at 250Hz or higher if EGF)
%
% smthKern  [1x2] defines gaus smoothing kern for power spec. 1st value is 
%           width in Hz [2] and 2nd value is sigma in Hz [0.5]
% 
% thetaRange[1x2] range to look for peak in 1st val is lower bound, 2nd is 
%           upper [6,12]

% upperFrq  [1] top limit of freq range returned. typically 25 or 125
%
% s2nWdth   [1] when calculating signal to noise for peak vs everything 
%           else this specifies how wide the band is - take half of this on 
%           each side [1]
%
% eegFreq   [not required] defaults to 250 if not specified
%
%
% --- RETURNS ----
% power - smoothed power spectrum
% freqs - frequencies that match power spectrum (truncated at 25Hz)
% maxPower - peak power from powerspectrum in theta range
% maxFreq - frequency at which peak power occured
% s2n - signal to noise - power in band around theta peak divided by power in rest of freqs
%
%
% --- EXAMPLE ---
%[power,freqs,maxPower,maxFreq,s2n] =eeg_power_spectra(eeg, [], [2,0.5], [6,12], 250, 1, Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HOUSE KEEPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 6 %eegFreq not specified default to 250
    eegFreq     =250; %Standard eeg
elseif nargin == 7 %eegFreq speicifed
    eegFreq     =varargin{1};
else 
    error('Specify either 6 or 7 arguments.');
end


lowerFrq         =2; %When calculating s2n don't count power below this value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



% Select out the eeg points to use - if pos2use is empty then we use all
% points - which is easy. If particular values are specified then need to
% select these but because pos and eeg are sampled at different rates need
% to convert between these.

if ~isempty(pos2use) %Select specific values - otherwise do nothing
    
    %N2L. Required because EEG is sampled at 5* the freq of position -
    %hence have to extrapolate from position index to eeg index
    eegVPosRate =eegFreq/50; %Dif in sampling rates
    N           =length(pos2use);
    offsets     =reshape(repmat([1-eegVPosRate:0],[N,1]),[],1);
    eeg2use     =sort( repmat(eegVPosRate.*pos2use(:), eegVPosRate, 1) + offsets );
    eeg         =eeg(eeg2use);
    clear N offsets eegVPosRate
end



%Do ffts
%fft - second value should be a power of 2 greater in length than the
%eeg series (power of two speeds up computation).

% fftLength       =2^18; % 2^16 Standardise length and resolutoin accross eegs and with intrinisc analysis
fftLength   =2^nextpow2(length(eeg));
fftRes      =fft(eeg, fftLength); %result is complex

freqs       =([0:(length(fftRes)-1)]*(eegFreq/length(fftRes)))';
freqs       =freqs(freqs<=upperFrq);

%NL. eqiv of taking the square modulus
power       =(fftRes.*conj(fftRes))/length(fftRes);

%For speed lose seciton of power above top freq
power       =power(1:length(freqs));

%Smooth to make more reliable - kernel size is specified by default vars
kernelLength =length(find(freqs<=smthKern(1)));
kern        =fspecial('gaussian',[1, kernelLength],...
    kernelLength*(smthKern(2)/smthKern(1) ));
power       =imfilter(power, kern', 'symmetric');

%Find max in theta range and return power plus freq
maxPower    =max(power(((freqs>thetaRange(1))&(freqs<thetaRange(2)))));
maxFreq     =mean(freqs(ismember(power, maxPower)));

%Finally get signal to noise for area around theta peak vs whole range excluding peak
%(this is what tom does and he uses +_ 0.5hz either side of peak)
peakFreqBand=freqs>maxFreq-(s2nWdth/2) & freqs < maxFreq+(s2nWdth/2); %Logical
minFreqBand         =freqs>=lowerFrq;
s2n         =nanmean(power(peakFreqBand))/ nanmean(power(~peakFreqBand & minFreqBand));
end

