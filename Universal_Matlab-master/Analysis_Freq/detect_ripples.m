function [ripInd, power, pkPower, meanFreq, duration] ...
    =detect_ripples(eeg, eegFreq, ripFreq, ripThresh)
% Return mask size of eeg indicating where SWRs are present

%Detects ripples using the Csisvari method which basically finds regions of
%the EEG with power in the ripple band above the mean plus std times some
%factor. Then exapands these out until power is back to the mean plus 0.5
%std. If these areas are large enough they are SWRs.
%
% NB. Must have eeg of sufficiently high freq to avoid problems with Niquist
% i.e. if ripple freq is 150:250 then you need upto 500Hz eeg
%
% NB2. Ripple duration and minimum power cutoff are hard coded variables -
% while exploring Stephs data I found that a upper limit threshold of 5SD
% worked quite well paried with lower threshold of 1SD. Not this is very
% stringent and might not be so good if these periods were being used to
% for replay decoding. I also found a minum duration of 80ms worked well.
% The resulting ripples looked plausible.
    
% -ARGS-
% eeg      Unfiltered eeg - doesn't need to be converted to volts as
%          detection method is relative to mean
%
% eegFreq  Sample frequency of eeg
%
% ripFreq  [Not manditory]. Range that we use to define ripples - would 
%          typically be [150,250] if not supplied default to [150,250]
%
% ripThesh [Not manditory]The detection threshold is defined as 
%           mean + std * ripThresh. If not supplied default is 4.
%
%
% -RETURNS-
% ripInd    [nSWRs x 2] For each ripple the index (into eeg) of the start
%           (first column) and end (second column)
%
% power     [lengthEEG,1] Power in defined ripple band for the entire eeg
%           after smoothing with 0.1s box car (used for SWR detection)
%
% pkPower   [nSWRs x 1] Peak power in ripple band seen in each detected
%           ripple
%
% meanFreq  [nSWRs x 1] For each ripple the mean frequency. Determined by
%           taking the mean freq of the filtered eeg over the extent of
%           each ripple
%
% duration  [nSWRs x 1] Length of each ripple in seconds.
%
%
% -EXAMPLE-
% [rips,power,pkPower,meanFreq,duration]   =detect_ripples(eeg, 1200, [150,250], 4)

% Alternativly the folliwing line with use default values for ripFreq and
% ripThresh (and will be equivalent to the line above.
% rips      =detect_ripples(eeg, 1200); 





% --- VARS ----------------------------------------------------------------

%Defult ripple freq range (i.e. if not defined as argument) and ripple
%threshold (i.e. how many times STD above mean)
dfRipFreq   =[150,250];
dfRipThresh =5; %Freyja uses 3 here - Csisvari uses 4

% Also set the minmum threshold for ripples - having detected an area of
% LFP with power above the threshold defined by ripThresh then work out to
% find a contiguous area above this threshold (i.e. mean + std*ripLowThresh
ripLowThresh =0.5; %Freyja uses 0 here - Csisvari uses 0.5

% Minmum duration also - ripples must be longer than this to be accepted.
% Freyja uses 20 or 30ms of LFP, 40 for MUA.
ripMinDur   =60/1000; %Minimum length in seconds (i.e. 40ms = 0.004)



% --- HOUSEKEEPING --------------------------------------------------------
%If all args haven't been set use defult
if nargin   <3
    ripFreq     =dfRipFreq;
end

if nargin   <4
    ripThresh   =dfRipThresh;
end





% --- MAIN ----------------------------------------------------------------
%Filter and get power
[ ~, ~, ~, ampEEG, freqEEG]=eeg_instant_freq_power( eeg, ripFreq, eegFreq);

power       =ampEEG.^2; % Power of filtered eeg
power(isnan(power))=0; %Set nans to 0

% Smooth with kern aprox half width of ripples i.e. 0.05s so eegFreq/20
kern        =ones(round(eegFreq/20),1);
kern        =kern./sum(kern);
power       =filter2(kern,power);

%Now detect high power epochs
meanPower   =nanmean(power);
stdPower    =nanstd(power);
ripsFull    =power>(meanPower + ripThresh*stdPower);
ripsLow     =power>(meanPower + ripLowThresh*stdPower);


%Now validate the putative ripples - find the start and end of each of
%the low treshold ripples
ripStr      =find(diff(ripsLow)==1);
ripEnd      =find(diff(ripsLow)==-1);

if ripStr(1)>ripEnd(1) %Deal with situations in which we start with an end
    ripEnd      =ripEnd(2:end);
end
ripStr      =ripStr(1:length(ripEnd));

duration    =(ripEnd-ripStr)./eegFreq; %Duration in seconds
%Eliminate the ones that are too short
ripStr      =ripStr(duration>ripMinDur);
ripEnd      =ripEnd(duration>ripMinDur);

%Finally check that each candidate ripple includes a period with activity
%above the higher threshold. ATM I have a loop here - might be possible to
%vectorise and improve if time. At the same time return some information
%about each event
[valid, pkPower, meanFreq]  =deal(zeros(length(ripStr)));
for n       =1:length(ripStr)
    valid(n)    =sum(ripsFull(ripStr(n):ripEnd(n))); %If sum>0 then  contains high power
    pkPower(n)  =max(power(ripStr(n):ripEnd(n)));
    meanFreq(n) =nanmean(freqEEG(ripStr(n):ripEnd(n)));
end

% Only return info about valid rips
valid       =logical(valid);
ripStr      =ripStr(valid);
ripEnd      =ripEnd(valid);
ripInd      =[ripStr,ripEnd];
pkPower     =pkPower(valid);
meanFreq    =meanFreq(valid);
duration    =(ripEnd-ripStr)./eegFreq;



    





