function [peaksDirection, peaksAvPower] = FindPeakDirectionsInFFT2(fftPower, fftDirection, fftFreq, MAX_FREQ, thresh_coeff)
% Code originally written by JK.
% 
% Takes an FFT2, computes the average power for each direction and then
% locates the directions with peak power.
%
% DM, May 2013. Cleaned up.
%

MIN_PEAK_DISTANCE_DEGREES = 10; %WARNING: this assumes each direction bin is 1 degree

% power_spect was a square/rectangle, but we wanted a circle, so
% mask out frequencies that are too high.
isGood = fftFreq <= MAX_FREQ;%
fftPower = fftPower(isGood);
fftDirection = fftDirection(isGood);

% Get (and smooth) the average power in each direction
[direction, avPower] = GetFFTPowerByDirection(fftPower,fftDirection);
maxPower = max(avPower);

if maxPower == 0
    peaksDirection = [];
    peaksAvPower = [];
    return;
end

% pad avPower so we can use find peaks with periodic boundaries
originalLen = numel(avPower);
avPower = [avPower avPower avPower]; 
direction = [direction direction direction];

% find the peaks
[pks,locs] = findpeaks(avPower, 'minpeakheight', maxPower*thresh_coeff, ...
                                'minpeakdistance', MIN_PEAK_DISTANCE_DEGREES);
                            
% account for padding
good = locs > originalLen & locs <= originalLen*2;
locs = locs(good);
pks = pks(good);

% findpeaks just returned the index of the peak, lets lookup the direction
locs = direction(locs); 

% restrict results to the interval [1 180]
good = locs >=1 & locs <= 180;
peaksDirection = locs(good);
peaksAvPower = pks(good);

end