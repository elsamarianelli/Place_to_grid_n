function ret = getThetaProps(varargin)
% Generates a structure, theta, with the following vectors, each of the same
% length as the input phase and amplitude vectors:
% oldPhase: same as input phase
% filteredEEG: same as inputAmplitude
% spikeCount: number of spikes in the given time bin
% cycleLabel: increasing integers, 0 for bad points
% amplitude: same as filteredEEG bud with NaN for bad points
% phase: similar to oldPhase, but with a shift and NaN for bad points
% instantaneousFiringRate: smoothed version of spikeCount
%
%%    Copyright (C) <2013>  <Daniel Manson> <d.manson@ucl.ac.uk> &
%%                          <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Boilerplate for standardised input/output behaviour
if nargin == 0, ret = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

ret.oldPhase = in.phase;
ret.oldAmplitude = in.filteredEEG;

% Get indices of spikes in to eeg
spkEegInd = ceil(in.spikeTS*in.esr);
spkEegInd(spkEegInd > numel(in.phase)) = numel(in.phase);

ret.spikeCount = accumarray(spkEegInd,1,size(in.phase));

% Adjust phase so that the cell's minimum spiking phase is zero
spikePhase = in.phase(spkEegInd);
minimumSpikingPhase = GetPhaseOfMinimumSpiking(spikePhase);
phaseAdjusted = fixAngle(in.phase - minimumSpikingPhase*(pi/180) + in.CANONICAL_MINIMUM_SPIKING_PHASE);

% Find points with negative frequency
isNegFreq = diff(unwrap(phaseAdjusted))<0;
isNegFreq(end+1) = isNegFreq(end);

% Define starts of theta cycles as the point where phase difference between two
% points is greater than pi.
phaseDiff = diff(phaseAdjusted);
isCycleStart = [true ; phaseDiff(2:end)<-pi ; true];
isCycleStart(isNegFreq) = false;
cycleLabel = cumsum(isCycleStart) + 1;

%Calculate power and find low power cycles 
power = in.filteredEEG.^2;
cycleTotalValidPower = accumarray(cycleLabel(~isNegFreq),power(~isNegFreq));
cycleValidBinCount = accumarray(cycleLabel(~isNegFreq),1);
cycleValidMeanPower = cycleTotalValidPower ./ cycleValidBinCount;
powerRejectionThreshold = prctile(cycleValidMeanPower,in.MEAN_POWER_PERCENTILE_THRESHOLD);
cycleHasBadPower = cycleValidMeanPower < powerRejectionThreshold;

% Find cycles that are too long or too short
cycleTotalBinCount = accumarray(cycleLabel,1); %i.e. unlike cycleValidBinCount this includes negFreq bins
cycleHasBadLength = cycleTotalBinCount>in.ALLOWED_THETA_CYCLE_LENGTH(2) | cycleTotalBinCount<in.ALLOWED_THETA_CYCLE_LENGTH(1);

% Remove bad stuff
isBadCycle = cycleHasBadLength | cycleHasBadPower;
isInBadCycle = isBadCycle(cycleLabel);
isBad = isInBadCycle | isNegFreq;
phaseAdjusted(isBad) = NaN;
amplitudeAdjusted = in.filteredEEG;
amplitudeAdjusted(isBad) = NaN;
cycleLabel(isBad) = 0;
ret.phase = phaseAdjusted;
ret.amplitude = amplitudeAdjusted; 
ret.cycleLabel = cycleLabel;
end


function [phaseMinimum] = GetPhaseOfMinimumSpiking(spikePhase)

kernelLength = 180;
sigma = kernelLength/4;
kernel = fspecial('gaussian',[kernelLength 1],sigma);

bins = -179.5:1:180;
phaseDistr = hist(spikePhase/pi*180, bins);
% figure(2), subplot(1,2,1), bar(bins, phaseDistr); axis tight
phaseDistr = imfilter(phaseDistr', kernel, 'circular');
% figure(2), subplot(1,2,2), bar(bins, phaseDistr); axis tight
phaseMinimum = bins(ceil(mean(find(phaseDistr==min(phaseDistr)))));
end

function [x] = fixAngle(x)
% Ensure angles are between -pi and pi.
% x must be in radians

x = mod(x+pi,2*pi)-pi;

end

function in = DealWithInputs(varargin)

    defaults.spikeTS = [];
    defaults.phase = [];
    defaults.filteredEEG = [];
    defaults.esr = [];
    defaults.CANONICAL_MINIMUM_SPIKING_PHASE = 0;
    defaults.MEAN_POWER_PERCENTILE_THRESHOLD = [20 42];
    defaults.ALLOWED_THETA_CYCLE_LENGTH = pi;
    
    VERSION = 1;
    
    % Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0
       in = defaults;
       return;
    end
    if ischar(varargin{1})
        if mod(nargin,2) 
            error(['%s called with %d inputs, the first of which is a string. ' ...
                ' Either the first input should be a structure or there should be an even' ...
                ' number of inputs specifying structure field names and values.'],mfilename,nargin);
                
        end
        user_in = cell2struct(varargin(2:2:end)',varargin(1:2:end));
    else
        user_in = varargin{1};
        if nargin > 1 && VERSION ~= varargin{2}
            error(['%s called with version number %g, but code is version %g.\n' ...
                    'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
                    'the function without a version number.'],mfilename,varargin{2},VERSION);
        end
    end
    in = setstructfields(defaults,user_in);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end