function [ret] = getSpikeProps(varargin)
%% This function builds the spk strucutre and augments runs structure
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

spkEegInd = ceil(in.spikeTS*in.eegSampleRate);
spkEegInd(spkEegInd > numel(in.theta.oldPhase)) = numel(in.theta.oldPhase);
spkPosInd = ceil(in.spikeTS*in.posSampleRate);
spkPosInd(spkPosInd > length(in.pos.xyDir_old)) = numel(in.pos.xyDir_old);

spk.time = in.spikeTS;
spk.posInd = spkPosInd; 
spk.eegInd = spkEegInd; 

spk.runLabel = in.pos.runLabel(spkPosInd);
spk.thetaCycleLabel = in.theta.cycleLabel(spkEegInd);

%  build mask that is true for spikes that are the first in their theta cycle
spk.isFirstInTheta = [true ; spk.thetaCycleLabel(1:end-1) ~= spk.thetaCycleLabel(2:end)];
spk.isLastInTheta = spk.isFirstInTheta([2:end 1]);

% calculate two kinds of numbering for spikes within a run
spk.numWithinRun = labeledCumsum(ones(size(spkPosInd)),spk.runLabel);
spk.thetaBatchLabelWithinRun = labeledCumsum(spk.isFirstInTheta,spk.runLabel);

% add some stuff to runs
in.runs.spikeCount = accumarray(spk.runLabel(spk.runLabel>0),1,size(in.runs.meanDir));
in.runs.rateInPosBins = in.runs.spikeCount./in.runs.durationInPosBins;

ret.spk = spk;
ret.runs = in.runs;
% PlotForDebug()

end

% function PlotForDebug()
%  %TODO: implement something here
% end


function S = labeledCumsum(X,L)
% e.g. X=[3   5 8  9 1  2    5    8 5  4  9   2]
%      L=[0   1 1  2 2  0    0    1 1  1  0   5]
% result=[NaN 5 13 9 10 NaN  NaN  8 13 17 NaN 2]

L = L(:);
X = X(:);

if numel(L) ~= numel(X)
    error('The two inputs should be vectors of the same length.')
end

% Do the full cumulative sum
X(isnan(X)) = 0;
S = cumsum(X);

mask = logical(L);

% Lookup the cumulative value just before the start of each segment
isStart = mask & [true ; L(1:end-1) ~= L(2:end)];
startInds = find(isStart);
if startInds(1) == 1
    S_starts = [0 ; S(startInds(2:end)-1)];
else
    S_starts = S(startInds-1);
end

% Subtract off the excess values (i.e. the ones obtained above)
L_safe = cumsum(isStart); % we do this to main case the labels in L were not sorted integers
S(mask) = S(mask) - S_starts(L_safe(mask));

% Put NaNs in the false blocks
S(L==0) = NaN; 
end

function in = DealWithInputs(varargin)

    defaults.pos = [];
    defaults.theta = [];
    defaults.runs = [];
    defaults.spikeTS = [];
    defaults.posSampleRate = 0;
    defaults.eegSampleRate = 0;

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