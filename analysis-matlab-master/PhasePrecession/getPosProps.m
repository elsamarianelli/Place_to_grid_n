function [ret] = getPosProps(varargin)
% This function uses the output of partionFields and returns the following vectors 
% of the same length as pos:
%
% pos.fieldLabel: the field label for each pos sample (or zero outside the fields)
% pos.runLabel: the run label for each pos smaple (or zero outside the fields)
% pos.r: the ratio of [distance from pos to centre] 
%                         to [distance from perimiter to centre]
%         or NaN outside the field
% pos.phi: the difference  between the pos's angle to the peak and the run's mean dir or Nan
%               outside the field.
% The values of r and phi have been smoothed in cartesian space
% pos.xy: the r and phi in cartesian coordinates 
% pos.xyDir: the direction of travel according to the xy coordiantes
% pos.d_meandir is r projected onto the run's mean direction
% pos.d_currentdir is r projected onto the current direction
% pos.xy_old: the basic xy values
% pos.xyDir_old: the direction arcording to the basic xy values
%
% and the following vectors of length nRuns:
% runs.
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

spkPosInd = ceil(in.spikeTS*in.psr);
spkPosInd(spkPosInd > length(in.xy)) = length(in.xy);

nPos = size(in.xy,1);

% store these in pos structure, to be extra clear what is what
pos.xy_old = in.xy;
pos.xyDir_old = in.dir;

% Convert data to bin units, this must be identical to the way we did it in partionFields
xyBinUnits = toBinUnits(in.xy, [min(in.xy); max(in.xy)], in.binsize);
xyBinInds = ceil(xyBinUnits);

% Pick out the field id for each pos
pos.fieldLabel = in.rmFieldLabel(sub2ind(size(in.rmFieldLabel),xyBinInds(:,1), xyBinInds(:,2)));

% Get the xy coordinates for all the points on the field perimiters
fieldPerimMask = bwperim(logical(in.rmFieldLabel));
fieldPerimMask([1 end],:) = false;
fieldPerimMask(:,[1 end]) = false;
[fieldPerimX, fieldPerimY] = find(fieldPerimMask);
fieldPerimInBinUnits = [fieldPerimX fieldPerimY] - 0.5; % we want the centre of the bin, kind of
fieldPerimXY = bsxfun(@times,fieldPerimInBinUnits,in.binsize);
fieldPerimXY = bsxfun(@plus,fieldPerimXY,min(in.xy));

posRUnsmoothed = nan(nPos,1);
posAngleFromPeak = nan(nPos,1);
perimAngleFromPeak = nan(length(fieldPerimXY),1);
for f = 1:size(in.peaksXY,1) % for each field

    %pick out this field's peak, the perimiter coordinates and the pos samples 
    %that lie within the field
    peaksXY_f = in.peaksXY(f,:);
    fieldPerimXY_f = fieldPerimXY(in.rmFieldLabel(fieldPerimMask)==f,:);
    xy_f = in.xy(pos.fieldLabel == f,:);
        
    %for each point on the perimiter and each pos sample within the field,
    %calculate the angle from the field peak
    perimAngleFromPeak_f = atan2(fieldPerimXY_f(:,2)-peaksXY_f(2),fieldPerimXY_f(:,1)-peaksXY_f(1));
    posAngleFromPeak_f = atan2(xy_f(:,2)-peaksXY_f(2),xy_f(:,1)-peaksXY_f(1));
    posAngleFromPeak(pos.fieldLabel == f) = posAngleFromPeak_f; % we need this again lower down
    perimAngleFromPeak(in.rmFieldLabel(fieldPerimMask)==f) = perimAngleFromPeak_f; %just used for plotting
    
    %for each pos sample work out which point on the perimiter is most
    %colinear with the field centre. +/- 2*pi complication with circular
    %distance between perimAngleFromPeak and posAngleFromPeak
    angleDiff_f = circ_abs(bsxfun(@minus,perimAngleFromPeak_f(:),posAngleFromPeak_f(:)'));
    [ignore,chosenPerimInd_f] =  min(angleDiff_f,[],1);
    
    %Calculate the distance to the peak from pos and from the chosen perim
    %point, and then calculate the ratio.
    distFromPosToPeak_f = sqrt(sum(bsxfun(@minus,xy_f,peaksXY_f).^2,2));
    distFromPerimToPeak_f = sqrt(sum(bsxfun(@minus,fieldPerimXY_f(chosenPerimInd_f,:),peaksXY_f).^2,2));
    posRUnsmoothed(pos.fieldLabel == f) = distFromPosToPeak_f./distFromPerimToPeak_f;
end    

% Label each non-zero contigouous section of posLabels (a "run") with a unique id 
pos.runLabel = LabelContiguousNonZeroSections(pos.fieldLabel);
posIsRun = pos.runLabel > 0;
runs.posStartInd = GetLabelStarts(pos.runLabel);
runs.posEndInd = GetLabelEnds(pos.runLabel);

% Find runs that are too short or speed drops too low or have no spikes
runsWithoutSpikes = true(size(runs.posStartInd));
spkRunLabels = pos.runLabel(spkPosInd);
runsWithoutSpikes(spkRunLabels(spkRunLabels>0)) = false;

xySpeedSmooth = smooth(in.spd,in.smoothingWidthInPosBins);
runs.durationInPosBins = (runs.posEndInd-runs.posStartInd + 1);
runs.minSpeed = accumarray(pos.runLabel(posIsRun),xySpeedSmooth(posIsRun),[],@min);
isBad = runs.minSpeed < in.RUN_MINIMUM_SPEED | runs.durationInPosBins < in.RUN_MINIMUM_DURATION | runsWithoutSpikes;
% isBad = runs.durationInPosBins < in.RUN_MINIMUM_DURATION | runsWithoutSpikes;

% Relabel runs from 1:n after rejecting isBad runs.
pos.runLabel = applyFilterToLabels( ~isBad, pos.runLabel );

runs.posStartInd = runs.posStartInd(~isBad);
runs.posEndInd = runs.posEndInd(~isBad);
runs.minSpeed = runs.minSpeed(~isBad);
runs.durationInPosBins = runs.durationInPosBins(~isBad);
posIsRun = pos.runLabel > 0;

% caclulate the mean and std direction for each run (inlined from circ_stats toolbox)
runComplexMeanDir = accumarray(pos.runLabel(posIsRun),exp(1i*in.dir(posIsRun)*(pi/180)),size(runs.posStartInd));
runs.meanDir = angle(runComplexMeanDir); %circ_mean
% runs.tortuosity = sqrt(2 - 2*abs(runComplexMeanDir)./runs.durationInPosBins); %circ_std
runs.tortuosity = 1 - abs(runComplexMeanDir)./runs.durationInPosBins; % circ_var

% calculate the angular distance between the run's mean direction and the pos's diretion to the peak centre 
posPhiUnsmoothed = nan(size(pos.fieldLabel));
posPhiUnsmoothed(posIsRun) = posAngleFromPeak(posIsRun) - runs.meanDir(pos.runLabel(posIsRun));

% Smooth (r,phi) coordinates in cartesian space: steps (a)-(c)
% (a) convert to cartestian space
[posXUnsmoothed,posYUnsmoothed] = pol2cart(posPhiUnsmoothed(:),posRUnsmoothed(:));
posXYUnsmoothed = [posXUnsmoothed posYUnsmoothed];

% (b) loop through runs, smoothing them with the correctly sized filter.
filtLength = floor((runs.posEndInd-runs.posStartInd+1)*in.RUN_SMOOTH_WINDOW_FRACTION);
pos.xy = nan(nPos,2);
for r=1:numel(runs.posStartInd) % for each run
    if filtLength(r) > 1
        win_r = blackman(filtLength(r));
        filt_r = fir1(filtLength(r)-1,in.SPATIAL_LOW_PASS_CUTOFF/in.psr*2,win_r);
        pos.xy(runs.posStartInd(r):runs.posEndInd(r),:) = filtfilt(filt_r,1,posXYUnsmoothed(runs.posStartInd(r):runs.posEndInd(r),:));
    end
end

% (c) convert back to polar space
[pos.phi, pos.r] = cart2pol(pos.xy(:,1),pos.xy(:,2));
pos.r(pos.r>1) = 1;

% Calculate direction of the smoothed data 
pos.xyDir = atan2(diff(pos.xy(:,2)),diff(pos.xy(:,1)));
pos.xyDir(end+1) = pos.xyDir(end);
pos.xyDir(runs.posEndInd) = pos.xyDir(runs.posEndInd-1);

% Project the distance value onto the current direction
pos.d_currentdir = pos.r.*cos(pos.xyDir-pos.phi);

% calculate the cumulative distance traveled on each run
dr = sqrt(sum(diff(pos.r,1).^2,2));
pos.d_cumulative = labeledCumsum([0 ; dr],pos.runLabel);

% calculate the cumulative sum of the expected normailiesd firing rate
pos.expectedRate_cumulative = labeledCumsum(1-pos.r,pos.runLabel);

% the direction projected onto the run mean direction is simply the x coordiante
pos.d_meandir = pos.xy(:,1); 

%smooth binned spikes to get an instantaneous firing rate
spikeCount = accumarray(spkPosInd,1,size(in.dir));
pos.instantaneousFiringRate = filter(in.ifr_kernel,1,spikeCount)*in.psr; 
pos.instantaneousFiringRate(~posIsRun) = NaN;

%find time spent within run
time = ones(size(pos.xyDir));
time = labeledCumsum(time,pos.runLabel);
pos.timeSpentWithinRun = time/in.psr;

% add some other bits to runs
runs.fieldNum = pos.fieldLabel(runs.posStartInd);
runs.meanSpeed = accumarray(pos.runLabel(posIsRun),in.spd(posIsRun),size(runs.fieldNum),@mean);
runs.centralPeripheral = accumarray(pos.runLabel(posIsRun),abs(pos.xy(posIsRun,2)),size(runs.fieldNum),@mean);

if in.PLOT_ON
    PlotForDebug(fieldPerimMask,perimAngleFromPeak,xyBinInds,pos.r,pos.phi,toBinUnits(in.peaksXY, [min(in.xy); max(in.xy)], in.binsize))
end

ret.pos = pos;
ret.runs = runs;
end


function PlotForDebug(fieldPerimMask,fieldPerimAngle,xyBinInds,posR,posPhi,peakBinUnits)
figure(1)
subplot(4,3,8)

angleCMInd = round(fieldPerimAngle/pi * 180) + 180;
angleCMInd(angleCMInd == 0) = 360;

im = single(~fieldPerimMask);
im(fieldPerimMask) = angleCMInd;
im = ind2rgb(im',hsv(360));
im(repmat(~fieldPerimMask',[1 1 3])) = 1;

isRun = ~isnan(posR);
runBinInds = xyBinInds(isRun,:);

runVals = zeros(size(fieldPerimMask'));
runVals(sub2ind(size(runVals),runBinInds(:,2),runBinInds(:,1))) = posR(isRun);
im(repmat(logical(runVals(:)),[3 1 1])) = repmat(runVals(logical(runVals)),[3 1 1]);

image(im); axis image xy;
hold on;
plot(peakBinUnits(:,1),peakBinUnits(:,2),'rx');
title('runs shaded by pdcd', 'FontSize', 8)
xlabel('bins', 'FontSize', 8); ylabel('bins', 'FontSize', 8); 
set(gca, 'FontSize', 8)
hold off;

figure(1)
subplot(4,3,9)

runMask = false(size(fieldPerimMask'));
runMask(sub2ind(size(runVals),runBinInds(:,2),runBinInds(:,1))) = true;

runVals = zeros(size(fieldPerimMask'));
runVals(sub2ind(size(runVals),runBinInds(:,2),runBinInds(:,1))) = posPhi(isRun);
runInds=round(runVals/pi*180) + 180;
runInds(runInds==0) = 360;
im = ind2rgb(runInds,hsv(360));
im(repmat(~runMask,[1 1 3])) = 1;
im(repmat(fieldPerimMask',[1 1 3])) = 0;

image(im); axis image xy;
hold on;
plot(peakBinUnits(:,1),peakBinUnits(:,2),'rx');
title('runs coloured by direction (hsv)', 'FontSize', 8)
set(gca, 'FontSize', 8)
xlabel('bins', 'FontSize', 8); ylabel('bins', 'FontSize', 8); 
hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following sub functions encapsulate some vectorisation tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = LabelContiguousNonZeroSections(X)
    % e.g. [0 0 6 6 6 10 10 0 0 0 6 6 0 0 0 1 1] gives...
    %      [0 0 1 1 1  2  2 0 0 0 3 3 0 0 0 4 4]
    X = X(:);
    L = cumsum(X &  [true ; X(1:end-1) ~= X(2:end)]);
    L(~X) = 0;
end

function starts = GetLabelStarts(X)
    % e.g. [0 0 6 6 6 10 10 0 0 0 6 0 0 0 0 1 1] gives...
    %      [    3     6           11        16 ] (white space for clarity here)
    starts = find(X &  [true ; X(1:end-1) ~= X(2:end)]);
end

function ends = GetLabelEnds(X)
    % e.g. [0 0 6 6 6 10 10 0 0 0 6 0 0 0 0 1 1] gives...
    %      [        5    7        11         17] (white space for clarity here)
    ends = find(X & [X(1:end-1) ~= X(2:end) ; true]);
end


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

function [x] = circ_abs(x)

x = abs(mod(x+pi,2*pi)-pi);

end

function in = DealWithInputs(varargin)

    defaults.xy = [];
    defaults.dir = [];
    defaults.spd = [];
    defaults.rmFieldLabel = [];
    defaults.spikeTS = 0;
    defaults.psr = 0;
    defaults.SPATIAL_LOW_PASS_CUTOFF = 3;
    defaults.binsize = [0.5 0.5];
    defaults.smoothingWidthInPosBins = 15;
    defaults.RUN_MINIMUM_SPEED = 2.5;
    defaults.RUN_MINIMUM_DURATION = 2;
    defaults.RUN_SMOOTH_WINDOW_FRACTION = 1/3;
    defaults.ifr_kernel = 0;
    defaults.PLOT_ON = false;

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