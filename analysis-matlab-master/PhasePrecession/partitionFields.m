function [ret] = partitionFields(varargin)
% Partition fields
% To partition spikes into fields. We find watersheds around the peaks of
% a super smooth ratemap
%
% Returns:
% peaksXY:  matrix of size px2, giving the (x,y) coordinates of the field
%           peaks, using
%           the same units as the input xy matrix.
% peaksRate: vector giving the firing rate at each of the peaks
% rmFieldLabel: matrix labeling each point in the ratemap with a field
%           number or zero outside fields
%
%%    Copyright (C) <2013>  <Daniel Manson> <d.manson@ucl.ac.uk> &
%%                           <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
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

%% Boilerplate for standardised input/output behaviour
if nargin == 0, ret = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

% Get indices of spikes in to pos
spikePosInd = ceil(in.spikeTS*in.psr);
spikePosInd(spikePosInd > length(in.xy)) = length(in.xy);

% Convert data to bin units
xyBinUnits = toBinUnits(in.xy, [min(in.xy); max(in.xy)], in.binsize);

% Create ratemap and smooth
xyBinInds = ceil(xyBinUnits);
dwell = accumarray(xyBinInds ,1)/in.psr;
spike = accumarray(xyBinInds(spikePosInd,:),1,size(dwell));
map = spike./dwell;
map(dwell == 0) = 0; % in order to do smoothing use 0/0=0 rather than NaN
map = imfilter(map, in.smoothingKernelForPartiton, 'symmetric');

% Watershed the map and label each region
watershedMask = watershed(-map);
rmFieldLabel = bwlabel(watershedMask );
spikesInField = accumarray(rmFieldLabel(rmFieldLabel>0), spike((rmFieldLabel>0)), [], @sum);

% Find the locations of the peaks in proper cm units
rmMaxMask = imregionalmax(map); % annoyingly watershed has already done this, but we don't get to see the results
rmMaxMask = bwmorph(rmMaxMask,'shrink', Inf);
[xi, yi] = find(rmMaxMask);
peaksInBinUnits = [xi yi] - 0.5; % we want the centre of the bin
peaksXY = bsxfun(@times,peaksInBinUnits,in.binsize);
peaksXY = bsxfun(@plus,peaksXY,min(in.xy));

% Threshold each field at some fraction of the relevant field's peak value
peaksRate(rmFieldLabel(rmMaxMask)) = map(rmMaxMask);
fieldThreshold = peaksRate' * in.FIELD_THRESHOLD_FRAC;
rmFieldMask = rmFieldLabel > 0;
rmFieldMask(rmFieldMask) = map(rmFieldLabel >0) > fieldThreshold(rmFieldLabel(rmFieldLabel >0));
rmFieldLabel(~rmFieldMask) = 0;

% Reorder peaksXY and peaksRate to use same labeling as rmFieldLabel
peakBinInds = ceil(peaksInBinUnits);
peakLabels = rmFieldLabel(sub2ind(size(rmFieldLabel),peakBinInds(:,1),peakBinInds(:,2)));
peaksXY(peakLabels,:) = peaksXY;
peaksRate(peakLabels) = peaksRate;
peakBinInds(peakLabels,:) = peakBinInds; %only used in plotting
peaksInBinUnits(peakLabels,:) = peaksInBinUnits;

% We may have been asked to discard small edge fields
if ~isnan(in.AREA_THRESHOLD_FRAC)
    
    %     % Work out field areas and whether the field peak is at the edge of the map
    %     fieldAreas = accumarray(rmFieldLabel(rmFieldLabel>0),1); %note to get area in cm^2 need to times by prod(in.binsize)
    %     isEdgeField = peaksInBinUnits(:,1) < 1 ...  %for this statement remember that peaksInBinUnits is the bin centre
    %         | peaksInBinUnits(:,2) < 1 ...
    %         | peaksInBinUnits(:,1) + 1 > size(map,1) ...
    %         | peaksInBinUnits(:,2) + 1 > size(map,2);
    
    % Work out field areas and whether the field peak is at the edge of the map
    fieldAreas = accumarray(rmFieldLabel(rmFieldLabel>0),1); %note to get area in cm^2 need to times by prod(in.binsize)
    dis2edge = [peaksInBinUnits(:,1), peaksInBinUnits(:,2),...
        size(map,1) - peaksInBinUnits(:,1),...
        size(map,2) - peaksInBinUnits(:,2)];
    minDis2Edge = min(dis2edge,[],2)*in.binsize(1);
    in.AREA_THRESHOLD_FRAC = max([0.5*ones(size(minDis2Edge)),0.75 - minDis2Edge*((0.75-0.5)/10)], [], 2);
    
    % Remove any edge fields that are too small
    if strcmp(in.type,'p')
        discardField = true(size(spikesInField));
        [ignore,maxSpikesInd] = max(spikesInField);
        discardField(maxSpikesInd) = false;
    else
        %         bigEnoughField = fieldAreas > in.AREA_THRESHOLD_FRAC * mean(fieldAreas);
        %         discardField = (isEdgeField & ~bigEnoughField) ;
        %         discardField = false(size(spikesInField));
        discardField = fieldAreas < in.AREA_THRESHOLD_FRAC * mean(fieldAreas);
    end
    rmFieldLabel = applyFilterToLabels( ~discardField, rmFieldLabel );
    peaksXY(discardField,:) = [];
    peaksRate(discardField) = [];
    
    ret.peaksXY = peaksXY;
    ret.peaksRate = peaksRate;
    ret.rmFieldLabel = rmFieldLabel;

end

if in.PLOT_ON
    PlotForDebug(map,rmFieldLabel,peakBinInds)
end

end

function PlotForDebug(map,fieldLabel,peakBinInds)
% clf;
figure(1)
subplot(4,3,4)
imagesc(map'), axis image xy
hold on;
plot(peakBinInds(:,1),peakBinInds(:,2),'kx');
title('2D firing rate map', 'FontSize', 8)
xlabel('bins', 'FontSize', 8); ylabel('bins', 'FontSize', 8); 
set(gca, 'FontSize', 8)
hold off;

subplot(4,3,7);
imagesc(fieldLabel'), axis image xy
hold on;
plot(peakBinInds(:,1),peakBinInds(:,2),'kx');
title('watershed field isolation', 'FontSize', 8)
xlabel('bins', 'FontSize', 8); ylabel('bins', 'FontSize', 8); 
set(gca, 'FontSize', 8)
hold off;

end

function in = DealWithInputs(varargin)
 
    defaults.xy = [];
    defaults.psr = [];
    defaults.spikeTS = [];
    defaults.type = [];
    defaults.binsize = 0;
    defaults.smoothingKernelForPartiton = 0;
    defaults.FIELD_THRESHOLD_FRAC = 0;
    defaults.AREA_THRESHOLD_FRAC = 0;
    defaults.PLOT_ON = 0;
    
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