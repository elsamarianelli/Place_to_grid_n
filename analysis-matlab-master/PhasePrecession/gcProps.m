function [ret] = gcProps(varargin)
%% Calculate spatial properties of grid cell firing
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
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

%% Get cell properties. Compute Gridness, orientation and scale
spkPosInd = ceil(in.spikeTS*in.psr);
spkPosInd(spkPosInd > length(in.xy)) = size(in.xy,1);

% Convert data to bin units
xyBinUnits = toBinUnits(in.xy, [min(in.xy); max(in.xy)], in.binsize);

% Create ratemap and smooth
xyBinInds = ceil(xyBinUnits);
dwell = accumarray(xyBinInds ,1)/in.psr;
spike = accumarray(xyBinInds(spkPosInd,:),1,size(dwell));

% Remember unvisited bins
nodwell = dwell==0;
dwell = filter2(ones(5)/25, dwell);

%Calculate bin-spike-values
spike = filter2(ones(5)/25, spike);

% Make ratemap
map = spike./dwell;
map(nodwell) = NaN;

% Compute sac
in.ac2D_in.x = map;
in.ac2D_in.nodwell = nodwell;
if in.ac2D_in.PLOT_ON
    figure(1)
    subplot(4,3,5);
    title('2D autocorrelation', 'FontSize', 8)
    set(gca, 'FontSize', 8)
    hold on
end
ac2D_out = autoCorr2D(in.ac2D_in);

% Get gridness, orientation, scale
in.acProps_in.sac = ac2D_out.autocorrelogram;
if in.acProps_in.PLOT_ON
    figure(1)
    subplot(4,3,6);
    set(gca, 'FontSize', 8)
    hold on
end
acProps_out = autoCorrProps(in.acProps_in);
if in.acProps_in.PLOT_ON
    title(sprintf('Grid cell metrics \n G %s | S %s | D %s',...
        num2str(acProps_out.gridness,2), num2str(acProps_out.scale,2),...
        num2str(acProps_out.orientation,2)), 'FontSize', 8)
end
ret.gridness = acProps_out.gridness;
ret.orientation = acProps_out.orientation;
ret.scale = acProps_out.scale*in.binsize(1);

end

function in = DealWithInputs(varargin)

defaults.xy = [];
defaults.psr = [];
defaults.spikeTS = [];
defaults.ac2D_in = autoCorr2D();
defaults.acProps_in = autoCorrProps();

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
