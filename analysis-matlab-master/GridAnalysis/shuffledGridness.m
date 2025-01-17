function [ret] = shuffledGridness(varargin)
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk> &
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

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   

    % Get size of inputs
    nPos = size(in.posXY,1);

    % Get random values, shuffle spikes, and define position indices of spikes
    shiftRange = (nPos/in.psr) - 2*in.MIN_SHIFT_SECONDS;
    shift = rand(1,in.SHUFFLES)*shiftRange + in.MIN_SHIFT_SECONDS;
    shift = [0 shift]; % the first element is going to be unshifted to give us the true gridness
    
    spikeTS = bsxfun(@plus,in.spikeTS,shift);
    spkPosInd = ceil(spikeTS*in.psr);

    spkPosInd = mod(spkPosInd - 1, nPos) + 1;

    % Create Ratemaps
    % Convert data to bin units
    xyBinUnits = toBinUnits(in.posXY, [min(in.posXY); max(in.posXY)], in.BINSIZE);

    % Define position indices of path and create occupancy map
    xyBinInds = ceil(xyBinUnits);
    dwell = accumarray(xyBinInds ,1)/in.psr;

    % Remember unvisited bins
    nodwell = dwell==0;
    dwell = imfilter(dwell,ones(in.SMOOTHING)/in.SMOOTHING^2);

    % Create a stack of spike maps for each timeshift.
    % Index to indicate which ratemap to bin into
    shiftInds = repmat(1:numel(shift), size(in.spikeTS,1), 1);
    % Indices for xy position of each spike x timestamp.
    xInds = xyBinInds(spkPosInd,1);
    yInds = xyBinInds(spkPosInd,2);

    % Bin spikes
    spike = accumarray([xInds(:) yInds(:) shiftInds(:)],1,[size(dwell),numel(shift)]);

    % Smooth
    spike = imfilter(spike,ones(in.SMOOTHING)/in.SMOOTHING^2);

    % Make ratemap
    map = bsxfun(@rdivide, spike, dwell);

    % Set ratemaps to NaN where unvisited
    % map = reshape(map,[size(map,1)*size(map,2), size(map,3)]);
    % map(nodwell(:),:) = NaN;
    % map = reshape(map,size(spike));

    autCorr2D_in.x = map;
    autCorr2D_in.nodwell = nodwell;
    autCorr2D_in.tol = in.TOL;
    autoCorr2D_out = autoCorr2D(autCorr2D_in);
    shuffledAutoCorr = autoCorr2D_out.autocorrelogram;
    
    % Compute the max power for shifted versions of the data
    gridnessShuffles = zeros(in.SHUFFLES + 1,1);
    autocorr_in = in.autocorr_in;
    for jj = 1:in.SHUFFLES + 1
        autocorr_in.sac = shuffledAutoCorr(:,:,jj);
        autocorr_out = autoCorrProps(autocorr_in);
        gridnessShuffles(jj) = autocorr_out.gridness;
    end

    ret.gridness = gridnessShuffles(1); %remember we have 0-shift for the first element
    ret.gridnessShuffles = gridnessShuffles(2:end);
    
    if in.PLOT_ON, PlotForDebug(ret.gridnessShuffles,ret.gridness); end;
end

function PlotForDebug(gridnessShuffles,gridness)
    cdfplot(gridnessShuffles); hold on;
    plot([gridness gridness],[0 1],'r-'); hold off;
    xlabel('gridness');
    ylabel('cumulative fraction');
    title('cumulative distribution of gridness for shuffle data');
end

function in = DealWithInputs(varargin)
    defaults.posXY = [];
    defaults.psr = [];
    defaults.spikeTS = [];
    defaults.BINSIZE = [2 2]; %cm
    defaults.SMOOTHING = 5;
    defaults.MIN_SHIFT_SECONDS = 1;
    defaults.SHUFFLES = 1000;
    defaults.TOL = 1e-10;
    defaults.PLOT_ON = true;

    defaults.autocorr_in.GET_MEAN_R_AND_SCALE = false;
    defaults.autocorr_in.GET_ORIENTATION = false;
    defaults.autocorr_in.FIELD_EXTENT_METHOD = 2;
    defaults.autocorr_in.PLOT_ON = false;
    
    VERSION = 1;
    
    % Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0
        in = defaults;
        return;
    end

    if isstruct(varargin{1})
        if nargin > 1 && VERSION ~= varargin{2}
            error(['%s called with version number %g, but code is version %g.\n' ...
                'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
                'the function without a version number.'],mfilename,varargin{2},VERSION);
        end
        in = ModifyExistingFields(defaults,varargin{1});
    else
        in = ModifyExistingFields(defaults,varargin{:});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
