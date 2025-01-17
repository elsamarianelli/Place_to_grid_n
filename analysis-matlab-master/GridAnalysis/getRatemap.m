function [ out ] = getRatemap ( varargin )
% See https://github.com/UCL/mTint/tree/master/LoadData for info

    if nargin == 0, out = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});
    
    if strcmpi(in.MAP_TYPE, 'rate')
        pos = in.posXY;
        binSize = [in.BIN_SIZE_CM in.BIN_SIZE_CM];
        maxPos = max(pos,[],1);
        minPos = min(pos,[],1);
    elseif strcmpi(in.MAP_TYPE, 'dir')
        pos = in.posDir;
        binSize = in.BIN_SIZE_DEG;
        maxPos = 360;
        minPos = 0;
    end
    
    spikePosInd = in.spikePosInd;
    
    % Get rid of pos samples with a nan for x or y value
    bad = any(isnan(pos),2);
    mapping = cumsum(~bad);
    pos(bad,:) = [];
    spikePosInd(bad(spikePosInd)) = [];
    spikePosInd = mapping(spikePosInd);
    
    % Convert degrees to bin units
    posBinInd = ceil(toBinUnits(pos,[minPos; maxPos],binSize));
    spikeBinInd = posBinInd(spikePosInd,:);
    
    % Calculate bin-dwell-values
    dwell = accumarray(posBinInd,1);
    dwell = dwell/in.pos_samprate;
    
    % Remember unvisted bins and then smooth
    nodwell = (dwell==0); 
    dwell = filter2(ones(in.SMOOTHING),dwell);
    
    % Calculate bin-spike-values and smooth
    spike = accumarray(spikeBinInd,1,size(dwell));
    spike = filter2(ones(in.SMOOTHING),spike);
    
    % Do division and set unvisited bins to NaN
    map = spike./dwell;
    map(nodwell) = NaN;
    pkfr = floor(max(map(find(~isnan(map))))*10)/10;
    out.pkfr = pkfr;
    out.map = map;
    out.nodwell = nodwell;
    if in.PLOT_ON, PlotForDebug(map,pkfr,binSize(1),in); end

end

function PlotForDebug(map,pkfr,binSize,in)
    firstBins = [0;0] + binSize/2;
    lastBins = (fliplr(size(map))'-0.5) * binSize;
    tintplot(struct('im',map,'AXIS_LIMS',[firstBins lastBins],'NUM_COLORS',in.NUM_COLORS));
    text(lastBins(1), lastBins(2),num2str(pkfr));
    axis off;
end

function in = DealWithInputs(varargin)
    defaults.pos_samprate = 50;
    defaults.spikePosInd = [];
    defaults.posXY = zeros(0,2);
    defaults.posDir = zeros(0,1);
    defaults.SMOOTHING = 5;
    defaults.BIN_SIZE_CM = 2.5; 
    defaults.BIN_SIZE_DEG = 6;
    defaults.PLOT_ON = true;
    defaults.MAP_TYPE = 'rate';
    defaults.NUM_COLORS = 5;
    
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