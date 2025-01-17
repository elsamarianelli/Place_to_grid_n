function ret = disc_DAC(varargin)
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   

    %1a. For the sake of efficiency and reducing memory footprint, cast posSamples
    %    from 8-byte doubles to 4-byte singles.  
    posSamp = single(in.posXY);
    nPos = length(posSamp);

    %1b. Produce displacement histogram for pos at spike events
    nSpikePos = length(in.spikeInd);
    dispHistSpike = BatchedDisplacementHist(posSamp(in.spikeInd,:), in.BATCH_SIZE, in.binSize);

    %1c. Produce displacement histogram for all of pos, or rather for a downsampled
    %    version of all of pos.
    nDownSampPos = round(nPos*in.downSampPos);
    downSampInds = round(linspace(1, nPos, nDownSampPos));
    dispHistDwell = BatchedDisplacementHist(posSamp(downSampInds,:),  in.BATCH_SIZE, in.binSize);

    %2a. Convert both histograms to densities
    dispHistDwell = dispHistDwell/in.binSize;
    dispHistSpike = dispHistSpike/in.binSize;


    if in.smthKernSigma ~= 0
        %2b. If requested, smooth the two histograms
        sigmaInBins = in.smthKernSigma/in.binSize; %convert from cm to bins
        kern = fspecial('gaussian',[1, sigmaInBins*in.smthKernWidthNSigmas], sigmaInBins);
        dispHistDwell =  imfilter(dispHistDwell, kern', 'symmetric', 'conv');
        dispHistSpike =  imfilter(dispHistSpike, kern', 'symmetric', 'conv');
    end

    %2c. Trim both histograms down to the region of meaningful data
    histLen = ceil(in.DISPLACEMENT_LIMIT_CM/in.binSize);
    dispHistDwell = dispHistDwell(1:histLen);
    dispHistSpike = dispHistSpike(1:histLen);
    ret.binCentres = (0:histLen-1)'*in.binSize;

    %3a. Subtract one probability distribution from the other 
    dispHist = dispHistSpike ./ dispHistDwell;

    ret.dispHist = dispHist;

    if in.PLOT_ON, PlotForDebug(ret.binCentres,ret.dispHist,in.DISPLACEMENT_LIMIT_CM); end;

end

function PlotForDebug(binCentres,dispHist,DISPLACEMENT_LIMIT_CM)
    plot(binCentres,dispHist,'r','LineWidth',2);
    xlim([0 DISPLACEMENT_LIMIT_CM]);
    xlabel('diplacement (cm)');
    ylabel('rate density (per cm)');
end

function dispHist = BatchedDisplacementHist(pos, BATCH_SIZE, binSize)

%1. Work out batch start and end indicies, the only complication being with the last batch
nPos = length(pos);
batchStarts = 1:BATCH_SIZE:nPos;
batchEnds = BATCH_SIZE:BATCH_SIZE:nPos;
nBatches = numel(batchStarts);
if numel(batchEnds) ~= nBatches, batchEnds(end+1) = nPos; end

%2. Work out sqrt((x_max-x_min)^2 + (y_max-y_min)^2), which is upper limit on distance
maxVals = max(pos,[],1);
minVals = min(pos,[],1);
maxDisp = round(sqrt(sum((maxVals-minVals).^2))/binSize) + 1;

% Prepare empty histogram, bin centres represent: 0, binSize, 2binSize, ...
dispHist = zeros(maxDisp,1); 
for aa = 1:nBatches
    xVals_a = pos(batchStarts(aa):batchEnds(aa),1);
    yVals_a = pos(batchStarts(aa):batchEnds(aa),2);
    
    %3. Process all pairs of points in this block
    
    %3a. Find distances
    dist_ab = sqrt( bsxfun(@minus,xVals_a,xVals_a').^2 + ...
                    bsxfun(@minus,yVals_a,yVals_a').^2   ); 
                
    %3b. Convert to bin indicies and add to histogram       
    dist_ab = round(dist_ab/binSize) + 1;
    newHistData = accumarray(dist_ab(:),1,[maxDisp 1]); % accumulate 1 here, but 2 in the loop below
    dispHist = dispHist + newHistData; 

    %4. Process all pairs consisting of a point in this block and a subsequent point
    for bb = aa+1:nBatches
        xVals_b = pos(batchStarts(bb):batchEnds(bb),1);
        yVals_b = pos(batchStarts(bb):batchEnds(bb),2);        
        
        %4a. Find distances
        dist_ab = sqrt( bsxfun(@minus,xVals_a,xVals_b').^2 + ...
                        bsxfun(@minus,yVals_a,yVals_b').^2   );
                    
        %4b. Convert to bin indicies and add to histogram   
        dist_ab = round(dist_ab/binSize) + 1;
        newHistData = accumarray(dist_ab(:),2,[maxDisp 1]); % accumulate 2 due to symmetry
        dispHist = dispHist + newHistData;
    end
end

end

function in = DealWithInputs(varargin)
    defaults.posXY = [];
    defaults.spikeInd = [];
    defaults.binSize = 1;
    defaults.smthKernSigma = 10;
    defaults.smthKernWidthNSigmas = 4;
    defaults.BATCH_SIZE = 2^8;
    defaults.downSampPos = 0.05;
    defaults.DISPLACEMENT_LIMIT_CM = 100;
    defaults.PLOT_ON =true;

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