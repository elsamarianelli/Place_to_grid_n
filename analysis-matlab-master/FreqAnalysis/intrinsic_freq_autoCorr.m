function ret = intrinsic_freq_autoCorr(varargin)
% See https://github.com/UCL/mTint/tree/master/FreqAnalysis for info

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   
    
    % #######################################################################
    % [Stage 1] Construct a cell array of spike-time histograms, one for each
    % contiguous block of trues in the posMask.

    %1a. Work out number of hist bins per pos sample
    acBinsPerPos = 1/in.posSampFreq/in.acBinSize; %e.g. 1/50Hz/0.002s = 10 bins per pos sample

    %1b. Work out full-length of auto-corr window, eausred in bins
    acWindowSizeBins = round(in.acWindow/in.acBinSize);

    %1c. Construct a single spike train histogram using all spike data 
    binCentres = (0.5:numel(in.posMask)*acBinsPerPos)*in.acBinSize;
    spkTrHist = hist(in.spikeTimes, binCentres(:));

    %1d. Find the start, end and length of each contigous block of true in posMask
    [starts,ends] = FindMaskCaps(in.posMask);
    chunkLengths = ends-starts+1;
    nChunks = numel(starts);

    %1e. Split the single histogram up into individual chunks, stored in a cell-array
    starts = (starts-1)*acBinsPerPos + 1;
    ends = ends*acBinsPerPos;
    histChunks = cell(nChunks,1);
    for ii=1:nChunks
        histChunks{ii} = spkTrHist(starts(ii):ends(ii))';
    end    


    % #######################################################################
    % [Stage 2] Get a weighted-average auto-correlogram 

    %2a. Remove chunks with insufficient spikes
    spikesInChunk = cellfun(@sum,histChunks);
    boringChunks = spikesInChunk<2; 
    histChunks(boringChunks) = [];
    chunkLengths(boringChunks) = [];

    %2b. Do autocorr for all remaining histograms
    autoCorrFunction = @(a) xcorr(a,acWindowSizeBins, 'unbiased');
    autoCorrGrid = cellfun(autoCorrFunction, histChunks,'UniformOutput',false); % Matlab doesn't actually execute this very efficiently
    autoCorrGrid = cell2mat(autoCorrGrid');

    %2c. Discard left half of autoCorrs, leaving centre and right half
    autoCorrGrid(1:acWindowSizeBins,:) = []; 

    %2d. Weight autocorr in each chunk by the chunk length
    autoCorrGrid = bsxfun(@times,autoCorrGrid,chunkLengths');

    %2e. Normalise the total mean by the total valid length
    totalValidLength = sum(chunkLengths);
    autoCorrSum = nansum(autoCorrGrid,2)./totalValidLength; 

    %2f. Produce a vector giving the autoCorr lags in seconds
    lags = (0:acWindowSizeBins)*in.acBinSize;

    % #######################################################################
    % [Stage 3] Send the autoCorr off for power spectrum analysis

    %3a. Mean normalise
    meanNormAc = autoCorrSum(2:end)-mean(autoCorrSum(2:end)); % don't use peak at bin 1

    %3b. Send to power_spectrum function for some serious work
    ret = power_spectrum(...
        meanNormAc,in.padToPow2,in.acBinSize,in.thetaRange,in.maxFreq,...
        in.smthKernelWidth,in.smthKernelSigma,in.s2nWdth,in.PLOT_ON,in.ymax,in.xmax);

    %3c. Add the autoCorr and lags to the output structure
    ret.autoCorrSum = autoCorrSum;
    ret.lags = lags;


    if in.PLOT_ON, PlotForDebug(in.maxFreq,ret.power,autoCorrGrid); end;

end

function PlotForDebug(maxFreq,power,autoCorrGrid)
    % Augment the power_spectrum plot
    hold on
    imagesc([0.6 1]*maxFreq, [0.6 1]*max(power),autoCorrGrid')
    hold off
    
end

function [starts,ends] = FindMaskCaps( mask )
    % mask is a logical array.
    % starts and ends are the indicies of the start and end of each contigous
    % block of true section in the mask.
    %
    d = diff(mask);
    starts = find(d==1); %d==1, means one element is true and previous element was false
    ends = find(d==-1)-1;

    %if first or last element was true, we have to add in the extra start and/or end
    if mask(1), starts = [1;starts]; end
    if mask(end), ends = [ends;numel(mask)]; end

end


function in = DealWithInputs(varargin)
    defaults.spikeTimes = [];
    defaults.posMask = [];
    defaults.thetaRange = [7 11];
    defaults.padToPow2 = 16;
    defaults.posSampFreq = 50;
    defaults.acBinSize = 0.002;
    defaults.acWindow = 0.5;
    defaults.smthKernelWidth = 2;
    defaults.smthKernelSigma = 0.1875;
    defaults.s2nWdth = 2;
    defaults.maxFreq = 25;
    defaults.ymax = [];
    defaults.xmax = 25;
    defaults.PLOT_ON = true;
    
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
        in = modifyExistingFields(defaults,varargin{1});
    else
        in = modifyExistingFields(defaults,varargin{:});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end