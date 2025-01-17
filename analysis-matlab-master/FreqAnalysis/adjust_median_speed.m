function ret = adjust_median_speed(varargin)
% See https://github.com/UCL/mTint/tree/master/FreqAnalysis for info

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   

% Sort speeds and get sorting order, so we know which pos sample each speed refers to
[orderedSpeeds,orderedSpeedIndex] = sort(in.speed(:));

% Working with the orderedSpeeds, find the index of:
%  a. the closest value to the grand median speed
%  b. the first speed over the minimum threshold
%  c. the fastest speed (everything after that in orderedSpeeds is nan)
indMedian = find(orderedSpeeds >= in.grandMedian, 1,'first'); %more or less correct
indFirstOverThreshold = find(orderedSpeeds >= in.minSpeed, 1,'first');
indLastNonNan = find(~isnan(orderedSpeeds), 1,'last');

% work out which side of the median has fewer data points
halfWidth = min(indMedian-indFirstOverThreshold,indLastNonNan-indMedian);

% return the pos samples corresponding to speeds in the valid range
ret.allowedPosSamples = orderedSpeedIndex(indMedian-halfWidth:indMedian+halfWidth);



if in.PLOT_ON
    PlotForDebug(indLastNonNan,orderedSpeeds,indMedian,halfWidth,in.grandMedian,in.minSpeed); 
end



end

function PlotForDebug(indLastNonNan,orderedSpeeds,indMedian,halfWidth,grandMedian,minSpeed)
    maxSpeed = orderedSpeeds(indLastNonNan);
    L = length(orderedSpeeds);
    hArea = area([0 maxSpeed],[1 1]*(indMedian+halfWidth),indMedian-halfWidth);
    set(hArea,'FaceColor',[0.9 0.9 1]);
    hold on;
    plot(orderedSpeeds(1:indLastNonNan),1:indLastNonNan,'k','LineWidth',2);
    xlabel('speed (cm/s)');
    ylabel('number of samples (cumulative)');
    if indLastNonNan ~= L
        plot([0 maxSpeed],[1 1]*(indLastNonNan+1),'--r');
        plot([0 maxSpeed],[L L],'--r');
    end
    xlim([0 maxSpeed]);
    ylim([0 L]);
    plot([0 maxSpeed],[1 1]*indMedian,'b','LineWidth',1)
    stem(grandMedian,indMedian,'b');
    plot([0 maxSpeed],[1 1]*indFirstOverThreshold,':b')
    stem(minSpeed,indFirstOverThreshold,':b');
    hold off;
end


function in = DealWithInputs(varargin)
    defaults.pos2use = [];
    defaults.eeg = [];
    defaults.header = [];
    defaults.thetaRange = [7 11];
    defaults.padToPow2 = NaN;
    defaults.sampFreq = 250;
    defaults.smthKernelWidth = 2;
    defaults.smthKernelSigma = 0.1875;
    defaults.s2nWdth = 2;
    defaults.maxFreq = 25;
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
        in = ModifyExistingFields(defaults,varargin{1});
    else
        in = ModifyExistingFields(defaults,varargin{:});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

