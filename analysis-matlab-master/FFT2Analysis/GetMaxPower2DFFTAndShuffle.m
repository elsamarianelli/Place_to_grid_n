function out = GetMaxPower2DFFTAndShuffle(varargin)
% See https://github.com/UCL/mTint/tree/master/FFT2Analysis for info

    if nargin == 0, out = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   

    nPos = size(in.posXY,1);
    minShiftPos = in.MIN_SHIFT_SECONDS * in.posFreq;

    % subtract minimum in x and y dimension and convert to bin indices
    posXY = bsxfun(@minus,in.posXY,min(in.posXY,[],1));
    posXYbin = ceil(posXY/in.BIN_SIZE_CM + eps);

    % Generate 2D-maps giving dwell count and spike count
    dwell_map = accumarray(posXYbin,single(1),[],@sum);
    spike_map = accumarray(posXYbin(in.spikePosInd,:),single(1),size(dwell_map),@sum);

    % Dividing one map by the other gives rate (in units of posFreq)
    rate_map = spike_map ./ dwell_map;

    % Get the FFT2 and its peak
    power_spect = GetMasked2DFFTfromRMap(rate_map, in.FFT_SIZE, in.MAX_FREQ);
    maxPower = max(max(power_spect));

    % Compute the max power for shifted versions of the data
    % Keep going until we get converence of the 95th percentile
    epsilion = Inf;
    perc95 = Inf;
    maxPowerShuffles = zeros(in.MAX_ITERATIONS,1);
    jj = 1;
    while jj <= in.MAX_ITERATIONS && epsilion > in.CONVERGENCE_THRESHOLD

        % shift the spikePosInd by a random time interval
        shifted_spikePosInd = ShiftSpikePosInd(in.spikePosInd, nPos, minShiftPos);

        % compute the new spike_map and the resulting rate_map
        spike_map = accumarray(posXYbin(shifted_spikePosInd,:),single(1),size(dwell_map),@sum);
        rate_map = spike_map ./ dwell_map;

        % Get the FFT2 and its peak
        power_spect = GetMasked2DFFTfromRMap(rate_map, in.FFT_SIZE, in.MAX_FREQ);
        maxPowerShuffles(jj) = max(max(power_spect));

        % Every 100th iteration check for convergence of 95th percentile
        if mod(jj, 100) == 0
            oldPerc95 = perc95;
            perc95 = prctile(maxPowerShuffles(~isnan( maxPowerShuffles(1:jj))), 95);
            disp([num2str(jj) '. perc95=' num2str(perc95)]);
            epsilion = abs(perc95 - oldPerc95);
        end

        jj = jj+1;
    end

    maxPowerShuffles(jj:end) = []; %remove unused elements (if there are any)

    out.maxPower = maxPower;
    out.maxPowerShuffles = maxPowerShuffles;
end


function shifted_spikePosInd = ShiftSpikePosInd(spikePosInd, nPos, minShiftPos)
% TAKES:
% spikePosInd - the index of pos sample corresponding to spike time
% nPos - the number of pos samples
% minShiftPos - the minium permitted temporal shift expresseed in units
%               of pos sample period (e.g. how many 50ths of a second).
%
% RETURNS:
% shifted_spikePosInd - the spikePosInd temporaly shifted, wrapping the
% end back round to the start.
%
% Original code by JK.
% DM, May 2013.  Rewritten and correted minor error.
%

shiftRange = nPos - 2*minShiftPos;
shift = rand(1)*shiftRange + minShiftPos;
shifted_spikePosInd = mod(floor(spikePosInd + shift) - 1, nPos) + 1;

end


function in = DealWithInputs(varargin)
    defaults.BIN_SIZE_CM = 2.5;
    defaults.MIN_SHIFT_SECONDS = 20;
    defaults.MAX_ITERATIONS = 1000;
    defaults.CONVERGENCE_THRESHOLD = 0.05;
    defaults.FFT_SIZE = 256;
    defaults.MIN_FREQ = 0;
    defaults.MAX_FREQ = 10;
    defaults.posFreq = 50;
    defaults.posXY = [];
    defaults.spikePosInd = [];
    
    VERSION = 1.01;
    
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