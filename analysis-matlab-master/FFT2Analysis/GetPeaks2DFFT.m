function out = GetPeaks2DFFT(varargin)
% See https://github.com/UCL/mTint/tree/master/FFT2Analysis for info

    if nargin == 0, out = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   

    S = in.FFT_SIZE;

    rate_map = in.rate_map;

    if isempty(rate_map)
        % subtract minimum in x and y dimension and convert to bin indices
        posXY = bsxfun(@minus,in.posXY,min(in.posXY,[],1)); 
        posXYbin = ceil(posXY/in.BIN_SIZE_CM + eps);

        % Generate 2D-maps giving dwell count and spike count
        dwell_map = accumarray(posXYbin,1,[],@sum,0);
        spike_map = accumarray(posXYbin(in.spikePosInd,:),1,size(dwell_map),@sum,0);

        % Dividing one map by the other gives rate (in units of posFreq)
        rate_map = spike_map ./ dwell_map;    
    end

    % Get the FFT2
    [fftPower, fftDirection, fftFreq] = GetMasked2DFFTfromRMap(rate_map, S, in.MAX_FREQ);
    fftPower = fftPower-in.powerThreshold;
    fftPower(fftPower<0)=0;

    % get the peak orientations in degrees
    [pkDirections] = FindPeakDirectionsInFFT2(fftPower, fftDirection, fftFreq, in.MAX_FREQ, in.PEAK_MAX_FRACTION);

    if isempty(pkDirections) %can be true if power_spect is all zeros following thresholding
        out.half_lambda = [];
        out.pkDirections = [];
        return;
    end

    % Produce a matrix of FFT2 powers, in which columns correspond
    % to the peak orientations and rows correspond to all wavelengths
    r_outer = min(size(fftPower)-1)/2;
    k_vals = (1:r_outer-1)';
    xInd = r_outer+1 + bsxfun(@times, k_vals, cosd(pkDirections));
    yInd = r_outer+1 + bsxfun(@times, k_vals, sind(pkDirections));
    xInd = ceil(xInd);   yInd = ceil(yInd);
    allFFTValsAtPeakDirs = fftPower(sub2ind(size(fftPower),yInd,xInd));

    % find the wavelength with maximum power at each of the peak orientaitons
    [~,kPeakInds] = max(allFFTValsAtPeakDirs,[],1);
    kPeaks = k_vals(kPeakInds) -0.5;
    half_lambda = in.BIN_SIZE_CM*S/2 *kPeaks.^-1; %DM, not sure what JK comment means: "divide by 2 to get the size of the region with firing only as opposed to firing +silent."

    % there shouldn't be any nans at this point, but JK does (something
    % similar to) the following...
    bad = isnan(half_lambda);
    half_lambda(bad) = [];
    pkDirections(bad) = [];

    out.half_lambda = half_lambda;
    out.pkDirections = pkDirections;
end


function in = DealWithInputs(varargin)
    defaults.BIN_SIZE_CM = 2.5;
    defaults.FFT_SIZE = 256;
    defaults.MIN_FREQ = 0;
    defaults.MAX_FREQ = 10;
    defaults.PEAK_MAX_FRACTION = 0.1;
    defaults.posXY = [];
    defaults.spikePosInd = [];
    defaults.powerThreshold = 0;
    defaults.rate_map = [];
    
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