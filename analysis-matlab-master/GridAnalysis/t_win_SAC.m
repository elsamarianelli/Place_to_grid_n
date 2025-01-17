function ret = t_win_SAC(varargin)
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   
    
    % #######################################################################
    % [Stage 0] Get some number
    n_samps = length(in.posSamp);
    n_spks = length(in.spikeInd);
    winSizeBins = min(in.winSizeSec*in.pos_sample_rate, n_samps); 
    downsample = ceil(in.pos_sample_rate/in.downsampFreq); % factor by which positions are downsampled.
    Pthresh = in.Pthresh/downsample; % take account of downsampling


    % #######################################################################
    % [Stage 1] Prepare for main loop

    %1a. Calculate number of spikes in the window for each spikeInd (ignoring spike itself)
    [~,endOfWindowSpkInd] = histc(in.spikeInd + winSizeBins,in.spikeInd);
    endOfWindowSpkInd(~endOfWindowSpkInd) = n_spks;
    nSpikesInWin = endOfWindowSpkInd - (1:n_spks)';
    
    %1b. Work out offset inidices to be used when storing spike data
    off_spike = cumsum([0 ; nSpikesInWin]);

    %1c. Work out number of downsampled pos bins in window and offset indicies for storing data
    nPosInWindow = min(winSizeBins, n_samps - in.spikeInd);
    nDownsampInWin = floor((nPosInWindow-1)/downsample)+1;
    off_dwell = cumsum([0 ; nDownsampInWin]);

    %1d. Pre-allocate dwell and spike arrays, singles for speed
    dwell = nan(off_dwell(end), 2,'single');
    spike = nan(off_spike(end), 2,'single');

    %1e. For speed, we convert posSamp to singles
    posSamp = single(in.posSamp);

    %1f. To reduce indexing inside the loop we get the posSamples for each spikeInd 
    spikePosSamp = posSamp(in.spikeInd,:);

    % #######################################################################
    % [Stage 2] Get displacement data for window starting at each spikeInd
    for ii = 1:n_spks;
        iiPos = spikePosSamp(ii,:);

        %2a. calculate dwell displacements
        winInd_dwell = in.spikeInd(ii)+1 : downsample : in.spikeInd(ii)+nPosInWindow(ii);
        dwell((1:nDownsampInWin(ii)) + off_dwell(ii),:) = bsxfun(@minus,posSamp(winInd_dwell,:),iiPos);

        %2b. calculate spike displacements
        winInd_spike = ii+1 : ii+nSpikesInWin(ii);
        spike((1:nSpikesInWin(ii)) + off_spike(ii),:) =  bsxfun(@minus,spikePosSamp(winInd_spike,:),iiPos);
    end

    % #######################################################################
    % [Stage 3] Bin displacement data in half-plane

    %3a. Any vectors with negative first coordinate should be rotated 180 degrees
    dwell = bsxfun(@times,dwell,sign(dwell(:,1)));
    spike = bsxfun(@times,spike,sign(spike(:,1)));

    %3b. Calculate limits and bin sizes
    dwell_abs_max = max(abs(dwell));
    p_range = [0 -dwell_abs_max(:,2); dwell_abs_max]; %we know the first cordiante is >=0
    binsize = (p_range(2,:)-p_range(1,:))/in.nbins;
    binsize = max(binsize); % Keep same binsize for x and y

    %3c. Bin dispalcement data for dwell and spikes
    siz = ceil(diff(p_range,[],1)/binsize);
    inds_dwell = ceil(toBinUnits(dwell, p_range, [binsize binsize]));
    H_dwell = accumarray(inds_dwell, 1, siz);
    inds_spike = ceil(toBinUnits(spike, p_range, [binsize binsize]));
    H_spike = accumarray(inds_spike, 1, siz);

    %3d. Swap (y,x) for (x,y)
    H_dwell = permute(H_dwell, [2 1]);
    H_spike = permute(H_spike, [2 1]);


    % #######################################################################
    % [Stage 4] Normalise H_spike using H_dwell

    %4a. Append rotated versions so that two half-planes become a symmetric full plane
    H_dwell = [rot90(H_dwell,2) H_dwell];
    H_spike = [rot90(H_spike,2) H_spike];

    %4b. Apply place field-like smoothing to H_dwell and H_spike
    if in.boxcarWidth>1
        b = ones(in.boxcarWidth);
        c = double(H_dwell>0);
        denom = filter2(b, c);
        denom( denom==0 ) = NaN;
        fH_dwell = filter2(b, H_dwell);
        fH_dwell = fH_dwell./denom;
        fH_spike = filter2(b, H_spike);
        fH_spike = fH_spike./denom;
    else
        fH_dwell = H_dwell;
        fH_dwell(H_dwell==0) = NaN;
        fH_spike = H_spike;
    end

    %4c. Do division
    H = fH_spike./fH_dwell;

    %4d. Set unvisted to NaN
    H(H_dwell < Pthresh) = NaN;

    %4e. provide output
    ret.H = H;
    ret.Hs = H_spike;
    ret.Hp = H_dwell;


    if in.PLOT_ON, PlotForDebug(H,in); end;

end

function PlotForDebug(H,in)
    if in.CLIM_AT_ZERO 
        im = max(H,0); %clip lower limit to zero
    else
        im = H-min(min(H)); % shift lower limit to become zero
    end
    im = im / max(im(:));
    cm = jet(60);
    im = round(im*60.999 + 0.5);
    bad = repmat(isnan(im),[1 1 3]);
    im = ind2rgb(im,cm);
    im(bad) = 1;
    image(im);
    axis image
    xlabel('bins'); ylabel('bins');
    title('time windowed SAC');
end

function in = DealWithInputs(varargin)
    defaults.posSamp = [];
    defaults.spikeInd = [];
    defaults.winSizeSec = 10;
    defaults.pos_sample_rate = 50;
    defaults.nbins = 71;
    defaults.boxcarWidth = 5;
    defaults.Pthresh = 100;
    defaults.downsampFreq = 50;
    defaults.PLOT_ON = true;
    defaults.CLIM_AT_ZERO = false;
    
    VERSION = 1.0;
    
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
