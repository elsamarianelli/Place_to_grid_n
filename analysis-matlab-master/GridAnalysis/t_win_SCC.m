function ret = t_win_SCC(varargin)
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info
% or readme.md in this folder

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   
    
    % #######################################################################
    % [Stage 0] Get some number
    n_samps = length(in.posSamp);
    n_spks_A = length(in.spikeIndA);
    n_spks_B = length(in.spikeIndB);
    winSizeBins = min(in.winSizeSec*in.pos_sample_rate, n_samps); 
    downsample = ceil(in.pos_sample_rate/in.downsampFreq); % factor by which positions are downsampled.
    Pthresh = in.Pthresh/downsample; % take account of downsampling

    % If we have a known dirft of orientation, we should correct that now
    if isempty(in.drift_orientation)
        posSamp = in.posSamp;
    else
        posSamp = cumrot(in.posSamp, in.drift_orientation);
    end
    
    % #######################################################################
    % [Stage 1] Prepare for main loop

    %1a. For each A spikeInd, calculate first and last B spikes in the window 
    %    and thus the number of B spikes.
    [~,startOfWindowBSpkInd] = histc(in.spikeIndA - winSizeBins, in.spikeIndB);
    startOfWindowBSpkInd = startOfWindowBSpkInd + 1;
    [~,endOfWindowBSpkInd] = histc(in.spikeIndA + winSizeBins,in.spikeIndB);
    startOfWindowBSpkInd(~startOfWindowBSpkInd) = 1;
    endOfWindowBSpkInd(~endOfWindowBSpkInd) = n_spks_B;
    nBSpikesInWin = endOfWindowBSpkInd - startOfWindowBSpkInd + 1;
    
    %1b. Work out offset inidices to be used when storing spike data
    off_spike = cumsum([0 ; nBSpikesInWin]);

    %1c. For each spike A, "clip" the start/end window index to stay within
    %    the 1:n_samps range of pos. And work out number of downsampled pos 
    %    bins in window and offset indicies for storing data.
    startOfWindowPosInd = max(in.spikeIndA -winSizeBins, 1);
    endOfWindowPosInd = min(in.spikeIndA + winSizeBins, n_samps);
    nPosInWindow =  endOfWindowPosInd - startOfWindowPosInd + 1;
    nDownsampInWin = double(floor((nPosInWindow-1)/downsample)+1); % TODO: check this !!!
    off_dwell = double(cumsum([0 ; nDownsampInWin])); 
    % Note dwell indices cannot be "single" because we run out of integer
    % -level precision when using large windows and lots of spikes.
    
    %PlotForDebug2(startOfWindowPosInd,endOfWindowPosInd,in.spikeIndA,in.spikeIndB,...
    %                    startOfWindowBSpkInd,endOfWindowBSpkInd)

    %1d. Pre-allocate dwell and spike arrays, singles for speed
    dwell = nan(off_dwell(end), 2,'single');
    spike = nan(off_spike(end), 2,'single');

    %1e. For speed, we convert posSamp to singles
    posSamp = single(posSamp);

    %1f. To reuce indexing inside the loop we get the posSamples for each
    %    spike in A and B.
    spikeAPosSamp = posSamp(in.spikeIndA,:);
    spikeBPosSamp = posSamp(in.spikeIndB,:);
    
    % #######################################################################
    % [Stage 2] Get displacement data for window starting at each spikeInd
    for ii = 1:n_spks_A
        iiPos = spikeAPosSamp(ii,:);

        %2a. calculate dwell displacements
        winInd_dwell = startOfWindowPosInd(ii) : downsample : endOfWindowPosInd(ii);
        dwell((1:nDownsampInWin(ii)) + off_dwell(ii),:) = bsxfun(@minus,posSamp(winInd_dwell,:),iiPos);

        %2b. calculate spike displacements
        winInd_spike = startOfWindowBSpkInd(ii): endOfWindowBSpkInd(ii);
        spike((1:nBSpikesInWin(ii)) + off_spike(ii),:) =  bsxfun(@minus,spikeBPosSamp(winInd_spike,:),iiPos);
    end


    if in.as_1d
        % #######################################################################
        % [Stage 3] Bin displacement data
        
        %3a. Convert 2d to 1d
        dwell = sqrt(sum(dwell.^2,2));
        spike = sqrt(sum(spike.^2,2));
        
        %3b. Calculate limits and bin sizes
        dwell_abs_max = max(dwell);
        binsize = dwell_abs_max/in.nbins;
        
        %3c. Bin dispalcement data for dwell and spikes
        siz = [ceil(dwell_abs_max/binsize) 1];
        inds_dwell = ceil(toBinUnits(dwell, [0 ; dwell_abs_max], binsize));
        H_dwell = accumarray(inds_dwell, 1, siz);
        inds_spike = ceil(toBinUnits(spike, [0 ; dwell_abs_max], binsize));
        H_spike = accumarray(inds_spike, 1, siz);
        
        % #######################################################################
        % [Stage 4] Normalise H_spike using H_dwell

        %4a. Apply place field-like smoothing to H_dwell and H_spike
        if in.boxcarWidth>1         
            denom = smooth(double(H_dwell>0),in.boxcarWidth);
            denom( H_dwell==0 ) = NaN;
            fH_dwell = smooth(H_dwell,in.boxcarWidth);
            fH_dwell = fH_dwell./denom;
            fH_spike = smooth(H_spike,in.boxcarWidth);
            fH_spike = fH_spike./denom;
        else
            fH_dwell = H_dwell;
            fH_dwell(H_dwell==0) = NaN;
            fH_spike = H_spike;
        end
        
    else
        %3b. Calculate limits and bin sizes
        dwell_abs_max = max(abs(dwell));
        p_range = [-dwell_abs_max; dwell_abs_max];
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

        %4a. Apply place field-like smoothing to H_dwell and H_spike
        if in.boxcarWidth>1
            b = ones(in.boxcarWidth);
            denom = double(H_dwell>0);
            denom = filter2(b, denom);
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
    end
    

    %4b. Do division
    H = fH_spike./fH_dwell;

    %4c. Set unvisted to NaN
    H(H_dwell < Pthresh) = NaN;

    %4d. provide output
    ret.H = H;
    ret.Hs = H_spike;
    ret.Hp = H_dwell;


    if in.PLOT_ON, PlotForDebug(H,in); end;

end

function PlotForDebug(H,in)
    if in.as_1d
        plot(H);
        xlabel('distance from spike')
        ylabel('rate (?)')
    else
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
        title('time windowed SCC');
    end
end

function PlotForDebug2(startOfWindowPosInd,endOfWindowPosInd,spikeIndA,spikeIndB,...
                    startOfWindowBSpkInd,endOfWindowBSpkInd)    
    old_ax = gca();
    figure(99)
    cla;
    plot(startOfWindowPosInd,'k')
    hold on;
    p1 = plot(endOfWindowPosInd,'k');
    p2 = plot(spikeIndA,'r');
    p3 = plot(spikeIndB(startOfWindowBSpkInd),'b');
    plot(spikeIndB(endOfWindowBSpkInd),'b')
    legend([p1,p2,p3],'first/last pos sample','A spike pos','first/last B spike pos',...
        'Location','NorthWest')
    xlabel('window (i.e. spike A) #')
    ylabel('pos sample index')
    axes(old_ax)
end

function in = DealWithInputs(varargin)
    defaults.posSamp = [];
    defaults.spikeIndA = [];
    defaults.spikeIndB = [];
    defaults.winSizeSec = 10;
    defaults.pos_sample_rate = 50;
    defaults.nbins = 71;
    defaults.boxcarWidth = 5;
    defaults.Pthresh = 100;
    defaults.downsampFreq = 50;
    defaults.PLOT_ON = true;
    defaults.as_1d = false;
    defaults.CLIM_AT_ZERO = false;
    defaults.drift_orientation = [];
    
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
