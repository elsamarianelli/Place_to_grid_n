function ret = eeg_power_spectra(varargin)
% See https://github.com/UCL/mTint/tree/master/FreqAnalysis for info

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   
    
    % Select the relevant eeg samples, according to pos2use, and mean normalise

    if ~isempty(in.pos2use)
        %For eeg at 250Hz, pos values refer to every 5th eeg value, thus for every 
        %pos value we need to pick out 5 eeg values.  At higher eeg freq this may
        %be different.
        eeg = reshape(in.eeg,in.sampsPerPos,[]);
        eeg = eeg(:,in.pos2use);
        eeg = eeg(:);
    else
        eeg = in.eeg;
    end
    eeg = eeg - nanmean(eeg); %Mean normalise for given selection

    % Then send all the inputs to power_spectrum which does the real work...
    ret = power_spectrum(...
        eeg,in.padToPow2,1/in.sampFreq,in.thetaRange,in.maxFreq,...
        in.smthKernelWidth,in.smthKernelSigma,in.s2nWdth,in.PLOT_ON,in.ymax, in.xmax);
    
end


function in = DealWithInputs(varargin)
    defaults.pos2use = [];
    defaults.eeg = [];
    defaults.sampsPerPos = 5;
    defaults.thetaRange = [7 11];
    defaults.padToPow2 = NaN;
    defaults.sampFreq = 250;
    defaults.smthKernelWidth = 2;
    defaults.smthKernelSigma = 0.1875;
    defaults.s2nWdth = 2;
    defaults.maxFreq = 125;
    defaults.ymax = [];
    defaults.xmax = 25;
    defaults.PLOT_ON = true;
    
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
        in = modifyExistingFields(defaults,varargin{1});
    else
        in = modifyExistingFields(defaults,varargin{:});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end