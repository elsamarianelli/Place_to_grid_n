function ret = eeg_instant_freq_power(varargin)
% See https://github.com/UCL/mTint/tree/master/FreqAnalysis for info
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk> &
%%                          <Daniel Manson> <d.manson@ucl.ac.uk>
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

    %1. Get a couple of values from eeg header
    [eegSampFreq,eegPerPos] = in.header.KeyValueMultiple('sample_rate','EEG_samples_per_position','num');

    %2. Mean normalise eeg
    eeg = double(in.eeg-mean(in.eeg));

    %3. Filter the EEG (bandpass on theta range)
    window = blackman(eegSampFreq+1); %Define window of interest - 1 second winowd here
    nqLim = eegSampFreq/2;
    EEGFilter = fir1(eegSampFreq,in.filterRange/nqLim,window); % Denominator - 125 is half the sampling rate - Nyquist frequency.
    filteredEEG = filtfilt(EEGFilter, 1, eeg); % The numerators are then the real actual frequencies you want to band pass.

    %4. Find complex-valued analytic function using Hilbert transform
    analyticEEG = hilbert(filteredEEG);
    clear eeg

    %5. Process analytic function to get phase, frequency and power
    %   (actually using true amp instead of power - power would just be that squared)
    phaseEEG = angle(analyticEEG); %Phase (radians) is the angle of the complex variable at each time point
    phaseEEG = unwrap(phaseEEG); %Adds 2pi if necessary to render plot of phase smooth
    ampEEG = abs(analyticEEG); %Modulus of analytic function i.e. instantaneous amplitude

    %6. Calculate instantaneous freq - note this vector is length of EEG-1 to get
    %  round this add a nan to the end. Will result in very limited timeslipage
    freqEEG = [diff(phaseEEG) ; nan]*eegSampFreq/(2*pi);

    %7. To compare against position need to down sample - usually eeg at 250Hz
    %   and pos at 50Hz, so eegPerPos is 5.
    ret.freqEEGposHz = nanmean(reshape(freqEEG, eegPerPos, []))';
    ret.ampEEGposHz = nanmean(reshape(ampEEG, eegPerPos, []))';
    ret.phaseEEGposHz = nanmean(reshape(phaseEEG, eegPerPos, []))';

    %8. Provide the return values
    ret.filteredEEG = filteredEEG;
    ret.phaseEEG = phaseEEG;
    ret.ampEEG = ampEEG;


    if in.PLOT_ON
        PlotForDebug(in.eeg,eegSampFreq,filteredEEG,analyticEEG,ampEEG,phaseEEG); 
    end

end


function PlotForDebug(eeg,eegSampFreq,filteredEEG,analyticEEG,ampEEG,phaseEEG)
    times = (1:numel(eeg))'/eegSampFreq; % recreating timestamps here is perhaps a bit confusing
    plot(times,eeg,'Color',[1 1 1]*0.8);
    xlim([0 times(end)]);
    xlabel('time (s)');
    ylabel('eeg (probably V)');
    hold all;
    plot(times,filteredEEG,'Color',[1 1 1]*0.5,'LineWidth',3);
    plot(times,real(analyticEEG),'r');
    plot(times,ampEEG,'b');
    plot(times,-ampEEG,'b');
    s = median(ampEEG)/50;
    phaseInd = round(mod(phaseEEG,2*pi)/(2*pi) * 100);
    phaseInd(phaseInd==0) = 100;
    phaseImg = ind2rgb(phaseInd',hsv(100));
    image([0 times(end)],[-s s],phaseImg); %This seems to be the best way to do a simple multicolor plot (scatter and patch are slow to render)
    hold off;
end

function in = DealWithInputs(varargin)
    defaults.eeg = [];
    defaults.filterRange = [7 11];
    defaults.header = [];
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
