function [ret] = thetaMod(varargin)
%% Calculates theta modulation properties of cells and EEG
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
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

%% Boilerplate for standardised input/output behaviour
if nargin == 0, ret = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

%% Compute theta modulation of the the LFP and the cell spiking
sw = in.spd>in.minSpeed;
sw = bwlabel(sw);
co = hist(sw,0:max(sw));
co = co(2:end);
co = find(co>in.minDuration*in.psr);
sw = ismember(sw,co);

in.eeg_in.pos2use = find(sw);
in.eeg_in.eeg = in.eeg;
in.eeg_in.sampFreq = in.esr;
in.eeg_in.sampsPerPos = round(in.esr/in.psr);
if in.eeg_in.PLOT_ON
    figure(1)
    subplot(4,3,2);
    set(gca, 'FontSize', 8)
    axis square
    hold on
end
eeg_out = eeg_power_spectra(in.eeg_in);
if in.eeg_in.PLOT_ON
    title(sprintf('EEG power spectrum | %s', num2str(eeg_out.s2n,2)), 'FontSize', 8)
end
ret.LFPThetaModulation = eeg_out.s2n;
ret.thetaFreq = eeg_out.maxFreq;

% And of the cell
in.ac_in.spikeTimes = in.spikeTS;
in.ac_in.posMask = sw;
in.ac_in.posSampFreq = in.psr;

if in.ac_in.PLOT_ON
    figure(1)
    subplot(4,3,3);
    set(gca, 'FontSize', 8)
    axis square
    hold on
end
ac_out = intrinsic_freq_autoCorr(in.ac_in);
if in.ac_in.PLOT_ON
    title(sprintf('spike train power spectrum | %s', num2str(ac_out.s2n,2)), 'FontSize', 8)
end
ret.thetaModulation = ac_out.s2n;
ret.intrFreq = ac_out.maxFreq;

end

function in = DealWithInputs(varargin)
defaults.eeg = [];
defaults.spikeTS = [];
defaults.psr = [];
defaults.esr = [];
defaults.eeg_in = eeg_power_spectra();
defaults.ac_in = intrinsic_freq_autoCorr();

VERSION = 1;

% Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    in = defaults;
    return;
end
if ischar(varargin{1})
    if mod(nargin,2)
        error(['%s called with %d inputs, the first of which is a string. ' ...
            ' Either the first input should be a structure or there should be an even' ...
            ' number of inputs specifying structure field names and values.'],mfilename,nargin);
        
    end
    user_in = cell2struct(varargin(2:2:end)',varargin(1:2:end));
else
    user_in = varargin{1};
    if nargin > 1 && VERSION ~= varargin{2}
        error(['%s called with version number %g, but code is version %g.\n' ...
            'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
            'the function without a version number.'],mfilename,varargin{2},VERSION);
    end
end
in = setstructfields(defaults,user_in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end