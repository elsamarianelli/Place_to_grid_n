function [trial] = loadPlace(in,loadOptions)
%% Function to load DACQ data
% Structure 'in' contains path and filename of data structure and tetrode
% and cluster number of cell, whose timestamps are to be loaded.
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

% persistent trial_cached filename_cached spikeData_cached
% % Reload trial data if trial has changed
% if ~isempty(filename_cached) && strcmp(in.filename,filename_cached)
%     trial = trial_cached;
%     spikeData = spikeData_cached;
% else
    [trial, spikeData] = ppTrialPlace(in.path,in.filename,loadOptions.EEGFilter,in.eegchannel);
%     spikeData_cached = spikeData;
%     trial_cached = trial;
%     filename_cached = in.filename;
% end

% Load Cell timestamps
trial.type = in.type;
trial.spikeTS = spikeData(in.tetrode).timestamp(spikeData(in.tetrode).cut == in.cluster);
trial.spikeTS(trial.spikeTS > (size(trial.xy,1)/trial.psr)) = []; 
end

function [trial, spikeData] = ppTrialPlace(path,filename,EEGFilter,ch)

% Setup loading presets
load_in = load_trial();
load_in.GET_TRODES = 1:8; %load all available tetrodes and cuts.
load_in.GET_EEG = 1;
load_in.data_path = path;
load_in.flnmroot = filename;
load_in.POS_PROCESS_BOX_CAR = 0.4;
T = [15000 30000 45000 60000]; %Trial length in bins
max_xy = 60;

% Load data
data = load_trial(load_in);

% Set values of trial structure to those needed for phase precession.
% Get sample rates of pos and eeg
trial.psr = data.posData.header.KeyValueExact('sample_rate','num');
trial.esr = data.eegData(ch).header.KeyValueExact('sample_rate','num');

% Get path and correct weird errors with trial length
trial.xy = data.posData.xy;
[ignore, TT] = min(abs(T-size(trial.xy,1)));
T = T(TT);
trial.xy = trial.xy(1:T,:);
trial.xy = bsxfun(@minus,trial.xy,min(trial.xy));
trial.xy(trial.xy>max_xy) = max_xy;

% Get direction and speed
trial.dir = data.posData.dir_disp(1:T);
trial.spd = data.posData.speed(1:T);

% Get EEG and correct weird trial length error.
trial.eeg = data.eegData(ch).eeg_and_time(1:((trial.esr/trial.psr)*T),1);
trial.filteredEEG = filtfilt(EEGFilter, 1, trial.eeg);
trial.phase = angle(hilbert(trial.filteredEEG));

spikeData = data.spikeData;
end