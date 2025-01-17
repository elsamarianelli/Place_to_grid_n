function [trial] = loadSargolini(in,loadOptions)
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

persistent trial_cached filename_cached
% Reload trial data if trial has changed
if ~isempty(filename_cached) && strcmp(in.filename,filename_cached)
    trial = trial_cached;    
else
    [trial] = ppTrial(in.path,in.filename,loadOptions.EEGFilter);
    trial_cached = trial;
    filename_cached = in.filename;
end

% Load Cell timestamps
trial.type = in.type;
spikeTS = load([in.path in.filename '_T' num2str(in.tetrode),...
    'C' num2str(in.cluster) '.mat'], '-mat', 'cellTS');
trial.spikeTS = spikeTS.cellTS;

%Remove spikes that occurred outside the trial (necessary for concatenated
%moser trials)
trial.spikeTS(trial.spikeTS > (size(trial.xy,1)/trial.psr)) = []; 

end

function [trial] = ppTrial(path,filename,EEGFilter)

% Load trial and cell data
% Load pos
data = load([path filename '_pos.mat'],'-mat','posx', 'posy', 'post');

% Correct for weird trial lengths in moser data.
T = [600 1200]; % Trials can be either 600 or 1200 seconds
[ignore, b] = min(abs(max(data.post)-T));
T = T(b);

% Define position sample rate from pos times
psr = round(1/mean(diff(data.post)));

% Number of total pos samples
npos = T*psr;

% Path
xy = [data.posx(1:npos) data.posy(1:npos)];
% Input data must be non nan.
if any(isnan(xy))
    all_pos = 1:size(xy,1);
    good_pos = find(~isnan(xy(:,1)) & ~isnan(xy(:,2)));
    xy = interp1(good_pos,xy(good_pos,:),all_pos,'linear', 'extrap');
end

% Load EEG data
data = load([path filename '_eeg.mat'],'-mat','EEG', 'Fs');
esr = data.Fs;

% Account for moser concatentating pos and spike data but not eeg.
if numel(data.EEG) < (esr/psr)*size(xy,1)
    xy = xy(1:((psr/esr)*numel(data.EEG)),:);
    eeg = data.EEG(1:((esr/psr)*size(xy,1)));
else
    eeg = data.EEG(1:(T*esr));
end
    
% Calculate speed from path
spd = diff(xy);
spd = sqrt(spd(:,1).^2 + spd(:,2).^2);
spd = [spd(:); 0]*psr;

% Calculate direction *from path*
dir = mod( atan2(diff(xy(:,2)), diff(xy(:,1)))*(180/pi), 360);
dir = [dir(:);dir(end)];

% filter EEG and construct an analytic function using the Hilbert transform.
filteredEEG = filtfilt(EEGFilter, 1, eeg);
analyticEEG = hilbert(filteredEEG);

% Instantaneous phase is angle of analytic function
EEGPhase = angle(analyticEEG);

% Create output
trial = struct('xy',xy,'dir',dir,'spd',spd,'eeg',eeg,'phase',EEGPhase,...
    'filteredEEG', filteredEEG,'psr',psr, 'esr',esr);

end