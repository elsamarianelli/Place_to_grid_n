function [out, summary] = ppMain(inArray, SETTINGS, loadfun)
%% Help
% inArray is an array of structures passed to @loadfun. Creation of inArray
% and design of @loadfun is left to the user's discretion, (see example below).
% @loadfun should accept inArray(i) and loadfun must output a structure
% 'trial' with following fields:
%
% xy            - [Nx2] double: xy co-ordinates of rat path
% dir           - [Nx1] double: directions in degrees of rat's movement-
%                   based (one-spot) direction.
% spd           - [Nx1] double: rat's speed.
% eeg           - [(esr/psr)*Nx1] double: raw EEG vector
% phase         - [(esr/psr)*Nx1] double: phase of filtered EEG
% filteredEEG   - [(esr/psr)*Nx1] double: filtered EEG
% psr           - [1x1] double: position sample rate
% esr           - [1x1] double: eeg sampling rate
% type          - [1x1] char: letter defining cell type eg ('g' for grid
%                   cell or 'p' for place)
% spikeTS       - [Sx1] double: times at which spikes in cell fired
%
% SETTINGS      - structure created by ppSettings.m containing all the
%                   constants required by ppMain and all functions called
%                   therein. Called as ppSettings(outFile)
%
% outFile       - filename for output variables made unique by
%                   concatenation with datestr
%
% EG for Sargolini 2006 data available from Moser website the function :
%
% loadfun = @loadSargolini;
% inArray = createInputArrayForSargoliniDataSet;% is used to create inArray
% outFile = 'SargoliniGridPP';
% [out, summary] = ppMain(createInputArrayForSargoliniDataSet,...
%                           ppSettings(outFile), @loadSargolini);
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

%% Main Loop
for t=1:numel(inArray)  % For each trial/cell in the list
    disp(t)
    %try
    %% Load trial and cell data
    [trial] = loadfun(inArray(t),SETTINGS.PPLOAD);
    
    %% Now get grid properties (if grid cell)
    if strcmp(trial.type,'g')
        SETTINGS.GCPROPS = modifyExistingFields(SETTINGS.GCPROPS, trial);
        [gcp] = gcProps(SETTINGS.GCPROPS);
    else
        gcp.gridness = NaN;
        gcp.orientation = NaN;
        gcp.scale = NaN;
    end
    
    %% Calculate LFP and spiking theta mnodulation
    SETTINGS.THETAMOD = modifyExistingFields(SETTINGS.THETAMOD, trial);
    eegSpkMod = thetaMod(SETTINGS.THETAMOD);
    
    %% Calculate directionality and skaggs info
    SETTINGS.DIRECTIONALITY = modifyExistingFields(SETTINGS.DIRECTIONALITY, trial);
    kl = directionality(SETTINGS.DIRECTIONALITY);
    
    if (~strcmp(trial.type, 'g') || gcp.gridness > SETTINGS.gridnessThreshold) &&...
            eegSpkMod.LFPThetaModulation > SETTINGS.LFPThetaModulationThreshold &&...
            eegSpkMod.thetaModulation > SETTINGS.intrThetaModulationThreshold
        
        %% Get theta properties for the cell/eeg
        SETTINGS.GETTHETAPROPS = modifyExistingFields(SETTINGS.GETTHETAPROPS, trial);
        theta = getThetaProps(SETTINGS.GETTHETAPROPS);
        
        %% Partition the fields in the ratemap
        SETTINGS.PARTITIONFIELDS = modifyExistingFields(SETTINGS.PARTITIONFIELDS, trial);
        fields = partitionFields(SETTINGS.PARTITIONFIELDS);
        
        %% Use partitioned fields to get the various kinds of transformed
        %  coordinates and some of the per-run info
        SETTINGS.GETPOSPROPS = modifyExistingFields(SETTINGS.GETPOSPROPS, trial);
        SETTINGS.GETPOSPROPS = modifyExistingFields(SETTINGS.GETPOSPROPS, fields);
        gpp_out = getPosProps(SETTINGS.GETPOSPROPS);
        
        %% Use theta pos properties to get properties for each spike and
        %  more per-run info
        getSpikeProps_in = struct('pos', gpp_out.pos, 'theta',theta,...
            'runs',gpp_out.runs, 'spikeTS',trial.spikeTS,...
            'posSampleRate',trial.psr, 'eegSampleRate',trial.esr);
        [gsp_out] = getSpikeProps(getSpikeProps_in);
        
        %% Define values/vectors on which to split the data.
        [runs] = runSplit(gsp_out.runs);
        
        %% Perform circular-linear regression, correlation, calculate p-values
        %  and confidence intervals.
        pos = gpp_out.pos;
        spk = gsp_out.spk;
        [regress] = ppRegression(theta,pos,spk,runs,SETTINGS.PPREGRESSION);
        
        %% Save variables in a struct to a mat file
        createOutputStructure(inArray(t), trial, theta, pos, spk, runs,...
            regress, gcp.gridness, gcp.orientation, gcp.scale, kl.divergence, eegSpkMod.thetaModulation,...
            eegSpkMod.LFPThetaModulation, kl.converged, t, SETTINGS.outPath, SETTINGS.outFile);
    else
        createOutputStructure(inArray(t), trial, NaN, NaN, NaN, NaN, NaN,...
            gcp.gridness, gcp.orientation, gcp.scale, kl.divergence, eegSpkMod.thetaModulation,...
            eegSpkMod.LFPThetaModulation, kl.converged, t, SETTINGS.outPath, SETTINGS.outFile);
    end
    
    % Save plots and close to stop them accumulating and hogging memory
    if SETTINGS.DIRECTIONALITY.PLOT_ON || SETTINGS.PARTITIONFIELDS.PLOT_ON || ...
            SETTINGS.GETPOSPROPS.PLOT_ON || SETTINGS.PPREGRESSION.PLOT_ON
        
        savePlots(SETTINGS,1,t)
        close(1)
    end
    
    %As it stands, need to rename output structure manually later.
    %catch err
    %    warning(['skipping main loop iteration ' num2str(t) ': ' getReport(err,'extended')]);
    %end
    
end

%% Create final output variable with data from all cells and a summary
[summary, out] = resultsSummary(SETTINGS.outPath,SETTINGS.outFile,t,SETTINGS.RESULTSSUMMARY);
save([SETTINGS.outPath,SETTINGS.outFile], 'out', 'summary')

%% Produce a PDF of all cells if plotted
if (SETTINGS.DIRECTIONALITY.PLOT_ON || SETTINGS.PARTITIONFIELDS.PLOT_ON || ...
        SETTINGS.GETPOSPROPS.PLOT_ON || SETTINGS.PPREGRESSION.PLOT_ON) && ...
        SETTINGS.outputCollatedPDF
    
    system(['"D:\Program Files (x86)\Adobe\Acrobat 9.0\Acrobat\AcroDist.exe"',...
        ' /F /J "D:\Program Files (x86)\Adobe\Acrobat 9.0\Acrobat\Settings\',...
        'Smallest File Size.joboptions" /Q /N ', '"', [SETTINGS.outPath SETTINGS.outFile], '.ps"']);
end

end