function regress_pdcd = pp2DFigures(out)
%% Just points to cells with the some chosen characteristics and plots 2D
%% phase plots (UNFINISHED)
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
%                           
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

gridness = [out(:).gridness]';
thetaModulation = [out(:).thetaModulation]';
lfpThetaModulation = [out(:).LFPThetaModulation]';
directionality = horzcat(out(:).directionality)';
region = {out(:).region}';

regress = [out(:).regress];
slope = vertcat(regress(:).slope);
slope = reshape(slope,[8,numel(regress)])';
shuffledPVal = vertcat(regress(:).p_shuffled);
shuffledPVal = reshape(shuffledPVal,[8,numel(regress)])';

gridness = repmat(gridness(:),[1,8]);
thetaModulation = repmat(thetaModulation(:),[1,8]);
lfpThetaModulation = repmat(lfpThetaModulation(:),[1,8]);
directionality = repmat(directionality(:),[1,8]);
region = repmat(region(:),[1,8]);

%% Grid Cells
% highThetaModNonDirGrids = gridness>=0 & thetaModulation>=5 &...
%     lfpThetaModulation>=10 & directionality<=0.1;
highThetaModNonDirSigNegGrids = slope<0 & shuffledPVal<=0.025 & gridness>=0 &...
    thetaModulation>=5 & lfpThetaModulation>=10 & directionality<=0.1;

% Numbers and indices into layer II grid cells
% layerII = highThetaModNonDirGrids & strcmpi(region,'MECLII');
layer_II_sig_neg = highThetaModNonDirSigNegGrids & strcmpi(region,'MECLII');

% Numbers and indices into layer III grid cells
% layerIII = highThetaModNonDirGrids & strcmpi(region,'MECLIII');
layer_III_sig_neg = highThetaModNonDirSigNegGrids & strcmpi(region,'MECLIII');

% Numbers and indices into layer V grid cells
% layerV = highThetaModNonDirGrids & strcmpi(region,'MECLV');
layer_V_sig_neg = highThetaModNonDirSigNegGrids & strcmpi(region,'MECLV');

% Numbers and indices into layer VI grid cells
% layerVI = highThetaModNonDirGrids & strcmpi(region,'MECLVI');
layer_VI_sig_neg = highThetaModNonDirSigNegGrids & strcmpi(region,'MECLVI');

%% Plot 2D and linearised phase precession for grouped data
% LAYER II
in.inStruct = out(highThetaModNonDirSigNegGrids(:,7));
regress_pdcd(1) = pp2DPhase(in);

% LAYER III
in.inStruct = out(layer_III_sig_neg(:,7));
regress_pdcd(2) = pp2DPhase2(in);

% LAYER V
in.inStruct = out(layer_V_sig_neg(:,7));
regress_pdcd(3) = pp2DPhase2(in);

% LAYER VI
in.inStruct = out(layer_VI_sig_neg(:,7));
regress_pdcd(4) = pp2DPhase2(in);