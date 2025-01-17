function [summary, out] = resultsSummary(varargin)
%% Produce a summary of results from ppMain, (categorising cells numbers)
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


% regressors = {'spks'; 'cfr'; 'ifr'; 'tim'; 'd'; 'pdmd'; 'pdcd';'tcn'};
% Example constants
% summarySet.SHALLOW_THRESHOLD = 90*(pi/180);
% summarySet.DIRECTIONALITY_THRESHOLD = 0.1;
% summarySet.GRIDNESS_THRESHOLD = 0;

%If resultsSummary is used inside ppMain then it will produce the
%concatenated array of structures with one element per cell for the whole
%data set. Otherwise it will take the concatenated array and produce a
%summary.
if nargin==4
    outPath = varargin{1};
    outFile = varargin{2};
    t = varargin{3};
    summarySet = varargin{4};
    x = load([outPath,outFile]);
    for n=(1+t-numel(fieldnames(x))):t;
        out(n) = x.(['out',num2str(n)]);
    end
    clear x
elseif nargin==2
    out = varargin{1};
    summarySet = varargin{2};
    t = numel(out);
end

%Get values on which to make selections
gridness = vertcat(out(:).gridness);
gridness = repmat(gridness,[1,8]);
%If gridness is nan - it wasn't a grid cell
gridness(isnan(gridness)) = 1;
directionality = vertcat(out(:).directionality);
directionality = repmat(directionality,[1,8]);
regress = vertcat(out(:).regress);
slope = vertcat(regress(:).slope);
slope = reshape(slope,[8,t])';
p_value = vertcat(regress(:).p_shuffled);
p_value = reshape(p_value,[8,t])';

% produce vectors of viable cells in each category, for each regressor.
directional = directionality>=summarySet.DIRECTIONALITY_THRESHOLD;
gridlike = gridness>=summarySet.GRIDNESS_THRESHOLD;
steep = p_value<=0.05 & (slope<=-summarySet.SHALLOW_THRESHOLD | slope>=summarySet.SHALLOW_THRESHOLD);
shallow = p_value<=0.05 & abs(slope)<summarySet.SHALLOW_THRESHOLD;

% Produce some summary values. These can be played for user needs.
summary.Ngridlike = sum(gridlike(:,1));
summary.Ndirectional = sum(gridlike(:,1) & directional(:,1));
summary.NIND = (gridlike & slope<0 & steep);
summary.PIND = (gridlike & slope>0 & steep);
summary.NWeakIND = (gridlike & slope<0 & shallow);
summary.PWeakIND = (gridlike & slope>0 & shallow);
summary.Uncategorised = (gridlike & ~directional & p_value>0.05);
summary.significant = p_value<=0.05;

end
% summary.Ngridlike = sum(gridlike);
% summary.Ndirectional = sum(gridlike & directional);
% summary.NIND = find(gridlike & ~directional & steep);
% summary.NDirIND = find(gridlike & directional & steep);
% summary.NWeakIND = find(gridlike & ~directional & shallow);
% summary.NDirWeakIND = find(gridlike & directional & shallow);
% summary.Uncategorised = find(gridlike & ~directional & ~p_value<=0.05);
% summary.DirUncategorised = find(gridlike & directional & ~p_value<=0.05);