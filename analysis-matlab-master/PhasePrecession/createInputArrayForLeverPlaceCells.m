function [inArray] = createInputArrayForLeverPlaceCells
% This function is associated with the construction of an input array
% necessary to run ppMain specifically for the Lever placec cell data set.
%
% This function must define the minimum number of inputs required to run
% the phase precession scripts.
% 
%   txtFile - full path and filename of location of text file containing
%       all the necessary input variables, arranged in columns, one for
%       each cell:
%
%   type - (REQUIRED by ppMain and functions therein) a variable indicating
%       'p' for place cell or 'g' for grids. Default behavior (if empty) is
%       for grid cells but the field must be present.
%
%   In addition any information needed for your loadfun to load the 
%   necessary trial/cell data (eg, filenames, paths, tetrode and cell
%   numbers. See txtFile.
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

txtFile = 'E:\Science\MetaData\Colin\Unique Place Cells.txt';

% First load relevant trial metadata
fid = fopen(txtFile,'r');
C = textscan(fid, '%s %d %s %d %d %s %d %d %d %d %s %s %s %s %s %d');
fclose(fid);
C = C';
for i = 1:size(C,1)
    if ~iscell(C{i})
        C(i,1:numel(C{i})) = num2cell(C{i});
    else
        C(i,1:numel(C{i})) = C{i};
    end
end
        
inArray = cell2struct(C,...
    {'lab', 'ratNo', 'date', 'day', 'room', 'env', 'session',...
    'trialNo', 'tetrode', 'cluster', 'type', 'region', 'filename',...
    'path', 'hemi', 'eegchannel'}, 1);

end