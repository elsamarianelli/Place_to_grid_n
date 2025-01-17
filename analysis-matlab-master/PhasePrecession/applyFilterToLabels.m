function varargout = applyFilterToLabels( M, varargin )
% [L_1new, L_2new,...] = applyFilterToLabels( M, L_1, L_2, ... )
%
% M is a logical mask specifying which label numbers to keep
% L_i is an array of positive integer labels.
%
% The function sets the undesired labels to 0 and renumbers the remaining 
% labels 1 to n, where n is the number of trues in M.
%
% e.g. L = [0 0 0 2 2 1 1 0 0 3 3 3 0 0 4 4 4 0 0 5 5 5 3]
% and M = [0 1 0 1 1];
% gives:   [0 0 0 1 1 0 0 0 0 0 0 0 0 0 2 2 2 0 0 3 3 3 0]
%
% If multiple L arrays are provided the relabeling is applied to each of
% them.
%
%%    Copyright (C) <2013>  <Daniel Manson> <d.manson@ucl.ac.uk>
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

newVals = M .* cumsum(M);

varargout = cell(size(varargin));

for ii=1:numel(varargin)
    L = varargin{ii};
    L(L>0) = newVals(L(L>0));
    varargout{ii} = L;
end

end