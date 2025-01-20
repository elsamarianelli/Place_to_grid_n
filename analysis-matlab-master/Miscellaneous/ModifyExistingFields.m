function structout = ModifyExistingFields(structin, varargin)
% ModifyExistingFields is based on Matlab's setstructfields.
%
% structout is the same as structin, but some field values have been
% modified.  Unlike in setstructfield no new fields are added.
%
% There are two ways to specify the new values, either the second input is 
% a structure, or the second and subsequent inputs are name-value pairs
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

structout = structin;

if nargin == 1 || (nargin == 2 && isempty(varargin{1})), return; end;

if ~isstruct(varargin{1})
    if ~ischar(varargin{1}) || mod(nargin-1,2)
        error(['ModifyExistingFields called with %d inputs, the second of which is not a struct. ' ...
            ' Either the secong input should be a structure or there should be 1+2n ' ...
            ' inputs specifying one starting structure and n field names & values to modify.'],nargin);
    end
    modifyFieldNames = varargin(1:2:end);
    modifyFieldValues = varargin(2:2:end);
else
    modifyStruct = varargin{1};
    modifyFieldNames = fieldnames(modifyStruct);
    modifyFieldValues = struct2cell(modifyStruct);
end

modifyFieldNameIsOk = isfield(structin,modifyFieldNames);

for ii = 1:length(modifyFieldNames)
    if modifyFieldNameIsOk(ii)
        value = modifyFieldValues{ii};
        name = modifyFieldNames{ii};
        if isstruct(value) && isstruct(structout.(name))
            structout.(name) = ModifyExistingFields(structout.(name), value);
        else
            structout.(name) = value;
        end
    end
end


end