classdef mTintHeader  
% Create an mTintHeader object from a file by finding the header text and then doing:
%  myHeader = mTintHeader(header_text_from_file);
% You can then use myHeader.KeyValue(key_name,[format]) to get a key value.
% To match a key name exactly, rather than just matching on the start of
% the key name, use myHeader.KeyValueExact(key_name,[format]).
% [format] is optional, by default a cell array of string is returned. Alternatively:
%     'num' convert key value string(s) to a double
%     'string' return key(s) as a single space-padded char array rather than a cell array
%     'int8', 'uint32', 'single', ... cast numeric output to the given Matlab type
%     'echo' displays the values on command line.
%
% When you create the header you can pass in a second input of the following
% form if you want to modify keys:
% myHeader = header(text, {'some_key_name1','hello','some_key_name2','good morning',...})
% (Matching is exact, if the key name is not found a new key is added.)
%
% To view the whole header do:
%   myHeader.ToString
%
% You can combine these two features to make a new header with modifications:
%   newHeader = mTintHeader(oldHeader.ToString,modifications)
%
% You could also potentially write out the modified header string to file.
% However it is adviseable to write to a new copy of the file rather than 
% overwritting the original.  (TODO: check newline behaviour).
%
% DM, Jan 2013
%
% TODO: since the header is immutable we could store it as a strucutre for
% use with exact queries ...need to keep the char array for inexact
% matches.
    properties (GetAccess = 'private', SetAccess = 'private')
        key_names = [];
        key_values = {};
        key_names_cellstr = {};
        key_names_width = 0;
    end
    
    methods
        %constructor
        function ob = mTintHeader(str,modifications)
            cell_array = textscan(str,'%s %[^\n\r]');
            ob.key_names = char(cell_array{1}{:});
            ob.key_values = cell_array{2};
            ob.key_names_cellstr = cell_array{1};
            ob.key_names_width = size(ob.key_names,2);
            
            if exist('modifications','var') && ~isempty(modifications)
                for ii=1:2:length(modifications)
                    ob = ob.ModifyKey(modifications{ii},modifications{ii+1},true);
                end
            end
        end
        
        % the five public methods
        function [values,keys] = KeyValue(ob,varargin)
            [values,keys] = ob.GetKeyValues(false,varargin{:});
        end
        function [values,keys] = KeyValueExact(ob,varargin)
            [values,keys] = ob.GetKeyValues(true,varargin{:});
        end    
        function varargout = KeyValueMultiple(ob,varargin) %exact, list of key names, and a format string (required)
            varargout = cell(nargout,1);
            for ii=1:numel(varargin)-1
                varargout{ii} = ob.GetKeyValues(true,varargin{ii},varargin{end});
            end
        end
        function ret = KeysPresent(ob,key_names) 
            %takes a cell array of exact key name strings and returns a logical array
            ret = false(size(key_names));
            for ii=1:numel(key_names)
                ret(ii) = ~isempty(ob.FindKeys(key_names{ii},true));
            end
        end
        function str = ToString(ob)
            data(2,:) = ob.key_values;
            data(1,:) = ob.key_names_cellstr;
            data(2,cellfun(@numel,data(2,:))==0) = {' '}; %must have a space if no value            
            str = sprintf('%s %s\r\n',data{:});
        end
        
    end
    
    %below are the private methods which are used by the public methods above
    methods(Access = 'private')
        function [values, keys] = GetKeyValues(ob,exact,target_name,format)
            values={};
            keys={};
            inds = ob.FindKeys(target_name,exact);
            if isempty(inds), return; end;
                
            keys = ob.key_names_cellstr(inds);
            values = ob.key_values(inds);
            if exist('format','var')
                values = ob.FormatKeys(values,format);
            end 
        end
        
        function ob = ModifyKey(ob,target_name,target_value,exact)
            target_name = char(target_name);
            target_value = char(target_value);
            target_width = size(target_name,2);
            
            if target_width > ob.key_names_width
                fprintf('Key name [%s] is longer than all existing key names.  Added a new key.\n',target_name);
                ob.key_names(end+1,1:target_width) = target_name;
                ob.key_names(1:end-1,ob.key_names_width+1:end) = ' '; %got 0-padding, wanted space-padding
                ob.key_names_cellstr(end+1) = {target_name};
                ob.key_values(end+1) = {target_value};
                ob.key_names_width = target_width;
                return;
            end
    
            inds = ob.FindKeys(target_name,exact);
           
           if isempty(inds)
                ob.key_names(end+1,1:target_width) = target_name;
                ob.key_names(end,target_width+1:end) = ' '; %got 0-padding, wanted space-padding
                ob.key_names_cellstr{end+1} = target_name;
                ob.key_values(end+1) = {target_value};
                fprintf('Key [%s] not found. Added a new key-value pair.\n',target_name);
                return;
            else
                ob.key_values(inds) = {target_value};
                fprintf('Modified value of existing key [%s].\n',target_name);
            end
           
        end
        
        function inds = FindKeys(ob,target,exact)
            inds = [];
            target = char(target);
            target_width = size(target,2);
            if target_width > ob.key_names_width, return; end;
            
            %cut all names down to the appropriate size to match target
            all_names = ob.key_names;
            if exact && target_width ~= ob.key_names_width 
                all_names = all_names(:,1:target_width+1);
                target(end+1) = ' '; %char pads with spaces.  We assume that keys do not themselves contain spaces.
            else % if exact && target_width = key_names_width, then don't need to change anything
                all_names = all_names(:,1:target_width);
            end
            
            %Find where the row is equal to the target string
            d = bsxfun(@eq,all_names,target);
            d = all(d,2);
            if exact
                inds = find(d,1,'first'); %if there are multiple exact matches, return only the first
            else
                inds = find(d);
            end
        end
    end
    
    % this method is static, meaning it's basically a normal function 
    % and doesn't need to know anything about the mTintHeader object.
     methods(Static = true,Access = 'private')
        function vals = FormatKeys(vals,format)
            switch format
            case 'num'
                vals = cellfun(@(v)sscanf(['0' v],'%f',1),vals); %we pad with zero, to ensure a match
            case 'string'
                if length(vals)>1
                    warning_dbstack(3,'more than one answer returned: strings may be padded with spaces');
                end
                vals=char(vals);
            case {'uint8','uint16','uint32','int8','int16','int32','single','double'}
                vals = cellfun(@(v)sscanf(['0' v],'%f',1),vals); %we pad with zero, to ensure a match
                vals = cast(vals,format);
            case 'echo'
                disp(char(vals));
            case 'default'
                %nothing to do
            otherwise
                error('unrecognized cast');
            end
        end
        

    end
    
end

function warning_dbstack(nStack,varargin)
    old_backtrace = warning('backtrace','off'); %we'll manually add the backtrace if it's wanted
    warning(varargin{:});
    if strcmp(old_backtrace.state,'on')
        dbstack(nStack+1);%+1 is for this function
        warning('backtrace','on'); % reset it back to what the user wanted
    end
end
