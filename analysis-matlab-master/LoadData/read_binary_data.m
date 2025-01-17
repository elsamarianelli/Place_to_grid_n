function [header, RawBinaryData] = read_binary_data(filename,type)
% Read header text and binary data from tint data files.
%
% Function [header, RawBinaryData] = read_binary_data(filename,type)
%
% Type can be 'set', 'pos', 'tet', 'eeg', 'tet_times'
%
% returns:
%    [header] an mTintHeader object
%    [RawBinaryData] a vector of bytes (i.e. agnostic of endian-ness)
%

BATCH_SIZE = 1024;
DATA_START_TOKEN = sprintf('\r\ndata_start');
DATA_END_TOKEN = sprintf('\r\ndata_end');

switch type
    case 'set'
        hasData = false;
    case {'pos'}
        hasData = true;
        binaryDataClass = '*uint8';
        skip = [];
    case 'tet_times'
        hasData = true;
        binaryDataClass = '4*int8=>int8';
        skip = 212;
    case {'eeg' 'tet' 'inp'}
        hasData = true;
        binaryDataClass = '*int8';
        skip = [];
    otherwise
        warning('unrecognized filetype');
        return;
end

f = fopen(filename,'r');
if f == -1
    warning('No such file:\n%s',filename);
    header = []; RawBinaryData = [];
    return;
end


if hasData
    %find start token
    headerLength = [];
    p = -BATCH_SIZE;
    while isempty(headerLength) && ~feof(f)
        p = p+BATCH_SIZE;
        headerLength = strfind(fread(f,BATCH_SIZE,'*uint8')',uint8(DATA_START_TOKEN));
    end
    headerLength = headerLength + p;
    if isempty(headerLength)
        fclose(f); error('could not find data start token in file:\n%s',filename);
    end
    
    %find end token
    dataLength = [];
    fseek(f,0,'eof');
    p = ftell(f);
    while isempty(dataLength) && p>0
        p = p-BATCH_SIZE;
        if p<0, p=0; end;
        fseek(f,p,'bof');
        dataLength = strfind(fread(f,BATCH_SIZE,'*uint8')',uint8(DATA_END_TOKEN));
    end
    if isempty(dataLength)
        fclose(f); error('could not find data end token in file:\n%s',filename);
    end
    dataLength = dataLength + p -headerLength  - length(DATA_START_TOKEN);
    
    frewind(f); %back to start of file
    headerStr = fread(f,headerLength,'uint8=>char');
    fseek(f,length(DATA_START_TOKEN)-1,'cof'); %seek through length of token
    if isempty(skip)
        RawBinaryData = fread(f,dataLength,binaryDataClass);
    else
        RawBinaryData = fread(f,dataLength,binaryDataClass,skip);
    end
    
else
    %whole file is a header
    headerStr = fread(f,'uint8=>char');
    RawBinaryData = [];
end

fclose(f);

header = mTintHeader(headerStr);

end


% ----------------------------------------------------------------------------------------------------------------------
