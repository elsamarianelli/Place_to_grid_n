function [header, RawBinaryData] = fix_binary_data(filename,type,suffixForOld,duration) 
% Read header text and binary data from tint data files.
%
% If file is missing some header values and/or data end tokens it may be
% possible to correct the file.  The old file is saved with a suffix.
%
% Function [header, RawBinaryData] = fix_binary_data(filename,type,suffixForOld,duration)
%
% Currently supported types are 'pos' or 'tet' or 'set'.
% TODO: fix eeg files.
%
% the returned values should be the same as a call to read_binary_data 
% using the newly modified files:
%    [header] an mTintHeader object 
%    [RawBinaryData] a vector of bytes (i.e. agnostic of endian-ness)
%

    BATCH_SIZE = 1024;
    DATA_START_TOKEN = sprintf('\r\ndata_start');
    DATA_END_TOKEN = sprintf('\r\ndata_end');
    POS_NAN = 1023;
    
    switch type
        case {'set'}
            hasData = false;
        case {'tet' 'pos'}
            hasData = true;
            binaryDataClass = '*uint8';
            skip = [];
        case {'inp' 'eeg'}
            warning(['No code written to deal with "' type '" files yet.']);
        otherwise
            warning('unrecognized filetype');
            return;
    end

    needsWriting = false;
                
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
        fileLen = p;
        while isempty(dataLength) && p>0
            p = p-BATCH_SIZE;
            if p<0, p=0; end;
            fseek(f,p,'bof');
            dataLength = strfind(fread(f,BATCH_SIZE,'*uint8')',uint8(DATA_END_TOKEN));
        end
        if isempty(dataLength)
            fprintf('[%s]:\n\tcould not find data end token in file.\n',filename);
            needsWriting = true;
            dataLength = fileLen - headerLength - length(DATA_START_TOKEN) + 1;
        else
            dataLength = dataLength + p -headerLength  - length(DATA_START_TOKEN);
        end

        frewind(f); %back to start of file
        headerStr = fread(f,headerLength,'uint8=>char');
        fseek(f,length(DATA_START_TOKEN)-1,'cof'); %seek through length of token
        if isempty(skip)
            RawBinaryData = fread(f,dataLength,binaryDataClass);
        else
            RawBinaryData = fread(f,dataLength,binaryDataClass,skip);
        end

    else
        headerStr = fread(f,'uint8=>char');
        RawBinaryData = [];
    end

    fclose(f);

    header = mTintHeader(headerStr);

    if exist('duration','var') && duration ~= header.KeyValueExact('duration','num')
        header = mTintHeader(header.ToString,{'duration',num2str(duration)});
        needsWriting = true;
    end
    
    % At this point we have header and RawBinaryData as in read_binary_data function
    % we now inspecct the RawBinaryData to try and work out what certain
    % header values ought to be, if neccessary modifying the header object
    % In the case of pos file we may also modify the RawBinaryData.
    switch type
        case {'tet'}
            if mod(length(RawBinaryData),216) ~= 0
                RawBinaryData = RawBinaryData(1:floor(length(RawBinaryData)/216)*216);
            end
            RawBinaryData = reshape(RawBinaryData,216,[]);
            trueNSpikes = find(any(RawBinaryData,1),1,'last');
            headerNSpikes = header.KeyValueExact('num_spikes','num');
            if headerNSpikes ~= trueNSpikes 
                fprintf('[%s]:\n\tHeader says num_spikes=%d, but num_spikes seems to be %d.\n',filename,headerNSpikes,trueNSpikes);
                header = mTintHeader(header.ToString,{'num_spikes',num2str(trueNSpikes)});
                needsWriting = true;
            end
        case {'pos'}
            headerNPos = header.KeyValue('num_pos_samples','num');
            sampRate = header.KeyValueExact('sample_rate','num');
            
            num_data_pairs = header.KeyValue('num_colours','num');
            RawBinaryData = reshape(RawBinaryData,num_data_pairs*4 + 4,[]);
            trueNPos = find(any(RawBinaryData,1),1,'last');
            
            if exist('duration','var') && headerNPos ~= sampRate * duration
                fprintf('[%s]:\n\tHeader says num_pos_samples=%d, but num_pos_samples seems to be %d.\n\tHowever now forcing it to %d to match duration\n',filename,headerNPos,trueNPos,duration*sampRate);
                header = mTintHeader(header.ToString,{'num_pos_samples',num2str(duration*sampRate)});
                if headerNPos < sampRate*duration
                    RawBinaryData(1,duration*sampRate) = 0; %extend it 
                    % TODO: this works for ratemaps in Tint, but fails for pos possibly because it 
                    % needs timestamps and POS_NAN rather than simply 0s.
                else
                    RawBinaryData(:,duration*sampRate+1:end) = []; %trim it
                end
            elseif headerNPos ~= trueNPos
                fprintf('[%s]:\n\tHeader says num_pos_samples=%d, but num_pos_samples seems to be %d.\n',filename,headerNPos,trueNPos);
                header = mTintHeader(header.ToString,{'num_pos_samples',num2str(trueNPos)});
                needsWriting = true;
            end
    end

    % Now, if the header has been modified or if the data_end token was
    % missing we copy the old file and write a new one afresh using the
    % modified header and the (unmodified) RawBinaryData.
    if needsWriting
        if ~exist([filename suffixForOld],'file') %if the copied version already exists then do not overwirte it (as are probably reading an already modified version)
            copyfile(filename,[filename suffixForOld]);
        end
        delete(filename);
        
        f = fopen(filename,'w');
        fprintf(f,'%s',header.ToString);
        if hasData
            fprintf(f,DATA_START_TOKEN(3:end)); %lose initial \r\n
            fwrite(f,RawBinaryData(:));
            fprintf(f,DATA_END_TOKEN,'');
        end
        fclose(f);
    end

end

