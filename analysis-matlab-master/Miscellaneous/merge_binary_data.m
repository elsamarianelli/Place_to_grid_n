function [outHeader, outRawBinaryData] = merge_binary_data(destFilename,srcFilenames,type,suffixForOld) 
% Read header text and binary data from tint data files.
%
% Function [header, RawBinaryData] = merge_binary_data(destFilename,srcFilenames,type,suffixForOld)
%
% srcFilenames should be a cell array of two or more strings 
% e.g. {'C:\my trial part a.pos', 'C:\my trial part b.pos'}.  
%
% Type should be one of: 'tet' 'pos' 'set' 'cut' 'eeg'.
% Note that code may not yet have been written for all of the above, if so
% a warning will be issued when you try and use it.
%
% suffixForOld should be a string such as '.old' it is used if the
% desFilename already exists - the old file will be renamed with this
% suffix and the new file will be put in its place.
%
% the returned values should be the same as a call to read_binary_data 
% using the newly merged files:
%    [header] an mTintHeader object 
%    [RawBinaryData] a vector of bytes (i.e. agnostic of endian-ness)
%
% Note: merged pos data will contain a jump at the join, which will 
% end up getting smoothed out when read in using matlab/tint code.  Unless
% you are merging together lots of very short trials this is unlikely to be
% a real problem.  Here we add some extra fields to the pos header so that
% in principle you could customise code to not do the smoothing.
%
% For convience this function can merge cut files even though they aren't
% actually binary data in the same way as other files
%
% USE AT YOUR OWN RISK!!!
%
% WARNING: this seems to work more or less, but be careful with it.
%
% DM, Jun 2014.
%
    DATA_START_TOKEN = sprintf('\r\ndata_start');
    DATA_END_TOKEN = sprintf('\r\ndata_end');
    [~,~,endian] = computer;  % 'L' or 'B'
    
    nFiles = numel(srcFilenames);
    headers = cell(nFiles,1);
    binaryData = cell(nFiles,1);
    for ii=1:nFiles
        [headers{ii},binaryData{ii}] = read_binary_data(srcFilenames{ii},type);
        binaryData{ii} = typecast(binaryData{ii},'uint8');%ignore whatever choice was made in read_bianry_data
    end
    
    % We will use the first header as the final header, but apply some
    % modifications to it.    
    KEYS_TO_SUM_OVER = {}; %add up values for these keys and make new key with total value
    KEYS_TO_LIST = {}; %concatenate values for these keys, with " * " delimiter, use "_original" prefix for new key
    %TODO: may want a KEYS_TO_CHECK_EQUALITY and raise warning if inequality detected
    %..or possibly the other way round...ie. assume equality is required
    %for all keys except some whitelist that is allowed to vary.
     
    blockSizeBytes = NaN; %if data consits of blocks that include a timestamp, we need to know block size...
    timeAtByte = NaN; % ..and where the timestamp is within the block
    switch type
        case {'tet'}
            KEYS_TO_SUM_OVER = {'duration','num_spikes'};
            KEYS_TO_LIST = {};
            hasData = true;
            timeAtByte = 1;
            blockSizeBytes = 50+4;
        case {'pos'}
            KEYS_TO_SUM_OVER = {'duration','num_pos_samples'};
            KEYS_TO_LIST = {'num_pos_samples'};
            hasData = true;
            timeAtByte = 1;
            blockSizeBytes = 20;
        case {'set'}
            KEYS_TO_SUM_OVER = {'duration'};
            hasData = false;
        case {'eeg'}
            KEYS_TO_SUM_OVER = {'duration','num_EEG_samples'};
            KEYS_TO_LIST = {'num_EEG_samples'};
            hasData = true; % data, but no timestamps
        case {'inp' 'cut'}
            warning(['No code written to deal with "' type '" files yet.']);
            return;
        otherwise
            warning('unrecognized filetype');
            return;
    end
    

    nSumOver = numel(KEYS_TO_SUM_OVER);
    nList = numel(KEYS_TO_LIST);
    nModifications = nSumOver + nList;
    
    % go through all headers summing up values such as duration and counts
    totals = zeros(nSumOver,1);
    tmp = cell(numel(KEYS_TO_SUM_OVER),1);
    for ii=1:nFiles
        [tmp{:}] = headers{ii}.KeyValueMultiple(KEYS_TO_SUM_OVER{:},'num');
        totals = totals + cell2mat(tmp);
    end
    
    
    % go through all headers concatenating values such as counts into a long string.
    lists = cell(nList,1);
    for ii=1:nFiles
        vals = {headers{ii}.KeyValueMultiple(KEYS_TO_LIST{:},'string')};
        for jj=1:nList
            lists{jj} = [lists{jj} ' * ' vals{jj}];
        end
    end
    
    
    % produce the modifications cell array of the form: 'name',value,'name',value,...
    modifications = cell(nModifications*2,1);
    for jj=1:nSumOver
        modifications{jj*2-1} = KEYS_TO_SUM_OVER{jj};
        modifications{jj*2} = num2str(totals(jj));
    end
    for jj=1:nList
        modifications{(jj+nSumOver)*2-1} = [KEYS_TO_LIST{jj} '_original'];
        modifications{(jj+nSumOver)*2} = lists{jj};
    end
    modifications{end+1} = 'merged_nFiles';
    modifications{end+1} = num2str(numel(srcFilenames));
    
    % build the new header, using the firt file's header and the
    % modifications computed above
    outHeader = mTintHeader(headers{1}.ToString,modifications);
    
    if hasData
        %If the file has a data section we need to concatenate across files.
        
        if ~isnan(timeAtByte)
            %But where the data contains timestamps it gets more complciated...
            totalOffset = 0;
            for ii=1:nFiles
                %reshape the binary data to the block size
                binaryData{ii} = reshape(binaryData{ii},blockSizeBytes,[]);
                
                %pick out the 4 bytes corresponding to timestamps, and read as int32
                times = reshape(binaryData{ii}(timeAtByte:timeAtByte+3,:),[],1);
                times = typecast(times,'uint32');
                
                % add the temporal offset
                if endian == 'L', times = swapbytes(times); end
                times = times + totalOffset;
                if endian == 'L', times = swapbytes(times); end
                
                % put the new time data back
                binaryData{ii}(timeAtByte:timeAtByte+3,:) = reshape(typecast(times,'uint8'),4,[]);
                
                %rehsape the binary data back to a single vector
                binaryData{ii} = reshape(binaryData{ii},[],1);
                
                % add to the totalOffset value, ready for the next file
                [duration,timebase] = headers{ii}.KeyValueMultiple('duration','timebase','num');
                totalOffset = totalOffset + duration*timebase;
            end
        end
        
        % combine the binary data by concatenating all the chunks.
        outRawBinaryData = cat(1,binaryData{:});
    end
    
    
    %if the output filename already exists then copy it and add the suffix for old
    %though if that copy already exists don't do the copy because it's safer
    %to "let the dead lie" in these circumstances.
    %this is more useful in the fix_binary_data function but possibly still
    %a nice thing to have here.
    if exist(destFilename,'file') && ~exist([destFilename suffixForOld],'file') 
        copyfile(destFilename,[destFilename suffixForOld]);
    end

    
    % write out the header and binary data (and the data_start/_end tokens)
    f = fopen(destFilename,'w');
    fprintf(f,'%s',outHeader.ToString);
    if hasData
        fprintf(f,DATA_START_TOKEN(3:end)); %lose initial \r\n
        fwrite(f,outRawBinaryData(:));
        fprintf(f,DATA_END_TOKEN,'');
    end
    fclose(f);


end

