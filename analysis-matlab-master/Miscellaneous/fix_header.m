function [header, RawBinaryData] = fix_header(filename,type,suffixForOld,modifications) 
% Read header text and binary data from tint data files.
%
% Modify the header using the key names-values specified in the cell array
% [modifications].
% The [type] value should be one of: 'set', 'pos', 'tet', 'eeg'.
% [suffixForOld] is added to the original file to create a copy, though if
% such a file alreayd exists, a new copy is not made.
%
% See also fix_binary_data and merge_binary_data

    f = fopen(filename,'r');
    if f == -1
        warning('No such file:\n%s',filename);
        header = []; RawBinaryData = [];
        return;
    end
    fclose(f);
    
    [header,RawBinaryData] = read_binary_data(filename,type);
    
    
    DATA_START_TOKEN = sprintf('\r\ndata_start');
    DATA_END_TOKEN = sprintf('\r\ndata_end');
    
    header = mTintHeader(header.ToString,modifications);
    

    if ~exist([filename suffixForOld],'file') %if the copied version already exists then do not overwirte it (as are probably reading an already modified version)
        copyfile(filename,[filename suffixForOld]);
    end
    delete(filename);

    f = fopen(filename,'w');
    fprintf(f,'%s',header.ToString);
    if ~isempty(RawBinaryData)
        fprintf(f,DATA_START_TOKEN(3:end)); %lose initial \r\n
        fwrite(f,RawBinaryData(:));
        fprintf(f,DATA_END_TOKEN,'');
    end
    fclose(f);

end

