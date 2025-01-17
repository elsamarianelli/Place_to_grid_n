function [varargout] = read_tetrode_file(flnm)
%
%
% Two forms:
% [header, timestamp, waveforms] = read_tetrode_file(flnm) does read waveforms
% [header, timestamp] = read_tetrode_file(flnm) does not read waveforms
%
% TAKES
% flnm : fully qualified filename including path plus extension e.g. '.1', '.2' 
%
% RETURNS: varargout of either 2 [header, timestamps] or 3 [header, timestamps, spike_waveforms]
% vars
% header : header information as key-value pairs (cell array of strings)
% timestamp : timings for each spike. [num_spikes x num_chans]
% spike : spike waveforms for each channel [num_spikes x samples_per_spike x num_chans ]



%--------------------------------------------------------------------------
%Deterimine if spike waveforms are required if only two arguments out then
%they are not, if three they are
if nargout==3,    waveforms=1;
else    waveforms=0;
end

% Determine if computer is big or little endian (most windows machines are latter)
% we don't use this information until near the end of the function
[ignore,ignore,endian] = computer;  % 'L' or 'B'

if waveforms
    % get header and read in all binary data as individual int8s
    [header, data] = read_binary_data(flnm,'tet');
else
    % get header and read in the first four int8s of every 216 bytes
    [header, data] = read_binary_data(flnm,'tet_times');
end

% error checking
if isempty(data)
    fprintf('\nCould not load tetrode %s: invalid data file', flnm(end));
    varargout{1} = 'invalid tetrode file';
    varargout{2} = 'invalid tetrode file';
    varargout{3} = 'invalid tetrode file';
    return
end
if ~all(header.KeysPresent({'num_chans','timebase','samples_per_spike','num_spikes'}))
    fprintf('Could not load tetrode %s data: header invalid', flnm(end));
    varargout{1} = 'invalid tetrode header';
    varargout{2} = 'invalid tetrode header';
    varargout{3} = 'invalid tetrode header';
    return
end

varargout{1} = header;

% Assume 4 channels and each spike is 50samps at 8bits
timebase = header.KeyValue('timebase','num');
num_spikes = header.KeyValue('num_spikes','num');

%Now we convert the list of bytes into the required form
%
%Note structure of data is such that for every spikes on each channel a 4
%byte timestamp is logged followed by the 50 waveform points each being 1
%byte. So 54bytes per channel, 216 per spikes. Key point is that timestamps
%are always the same for each channel, so fully redundant (i.e. just read
%first ch). The timestamp is in bigendian ordering so if this computer lives
%in a little endian world we will need to swap bytes to get meaningful values.
%Waveform points are just int8s so aren't affect by byte ordering.
%Data also seems to be padded beyond the end of num_spikes.

switch waveforms
    case 1 %data is bytes corresponding to timestamps and waveform
        data = reshape(data,216,[]);
        timestamps = data(1:4,1:num_spikes);
        spikes= data([5:54, 59:108, 113:162, 167:216],1:num_spikes)';
        spikes=reshape(spikes,[num_spikes, 50,4]); % int8 [nSpike x 50 x 4]
        varargout{3} = spikes;
        
    case 0 %data is bytes corresponding only to timestamp
        data = reshape(data,4,[]);
        timestamps = data(:,1:num_spikes);
end

%convert timestamp bytes to int32 with same endian-ness as this computer
timestamps = typecast(timestamps(:),'int32');
if endian == 'L' 
    timestamps = swapbytes(timestamps);
end
%convert timestamps to seconds
timestamps = double(timestamps)/timebase;
%repmat for back compatibility
varargout{2} = repmat(timestamps,[1 4]);



