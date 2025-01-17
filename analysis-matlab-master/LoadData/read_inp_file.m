function [ts, type, key, chan_states, header] = read_inp_file(flnm)
% Read Axona .inp file
%
% ts gives timestamps in seconds
% type is 'K' 'I' or 'O' for key, input or output
% key gives ascii or special keycode, valid at indices where type=='K'
% chan_states is 16xn and is valid at inidices where type=='I' or 'O'
%     it gives the logical state of each of 16 channels.
%
%

% Extract data from file
[header, data] = read_binary_data(flnm,'inp');
[~, ~, endian] = computer;  % 'L' or 'B', windows is L
[timebase, num_events, bytes_per_sample] = header.KeyValueMultiple(...
                'timebase', 'num_inp_samples', 'bytes_per_sample', 'num');

% Reshape data by num_events columns
if mod(numel(data),bytes_per_sample) ~= 0
    error('data length is not multiple of bytesPerSmaple');
end
data = reshape(data,bytes_per_sample, []);
if size(data,2) > num_events
    data = data(:, 1:num_events);
end

% Prepare timestamps
ts = data(1:4,:);
ts = typecast(ts(:),'uint32');
if endian == 'L', ts = swapbytes(ts); end
ts = double(ts)/timebase; % convert ts into seconds

% Deal with 3-byte payload...

type = char(data(5,:)'); % I=Input, O=Output, K=ASCII key or F1, F2 key etc.

% extract ASCII or special keycode, putting 0 where type~='K'
% See the axona file formats pdf for details of the special key codes
key = char(zeros(num_events, 1));
key(type=='K' & data(6,:)'==0) = data(7, type=='K' & data(6,:)'==0);
key(type=='K' & data(6,:)'~=0) = data(6, type=='K' & data(6,:)'~=0);

% For type='I' or 'O', we need to read the 16 channel logical states
chan_states = typecast(reshape(data(6:7,:),[], 1),'uint16');
if endian == 'L', chan_states = swapbytes(chan_states); end;
chan_states = reshape(bwunpack(uint32(chan_states)), 32, []);
chan_states = chan_states(1:16, :);
chan_states(:, type=='K') = 0;

end