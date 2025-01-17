function [eeg_and_time, header, factor] = read_eeg_file(flnm,eegNum,setfile_header,CONVERT_TO_VOLTS)
% reads eeg file into header object and data array
% if CONVERT_TO_VOLTS is true it then looks up various key values to work
% out the conversion factor required to scale the raw data to volts.
% Finally it performs the conversion.

[header, eeg] = read_binary_data(flnm,'eeg');
[sample_rate,bytes_per_samp] = header.KeyValueMultiple('sample_rate','bytes_per_sample','num');


if bytes_per_samp == 1 % and sample_rate is 250Hz
    num_EEG_samples = header.KeyValueExact('num_EEG_samples','num');
    
elseif bytes_per_samp == 2 % and sample_rate is 4800Hz   
    num_EEG_samples = header.KeyValueExact('num_EGF_samples','num');
    [ignore,ignore,endian] = computer;  % 'L' or 'B' i.e. little or big
    eeg = typecast(eeg,'int16');
    if endian == 'B' %TODO: why is this B whereas in read_tetrode it's L? 
        eeg = swapbytes(eeg);
    end
    
else
    error('eeg bytes_per_samp is neither 1 nor 2.');
end
eeg = single(eeg); %hopefully no one needs it to be double

if num_EEG_samples ~= numel(eeg)
    error('eeg data is wrong size');
end

if nargin==1 || ~CONVERT_TO_VOLTS
    factor = 1;
    
else %do perform conversion to volts 
    
    %There are three separate factors which must be identified before we can
    %perform the conversion to volts.

    % 1) Find out which channel eeg was actually on and lookup gain for the
    % channel.
    chn = setfile_header.KeyValueExact(['EEG_ch_' num2str(eegNum)],'num');
    eegChanGain = setfile_header.KeyValueExact(['gain_ch_', num2str(chn-1)],'num');

    %2) Get conversion factor based on system type
    ADC_fullscale_mv = setfile_header.KeyValueExact('ADC_fullscale_mv','num');
    switch ADC_fullscale_mv
        case {1500,3680,3000} %1500 for USB system, other two are for older systems
            convFact = ADC_fullscale_mv/1000; 
        otherwise, error('System type is not matched in EEG header');
    end

    %3) Conversion factor based on bit rate.
    %   This factor will get eeg into range between -1 and 1 (though some values
    %   might exceed this range if it has been filtered
    maxIntVal = 2^(bytes_per_samp*8-1);

    %4) perform conversion
    factor = 1/maxIntVal*convFact/eegChanGain;
    eeg = eeg*factor;
end

% These final three lines are slow, and ultiamtely aren't actually
% required.  So maybe get rid of them?
times = (single(1):single(num_EEG_samples))'; 
times = times/single(sample_rate);
eeg_and_time = [eeg times];


end
