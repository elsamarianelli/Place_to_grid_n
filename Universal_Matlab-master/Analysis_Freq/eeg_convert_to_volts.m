function [ eeg ] = eeg_convert_to_volts( eeg, eegCh, header, varargin )
%EEG_CONVERT_TO_VOLTS Takes eeg from DACQ and converts to volts as a 'double'
% EEG is recorded as by DACQ as an int8 or int16 value. Need to convert that to 
% a value in volts using information about the hardware that was used for the recording 
% and also the gain applied to the eeg channel. Also conver to double for ease of use.
% Standard (250Hz) eeg is 8bit whereas high rate (4.8kHz egf files) are 16bit. User
% needs to specify 4th input variable to indicate which of these to use otherewise
% default is 8bit (for back compatability).
%
% ARGS - takes 3 plus 1 optional
% eeg       [nEEGpts x 2] the eeg reported by read_eeg should be int8,
%           first column is the eeg value and second is the timestamp 
% eegCh     channel that eeg was recorded on (usually a multiple of 4) [actual channel]
% header    main session header, typically returned in tint.header
% [optional] type - either '8bit' or '16bit', if not supplied defaults to '8bit'
%
% RETURNS
% eeg - eeg in volts
%
% Example
% eeg = eeg_convert_to_volts(eeg, myChannel, data.settings, '16bit')
% or
% eeg = eeg_convert_to_volts(eeg, myChannel, data.settings)

% --- HOUSE KEEPING ---------------------------------------------------------------
switch nargin
    case 3
        eegType='8bit'; %Standard eeg
        
    case 4
        eegType=varargin{1};
        clear varargin
        
    otherwise
        error('Invalid number of arguments - must specify 3 or 4');
end

% --- MAIN -------------------------------------------------------------------------

%1) Get conversion factor based on hardware
%Conversion actually depends on which hardware was used - usb uses
%lower intrinisic gain. The version of the system can be deterimined
%from DACQ software version
softwareV = key_value('ADC_fullscale_mv', header);
softwareV = str2double(softwareV{1});
switch softwareV
    case 1500,        convFact = 1.5; %for USB system
    case 3680,     convFact = 3.68; %for old system
    case 3000, convFact = 3.00; %another old system
    otherwise,       error('System type is not matched in EEG header');
end


%2) Get eeg into range between -1 and 1 (though some values might exceed this range if
%it has been filtered
%EEG as returned from read_eeg_file is a two column matrix, first column
%being the voltage data and the second being timestamps. This second column
%is removed by eeg_convert_to_volts.m but double check here that it is not
%present. First check eeg is column form, then just take first column.
if size(eeg,2)>size(eeg,1)
    eeg         =eeg';
end
eeg         =eeg(:,1);

eeg         =double(eeg); %Should be double anyway but just incase
switch eegType
    case '8bit' %Standard 250Hz eeg
        %First check that values are in expected range (-128 to 127)
        percAboveRange=sum(eeg>=128)/length(eeg); %Value 0 to 1
        if percAboveRange>0.01
            warning('More than 1% values above int8 range expected');
        end
        eeg=eeg/128; %Mainly values between -1 and 1
        
    case '16bit'
          eeg=eeg/32768; %Mainly values between -1 and 1
          
    otherwise
         error('Unrecognised eeg type specified - ''8bit'' or ''16bit'' allowed');
end

%3) Convert to volts using conversion and correct value from header - note
%that eeg channel is 0-indexed in the header, so to get from the value in
%data.eeg.actualCh (which is 1-indexed)to the appropriate header have to
%subtract 1
eegGain = key_value(['gain_ch_', num2str(eegCh-1)], header,'num', 'exact');  %Use the channel to get the gain
eeg=(convFact/eegGain)*eeg;




