function ret = load_trial(varargin)
% See https://github.com/UCL/mTint/tree/master/LoadData for info

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   
    
    nTrodes = max([in.GET_TRODES(:)' 4]); % i.e. 4 or more

    % work out the full filepath for all the files (may not use all of them)
    set_file = fullfile(in.data_path,[in.flnmroot '.set']);
    pos_file = fullfile(in.data_path ,[in.flnmroot '.pos']);
    eeg_file = build_eeg_file_list(in.data_path,in.flnmroot); %see subfunction, this is a bit complicated
    if in.USE_CLU_FILES
        cut_file = @(x)fullfile(in.data_path,[in.flnmroot in.cutSuffix '.clu.' num2str(x)]); %annoymous function to build cut filename string       
    else
        cut_file = @(x)fullfile(in.data_path,[in.flnmroot in.cutSuffix num2str(x) '.cut']); %annoymous function to build cut filename string
    end
    trode_file = @(x)fullfile(in.data_path,[in.flnmroot '.' num2str(x)]); %annoymous function to build tet filename string
    
    % Read settings file  
    setfile_header = read_binary_data(set_file,'set');

    % make an empty structure for all tetrodes
    spikeData.exists = 0;
    spikeData.header = [];
    spikeData.timestamp = [];
    spikeData.spike = [];
    spikeData.pos_sample = [];
    spikeData.cut = [];
    spikeData = repmat(spikeData,nTrodes,1);
    
    % For each tetrode requested, read tetrode and cut files

    for t = in.GET_TRODES
        spikeData(t).trode_id = t;
        if exist(trode_file(t),'file')
            spikeData(t).exists = 1;
            if in.GET_WAVEFORMS
                [spikeData(t).header, spikeData(t).timestamp, spikeData(t).spike] = read_tetrode_file(trode_file(t));
                spikeData(t).spike = double(spikeData(t).spike);
            else
                [spikeData(t).header, spikeData(t).timestamp] = read_tetrode_file(trode_file(t));
            end
        end
        if in.GET_CUT && exist(cut_file(t),'file')
            spikeData(t).cut = read_cut_file(cut_file(t));
        end
    end 
    
    % Read pos file
    if in.GET_POS && exist(pos_file,'file')
        posData.exists = 1;
        [posData.led_pos, posData.led_pix, posData.header] = read_pos_file(pos_file);
        if ~isempty(in.POS_HEADER_MODIFY)
            fprintf('[%s] Modifications requested for pos header\n',in.flnmroot)
            posData.header = mTintHeader(posData.header.ToString,in.POS_HEADER_MODIFY); %replace the origial header with a new header that has some values modified
        end
        
        [posData.xy, posData.dir, posData.dir_disp, posData.speed, posData.times] = postprocess_pos_data(posData, in.POS_PROCESS_MAX_SPEED, in.POS_PROCESS_BOX_CAR, setfile_header);
        if in.CONVERT_POS_TO_CM
            posData.xy = posData.xy / posData.header.KeyValue('pixels_per_metre',  'num') * 100;
        end
        posData.sample_rate = posData.header.KeyValue('sample_rate',  'num');
        posData.acc = makeAcceleration(posData.speed, posData.sample_rate,in.ACC_BOX_CAR);

        trialData.trial_start_time = posData.times(1) - (1 / posData.sample_rate);
        trialData.trial_end_time =  posData.header.KeyValue('duration', 'num');

        % Fill in .pos_sample field of the tetrode structure - the position (and eeg) sample co-occuring with a given spike. Timestamp is in seconds.
        for t = in.GET_TRODES
            if spikeData(t).exists == 1
                spikeData(t).pos_sample = ceil(mean(spikeData(t).timestamp,2) .* posData.sample_rate);
            end
        end
    else
        posData.exists = 0;
        posData.led_pos = [];
        posData.led_pix = [];
        posData.header = [];
        posData.xy = [];
        posData.dir = [];
        posData.speed = [];
        posData.times = [];
        posData.sample_rate = [];
        posData.acc = [];
        trialData.trial_start_time = [];
        trialData.trial_end_time = [];
    end

    % Read eeg data
    eegData.exists = 0;
    eegData.eeg_and_time = [];
    eegData.header = [];
    eegData.sample_rate = [];
    eegData.eeg_range = [];
    if in.GET_EEG 
        eegData = repmat(eegData,numel(eeg_file),1);
        for ii=1:numel(eeg_file)
            if exist(eeg_file{ii},'file')
                eegData(ii).exists = 1;
                [eegData(ii).eeg_and_time, eegData(ii).header, eegData(ii).conversionFactor] = read_eeg_file(eeg_file{ii},ii,setfile_header,in.CONVERT_EEG_TO_VOLTS);
                eegData(ii).sample_rate = eegData(ii).header.KeyValue('sample_rate',  'num');
                if in.GET_EEG_RANGE
                    mean_eeg = mean(eegData(ii).eeg_and_time(:,1));
                    std_eeg = std(eegData(ii).eeg_and_time(:,1));
                    eegData(ii).eeg_range = mean_eeg + (2 .* std_eeg);
                end
                if isempty(trialData.trial_start_time)
                    trialData.trial_start_time = eegData(ii).eeg_and_time(1,2) - (1 / eegData(ii).sample_rate);
                    trialData.trial_end_time = eegData(ii).eeg_and_time(end,2);
                end
            end
        end
    end

    if ~isempty(in.DISCARD_MINUTES)
        
        trialData.discarded_minutes = in.DISCARD_MINUTES;
        
        % for each of the data files, produce a vector of logicals specifying rows to be removed
        badPos = false(size(posData.times));
        badSpike = cellfun(@(x)false(size(x,1),1),{spikeData(:).timestamp},'UniformOutput',false);
        badEEG = cellfun(@(x)false(size(x,1),1),{eegData(:).eeg_and_time},'UniformOutput',false);
        
        for ii=1:size(in.DISCARD_MINUTES,1)
            dstart = in.DISCARD_MINUTES(ii,1) *60;
            dend = in.DISCARD_MINUTES(ii,2) *60;
            if isnan(dend), dend = dstart+1; end
            
            badPos = badPos | (posData.times > dstart & posData.times < dend);
            for t=in.GET_TRODES
                badSpike{t} = badSpike{t} | (spikeData(t).timestamp(:,1) > dstart & spikeData(t).timestamp(:,1) < dend);
            end
            for eg = 1:numel(eegData)
                if ~isempty(eegData(eg).eeg_and_time)
                    badEEG{eg} = badEEG{eg} | (eegData(eg).eeg_and_time(:,2) > dstart & eegData(eg).eeg_and_time(:,2) < dend); 
                end
            end
        end
        
        
        % remove the required rows in pos
        if ~isempty(posData.times)
            fprintf('[%s] Discarded %0.0f%% of posData\n',in.flnmroot,mean(badPos)*100);
            posData.led_pos(badPos,:) = [];
            posData.led_pix(badPos,:) = [];
            posData.xy(badPos,:) = [];
            posData.dir(badPos) = [];
            posData.speed(badPos) = [];
            posData.times(badPos) = [];
            posData.acc(badPos) = [];
        end
        
        % remove the required rows in each tetrode
        for t=in.GET_TRODES
            fprintf('[%s] Discarded %0.0f%% of tetrode-%d data\n',in.flnmroot,mean(badSpike{t})*100,t);
            spikeData(t).timestamp(badSpike{t},:) = [];
            spikeData(t).pos_sample(badSpike{t},:) = [];
            if ~isempty(spikeData(t).cut)
                spikeData(t).cut(badSpike{t},:) = [];
            end
            if ~isempty(spikeData(t).spike)
                spikeData(t).spike(badSpike{t},:,:) = [];
            end
        end
        
        for eg = 1:numel(eegData)
            if ~isempty(eegData(eg).eeg_and_time)
                fprintf('[%s] Discarded %0.0f%% of eeg-%d data\n',in.flnmroot,mean(badEEG{eg})*100,eg);
                eegData(eg).eeg_and_time(badEEG{eg},:) = []; %M-Lint is wrong here in R2012b
            end
        end
        
    else
        trialData.discarded_minutes = [];
    end
    
    %collect everything to return
    ret.setfile_header = setfile_header;
    ret.eegData = eegData;
    ret.trialData = trialData;
    ret.posData = posData;
    ret.spikeData = spikeData;

end

function acc = makeAcceleration(speed, pos_sample_rate,box_car)
dt = 1/pos_sample_rate;
acc = [diff(speed)/dt ; 0];
acc(end) = acc(end-1); %useful to have it the same length as speed

acc = filter(ones(box_car,1)/box_car,1,acc); %box-car running average
end


function eeg_files_list = build_eeg_file_list(data_path, flnmroot)
if data_path(end) ~= filesep, data_path(end+1) = filesep; end % TODO: tidy up this probelm more generally
file_list = dir([data_path flnmroot '.e*']);  %get a list of files that match the pattern

eeg_files_list = cell(4,1);
for ii=1:numel(file_list)
    
    % for each file, decide whether it is the 1st, 2nd, 3rd etc. eeg
    [ignore,ignore,ext] = fileparts(file_list(ii).name);
    if ~isempty(strfind(lower(ext),'f'))
        continue
    end
    number = regexp(ext,'[\d]','match');
    if isempty(number),
        number = 1;
    else
        number = str2double(number{1});
    end
    
    % store the filename in the cell array according to it's number
    eeg_files_list{number} = [data_path file_list(ii).name];
end

end


function in = DealWithInputs(varargin)
    defaults.GET_POS = true;
    defaults.GET_TRODES = [1 2 3 4];
    defaults.GET_WAVEFORMS = false;
    defaults.GET_CUT = true;
    defaults.GET_EEG = false;
    defaults.CONVERT_EEG_TO_VOLTS = true;
    defaults.CONVERT_POS_TO_CM = true;
    defaults.data_path  = 'C:\ (example)';
    defaults.flnmroot = 'R1884 2008-01-27 3 (example)';
    defaults.POS_PROCESS_MAX_SPEED = 1;
    defaults.POS_PROCESS_BOX_CAR = 0.4;
    defaults.ACC_BOX_CAR = 5;
    defaults.POS_HEADER_MODIFY = {};
    defaults.DISCARD_MINUTES = [];
    defaults.GET_EEG_RANGE = false;
    defaults.cutSuffix = '_';
    defaults.USE_CLU_FILES = false;
    
    VERSION = 2.06;
    
    
    % Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0
        in = defaults;
        return;
    end

    if isstruct(varargin{1})
        if nargin > 1 && VERSION ~= varargin{2}
            error(['%s called with version number %g, but code is version %g.\n' ...
                'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
                'the function without a version number.'],mfilename,varargin{2},VERSION);
        end
        in = ModifyExistingFields(defaults,varargin{1});
    else
        in = ModifyExistingFields(defaults,varargin{:});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
