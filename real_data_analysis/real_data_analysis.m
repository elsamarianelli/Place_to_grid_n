%% code to look at real grid data
%  looking at real single cell recordings from old experiments, to try and
%  find cells which surpass gridness thresholds and plot using
%  ucl-hippocampus github functions 

addpath("real_data")
addpath("analysis-matlab\LoadData\")
addpath("analysis-matlab\Miscellaneous\")
addpath("analysis-matlab\GridAnalysis\")

%% [1] load data

% Define the folder containing the data files
data_folder = 'C:\Users\Elsa Marianelli\Desktop\Rodent_Data';
% data_folder = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\real_data_analysis\real_data';
file_list = dir(fullfile(data_folder, '*')); % Get all files

% Extract filenames
file_names = {file_list.name};

% Remove '.' and '..' directories
file_names = file_names(~ismember(file_names, {'.', '..'}));

% Regular expression to extract the file code (everything before the second underscore or before the first dot)
file_codes = regexp(file_names, '^[^._]+_[^._]+', 'match', 'once');

% Get unique file codes
unique_file_codes = unique(file_codes);
disp('Unique file codes:')
disp(unique_file_codes)

%% [2] get gridness scores from ALL tetrodes and cells and plot those with 
% gridness above  over 95th percentile of shuffled gridness data

% trial cell with above threshold gridness for: 20151025_R2335, tetrode 5, cell 4 
                                               %20151024_R2335, tetrode 4, cell 6 

numTetrodes = 16; % max number of tetrodes in any file 
numCells = 23; % max number of cells in any file
in.GET_TRODES = 1:numTetrodes;
in.data_path = data_folder;
in.cutSuffix='_';
in.GET_POS=true;
in.GET_CUT=true;

gridnessScores = nan(size(unique_file_codes, 2), numTetrodes, numCells); 
threshold_score = nan(size(unique_file_codes, 2), numTetrodes, numCells);
% above stores 1s for cells which are both over .5 and shuffled threshold
above = nan(size(unique_file_codes, 2), numTetrodes, numCells); 

% loop through each unique field code, each tetrode and all it's cells, and
% find cells above gridness threshold (for both .5 threshold and shuffled
% 95th percentile threshold)

for f = 2:length(unique_file_codes)
    in.flnmroot = [unique_file_codes{f}, '_screening'];
    trial = load_trial(in);
    for t = 1:numTetrodes
         
        sp = trial.spikeData(t);
        uniqueCells = unique(sp.cut); % Get unique cell IDs
        
        if ~isempty(uniqueCells) & ~all(uniqueCells==0)
            for c = uniqueCells(2:end)'
                
                % Extract spike timestamps for the given cell
                in.spikeTS = sp.timestamp(sp.cut == c);
                in.spikePosInd = sp.pos_sample(sp.cut == c);
                
                % Position and shuffle parameters
                in.posXY = trial.posData.xy;
                in.MAP_TYPE = 'rate';
                in.SHUFFLES = 1000; 
                in.psr = trial.posData.sample_rate; % pos sample rate
                in.PLOT_ON = false;
                in.autocorr_in.PLOT_ON = false;
                in.SMOOTHING = 8; % boxcar width in bins for smoothing
    
                % Compute shuffled gridness score
                ret = shuffledGridness(in); % randomly permute spike times 
                                            % relative to behavior and
                                            % recalculate gridness
             
                % Compute the 95th percentile threshold
                shuffled_values = ret.gridnessShuffles; 
                threshold = prctile(shuffled_values, 90);
    
                % Store results in a structured format
                gridnessScores(f, t, c) = ret.gridness;
                threshold_score(f, t,c) = threshold;
                above_check = ((ret.gridness) > threshold) && (ret.gridness>.5);
                above(f, t, c) = above_check;
                
                % display tetrode and cell 
                if above_check
                    CHECK = ' ..... GRIDNESS ABOVE THRESHOLD';
                else
                    CHECK = '';
                end

                disp(['Checking gridness for: file ', unique_file_codes{f}, ...
                    ', tetrode ', num2str(t), ', cell ', num2str(c), CHECK]);
    
            end
        else 
            continue        
        end
    end
end

%% [3] plotting for cells with gridness above both thresholds 

% Create a logical mask where the cells contain 1
mask = (above == 1);
% Find row and column indices
[f,t, c] = ind2sub(size(mask), find(mask));
indices = [f, t, c];

% plotting all alongside each other 
t = tiledlayout(10, 6, 'Padding', 'tight', 'TileSpacing', 'none');
for idx = 1:size(indices)
    f = indices(idx, 1);
    in.flnmroot = [unique_file_codes{f}, '_screening'];
    trial = load_trial(in);
    t = indices(idx, 2); 
    c = indices(idx, 3);
    in.SMOOTHING = 7;

    % use this in a loop to generate individual plots for each cell
    % stats = get_figure_and_stats(trial, t, c, in);

    % use this ifyou just want to look at the metiric sac and the cells
    % side by side all together
    stats_simpler = get_figure_and_stats_for_tiled_layout(trial, t, c, in);
end

% %% saving record of which cells are over threshold
% save_folder = 'C:\Users\Elsa Marianelli\Desktop\Rodent_Data\Results';
% if ~exist(save_folder, 'dir')
%     mkdir(save_folder);
% end
% save_path = fullfile(save_folder, 'indices.mat');
% save(save_path, 'indices');
% 
% %%
% % plotting a single cell to try and fix grid metric methods
% t = 8; c = 2; % unclear one
% t = 8; c = 1; % clearer
% in.SMOOTHING = 8;
% stats = get_figure_and_stats(trial, t, c, in);


