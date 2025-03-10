function [trajectories, bounds] = get_real_traj(data_folder)
% code to extract foraging trajectories from other experiments to use in
% simulations

addpath("real_data")
addpath("analysis-matlab\LoadData\")
addpath('C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\real_data_analysis')
%% [1] load data

% Define the folder containing the data files
% data_folder = 'C:\Users\Elsa Marianelli\Desktop\Rodent_Data';
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
numTetrodes = 16; % max number of tetrodes in any file 
numCells = 23; % max number of cells in any file
in.GET_TRODES = 1:numTetrodes;
in.data_path = data_folder;

trajectories = {};
bounds = {};
for f = 2:length(unique_file_codes)
    in.flnmroot = [unique_file_codes{f}, '_screening'];
    trial = load_trial(in);
    pos = trial.posData.xy;
    % plot(pos(:,1), pos(:,2), '.');
    trajectories{f} = pos;
    bound = [max(pos(:,1)), max(pos(:,2))];
    bounds{f} = bound;
end


end