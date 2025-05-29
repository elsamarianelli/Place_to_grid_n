function traj = load_premade_traj(iter)

% Define trajectory folder path
traj_folder = fullfile(pwd, 'premade_traj_and_env'); % Relative to Place_to_grid_n
% traj_folder = fullfile(pwd, 'grids_data', 'premade_traj_and_env'); % Relative to Place_to_grid_n

% Get list of all saved trajectory files
traj_files = dir(fullfile(traj_folder, 'trajectory_*.mat'));

% Ensure we have valid files
if isempty(traj_files)
    error('No trajectory files found in %s', traj_folder);
end

% Select a trajectory based on iteration index
traj_index = mod(iter-1, length(traj_files)) + 1; % Loops through the 13 files

% Load the selected trajectory file
traj_data = load(fullfile(traj_folder, traj_files(traj_index).name));

% Extract the stored trajectory
traj = traj_data.traj;

fprintf('Using trajectory: %s for iteration %d\n', traj_files(traj_index).name, iter);

end