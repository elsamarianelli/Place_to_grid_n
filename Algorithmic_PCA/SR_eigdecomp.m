%% Generate and save Grid maps 
%  generated from doing PCA on SR Matrix learnt through TDLM during random foraging with 
%  place cell firing as successor features. 

addpath('/Users/elsamarianelli/Documents/GitHub')
addpath('/Users/elsamarianelli/Documents/GitHub/bound_warped_grids_new/Algorithmic_PCA/Functions/')
addpath('/Users/elsamarianelli/Documents/GitHub/bound_warped_grids_new/General_funcitons/')

% Define environment dimensions and the number of cells
dim_x = 351;
dim_y = 252;
n_polys = 1;
polys = cell(n_polys, 1);
polys{1} = [0 0, 349, 0, 349 250, 0 250, 0 0] + 2; 
n_cells = 200;
n_steps = 36000;
output_dir = 'results_grid_maps_rect'; 
n_iterations = 1; 

% Create the directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for iter =1:n_iterations

    % Create a subfolder for this iteration
    subfolder = fullfile(output_dir, sprintf('iteration_%.1f', iter));
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end

    fprintf('  Iteration %d/%d\n', iter, n_iterations);
    
    % [1] Generating Environment
    fprintf('  Step 1/8: Generating Environment\n');
    env = GenerateEnv(polys, dim_x, dim_y, 'trapezoid');
    
    % [2] Generating Trajectory
    fprintf('  Step 2/8: Generating Trajectory\n');
    traj = generate_trajectory(env, n_steps);
    
    % [3] Populating Place Cells
    fprintf('  Step 3/8: Populating Place Cells\n');
    [PlaceCellsUni, PlaceCellsTanni, env] = generate_place_cells(env, n_cells, dim_x, dim_y, 2, 'random');
    
    % Save place cells to the subfolder
    place_cells_file = fullfile(subfolder, sprintf('place_cells_iter_%d.mat', iter));
    save(place_cells_file, 'PlaceCellsUni', 'PlaceCellsTanni');
    
    % [4] Training SR Matrix
    fprintf('  Step 4/8: Training SR matrix using TDLR\n');
    % Initialize the model parameters for training
    M_U = ones(n_cells); % Initialize connection matrix
    R_U = zeros(n_cells, 1); % Initialize reward vector
    M_T = ones(n_cells); 
    R_T = zeros(n_cells, 1); 
    time_lag = 1; % Time lag of 100ms (a theta cycle)

    % Train the model using the generated trajectory
    [M_U, ~] = trainModel(PlaceCellsUni, M_U, R_U, traj, time_lag);
    [M_T, ~] = trainModel(PlaceCellsTanni, M_T, R_T, traj, time_lag);

    % [5] Get Grid and Place cells from trained SR model
    fprintf('  Step 5/8: Getting Place and Grid\n');

    % Recalculate place and grid cells based on the trained model
    cells_U = getPlace(PlaceCellsUni, M_U, env);
    cells_U = getGrid(cells_U, M_U, env, 'off'); % off/on for non-negativity constraint 
    cells_T = getPlace(PlaceCellsTanni, M_T, env);
    cells_T = getGrid(cells_T, M_T, env, 'off'); 
  
    % [6] Saving Grid and Place cells
    fprintf('Step 6/8: Save\n')
    % Save grid cells to the subfolder
    SR_cells_file = fullfile(subfolder, sprintf('grid_place_cells_SR_iter_%d.mat', iter));
    save(SR_cells_file, 'cells_U', 'cells_T');
  
end

