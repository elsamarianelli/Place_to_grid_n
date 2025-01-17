%% Generate and save gridness and scale data form uniform vs sanders place cell --> grid cells
%  saves gridness scales in binned environemnt for uniform and sanders
%  style PCs to GCs
%  no non negativity included yet. 

% Define environment dimensions and the number of cells
dim_x = 351;
dim_y = 252;
n_polys = 1;
polys = cell(n_polys, 1);
% polys{1} = [0 (250/2)-30, 349 0, 349 250, 0 (250/2)+30, 0 (250/2)-30] + 2;
polys{1} = [0 0, 349, 0, 349 250, 0 250, 0 0] + 2; 
n_cells = 200;
n_steps = 360000;
output_dir = 'gridness_results_circ'; 
n_iterations = 20; 

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
    [PlaceCellsUni, PlaceCellsTanni, env] = generate_place_cells(env, n_cells, dim_x, dim_y, 2);
    
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
    [M_U, R_U] = trainModel(PlaceCellsUni, M_U, R_U, traj, time_lag);
    [M_T, R_T] = trainModel(PlaceCellsTanni, M_T, R_T, traj, time_lag);

    % [5] Get Grid and Place cells from trained SR model
    fprintf('  Step 5/8: Getting Place and Grid\n');
    % Recalculate place and grid cells based on the trained model
    cells_U = getPlace(PlaceCellsUni, M_U, env);
    cells_U = getGrid(cells_U, M_U, env, 'off'); % off/on for non-negativity constraint 
    cells_T = getPlace(PlaceCellsTanni, M_T, env);
    cells_T = getGrid(cells_T, M_T, env, 'off'); % off/on for non-negativity constraint 
  
    % [6] Saving Grid and Place cells
    fprintf('  Step 6/8: Save\n')
    % Save grid cells to the subfolder
    SR_cells_file = fullfile(subfolder, sprintf('grid_place_cells_SR_iter_%d.mat', iter));
    save(SR_cells_file, 'cells_U', 'cells_T');
    
    % [7] Analyzing Gridness
    fprintf('  Step 7/8: Analyzing Gridness\n');
    [square_cells_U, gridness_square_U, scale_square_U] = analGrids(cells_U, 'square');
    [hexagon_cells_U, gridness_hex_U, scale_hex_U] = analGrids(cells_U, 'hexagon');

    [square_cells_T, gridness_square_T, scale_square_T] = analGrids(cells_T, 'square');
    [hexagon_cells_T, gridness_hex_T, scale_hex_T] = analGrids(cells_T, 'hexagon');
    
    % grid_scales_binned_U = get_binned_scales(cells_U, n_cells);
    % grid_scales_binned_T = get_binned_scales(cells_T, n_cells);

    % [8] Saving Results
    fprintf('  Step 8/8: Saving Results\n');
    results_file_U = fullfile(subfolder, sprintf('U_results%d.mat', iter));
    results_file_T = fullfile(subfolder, sprintf('T_results%d.mat', iter));

    save(results_file_U, 'square_cells_U', 'gridness_square_U', 'scale_square_U', 'hexagon_cells_U', 'gridness_hex_U', 'scale_hex_U');%, 'grid_scales_binned_U');
    save(results_file_T, 'square_cells_T', 'gridness_square_T', 'scale_square_T', 'hexagon_cells_T', 'gridness_hex_T', 'scale_hex_T');%, 'grid_scales_binned_T');

end
