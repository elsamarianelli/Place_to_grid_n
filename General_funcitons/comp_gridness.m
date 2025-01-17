%% Compare Gridness (hex vs square, standard vs expanded) for varying boundary effect parameters
%  Boundary_effects: gives the degree to which the place cells are
%  compressed closer to thr boundaries, 0.5 is very compressed, and 100 is
%  evenly distributed (sampling probability = distance to closest
%  boundary^(-1/boundary_effect)
%  Foraging trajectory --> from which a place cell activity X time
%  matrix is generated --> PCXPC activity covariance matrix -->
%  eigendecomposition eig(x) to take eigvectors which are projected back into 
%  2d space as a grid cell firing map (projected back by using as weights for 
%  Place cell inputs)
% Define environment dimensions and the number of cells
dim_x = 252;
dim_y = 252;
n_cells = 200;
n_steps = 200000;
boundary_effects = [0.5, 0.8, 1, 2, 100]; 
output_dir = 'gridness_results'; 
n_iterations = 10; % Number of iterations for each boundary effect

% Create the directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:1%length(boundary_effects)
    b_effect = boundary_effects(i);
    fprintf('Processing boundary effect: %.2f\n', b_effect);
    
    % Create a subfolder for this boundary effect
    subfolder = fullfile(output_dir, sprintf('boundary_effect_%.1f', b_effect));
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end

    subfolder = fullfile(output_dir, sprintf('uniform'));

    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end
    
    for iter = 1:n_iterations
        fprintf('  Iteration %d/%d\n', iter, n_iterations);
        
        % [1] Generating Environment
        fprintf('  Step 1/8: Generating Environment\n');
        env = generate_environment(dim_x, dim_y);
        
        % [2] Generating Trajectory
        fprintf('  Step 2/8: Generating Trajectory\n');
        traj = generate_trajectory(env, n_steps);
        
        % [3] Populating Place Cells
        fprintf('  Step 3/8: Populating Place Cells\n');
        [PlaceCellsUni, PlaceCellsTanni, env] = generate_place_cells(env, n_cells, dim_x, dim_y, b_effect);
        
        % Save place cells to the subfolder
        place_cells_file = fullfile(subfolder, sprintf('place_cells_iter_%d.mat', iter));
        save(place_cells_file, 'PlaceCellsUni', 'PlaceCellsTanni');
        % save(place_cells_file, 'PlaceCellsTanni');

        % [4] Generating PC Activity x Time Matrix
        fprintf('  Step 4/8: Generating PC Activity x Time Matrix\n');
        PCxTime_U = generate_pc_activity_matrix(PlaceCellsUni, traj, n_steps);
        % PCxTime_T = generate_pc_activity_matrix(PlaceCellsTanni, traj, n_steps);
        
        % [5] Generating Covariance Matrices
        fprintf('  Step 5/8: Generating Covariance Matrices\n');
        cov_M_U = cov(PCxTime_U');
        % cov_M_T = cov(PCxTime_T');
        
        % [6] Getting Grid Cells
        fprintf('  Step 6/8: Getting Grid Cells\n');
        cells_cov_U = get_grid(PlaceCellsUni, cov_M_U, env);
        % cells_cov_T = get_grid(PlaceCellsTanni, cov_M_T, env);
        
        % Save grid cells to the subfolder
        grid_cells_file = fullfile(subfolder, sprintf('grid_cells_iter_%d.mat', iter));
        save(grid_cells_file, 'cells_cov_U');
        % save(grid_cells_file, 'cells_cov_T');

        % [7] Analyzing Gridness
        fprintf('  Step 7/8: Analyzing Gridness\n');
        [square_cells_U, gridness_square_U, scale_square_U] = analGrids(cells_cov_U, 'square');
        [hexagon_cells_U, gridness_hex_U, scale_hex_U] = analGrids(cells_cov_U, 'hexagon');
        % [square_cells_T, gridness_square_T, scale_square_T] = analGrids(cells_cov_T, 'square');
        % [hexagon_cells_T, gridness_hex_T, scale_hex_T] = analGrids(cells_cov_T, 'hexagon');
        
        % [8] Saving Results
        fprintf('  Step 8/8: Saving Results\n');
        results_file = fullfile(subfolder, sprintf('results_iter_%d.mat', iter));
        save(results_file, 'square_cells_U', 'gridness_square_U', 'scale_square_U', 'hexagon_cells_U', 'gridness_hex_U', 'scale_hex_U');
        % save(results_file, 'square_cells_T', 'gridness_square_T', 'scale_square_T', 'hexagon_cells_T', 'gridness_hex_T', 'scale_hex_T');

    end
    
    fprintf('Completed boundary effect: %.2f\n\n', b_effect);
end



% Define the output directory and boundary effects
output_dir = 'gridness_results';
boundary_effects = [0.5, 0.8, 1, 2, 100]; 
n_iterations = 3;  % Number of iterations for each boundary effect

% Call the function to plot the results
plot_binned_results(output_dir, boundary_effects, n_iterations);

% generate environment
function env = generate_environment(dim_x, dim_y)
    n_polys = 1;
    polys = cell(n_polys, 1);
    polys{1} = [0 0, (dim_x-2) 0, (dim_x-2) (dim_y-2), 0 (dim_y-2), 0 0] + 2;
    env = GenerateEnv(polys, dim_x, dim_y);
end

% generate trajectory
function traj = generate_trajectory(env, n_steps)
    traj = HasselmoForage(env, n_steps);
    traj = round(traj);
end

% generate place cells
function [PlaceCellsUni, PlaceCellsTanni, env] = generate_place_cells(env, n_cells, dim_x, dim_y, boundary_effect)
    [xy_field, env] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, boundary_effect);
    av_bound_dist = nanmean(env.dwmap, 'all');
    fw = fieldWidth(av_bound_dist) / 3;
    PlaceCellsUni = generateUniformPCs(env, n_cells, xy_field, fw);
    PlaceCellsTanni = generateSanderPCs(env, n_cells, xy_field);
end

%  generate PC activity x time matrix
function PCxTime = generate_pc_activity_matrix(PlaceCells, traj, n_steps)
    x_indices = traj(:, 1);
    y_indices = traj(:, 2);
    fmaps = cellfun(@(PC) PC.fmap, PlaceCells, 'UniformOutput', false);
    PCxTime = zeros(length(fmaps), n_steps);
    for i = 1:n_steps
        x = x_indices(i);
        y = y_indices(i);
        for j = 1:length(fmaps)
            PC_fmap = fmaps{j};
            PCxTime(j, i) = PC_fmap(x, y);
        end
    end
end



% plot results
function plot_results(output_dir)
    % Load all .mat files from the directory and plot results
    files = dir(fullfile(output_dir, '*.mat'));
    figure;
    edges = linspace(-1.5, 1.5, 21); 
    for i = 1:length(files)
        data = load(fullfile(output_dir, files(i).name));
        histogram(data.gridness_square_U(:,1), edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'b'); hold on;
        histogram(data.gridness_square_T(:,1), edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'r');
        histogram(data.gridness_hex_U(:,1), edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'g');
        histogram(data.gridness_hex_T(:,1), edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'm');
    end
    xlabel('Gridness');
    ylabel('Probability Density');
    legend('SUstd', 'STstd', 'HUstd', 'HTstd');
    title('Distributions of Gridness for Different Boundary Effects');
end
