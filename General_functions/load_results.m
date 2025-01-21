function results = load_results(output_dir, n_iterations)
% Initialize the structure to hold all data across iterations
results = struct();

% Pre-allocate arrays to hold data from all iterations
results.PlaceCellsUni = cell(n_iterations, 1);
results.PlaceCellsTanni = cell(n_iterations, 1);
results.cells_U = cell(n_iterations, 1);
results.cells_T = cell(n_iterations, 1);
results.gridness_square_U = cell(n_iterations, 1);
results.scale_square_U = cell(n_iterations, 1);
results.gridness_hex_U = cell(n_iterations, 1);
results.scale_hex_U = cell(n_iterations, 1);
results.grid_scales_binned_U = cell(n_iterations, 1);
results.gridness_square_T = cell(n_iterations, 1);
results.scale_square_T = cell(n_iterations, 1);
results.gridness_hex_T = cell(n_iterations, 1);
results.scale_hex_T = cell(n_iterations, 1);
results.grid_scales_binned_T = cell(n_iterations, 1);

for iter = 1:n_iterations
    % Specify the subfolder for this iteration
    subfolder = fullfile(output_dir, sprintf('iteration_%.1f', iter));
    
    % Load Place Cells
    place_cells_file = fullfile(subfolder, sprintf('place_cells_iter_%d.mat', iter));
    load(place_cells_file, 'PlaceCellsUni', 'PlaceCellsTanni');
    
    % Store in the structure
    results.PlaceCellsUni{iter} = PlaceCellsUni;
    results.PlaceCellsTanni{iter} = PlaceCellsTanni;
    
    % Load Grid and Place cells trained by SR model
    SR_cells_file = fullfile(subfolder, sprintf('grid_place_cells_SR_iter_%d.mat', iter));
    load(SR_cells_file, 'cells_U', 'cells_T');
    
    % Store in the structure
    results.cells_U{iter} = cells_U;
    results.cells_T{iter} = cells_T;
    
    % Load Uniform model results
    results_file_U = fullfile(subfolder, sprintf('U_results%d.mat', iter));
    load(results_file_U, 'square_cells_U', 'gridness_square_U', 'scale_square_U', 'hexagon_cells_U', 'gridness_hex_U', 'scale_hex_U');%, 'grid_scales_binned_U');
    
    % Store in the structure
    results.gridness_square_U{iter} = gridness_square_U;
    results.scale_square_U{iter} = scale_square_U;
    results.gridness_hex_U{iter} = gridness_hex_U;
    results.scale_hex_U{iter} = scale_hex_U;
    % results.grid_scales_binned_U{iter} = grid_scales_binned_U;
    
    % Load Sanders model results
    results_file_T = fullfile(subfolder, sprintf('T_results%d.mat', iter));
    load(results_file_T, 'square_cells_T', 'gridness_square_T', 'scale_square_T', 'hexagon_cells_T', 'gridness_hex_T', 'scale_hex_T');%, 'grid_scales_binned_T');
    
    % Store in the structure
    results.gridness_square_T{iter} = gridness_square_T;
    results.scale_square_T{iter} = scale_square_T;
    results.gridness_hex_T{iter} = gridness_hex_T;
    results.scale_hex_T{iter} = scale_hex_T;
    % results.grid_scales_binned_T{iter} = grid_scales_binned_T;
    disp(iter)
end

end