%% Generate and save Grid maps 
%  generated from doing PCA on SR Matrix learnt through TDLM during random foraging with 
%  place cell firing as successor features. 
addpath('/Users/elsamarianelli/Documents/GitHub')
addpath('/Users/elsamarianelli/Documents/GitHub/bound_warped_grids_new/Algorithmic_PCA/Functions/')
addpath('/Users/elsamarianelli/Documents/GitHub/bound_warped_grids_new/General_functions/')
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Algorithmic_PCA\Functions'

% Define environment dimensions and the number of cells
dim_x = 351;
dim_y = 252;
n_polys = 1;
polys = cell(n_polys, 1);
polys{1} = [0 0, 349 0, 349 250, 0 250, 0 0] + 2; % rectangular, change for warped trapezoid
% warp = 20
% polys{1} = [0 0, 349 0+warp, 349 250-warp, 0 250, 0 0] + 2
n_cells = 200;
n_steps = 360000;
output_dir = 'results_grid_maps_rect_2'; 
n_iterations = 10; 

% alternative is to use real trajectories and boundaries from rodent behavioural data...
data_folder = 'C:\Users\Elsa Marianelli\Desktop\Rodent_Data';
[trajectories, bounds] = get_real_traj(data_folder);
% plot a sample trajectory 
traj = trajectories{2}; 
plot(traj(:,1), traj(:,2), '.');

% Create the directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for iter = 1:n_iterations

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
    fprintf('  Step 2/8: Generating Trajectory\n'); % hasselmo version
    traj = generate_trajectory(env, n_steps);
    
    % [] alternatively take real trajectory and use that environment ^^
    % [3] Populating Place Cells
    fprintf('  Step 3/8: Populating Place Cells\n');
    [PlaceCellsUni, PlaceCellsTanni, env] = generate_place_cells(env, n_cells, dim_x, dim_y, 2, 'random');
    for i = 1:5:250
        figure; 
        imagesc(PlaceCellsTanni{i}.fmap); title(i)
    end
    % Save place cells to the subfolder
    place_cells_file = fullfile(subfolder, sprintf('orig_place_cells'));
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
    cells_U = getGrid(cells_U, M_U, env, 'off'); % only have as off 
    cells_T = getPlace(PlaceCellsTanni, M_T, env);
    cells_T = getGrid(cells_T, M_T, env, 'off'); 

    % [6] Saving Grid and Place cells
    fprintf('Step 6/8: Save\n')
    % Save grid cells to the subfolder
    SR_cells_file = fullfile(subfolder, sprintf('grid_place_cells_SR'));
    save(SR_cells_file, 'cells_U', 'cells_T');
         

end

%% plotting figure 1 - Place cells in environment, trajectory, example place cells, and trained SR matrics 
% all with uniform place cells
map = env.map;
% Create a figure
figure;

% Define font sizes for better readability
titleFontSize = 14;
axisFontSize = 12;
scatterSize = 30;

% Subplot 1: Uniform Place Cells in the Environment with Boundaries
ax1 = subplot(2,2,1);
imagesc(~map);       % Display the environment, flipping black/white for visibility
colormap(ax1, gray); % Set grayscale colormap for this subplot only
hold on
[xy_field_u, env] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, 1000);
scatter(xy_field_u(:,1), xy_field_u(:,2), scatterSize, 'k', 'filled');
hold off
title('Place Cell Locations (n = 200)', 'FontSize', titleFontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', axisFontSize);

% Subplot 2: Trajectory in the Environment Boundaries
ax2 = subplot(2,2,2);
imagesc(~map);       
colormap(ax2, gray); 
hold on
plot(traj(:,1), traj(:,2), 'k', 'LineWidth', .2);
hold off;
title('Simulated Trajectory (2 Hours)', 'FontSize', titleFontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', axisFontSize);

% Subplot 3: A Single Place Cell Firing Field (Randomly Picked)
ax3 = subplot(2,2,3);
imagesc(PlaceCellsUni{50}.fmap);
colormap(ax3, jet); 
title('Example Place Cell Firing Field', 'FontSize', titleFontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', axisFontSize);

% Subplot 4: Trained SR Matrix for Uniform Cells
ax4 = subplot(2,2,4);
imagesc(M_U); 
colormap(ax4, parula); 
colorbar;
title('TDLR trained SR Matrix', 'FontSize', titleFontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', axisFontSize);

% Set global figure properties
set(gcf, 'Color', 'w'); 

%% plotting figure 2 - unoform vs boundary warped place fields, 
%% Plotting figure: Place cells, trajectory, example place cells, and SR matrix
map = env.map;
figure;

% Define font sizes for better readability
titleFontSize = 14;
axisFontSize = 12;
scatterSize = 30;

% --------- Row 1: Uniform Place Cells --------- %

% Subplot (1,1): Uniform Place Cells in the Environment
ax1 = subplot(2,5,1);
imagesc(~map);       
colormap(ax1, gray);
hold on;
[xy_field_u, env] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, 1000);
scatter(xy_field_u(:,2), xy_field_u(:,1), scatterSize, 'k', 'filled');hold off;

title('Uniform Place Cells (n=200)', 'FontSize', titleFontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', axisFontSize);

% Subplot (1,2-1,4): Three Example Uniform Place Cell Firing Fields
example_indices = [161, 136, 146, 76]; % Three example place cells
for i = 1:4
    ax = subplot(2,5,i+1);
    imagesc(PlaceCellsUni{example_indices(i)}.fmap);
    colormap(ax, jet); % Use jet (rainbow) for firing fields
    title(sprintf('Uniform Example %d', i), 'FontSize', titleFontSize, 'FontWeight', 'bold');
    set(gca, 'FontSize', axisFontSize);
end

% --------- Row 2: Tanni Place Cells --------- %

% Subplot (2,1): Tanni Place Cells in the Environment
ax6 = subplot(2,5,6);
imagesc(~map);
colormap(ax6, gray);
hold on;
[xy_field_t, env] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, 1.7);
scatter(xy_field_t(:,2), xy_field_t(:,1), scatterSize, 'k', 'filled');hold off;

title('BW Place Cells (n=200)', 'FontSize', titleFontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', axisFontSize);

% Subplot (2,2-2,4): Three Example Tanni Place Cell Firing Fields
for i = 1:4
    ax = subplot(2,5,i+6);
    imagesc(PlaceCellsTanni{example_indices(i)}.fmap);
    colormap(ax, jet);
    title(sprintf('BW Example %d', i), 'FontSize', titleFontSize, 'FontWeight', 'bold');
    set(gca, 'FontSize', axisFontSize);
end

% Set global figure properties
sgtitle('Comparison of Uniform and Tanni Place Cells', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'Color', 'w'); % Set background to white for better contrast
