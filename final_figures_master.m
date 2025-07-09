%% [1] Load Place and Grid Cell Data

% load data
folder = '/Users/Elsa Marianelli/Documents/';
sub = 'grids_data_square_env_uni_traj';
folder = [folder sub '/'];
basename = 'data_covar_uniform_densityVaried_sizeVaried_widthControled';
nIterations = 5;

[Info, Place_cells, Grid_cells] = load_SR_or_covar_data(folder, basename, nIterations);

% options for metrics and plotting 
env_shape = 'rectangle'; % 'rectangle' OR 'trapezoid' 
binning_mode = 'three';  % if rectangle use 'three' OR if trapezoid use 'trapezoid_lr'
pc_size = false; % if size Varied^

% Create subfolder for saving figures
subfolder = fullfile([sub '_figs_square_uniform_traj'], basename);
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% [2a] Visualise place cell centre distribution
% trapezoid/rectangle
plot_place_cell_centres(Info, Place_cells, subfolder, 'rectangle');

%% [2b] see how place cell field widths and distances vary in environmental regions
plot_field_width(Place_cells, Info, subfolder, pc_size)

%% [3] Visualise Example Place Cells
figure;
for i = 1:10
    map = Place_cells{i*12}.fmap;
    subplot(2, 5, i); imagesc(map); colormap("jet"); hold on;
    % plot(boundary_coords(:,1), boundary_coords(:,2), 'w.');
    axis off; axis image;
end
saveas(gcf, fullfile(subfolder, 'Figure_B_place.svg'));

%% [4] Visualise Example Grid Cells
figure;
for i = 1:20
    map = Grid_cells{1}{i*10}.map;  
    bg_val = mode(map(:));                % Most frequent value
    map(map == bg_val) = NaN;             % Set background to 0
    subplot(4, 5, i); imagesc(rot90(map)); colormap("jet"); hold on;
    % plot(boundary_coords(:,1), boundary_coords(:,2), 'k.');
    axis off; axis image;
end
saveas(gcf, fullfile(subfolder, 'Figure_B_grid.svg'));

%% [5] Compute Grid Cell Metrics in Environment Bins
% removes any grids where one of the environmental bin cannot fit an
% ellipse for comprison between bins to be fair, also skips girds
% atomatically if not enough peaks are found for speed
metrics_all = compute_grid_metrics_binned(Grid_cells, binning_mode);
save(fullfile(subfolder, 'metrics.mat'), 'metrics_all');

%% [6] Plots and saves grid metrics and respective variability of each environemntel region from mean 
show_variability = false;
plot_grid_metrics(metrics_all, binning_mode, subfolder, show_variability);         

%% [] add heat map plots!

%% [7] Visualise PCs ascribed to each gridness metrics
threshold = .5;
plot_gridness_by_region(metrics_all, threshold, subfolder)        