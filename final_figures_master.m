%% [1] Load Place and Grid Cell Data

% load data
folder = '/Users/Elsa Marianelli/Documents/';
sub = 'grids_data_trap_env_uni_traj';
folder = [folder sub '/'];
basename = 'data_covar_uniform_densityVaried_sizeConstant_WidthControled';
nIterations = 5;

[Info, Place_cells, Grid_cells] = load_SR_or_covar_data(folder, basename, nIterations);

% options for metrics and plotting 
env_shape = 'trapezoid'; % 'rectangle' OR 'trapezoid' 
binning_mode = 'trapezoid_lr';  % if rectangle use 'three' OR if trapezoid use 'trapezoid_lr'
pc_size = true; % if size Varied^

% Create subfolder for saving figures
subfolder = fullfile([sub '_trapezoid_figs'], basename);
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% [2a] Visualise place cell centre distribution
% trapezoid/rectangle
plot_place_cell_centres(Info, Place_cells, subfolder, env_shape);

%% [2b] see how place cell field widths and distances vary in environmental regions
plot_field_width(Place_cells, Info, subfolder, pc_size)
plot_field_width_trap_lr(Place_cells_all, In, subfolder, pc_size, true)

%% [3] Visualise Example Place Cells

if strcmp(env_shape, 'trapezoid')
    % get region outside of trap for plotting
    [x, y] = find(~ (rot90(Info.env.dwmap,-1) >=0)); 
end

figure;
for i = 1:10
    map = Place_cells{i*17}.fmap;
    subplot(2, 5, i); imagesc(map); hold on;
    plot(x, y, 'w.');
    axis off; axis image;
end
saveas(gcf, fullfile(subfolder, 'Figure_B_place.svg'));

%% [4] Visualise Example Grid Cells
figure;
for i = 1:20
    map = Grid_cells{1}{i*10}.map;  
    % map = cells{i*10}.grid;
    bg_val = mode(map(:));                % Most frequent value
    map(map == bg_val) = NaN;             % Set background to 0
    subplot(4, 5, i); imagesc(rot90(map)); hold on;
    plot(x, y, 'w.');
    axis off; axis image;
end
saveas(gcf, fullfile(subfolder, 'Figure_B_grid.svg'));

%% [4B] Example of ellipse fitted to an example grid cell for each bin
for ind= 70:5:80
    map = Grid_cells{1}{ind}.map;
    visualise_ellipses(map, 'trapezoid')
    figure; imagesc(rot90(map));hold on;     plot(x, y, 'w.');
end
saveas(gcf, fullfile(subfolder, ['Example_binned_grid_with_ellipses' num2str(ind) '.svg']));
saveas(gcf, fullfile(subfolder, ['Example_binned_grid' num2str(ind) '.svg']));

close all

%% [5] Compute Grid Cell Metrics in Environment Bins
% removes any grids where one of the environmental bin cannot fit an
% ellipse for comprison between bins to be fair, also skips girds
% atomatically if not enough peaks are found for speed
binning_mode = 'trapezoid_lr';
metrics_all = compute_grid_metrics_binned(Grid_cells, binning_mode);
save(fullfile(subfolder, 'metrics.mat'), 'metrics_all');

%% [6] Plots and saves grid metrics and respective variability of each environemntel region from mean
%  bar plots, heatmaps, and metrics table with comparisons
plot_grid_metrics(metrics_all,subfolder)
plot_grid_metrics_halves(metrics_all, subfolder)

%% [7] Visualise PCs ascribed to each gridness metrics
threshold = .5;
plot_gridness_by_region(metrics_all, threshold, subfolder) 
