%% [1] Load Place and Grid Cell Data

% load data
folder = '/Users/elsamarianelli/Documents/data_grids/';
sub = 'grids_data_square_env_uni_traj';

folder = [folder sub '/'];
basename = 'data_covar_uniform_densityConstant_sizeConstant_WidthControled';
% basename = 'data_SR_generate_densityConstant_sizeConstant_WidthControled';

nIterations = 5;

[Info, Place_cells, Grid_cells] = load_SR_or_covar_data(folder, basename, nIterations);

% options for metrics and plotting 
env_shape = 'rectangle'; % 'rectangle' OR 'trapezoid' 
binning_mode = 'three';  % if rectangle use 'three' OR if trapezoid use 'trapezoid_lr'
pc_size = false; % if size Varied^

% Create subfolder for saving figures
subfolder = fullfile([sub '_new_data'], basename);
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% [2a] Visualise place cell centre distribution
% trapezoid/rectangle
plot_place_cell_centres(Info, Place_cells, subfolder, env_shape);

%% [2b] see how place cell field widths and distances vary in environmental regions
% plot_field_width(Place_cells, Info, subfolder, pc_size)
% plot_field_width_trap_lr(Place_cells, Info, subfolder, pc_size, true)

%% [3] Visualise Example Place Cells

if strcmp(env_shape, 'trapezoid')
    % get region outside of trap for plotting
    [x, y] = find(~ (rot90(Info.env.dwmap,-1) >=0)); 
end

figure;
for i = 1:10
    map = Place_cells{i*17}.fmap;
    subplot(2, 5, i); imagesc(map); hold on;
    % plot(x, y, 'w.');
    axis off; axis image;
end
saveas(gcf, fullfile(subfolder, 'Figure_B_place.svg'));

%% [4] Visualise Example Grid Cells
figure;
for i = 1:20
    map = Grid_cells{1}{i*12}.map;  
    % map = cells{i*10}.grid;
    bg_val = mode(map(:));                % Most frequent value
    map(map == bg_val) = NaN;             % Set background to 0
    subplot(4, 5, i); imagesc(rot90(map)); hold on;
    % plot(x, y, 'w.');
    axis off; axis image;
end
saveas(gcf, fullfile(subfolder, 'Figure_B_grid.svg'));

%% [4B] Example of ellipse fitted to an example grid cell for each bin
for ind= 24:1:250

    figure; imagesc(rot90(map));
    % hold on;  plot(x, y, 'w.');
    % figure; imagesc(map);
    map = Grid_cells{1}{ind}.map;
    figure; visualise_ellipses(map, env_shape);

end
saveas(gcf, fullfile(subfolder, ['Example_binned_grid_with_ellipses' num2str(ind) '.svg']));
saveas(gcf, fullfile(subfolder, ['Example_binned_grid' num2str(ind) '.svg']));

close all

%% [5] Compute Grid Cell Metrics in Environment Bins
% removes any grids where one of the environmental bin cannot fit an
% ellipse for comprison between bins to be fair, also skips girds
% atomatically if not enough peaks are found for speed
metrics_all = compute_grid_metrics_binned(Grid_cells, binning_mode);
save(fullfile(subfolder, 'metrics.mat'), 'metrics_all');

%% [6] Plots and saves grid metrics and respective variability of each environemntel region from mean
%  bar plots, heatmaps, and metrics table with comparisons
plot_grid_metrics(metrics_all,subfolder)
% plot_grid_metrics_halves(metrics_all, subfolder)

%% [7] Visualise PCs ascribed to each gridness metrics
threshold = .5;
plot_gridness_by_region(metrics_all, threshold, subfolder) 

%% [8] additional plots with reloaded metrics 
%% rewrite metric tables correctly 

% for each subfolder with basename from list...
basenames = {'data_covar_uniform_densityConstant_sizeConstant_WidthControled';...
'data_covar_uniform_densityVaried_sizeVaried_WidthControled';...
'data_covar_uniform_densityConstant_sizeVaried_WidthControled';...
'data_covar_uniform_densityVaried_sizeConstant_WidthControled'};

% basenames = {'data_SR_generate_densityConstant_sizeConstant_WidthControled';...
% 'data_SR_generate_densityVaried_sizeVaried_WidthControled';...
% 'data_SR_generate_densityConstant_sizeVaried_WidthControled';...
% 'data_SR_generate_densityVaried_sizeConstant_WidthControled'};
folder = '/Users/elsamarianelli/Documents/data_grids/grids_data_square_env_uni_traj_plots_and_comps';
% folder = '/Users/Elsa Marianelli/Documents/';
% sub = 'grids_data_square_env_uni_traj_plots_and_comps';
sub = 'grids_data_trap_env_uni_traj_trapezoid_figs'; 
folder = [folder sub '/'];
nIterations = 5;

for b = 1:length(basenames)
    filename = fullfile(folder, basenames{b}, 'metrics.mat');
    metrics_all = load(filename);
    metrics_all = metrics_all.metrics_all;
    plot_grid_metrics(metrics_all,fullfile(folder, basenames{b}))
    % corner_scales = nan(size(metrics_all));
    % for i = 1:numel(metrics_all)
    %     m = metrics_all{i};
    %     if ~isempty(m) && isstruct(m) && numel(m.three) >= 1
    %         corner_scales(i) = m.three(1).scale_h;
    %     end
    % end
    % el_fit  = find(~isnan(corner_scales));
    % disp(el_fit)
end

%% comb 
% === SETTINGS ===
sub = 'grids_data_trap_env_uni_traj_trapezoid_figs';

folder = [folder sub '/'];

basenames = {'data_covar_uniform_densityConstant_sizeConstant_WidthControled';...
             'data_covar_uniform_densityVaried_sizeVaried_WidthControled';...
             'data_covar_uniform_densityConstant_sizeVaried_WidthControled';...
             'data_covar_uniform_densityVaried_sizeConstant_WidthControled'};

outFile = fullfile(folder,'All_Metrics_comb.xlsx');

% List the metrics you want to keep (leave empty to keep all)
metricsToInclude = {};  % e.g. {'Gridness','Scale'}
metric_titles = {'Gridness hex','Gridness hex no el','Scale hex','Eccentricity','Ellipse Orientation','Peak Orientation'};

% === GATHER ALL CSVs ===
AllMeans = [];
AllComps = [];
allMetrics = strings(0,1);

for b = 1:numel(basenames)
    subfolder = fullfile(folder, basenames{b});
    meansPath = fullfile(subfolder,'Metric_Means_by_Region_f.csv');
    compsPath = fullfile(subfolder,'Metric_Region_Comparisons_paired_f.csv');

    if isfile(meansPath)
        Tm = readtable(meansPath,'TextType','string');
        Tm.Dataset = repmat(string(basenames{b}),height(Tm),1);
        AllMeans = [AllMeans; Tm];
        allMetrics = [allMetrics; Tm.Metric];
    end
    if isfile(compsPath)
        Tc = readtable(compsPath,'TextType','string');
        Tc.Dataset = repmat(string(basenames{b}),height(Tc),1);
        AllComps = [AllComps; Tc];
        allMetrics = [allMetrics; Tc.Metric];
    end
end

% === FILTER METRICS ===
if isempty(metricsToInclude)
    metricsSel = unique(allMetrics);
else
    metricsSel = intersect(string(metricsToInclude), unique(allMetrics));
end

if ~isempty(AllMeans)
    AllMeans = AllMeans(ismember(AllMeans.Metric,metricsSel),:);
end
if ~isempty(AllComps)
    AllComps = AllComps(ismember(AllComps.Metric,metricsSel),:);
end

% === WRITE TO EXCEL ===
if ~isempty(AllMeans)
    writetable(AllMeans, outFile, 'Sheet','Means_by_Region');
end
if ~isempty(AllComps)
    writetable(AllComps, outFile, 'Sheet','Region_Comparisons');
end

