% [1] load in data
% 1)	Method – example firing patterns, illustration of the different 
% methods of dimensionality reduction
% 2)	General parameter sweeps for the simplest version (covar) - 
% specifically changes in place field size with distance to the wall, 
% changes in place field density with distance to the wall, 
% changes in place field size and density with distance to the wall. 
% Plot grid scale and ellipticity in nine areas of the environment, 
% bar charts of those parameters in edge, corner, central areas
% 3)     Then look at the same outputs but for square vs trapezoidal environments - 
% in this case just comparing grid field statistics betwen the large and small end 
% of the trapezoid (as in Fig 4 of Krupic et al., 2015)

%% Figure 1 - How does changing density and size parameters of Place cells affect Grids?

%% First, load the data
folder = '/Users/elsamarianelli/Documents/grids_data/';
base = 'data_';
name =    'covar_hasselmo_densityVaried_sizeVaried';
basename = [base, name];
nIterations =4;

[Info, Place_cells, Grid_cells] = load_SR_or_covar_data(folder, basename, nIterations);

% initiate folder to save figuyres
subfolder = fullfile('grids_figures', name);  % or any folder name

% Create subfolder if it doesn't exist
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end
%% [A] Visualise how distribution of place cell centres vary
% Extract all fmaps into a cell array
fmaps = cellfun(@(pc) pc.fmap, Place_cells, 'UniformOutput', false);

% Use cellfun to get linear index of the max in each fmap
linear_idx = cellfun(@(f) find(f == max(f(:)), 1), fmaps);  % take first if multiple
[rows, cols] = cellfun(@(f, idx) ind2sub(size(f), idx), fmaps, num2cell(linear_idx));

figure; hold on; plot(cols, rows, 'o', 'MarkerSize', 5, 'Color', 'k', 'MarkerFaceColor','k');

% Get environemnt bin lines 
xlim tight; ylim tight; xLimits = xlim; yLimits = ylim;
x_thirds = xLimits(1) + diff(xLimits) * [1/3, 2/3];
y_thirds = yLimits(1) + diff(yLimits) * [1/3, 2/3];

% Add lines at 1/3 and 2/3
for x = x_thirds
    xline(x, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end
for y = y_thirds
    yline(y, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end

% Clean appearance
axis equal; set(gca, 'Box', 'on', 'XTick', [], 'YTick', []);

% Save the figure
saveas(gcf, fullfile(subfolder, 'Figure_A_place_centres.fig'));

%% [B] Visualise how size of place cells vary 
figure; 
for i = 1:10
    map = Place_cells{i*6}.fmap;
    subplot(2, 5, i); imagesc(map)
    axis off; axis image      
end

% Save the figure
saveas(gcf, fullfile(subfolder, 'Figure_B_place.fig'));

%% [b]2 Visualise some grids
figure; 
for i = 1:40
    map = Grid_cells{1}{i*4}.map;
    subplot(4, 10, i); imagesc(map); colormap("parula")
    axis off; axis image      
end

% Save the figure
saveas(gcf, fullfile(subfolder, 'Figure_B_grid.fig'));

%% [C] evaluate scale and elipse in different environement bins

% SETTINGS 
shape       = 'hexagon';
nIterations = numel(Grid_cells);
nPCs        = numel(Grid_cells{1});

% initialise storage
metrics_all = cell(nIterations, nPCs);  % Holds per-PC metric structs

% loop to get metrics for each bin of environment( and also in corner,
% side, centre configuration for later plotting)
for iter = 1:nIterations
    for pc = 1:nPCs
        try
            map = Grid_cells{iter}{pc}.map;
            [~, binned_three, binned_nine] = get_binned_metrics(map, shape);

            % Store all bin data (3-bin and 9-bin)
            metrics_all{iter, pc}.three = binned_three;
            metrics_all{iter, pc}.nine  = binned_nine;

        catch ME
            warning("Error iter %d, PC %d: %s", iter, pc, ME.message);
            metrics_all{iter, pc} = [];  % Still indexable
        end
    end
end

%%  Plot BAR CHARTS of mean ± std 
metrics = {'expGrd_h_el', 'scale_h_el', 'eccent', ...
    'orient', 'ellipicity'};

for m=1:numel(metrics)
    metric = metrics{m};
    nIterations = 4;        % Number of iterations
    nPCs =250;            % Number of cells per iteration
    region_count = 3;                        % Corners, Edges, Center
    
    % Initialize matrix: [nIterations x 3]
    mean_vals = nan(nIterations, region_count);
    
    for iter = 1:nIterations
        vals = nan(nPCs, region_count);
        for pc = 1:nPCs         
           vals(pc, :) =arrayfun(@(s) s.(metric), metrics_all{iter, pc}.three);
        end
        % Mean across cells for each region (omit NaNs)
        mean_vals(iter, :) = mean(vals, 1, 'omitnan');
    end
    
    % === Compute means and stds across iterations ===
    grand_mean = mean(mean_vals, 1, 'omitnan');
    grand_std  = std(mean_vals, 0, 1, 'omitnan');
    
    % === Plot ===
    figure; hold on;
    bar_x = 1:region_count;
    b = bar(bar_x, grand_mean, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    % Add error bars
    er = errorbar(bar_x, grand_mean, grand_std, 'k', ...
        'LineStyle', 'none', 'LineWidth', 1.8);
    
    % Formatting
    set(gca, 'XTick', bar_x, 'XTickLabel', {'Corners', 'Edges', 'Center'}, 'FontSize', 16);
    ylabel(strrep(metric, '_', '\_'), 'FontSize', 18);
    xlabel('Environment Region', 'FontSize', 18);
    title(['Mean ' strrep(metric, '_', '\_') ' Across Regions'], 'FontSize', 20);
    box off;
    % Optional: Save
    % saveas(gcf, fullfile(subfolder, ['Barchart_' metric '.fig']));

end

%% PLOT HEAT MAPS
% === Settings ===
metric_names = metrics;
nBins = 9;
nMetrics = numel(metric_names);
bin_dims = [3, 3];  % For reshaping heatmaps

% === Loop through metrics ===
for m = 1:nMetrics
    metric = metric_names{m};
    all_vals = nan(nIterations, nPCs, nBins);  % Preallocate

    for iter = 1:nIterations
        for pc = 1:nPCs
            try
                bin_data = metrics_all{iter, pc}.nine;
                for b = 1:nBins
                    val = bin_data(b).(metric);
                    if isnumeric(val) && isscalar(val) && ~isnan(val)
                        all_vals(iter, pc, b) = val;
                    end
                end
            catch
                continue;
            end
        end
    end

    % === Average across PCs and iterations ===
    mean_vals = squeeze(mean(mean(all_vals, 2, 'omitnan'), 1));  % 1 x 9
    heat_data = reshape(mean_vals, bin_dims);  % 3x3 grid

    % === Plot ===
    figure('Color','w');
    h = heatmap(heat_data, ...
        'CellLabelColor', 'none', ...
        'Colormap', parula, ...
        'ColorbarVisible', 'on');
    title(metric)
    h.GridVisible = 'off';
    h.FontSize = 14;

    % Optional: Save
    saveas(gcf, fullfile(subfolder, ['Heatmap_' metric '.fig']));

end