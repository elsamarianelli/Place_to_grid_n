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
name =    'covar_hasselmo_densityConstant_sizeVaried';
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
% saveas(gcf, fullfile(subfolder, 'Figure_A_place_centres.fig'));

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
for i = 1:10
    map = Grid_cells{1}{i*3}.map;
    subplot(2, 5, i); imagesc(map)
    axis off; axis image      
end

% Save the figure
% saveas(gcf, fullfile(subfolder, 'Figure_B_grid.fig'));

%% [C] evaluate scale and elipse in different environement bins
% SETTINGS 
shape       = 'hexagon';
nIterations = numel(Grid_cells);
nPCs        = numel(Grid_cells{1});
% initialise
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

%%
% Compute means & stds
mean_three = squeeze(mean(mean(allVals_three, 2, 'omitnan'), 1, 'omitnan'));
std_three = squeeze(std(mean(allVals_three, 2, 'omitnan'), 0, 1, 'omitnan'));

% Plot
figure;
bar(mean_three, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'k'); hold on;
errorbar(1:3, mean_three, std_three, 'k.', 'LineWidth', 1.5);

xticks(1:3); xticklabels({'Corners', 'Edges', 'Center'});
ylabel('Grid Scale');
title('Grid Scale by Region');

% Save the figure and vals
saveas(gcf, fullfile(subfolder, 'Figure_C_BarChart.fig'));
save(fullfile(subfolder, 'all_binned_metrics.mat'), 'allVals', 'means', 'stds', 'metric_name', 'shape', 'bin_struct');

%% [D]Heat map
% Compute mean over iterations and PCs
mean_3x3 = reshape(mean(allVals_nine, [1,2], 'omitnan'), 3, 3);

figure;
imagesc(mean_3x3);
colorbar;
title('Mean Grid Scale (3x3 Environment Bins)');
axis off;

saveas(gcf, fullfile(subfolder, 'Figure_D_HeatMap.fig'));
%% ellipse fitting



