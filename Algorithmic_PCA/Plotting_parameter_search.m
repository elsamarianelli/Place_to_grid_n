%% Plotting Parameter search output from covariance_eigendecomp.m 
% Shows how different numbers of trajectory steps affect different gridness metrics.
addpath param_sweeps
% (1) SETTINGS
baseNames = { ...
    'param_sweeps\SR_Covar_check_SR_Uniform', ...
    'param_sweeps\SR_Covar_check_SR_Tanni',...
    'param_sweeps\SR_Covar_check_covar_Tanni',...
    'param_sweeps\SR_Covar_check_covar_Uniform',...
    'param_sweeps\SR_Covar_check_covar_uni_trajectory_Tanni',...
    'param_sweeps\SR_Covar_check_covar_uni_trajectory_Uniform'};

labels = { 'SR Uniform','SR Tanni', 'Covar Tanni', 'Covar Uniform', 'Covar no Traj Tanni', 'Covar no Traj Uniform'};

nIterations  = 5;                                        % Number of iterations for each parameter
metricNames  = {'expGrd_h', 'expGrd_s', 'stGrd_h', 'stGrd_s'}; % List of metrics to plot
threshold    = 0.8;                                      % Threshold for "gridness"
nPCs         = 250;                                      % Number of principal components
parameter    = 'Place cell type, and PCA type';               % for plotting

% (2) SETUP
colors = lines(numel(metricNames));       

% Count base folders and metrics
nMetrics = numel(metricNames);
nBases   = numel(baseNames);

% Preallocate result arrays (per base)
meanProps    = nan(nMetrics, nBases);
stdProps     = nan(nMetrics, nBases);
meanScale_h  = nan(nBases, 1);
stdScale_h   = nan(nBases, 1);
meanScale_s  = nan(nBases, 1);
stdScale_s   = nan(nBases, 1);

% Loop over base directories first
for s = 1:nBases
    baseDir = baseNames{s};
    fprintf('Processing baseDir: %s\n', baseDir);

    allProps = nan(nMetrics, nIterations);
    allScale_h = nan(1, nIterations);
    allScale_s = nan(1, nIterations);

    for iter = 1:nIterations
        fprintf('  Iteration: %.1f\n', iter);

        % Load the data once per iteration
        fname = fullfile(baseDir, sprintf('iteration_%.1f', iter), ...
            'metrics_and_maps.mat');

        if isfile(fname)
            data = load(fname);
            % Remove NaNs from GC_metrics
            is_valid = ~cellfun(@(x) isequaln(x, NaN), data.GC_metrics);
            data.GC_metrics = data.GC_metrics(is_valid);

            for m = 1:nMetrics
                metricName = metricNames{m};
                fprintf('    Metric: %s\n', metricName);

                try
                    metric_vals = cellfun(@(x) x.(metricName), data.GC_metrics, 'UniformOutput', false);
                    metric_vals = cell2mat(metric_vals);
                    allProps(m, iter) = sum(metric_vals > threshold) / numel(metric_vals);
                catch
                    warning('    Failed to extract metric %s at %s (iter %.1f)', metricName, baseDir, iter);
                end
            end

            % Extract scale fields (optional: only if they're guaranteed to be present)
            try
                scale_h_vals = cell2mat(cellfun(@(x) x.scale_h, data.GC_metrics, 'UniformOutput', false));
                scale_s_vals = cell2mat(cellfun(@(x) x.scale_s, data.GC_metrics, 'UniformOutput', false));

                allScale_h(iter) = mean(scale_h_vals, 'omitnan');
                allScale_s(iter) = mean(scale_s_vals, 'omitnan');
            catch
                warning('    Failed to extract scale_h or scale_s at %s (iter %.1f)', baseDir, iter);
            end

        else
            warning('Missing file: %s', fname);
        end
    end

    % Compute stats over iterations
    meanProps(:, s) = mean(allProps, 2, 'omitnan');
    stdProps(:, s)  = std(allProps, 0, 2, 'omitnan');
    meanScale_h(s) = mean(allScale_h, 'omitnan');
    stdScale_h(s)  = std(allScale_h, 0, 'omitnan');
    meanScale_s(s) = mean(allScale_s, 'omitnan');
    stdScale_s(s)  = std(allScale_s, 0, 'omitnan');
end

%% (1) PLOT: Gridness and Field Scale Comparison
figure('Color', 'w', 'Position', [100, 100, 1300, 600]);

% === FORMATTING ===
colors = lines(numel(metricNames));
markers = {'o', 's', 'd', '^'};
xvals = 1:nBases;
xtick_labels = labels;

% === CLEAN LEGEND LABELS ===
% Converts: 'expGrd_h' -> 'Expanded Hexagonal'
cleanNames = @(str) regexprep(str, ...
    {'^expGrd_', '^stGrd_', '_h$', '_s$'}, ...
    {'Expanded ', 'Standard ', 'Hexagonal', 'Square'});

displayNames = cellfun(cleanNames, metricNames, 'UniformOutput', false);

% === SUBPLOT 1: Gridness Metrics ===
subplot(1, 2, 1); hold on;
for m = 1:nMetrics
    errorbar(xvals, meanProps(m, :), stdProps(m, :), ...
        '-o', ...
        'Color', colors(m, :), ...
        'Marker', markers{mod(m - 1, numel(markers)) + 1}, ...
        'MarkerFaceColor', colors(m, :), ...
        'LineWidth', 2, ...
        'MarkerSize', 10, ...
        'DisplayName', displayNames{m});
end

xlim([0.5, nBases + 0.5]);  % Add padding on edges
xticks(xvals);
xticklabels(xtick_labels);
xtickangle(20);

xlabel(parameter, 'FontSize', 18, 'FontWeight', 'bold');
ylabel(sprintf('%% of PCs with Gridness > %.2f', threshold), 'FontSize', 18, 'FontWeight', 'bold');
title('Gridness Metrics by Condition', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 14);
set(gca, 'FontSize', 16, 'LineWidth', 1.5, ...
    'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');
grid on;

% === SUBPLOT 2: Field Scale Comparison ===
subplot(1, 2, 2); hold on;

scaleNames = {'Hexagonal Field Scale', 'Square Field Scale'};
scaleMeans = [meanScale_h'; meanScale_s'];
scaleStds = [stdScale_h'; stdScale_s'];

for m = 1:2
    errorbar(xvals, scaleMeans(m, :), scaleStds(m, :), ...
        '-o', ...
        'Color', colors(m, :), ...
        'Marker', markers{m}, ...
        'MarkerFaceColor', colors(m, :), ...
        'LineWidth', 2, ...
        'MarkerSize', 10, ...
        'DisplayName', scaleNames{m});
end

xlim([0.5, nBases + 0.5]);
xticks(xvals);
xticklabels(xtick_labels);
xtickangle(20);

xlabel(parameter, 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Mean Field Scale (pixels)', 'FontSize', 18, 'FontWeight', 'bold');
title('Field Scale Comparison by Condition', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 14);
set(gca, 'FontSize', 16, 'LineWidth', 1.5, ...
    'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');
grid on;

%% (2) Visualise Gridness Scores > Threshold for Each PC
baseName = 'param_sweeps\Tanni_Covar_ED_boundaryEffect_2.5';
nIterations = 5;
nPCs = 250;
threshold = 0.8;

% Call external visualisation function
gridness_visualisation_allPC(baseName, nIterations, nPCs, threshold);

% (2) Combine .fig Panels into One Layout
figFiles = {'bound_effect.fig', 'number_place_cells.fig', ...
            'traj_length.fig', 'random_vs_arrayed.fig'};
figAxes = cell(size(figFiles));

% Open and collect axes
for i = 1:numel(figFiles)
    f = openfig(figFiles{i}, 'invisible');
    axs = flipud(findall(f, 'Type', 'axes'));  % Flip to preserve subplot order
    figAxes{i} = axs;
    close(f);  % Immediately close after use
end

% Create new combined figure
figure;
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Copy axes content
for row = 1:4
    for col = 1:2
        tileIdx = (row - 1) * 2 + col;
        nexttile(tileIdx);
        ax = figAxes{row}(col);
        copyobj(ax.Children, gca);
        title('');
        xlabel(ax.XLabel.String);
        ylabel(ax.YLabel.String);
        set(gca, 'FontSize', 12);
        if row == 1
            legend(findobj(gca, '-regexp', 'Type', 'line|errorbar'), 'Location', 'best');
        end
    end
end

% Add column titles
colTitles = {'Gridness', 'Scale'};
for i = 1:2
    annotation(gcf, 'textbox', [0.22 + 0.48*(i-1), 0.965, 0.2, 0.03], ...
        'String', colTitles{i}, 'EdgeColor', 'none', ...
        'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Consistent Y labels
customYLabels = {'% PCs > Threshold', 'Field Scale'};
for i = 1:4
    nexttile((i - 1) * 2 + 1); ylabel(customYLabels{1}, 'FontSize', 12);
    nexttile((i - 1) * 2 + 2); ylabel(customYLabels{2}, 'FontSize', 12);
end

% Custom ticks on last row
tick_positions = [0, 4.5];
tick_labels = {'Random', 'Arrayed'};
nexttile(7); xticks(tick_positions); xticklabels(tick_labels);
nexttile(8); xticks(tick_positions); xticklabels(tick_labels);

% Specific axis label overrides
nexttile(5); xlabel('Number of Steps');
nexttile(6); xlabel('Number of Steps');


%% (3) Binned Parameter Sweep Plot

% === Settings ===
settings = baseNames;                       % cell array of base folder names
nIterations = 5;
nPCs = 250;
bin_types = 3;
metrics_to_extract = {'expGrd_h', 'scale_h'};  % metrics to collect
nSettings = numel(settings);

% === Preallocate storage ===
all_means = struct(); all_stds = struct();
for m = 1:numel(metrics_to_extract)
    metric = metrics_to_extract{m};
    all_means.(metric) = nan(nSettings, bin_types);
    all_stds.(metric)  = nan(nSettings, bin_types);
end

% === Loop over parameter settings ===
for s = 1:nSettings
    settingLabel = settings{s};
    fprintf('\n=== Processing setting %d/%d: %s ===\n', s, nSettings, settingLabel);

    % Temporary containers for this setting
    allVals = struct();
    for m = 1:numel(metrics_to_extract)
        allVals.(metrics_to_extract{m}) = nan(nIterations, nPCs, bin_types);
    end

    for iter = 1:nIterations
        fprintf(' - Iteration %d/%d\n', iter, nIterations);
        fname = fullfile(settingLabel, sprintf('iteration_%.1f', iter), 'metrics_and_maps.mat');
        if ~isfile(fname)
            warning('Missing file: %s', fname);
            continue;
        end
        data = load(fname);

        for pc = 1:nPCs
            try
                map = data.GC_metrics{pc}.map;
                [~, binned] = get_binned_metrics(map, 'hexagon', 'three');

                % Store each selected metric
                for m = 1:numel(metrics_to_extract)
                    for b = 1:bin_types
                        val = binned(b).(metrics_to_extract{m});
                        allVals.(metrics_to_extract{m})(iter, pc, b) = val;
                    end
                end
            catch
                warning('PC %d failed at iter %d, setting %s', pc, iter, settingLabel);
                continue;
            end
        end
    end

    % === Compute and store statistics ===
    for m = 1:numel(metrics_to_extract)
        metric = metrics_to_extract{m};
        metric_vals = allVals.(metric);

        mean_iter = squeeze(mean(metric_vals, 2, 'omitnan'));  % [nIter x 3]
        all_means.(metric)(s, :) = mean(mean_iter, 1, 'omitnan');
        all_stds.(metric)(s, :)  = std(mean_iter, 0, 1, 'omitnan');
    end
end

% === Plotting
plot_metric = 'scale_h';  % Change as needed
figure; hold on;
colors = lines(3);
bin_labels = {'Corners', 'Edges', 'Center'};

for b = 1:bin_types

    errorbar((1:numel(settings)), all_means.(plot_metric)(:, b), ...
                      all_stds.(plot_metric)(:, b), ...
                      '-o', 'DisplayName', bin_labels{b}, ...
                      'LineWidth', 1.5, 'Color', colors(b,:));
end

xlabel('Condition');
ylabel('Scale', 'FontSize', 12);
title(['Binned ' plot_metric ' vs Boundary Effect'], 'FontSize', 13);
legend('Location', 'best');
set(gca, 'FontSize', 12);
grid on;
xticks(1:4); xticklabels(labels);

%% (4) Place + Grid Cell Comparison
suffixes = 1;                      % Folder index
pc_ids = 1:20:241;                 % PC rows (13 total)
iter = 1;
nRows = length(pc_ids) + 1;
nCols = length(suffixes);

figure;
t = tiledlayout(nRows, nCols, 'TileSpacing', 'none', 'Padding', 'none');

for c = 1:nCols
    suf = suffixes(c);
    baseName = ['param_sweeps\Tanni_Covar_ED_covar_Uni_std_' num2str(suf)];
    
    % Load place cell (row 1)
    place_fname = fullfile(baseName, sprintf('iteration_%.1f', iter), 'orig_place_cells.mat');
    place_map = load(place_fname).PlaceCellsUni{10}.fmap;

    nexttile((1 - 1) * nCols + c);
    imagesc(place_map); axis off; colormap(jet);
    title(['PCs\_' num2str(suf)], 'Interpreter', 'none');

    % Load and plot grid cells
    gc_fname = fullfile(baseName, sprintf('iteration_%.1f', iter), 'metrics_and_maps.mat');
    data = load(gc_fname);

    for r = 1:length(pc_ids)
        pc = pc_ids(r);
        tileIdx = r * nCols + c;
        nexttile(tileIdx);
        try
            imagesc(data.GC_metrics{pc}.sac);
        catch
            warning('Missing or bad SAC for PC %d in %s', pc, baseName);
            continue;
        end
        axis off; colormap(jet);
    end
end


