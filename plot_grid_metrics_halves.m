function plot_grid_metrics_halves(metrics_all, subfolder)
% Plots bar plots of grid metrics from trapezoid_lr binned maps (Left vs Right)
% Inputs:
%   - metrics_all: output from compute_grid_metrics_binned(Grid_cells, 'trapezoid_lr')
%   - subfolder: folder to save figures and summary CSVs

% === CONFIG ===
metrics = {'expGrd_h_el','expGrd_h','scale_h', ...
           'expGrd_s_el','expGrd_s','scale_s', ...
           'eccent','orient','orient_peak'};
metric_titles = {'Gridness hex','Gridness hex no el','Scale hex', ...
                 'Gridness sq','Gridness sq no el','Scale sq', ...
                 'Eccentricity','Ellipse Orientation','Peak Orientation'};
xtick_labels = {'Left', 'Right'};

[nIterations, nPCs] = size(metrics_all);
if ~exist(subfolder, 'dir'); mkdir(subfolder); end

% Optional y-axis limits
metric_ranges = struct( ...
    'expGrd_h_el', [.3 1.3], 'expGrd_h', [.3 1], 'scale_h', [30 70], ...
    'expGrd_s_el', [.3 1.3], 'expGrd_s', [.3 1], 'scale_s', [30 70], ...
    'eccent', [0.3 0.9], 'orient', [30 90], 'orient_peak', [0 60]);

% === Data containers ===
mean_std_data = {};
comparison_data = {};

bar_colors_main = {[.8 .8 .8], [.3 .3 .3]}; % Light and dark gray

% === LOOP THROUGH METRICS ===
for m = 1:numel(metrics)
    metric = metrics{m};
    title_str = metric_titles{m};
    vals_all = nan(nIterations, 2);  % [Left, Right]

    for iter = 1:nIterations
        vals = nan(nPCs, 2);  % Left and Right for this iteration

        for pc = 1:nPCs
            data = metrics_all{iter, pc};
            if isstruct(data) && isfield(data, 'trapezoid_lr')
                try
                    tbl = struct2table(data.trapezoid_lr);
                    if ismember(metric, tbl.Properties.VariableNames)
                        vals(pc, :) = tbl.(metric)';
                    end
                catch
                    % skip if structure invalid
                end
            end
        end

        vals_all(iter, :) = mean(vals, 1, 'omitnan');
    end

    grand_mean = mean(vals_all, 1, 'omitnan');
    grand_std = std(vals_all, 0, 1, 'omitnan');

    % Store for summary table
    mean_std_data(end+1,:) = {title_str, 'Left', grand_mean(1), grand_std(1)};
    mean_std_data(end+1,:) = {title_str, 'Right', grand_mean(2), grand_std(2)};

    % Paired t-test: Left vs Right
    try
        [~, p] = ttest(vals_all(:,1), vals_all(:,2));
    catch
        p = NaN;
    end
    comparison_data(end+1,:) = {title_str, 'Left vs Right', p, pval_star(p)};

    % === Plot bar ===
    fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 4 4]);
    hold on;
    for i = 1:2
        bar(i, grand_mean(i), 'FaceColor', bar_colors_main{i}, ...
            'EdgeColor', 'k', 'BarWidth', 0.8, 'LineWidth', 2);
    end
    errorbar(1:2, grand_mean, grand_std, 'k', 'LineStyle', 'none', 'LineWidth', 2);

    set(gca, 'XTick', 1:2, 'XTickLabel', xtick_labels, 'FontSize', 16, 'LineWidth', 2);
    title(title_str, 'FontSize', 16);
    % if isfield(metric_ranges, metric)
        % ylim(metric_ranges.(metric));
    % end
    box off;
    saveas(fig, fullfile(subfolder, ['Bar_', metric, '.svg']));
    saveas(fig, fullfile(subfolder, ['Bar_', metric, '.fig']));
    close(fig);
end

% === Save Summary Tables ===
mean_std_tbl = cell2table(mean_std_data, 'VariableNames', {'Metric', 'Region', 'Mean', 'Std'});
writetable(mean_std_tbl, fullfile(subfolder, 'Metric_Means_by_Region.csv'));

comparison_tbl = cell2table(comparison_data, 'VariableNames', {'Metric', 'Comparison', 'p_value', 'Significance'});
writetable(comparison_tbl, fullfile(subfolder, 'Metric_Region_Comparisons.csv'));

end

function stars = pval_star(p)
if isnan(p)
    stars = '';
elseif p < 0.001
    stars = '***';
elseif p < 0.01
    stars = '**';
elseif p < 0.05
    stars = '*';
else
    stars = 'n.s.';
end
end
