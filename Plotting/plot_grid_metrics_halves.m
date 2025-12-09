function plot_grid_metrics_halves(metrics_all, subfolder)
% Plots bar plots of grid metrics from trapezoid_lr binned maps (Left vs Right)
% Averages across PCs per iteration, compares Left vs Right with a
% STANDARD paired Student's t-test (no Welch), and saves figures + CSVs.
%
% INPUTS
%   metrics_all : nIterations x nPCs cell array. Each cell either [] or a
%                 struct with field .trapezoid_lr (1x2 struct array: Left, Right)
%                 that holds the metric fields listed below.
%   subfolder   : output folder for figures and tables

% === CONFIG ===
metrics = {'expGrd_h_el','expGrd_h','scale_h', ...
           'expGrd_s_el','expGrd_s','scale_s', ...
           'eccent','orient','orient_peak'};
metric_titles = {'Gridness hex','Gridness hex no el','Scale hex', ...
                 'Gridness sq','Gridness sq no el','Scale sq', ...
                 'Eccentricity','Ellipse Orientation','Peak Orientation'};
xtick_labels = {'Left','Right'};

[nIterations, nPCs] = size(metrics_all);
if ~exist(subfolder, 'dir'); mkdir(subfolder); end

% Optional y-axis limits (comment out to disable per metric)
metric_ranges = struct( ...
    'expGrd_h_el', [.3 1.3], 'expGrd_h', [.3 1], 'scale_h', [30 70], ...
    'expGrd_s_el', [.3 1.3], 'expGrd_s', [.3 1], 'scale_s', [30 70], ...
    'eccent', [0.3 0.9], 'orient', [30 90], 'orient_peak', [0 60]);

% === Data containers ===
mean_std_data = {};   % rows: {Metric, Region, Mean, Std}
comparison_data = {}; % rows: {Metric, Comparison, t_value, df, p_value, CI_lower, CI_upper, Cohen_d_paired, Significance}

bar_colors_main = {[.8 .8 .8], [.3 .3 .3]}; % Left, Right

% === LOOP THROUGH METRICS ===
for m = 1:numel(metrics)
    metric = metrics{m};
    title_str = metric_titles{m};
    vals_all = nan(nIterations, 2);  % columns: [Left, Right]

    % ----- aggregate PCs -> iteration means -----
    for iter = 1:nIterations
        vals = nan(nPCs, 2);  % Left and Right for this iteration

        for pc = 1:nPCs
            data = metrics_all{iter, pc};
            if isstruct(data) && isfield(data, 'trapezoid_lr') && ~isempty(data.trapezoid_lr)
                try
                    % Expect 1x2 struct array with metric fields
                    tbl = struct2table(data.trapezoid_lr);
                    if ismember(metric, tbl.Properties.VariableNames)
                        v = tbl.(metric);  % should be [2x1] or [1x2]
                        v = v(:).';        % force row [Left Right]
                        if numel(v) >= 2
                            vals(pc, :) = v(1:2);
                        end
                    end
                catch
                    % skip malformed entries
                end
            end
        end

        % mean across PCs for this iteration
        vals_all(iter, :) = mean(vals, 1, 'omitnan');
    end

    % ----- across-iteration summary -----
    grand_mean = mean(vals_all, 1, 'omitnan');
    grand_std  = std(vals_all, 0, 1, 'omitnan');

    % store means/SDs
    mean_std_data(end+1,:) = {title_str, 'Left',  grand_mean(1), grand_std(1)};
    mean_std_data(end+1,:) = {title_str, 'Right', grand_mean(2), grand_std(2)};

    % ----- STANDARD paired t-test: Left vs Right -----
    g1 = vals_all(:,1);  % Left
    g2 = vals_all(:,2);  % Right
    try
        [~, p, ci, stats] = ttest(g1, g2);  % paired Student's t-test
        t_stat = stats.tstat;
        df     = stats.df;
    catch
        p = NaN; ci = [NaN NaN]; t_stat = NaN; df = NaN;
    end

    % Effect size for paired data: Cohen's d_z = mean(diff) / sd(diff)
    dvec      = g1 - g2;
    cohen_d_p = (mean(dvec, 'omitnan')) ./ (std(dvec, 'omitnan'));

    comparison_data(end+1,:) = { ...
        title_str, 'Left vs Right', t_stat, df, p, ci(1), ci(2), cohen_d_p, pval_star(p)};

    % ----- Plot bar with error bars -----
    fig = figure('Visible','off','Units','inches','Position',[1 1 4 4]);
    hold on;
    for i = 1:2
        bar(i, grand_mean(i), 'FaceColor', bar_colors_main{i}, ...
            'EdgeColor', 'k', 'BarWidth', 0.8, 'LineWidth', 2);
    end
    errorbar(1:2, grand_mean, grand_std, 'k', 'LineStyle','none', 'LineWidth', 2);
    set(gca, 'XTick', 1:2, 'XTickLabel', xtick_labels, 'FontSize', 16, 'LineWidth', 2);
    title(title_str, 'FontSize', 16);
    if isfield(metric_ranges, metric), ylim(metric_ranges.(metric)); end
    box off;

    % save figures
    saveas(fig, fullfile(subfolder, ['Bar_', metric, '.svg']));
    saveas(fig, fullfile(subfolder, ['Bar_', metric, '.fig']));
    close(fig);
end

% === Save Summary Tables ===
mean_std_tbl = cell2table(mean_std_data, ...
    'VariableNames', {'Metric','Region','Mean','Std'});
writetable(mean_std_tbl, fullfile(subfolder, 'Metric_Means_by_Region_HALVES.csv'));

comparison_tbl = cell2table(comparison_data, ...
    'VariableNames', {'Metric','Comparison','t_value','df','p_value','CI_lower','CI_upper','Cohens_d_paired','Significance'});
writetable(comparison_tbl, fullfile(subfolder, 'Metric_Region_Comparisons_HALVES.csv'));

end


