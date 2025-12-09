function plot_grid_metrics_ANOVA(metrics_all, subfolder)
% Plots grid metrics as bar plots and heatmaps from binned metric data.
% Uses paired Student's t-tests for pairwise region comparisons
% (since each iteration contributes all 3 regions).
% Includes ANOVA 1x3  repeated measures (corners/edges/center) check 
% only runs paired t tests if ANOVA is significant
% N in the output table = sum over iterations of how many PC values
% actually contributed to each iteration's mean for each region.

% === CONFIG ===
metrics = {'expGrd_h_el','expGrd_h','scale_h','expGrd_s_el','expGrd_s','scale_s','eccent','orient','orient_peak'};
metric_titles = {'Gridness hex','Gridness hex no el','Scale hex','Gridness sq','Gridness sq no el','Scale sq','Eccentricity','Ellipse Orientation','Peak Orientation'};
xtick_labels = {'Corners','Edges','Center'};
region_count = numel(xtick_labels);
[nIterations, nPCs] = size(metrics_all);
if ~exist(subfolder, 'dir'); mkdir(subfolder); end

% Optional y-axis/colorbar ranges
metric_ranges = struct( ...
    'expGrd_h_el', [.3 1.3], 'expGrd_h', [.3 1], 'scale_h', [30 70], ...
    'expGrd_s_el', [.3 1.3], 'expGrd_s', [.3 1], 'scale_s', [30 70], ...
    'eccent', [0.3 0.9], 'orient', [30 90], 'orient_peak', [0 60]);

% Store data for tables
% Includes N = sum over iterations of per-iteration PC counts used
mean_std_data = {};   % {'Metric','Region','Mean','Std','N'}
comparison_data = {};

% === BAR PLOTS ===
bar_colors_main = {[.8 .8 .8], [0.5 0.5 .5], [.3 .3 .3]};

for m = 1:numel(metrics)
    metric = metrics{m};
    title_str = metric_titles{m};

    vals_all  = nan(nIterations, region_count);   % iteration means
    pc_counts = zeros(nIterations, region_count); % per-iteration contributing PC counts

    % ----- aggregate PCs -> iteration means per region -----
    for iter = 1:nIterations
        vals = nan(nPCs, region_count);
        for pc = 1:nPCs
            disp(pc)
            data = metrics_all{iter, pc};
            if isstruct(data) && isfield(data, 'three') && ~isempty(data.three)
                try
                    tbl = struct2table(data.three);
                    if ismember(metric, tbl.Properties.VariableNames)
                        v = tbl.(metric);   % 1x3 or 3x1
                        v = v(:).';         % [Corners Edges Center]
                        if numel(v) >= 3
                            vals(pc, :) = v(1:3);
                        end
                    end
                catch
                    % skip 
                end
            end
        end
        vals_all(iter, :)  = mean(vals, 1, 'omitnan');   % iteration mean
        pc_counts(iter, :) = sum(~isnan(vals), 1);       % how many PCs used per region this iteration
    end

    % ----- across-iteration summary -----
    grand_mean = mean(vals_all, 1, 'omitnan');
    grand_std  = std(vals_all, 0, 1, 'omitnan');
    N_regions  = sum(pc_counts, 1);  % total PCs contributing across all iterations

    % Store region stats (+ N)
    mean_std_data(end+1,:) = {title_str, 'Corners', grand_mean(1), grand_std(1), N_regions(1)};
    mean_std_data(end+1,:) = {title_str, 'Edges',   grand_mean(2), grand_std(2), N_regions(2)};
    mean_std_data(end+1,:) = {title_str, 'Center',  grand_mean(3), grand_std(3), N_regions(3)};

    % === Pairwise comparisons (PAIRED Student's t-test) ===
    comp_pairs = { ...
        {'Corners', vals_all(:,1)}, ...
        {'Edges',   vals_all(:,2)}, ...
        {'Center',  vals_all(:,3)} };

    pair_labels = {'Corners vs Edges', 'Corners vs Center', 'Edges vs Center'};
    pair_idx    = [1 2; 1 3; 2 3];

    for cp = 1:3
        g1_data = comp_pairs{pair_idx(cp,1)}{2};
        g2_data = comp_pairs{pair_idx(cp,2)}{2};

        % PAIRED Student's t-test
        try
            [~, p, ci, stats] = ttest(g1_data, g2_data);
            t_stat = stats.tstat;
            df     = stats.df;   % df = nPairs - 1
        catch
            p = NaN; ci = [NaN NaN]; t_stat = NaN; df = NaN;
        end

        % Paired effect size (Cohen's dz)
        dvec      = g1_data - g2_data;
        cohen_d_p = mean(dvec, 'omitnan') ./ std(dvec, 'omitnan');

        comparison_data(end+1,:) = { ...
            title_str, ...
            pair_labels{cp}, ...
            t_stat, ...
            df, ...
            p, ...
            ci(1), ...
            ci(2), ...
            cohen_d_p, ...
            pval_star(p)};
    end

    % ----- Plot bars -----
    fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 4 4]);
    hold on;
    for i = 1:region_count
        bar(i, grand_mean(i), 'FaceColor', bar_colors_main{i}, 'EdgeColor', 'k', 'BarWidth', 0.8, 'LineWidth', 2);
    end
    errorbar(1:region_count, grand_mean, grand_std, 'k', 'LineStyle', 'none', 'LineWidth', 2);
    set(gca, 'XTick', 1:region_count, 'XTickLabel', xtick_labels, 'FontSize', 16, 'LineWidth', 2);
    title(title_str, 'FontSize', 16);
    if isfield(metric_ranges, metric), ylim(metric_ranges.(metric)); end
    box off;
    saveas(fig, fullfile(subfolder, ['Bar_', metric, '.fig']));
    saveas(fig, fullfile(subfolder, ['Bar_', metric, '.svg']));
    close(fig);
end

% === HEATMAPS ===
for m = 1:numel(metrics)
    metric = metrics{m};
    title_str = metric_titles{m};
    heatmaps = {};

    for iter = 1:nIterations
        for pc = 1:nPCs
            data = metrics_all{iter, pc};
            if isstruct(data) && isfield(data, 'nine') && ~isempty(data.nine)
                try
                    mat = nan(3);
                    for i = 1:3
                        for j = 1:3
                            mat(i,j) = data.nine(i,j).(metric);
                        end
                    end
                    heatmaps{end+1} = mat; %#ok<AGROW>
                catch
                    % skip
                end
            end
        end
    end

    if ~isempty(heatmaps)
        avg_map = mean(cat(3, heatmaps{:}), 3, 'omitnan');
        fig = figure('Visible', 'off');
        imagesc(avg_map);
        colormap(gray);
        if isfield(metric_ranges, metric), clim(metric_ranges.(metric)); end
        axis image;
        set(gca, 'XTick', [], 'YTick', []);
        cb = colorbar; cb.FontSize = 16;
        title(['Heatmap: ', title_str]);
        saveas(fig, fullfile(subfolder, ['Heatmap_', metric, '.fig']));
        saveas(fig, fullfile(subfolder, ['Heatmap_', metric, '.svg']));
        close(fig);
    end
end

% === Save Summary Tables ===
mean_std_tbl = cell2table(mean_std_data, 'VariableNames', {'Metric', 'Region', 'Mean', 'Std', 'N'});
writetable(mean_std_tbl, fullfile(subfolder, 'Metric_Means_by_Region_f.csv'));

comparison_tbl = cell2table(comparison_data, ...
    'VariableNames', {'Metric', 'Comparison','t_value','df','p_value','CI_lower','CI_upper','Cohens_d_paired','Significance'});
writetable(comparison_tbl, fullfile(subfolder, 'Metric_Region_Comparisons_paired_f.csv'));
end

% ===== Helper =====
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
