function plot_grid_metrics_ANOVA(metrics_all, subfolder)
% plot_grid_metrics
% From a metrics_all cell array (iterations × PCs), this function:
%   1) Computes iteration-level means for each metric in each region
%      (Corners / Edges / Center).
%   2) Runs a 1×3 repeated-measures ANOVA across regions (per metric).
%   3) If the ANOVA is significant, runs paired t-tests between regions.
%   4) Plots bar charts with error bars (mean ± SD across iterations).
%   5) Plots 3×3 average heatmaps for each metric.
%   6) Saves:
%        - Metric_Means_by_Region_f.csv
%        - Metric_Region_ANOVA.csv
%        - Metric_Region_Comparisons_paired_f.csv
%
% N in the output table = total number of PCs that contributed to the
% iteration means for that region (summed across iterations).

%% config
% Metric names as stored in metric structs
metrics = { ...
    'expGrd_h_el', ...   % expanding hex gridness (ellipse-regularised)
    'expGrd_h',    ...   % expanding hex gridness (standard)
    'scale_h',     ...   % hex grid scale
    'expGrd_s_el', ...
    'expGrd_s',    ...
    'scale_s',     ...
    'eccent',      ...   % eccentricity
    'orient',      ...   % ellipse orientation (angle; circular)
    'orient_peak'};      % peak orientation (angle; circular)

% labels for plots / tables
metric_titles = { ...
    'Gridness hex', ...
    'Gridness hex no el', ...
    'Scale hex', ...
    'Gridness sq', ...
    'Gridness sq no el', ...
    'Scale sq', ...
    'Eccentricity', ...
    'Ellipse Orientation', ...
    'Peak Orientation'};

% Regions in the 3×3 bin layout
region_labels = {'Corners','Edges','Center'};
nRegions      = numel(region_labels);

% Count iterations of PCs
[nIterations, nPCs] = size(metrics_all);

% Make output folder if needed
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

% y-limits for bar plots / colour limits for heatmaps
metric_ranges = struct( ...
    'expGrd_h_el', [.3 1.3], ...
    'expGrd_h',    [.3 1],   ...
    'scale_h',     [30 70],  ...
    'expGrd_s_el', [.3 1.3], ...
    'expGrd_s',    [.3 1],   ...
    'scale_s',     [30 70],  ...
    'eccent',      [0.3 0.9], ...
    'orient',      [30 90],   ...
    'orient_peak', [0 60]);

% Colours for bars (light to dark)
bar_colors = {[.8 .8 .8], [0.5 0.5 0.5], [.3 .3 .3]};

% Significance threshold for ANOVA
alpha = 0.05;

%% STORAGE

% For mean ± SD and N per region
mean_std_rows = {};  % columns: {'Metric','Region','Mean','Std','N'}

% For ANOVA (1 row per metric)
anova_rows    = {};  % columns: {'Metric','F_value','df1','df2','p_value'}

% For pairwise comparisons (3 rows per metric)
comp_rows     = {};  % columns: {'Metric','Comparison','t','df','p', ...
                    %           'CI_low','CI_high','Cohens_d','Stars'}

%%  MAIN LOOP

for m = 1:numel(metrics)
    metric_name  = metrics{m};
    metric_label = metric_titles{m};

    % iteration_means(i, r): mean of this metric in region r for iteration i
    iteration_means = nan(nIterations, nRegions);

    % pc_counts(i, r): how many PCs contributed to region r in iteration i
    pc_counts = zeros(nIterations, nRegions);

    %% 1) GATHER PC VALUES AND COMPUTE ITERATION MEANS

    for it = 1:nIterations

        % Temporary storage for this iteration:
        % rows = PCs, columns = [Corners Edges Center]
        vals_this_iter = nan(nPCs, nRegions);

        for pc = 1:nPCs
            data = metrics_all{it, pc};

            % We only care about rectangular 3×3 binning here: data.three
            if isstruct(data) && isfield(data, 'three') && ~isempty(data.three)
                try
                    tbl = struct2table(data.three);  % one row per region

                    % Only proceed if this metric exists in the table
                    if ismember(metric_name, tbl.Properties.VariableNames)
                        v = tbl.(metric_name);   % can be 1×3 or 3×1
                        v = v(:).';             % force to [1×3] row

                        % Take the first 3 elements as [Corners Edges Center]
                        if numel(v) >= 3
                            vals_this_iter(pc, :) = v(1:3);
                        end
                    end
                catch
                    % Ignore any PC with malformed data; stays as NaN
                end
            end
        end

        % Mean across PCs for this iteration and region (omit NaNs)
        iteration_means(it, :) = mean(vals_this_iter, 1, 'omitnan');

        % Count how many PCs contributed non-NaN values to each region
        pc_counts(it, :) = sum(~isnan(vals_this_iter), 1);
    end

    %% 2) SUMMARY ACROSS ITERATIONS: MEAN, SD, N

    region_means = mean(iteration_means, 1, 'omitnan');
    region_stds  = std(iteration_means, 0, 1, 'omitnan');
    % Total contributing PCs per region across all iterations
    region_N     = sum(pc_counts, 1);

    % Store in summary table
    for r = 1:nRegions
        mean_std_rows(end+1, :) = { ...
            metric_label, ...
            region_labels{r}, ...
            region_means(r), ...
            region_stds(r), ...
            region_N(r)};
    end

    %% 3) 1×3 REPEATED-MEASURES ANOVA (Corners / Edges / Center)

    % iteration_means: rows = iterations, columns = [Corners Edges Center]
    % We treat each iteration as a "subject" measured in 3 conditions.
    data_tbl = array2table(iteration_means, ...
        'VariableNames', region_labels);

    % Define within-subject factor "Region" with 3 levels
    within_tbl = table(categorical(region_labels'), ...
        'VariableNames', {'Region'});

    % Default outputs in case something breaks
    F_value = NaN; p_value = NaN; df1 = NaN; df2 = NaN;

    try
        rm = fitrm(data_tbl, 'Corners-Center ~ 1', 'WithinDesign', within_tbl);
        ranovatbl = ranova(rm, 'WithinModel', 'Region');

        % We want the row corresponding to the Region effect
        idx_region = strcmp(ranovatbl.Properties.RowNames, 'Region');
        idx_error  = strcmp(ranovatbl.Properties.RowNames, 'Error(Region)');

        F_value = ranovatbl.F(idx_region);
        p_value = ranovatbl.pValue(idx_region);
        df1     = ranovatbl.DF(idx_region);      % numerator df
        if any(idx_error)
            df2 = ranovatbl.DF(idx_error);       % denominator df
        else
            df2 = NaN;
        end
    catch
        % If ANOVA fails (e.g., all NaNs), we leave NaNs and move on
    end

    % Store ANOVA result for this metric
    anova_rows(end+1, :) = { ...
        metric_label, F_value, df1, df2, p_value};

    %% 4) PAIRED T-TESTS BETWEEN REGIONS (ONLY IF ANOVA IS SIGNIFICANT)

    % We compare:
    %   1) Corners vs Edges
    %   2) Corners vs Center
    %   3) Edges   vs Center
    pair_defs = { ...
        [1 2], 'Corners vs Edges'; ...
        [1 3], 'Corners vs Center'; ...
        [2 3], 'Edges vs Center'};

    for k = 1:size(pair_defs, 1)
        idx_pair = pair_defs{k,1};
        label    = pair_defs{k,2};

        x = iteration_means(:, idx_pair(1));  % first region
        y = iteration_means(:, idx_pair(2));  % second region

        % If ANOVA is not significant or failed, we do not interpret the
        % pairwise t-tests. We simply record NaNs.
        if ~isnan(p_value) && (p_value < alpha)
            % Paired t-test across iterations
            try
                [~, p, ci, stats] = ttest(x, y);
                t_val = stats.tstat;
                df    = stats.df;
            catch
                p = NaN; ci = [NaN NaN]; t_val = NaN; df = NaN;
            end

            % Paired effect size (Cohen's dz)
            diff_xy   = x - y;
            dz        = mean(diff_xy, 'omitnan') ./ std(diff_xy, 'omitnan');
        else
            % ANOVA not sig → no meaningful pairwise stats
            t_val = NaN;
            df    = sum(~isnan(x - y)) - 1;
            p     = NaN;
            ci    = [NaN NaN];
            dz    = NaN;
        end

        comp_rows(end+1, :) = { ...
            metric_label, ...
            label, ...
            t_val, ...
            df, ...
            p, ...
            ci(1), ...
            ci(2), ...
            dz, ...
            pval_star(p)};
    end

    %% 5) BAR PLOT FOR THIS METRIC

    fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 4 4]);
    hold on;

    % Draw bars
    for r = 1:nRegions
        bar(r, region_means(r), ...
            'FaceColor', bar_colors{r}, ...
            'EdgeColor', 'k', ...
            'BarWidth', 0.8, ...
            'LineWidth', 2);
    end

    % Error bars (SD across iterations)
    errorbar(1:nRegions, region_means, region_stds, ...
        'k', 'LineStyle', 'none', 'LineWidth', 2);

    set(gca, 'XTick', 1:nRegions, ...
             'XTickLabel', region_labels, ...
             'FontSize', 14, ...
             'LineWidth', 1.5);
    title(metric_label, 'FontSize', 16);
    box off;

    % Apply y-limits if defined for this metric
    if isfield(metric_ranges, metric_name)
        ylim(metric_ranges.(metric_name));
    end

    saveas(fig, fullfile(subfolder, ['Bar_', metric_name, '.fig']));
    saveas(fig, fullfile(subfolder, ['Bar_', metric_name, '.svg']));
    close(fig);

end  % end loop over metrics

%% 6) HEATMAPS (3×3) FOR EACH METRIC

% This section is unchanged in logic: we average 3×3 maps (data.nine)
% across PCs and iterations, and then plot a simple heatmap.

[nIterations, nPCs] = size(metrics_all);  %#ok<ASGLU> (same as before)

for m = 1:numel(metrics)
    metric_name  = metrics{m};
    metric_label = metric_titles{m};
    maps_for_metric = {};  % collect all 3×3 maps

    for it = 1:nIterations
        for pc = 1:nPCs
            data = metrics_all{it, pc};

            if isstruct(data) && isfield(data, 'nine') && ~isempty(data.nine)
                try
                    mat = nan(3,3);
                    for i = 1:3
                        for j = 1:3
                            mat(i,j) = data.nine(i,j).(metric_name);
                        end
                    end
                    maps_for_metric{end+1} = mat; %#ok<AGROW>
                catch
                    % skip malformed entries
                end
            end
        end
    end

    if ~isempty(maps_for_metric)
        % Average over the third dimension (stack of 3×3 maps)
        avg_map = mean(cat(3, maps_for_metric{:}), 3, 'omitnan');

        fig = figure('Visible', 'off');
        imagesc(avg_map);
        colormap(gray);
        axis image;
        set(gca, 'XTick', [], 'YTick', []);
        cb = colorbar;
        cb.FontSize = 14;
        title(['Heatmap: ', metric_label], 'FontSize', 16);

        % Apply colour limits if defined
        if isfield(metric_ranges, metric_name)
            clim(metric_ranges.(metric_name));
        end

        saveas(fig, fullfile(subfolder, ['Heatmap_', metric_name, '.fig']));
        saveas(fig, fullfile(subfolder, ['Heatmap_', metric_name, '.svg']));
        close(fig);
    end
end

%% 7) WRITE OUT TABLES AS CSV

% Means and SDs and N
mean_std_tbl = cell2table(mean_std_rows, ...
    'VariableNames', {'Metric', 'Region', 'Mean', 'Std', 'N'});
writetable(mean_std_tbl, fullfile(subfolder, 'Metric_Means_by_Region_f.csv'));

% ANOVA summary
anova_tbl = cell2table(anova_rows, ...
    'VariableNames', {'Metric','F_value','df1','df2','p_value'});
writetable(anova_tbl, fullfile(subfolder, 'Metric_Region_ANOVA.csv'));

% Pairwise comparisons
comp_tbl = cell2table(comp_rows, ...
    'VariableNames', {'Metric', 'Comparison', 't_value', 'df', ...
                      'p_value', 'CI_lower', 'CI_upper', ...
                      'Cohens_d_paired', 'Significance'});
writetable(comp_tbl, fullfile(subfolder, 'Metric_Region_Comparisons_paired_f.csv'));

end  % end main function


%% === HELPER: CONVERT P-VALUE TO STARS ===
function stars = pval_star(p)
% Converts a p-value into a string of stars for quick lookup.
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
