function summary_table = generate_grid_metric_summary(metrics_all)
% Generate summary stats for 3-region binned grid metrics
% Inputs:
%   metrics_all - [nIterations x nPCs] cell array of structs with .three field
% Output:
%   summary_table - table with means, stds, p-values, and significance stars

    % Config
    region_names = {'Corners', 'Edges', 'Center'};
    region_ids = [1, 2, 3]; % order of regions in .three
    comparisons = {'Corners vs Edges', 'Corners vs Center', 'Edges vs Center'};
    compare_pairs = {[1 2], [1 3], [2 3]};

    % Detect all metrics
    example_entry = find(~cellfun(@isempty, metrics_all), 1, 'first');
    metric_list = fieldnames(metrics_all{example_entry}.three);

    % Preallocate result containers
    rows = {};
    means = [];
    stds = [];
    pvals = [];
    stars = {};

    for m = 1:numel(metric_list)
        metric = metric_list{m};

        % Extract values from metrics_all
        values_by_region = cell(1, 3);
        for r = 1:3
            vals = [];
            for i = 1:numel(metrics_all)
                entry = metrics_all{i};
                if isstruct(entry) && isfield(entry, 'three')
                    val = entry.three(region_ids(r)).(metric);
                    if isnumeric(val) && ~isnan(val)
                        vals(end+1) = val;
                    end
                end
            end
            values_by_region{r} = vals;
            rows{end+1,1} = sprintf('%s – %s', metric, region_names{r});
            means(end+1,1) = mean(vals, 'omitnan');
            stds(end+1,1) = std(vals, 'omitnan');
        end

        % Pairwise comparisons
        for c = 1:length(compare_pairs)
            idx1 = compare_pairs{c}(1);
            idx2 = compare_pairs{c}(2);
            group1 = values_by_region{idx1};
            group2 = values_by_region{idx2};
            [~, p] = ttest2(group1, group2);
            pvals(end+1,1) = p;
            stars{end+1,1} = significance_star(p);
            rows{end+1,1} = sprintf('%s – %s', metric, comparisons{c});
            means(end+1,1) = NaN;
            stds(end+1,1) = NaN;
        end
    end

    % Format table
    summary_table = table(rows, means, stds, pvals, stars, ...
        'VariableNames', {'Metric_Region', 'Mean', 'Std', 'p_value', 'Significance'});
end

function s = significance_star(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end
