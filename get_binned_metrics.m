function [binned_sac, binned_metrics_three, binned_metrics_nine, binned_metrics_variability] = get_binned_metrics(full_map, shape)
% GET_BINNED_METRICS: Calculates grid-related metrics in 3x3 spatial bins.
%
% Outputs:
% - binned_sac: raw SACs in each bin
% - binned_metrics_nine: metrics in each of 9 bins
% - binned_metrics_three: metrics aggregated into corner/edge/center
% - binned_metrics_variability: deviation from normalized mean (per metric)

% [1] Divide map into 3x3 blocks
rows = size(full_map,1); cols = size(full_map,2);
row_blocks = round([rows/3, rows/3, rows - 2*round(rows/3)]);
col_blocks = round([cols/3, cols/3, cols - 2*round(cols/3)]);
mini_maps = mat2cell(full_map, row_blocks, col_blocks);

% [2] Compute SACs for each bin
binned_sac = cellfun(@xPearson, mini_maps, 'UniformOutput', false);

% [3] Metric names to compute
metric_names = {'stGrd_h', 'expGrd_h', 'scale_h', ...
                'stGrd_h_el', 'expGrd_h_el', 'scale_h_el', ...
                'eccent', 'orient', 'xyScale', 'abScale', 'ellipicity'};

% Init storage
binned_metrics_nine = struct();
binned_metrics_three = struct();
binned_metrics_variability = struct();

% [4] Loop over each bin to compute metrics
for i = 1:3
    for j = 1:3
        sac = binned_sac{i,j};
        map = mini_maps{i,j};
        try
            [stGrd, expGrd, scale] = multiGridness(sac, shape, false, map);
            [xyScale, eccent, orient, abScale, ellip] = gridEllipse_fit(sac, true);
            regSac = gridEllipse_correct(sac, abScale, orient);
            [stGrd_el, expGrd_el, scale_el] = multiGridness(regSac, shape, true, map);

            binned_metrics_nine(i,j).stGrd_h     = stGrd;
            binned_metrics_nine(i,j).expGrd_h    = expGrd;
            binned_metrics_nine(i,j).scale_h     = scale;
            binned_metrics_nine(i,j).stGrd_h_el  = stGrd_el;
            binned_metrics_nine(i,j).expGrd_h_el = expGrd_el;
            binned_metrics_nine(i,j).scale_h_el  = scale_el;
            binned_metrics_nine(i,j).eccent      = eccent;
            binned_metrics_nine(i,j).orient      = rad2deg(orient); % convert to degrees
            binned_metrics_nine(i,j).xyScale     = xyScale;
            binned_metrics_nine(i,j).abScale     = abScale;
            binned_metrics_nine(i,j).ellipicity  = ellip;
        catch
            for m = 1:numel(metric_names)
                binned_metrics_nine(i,j).(metric_names{m}) = NaN;
            end
        end
    end
end

% [5] Indices for grouping into regions
corners = {[1,1], [1,3], [3,1], [3,3]};
edges   = {[1,2], [2,1], [2,3], [3,2]};
center  = {[2,2]};

% [6] Aggregate metrics into 3 regions (mean per region) + variability
for m = 1:length(metric_names)
    name = metric_names{m};
    metric_grid = nan(3,3);

    % Get 3x3 grid for this metric
    for i = 1:3
        for j = 1:3
            val = binned_metrics_nine(i,j).(name);
            if isnumeric(val) && isscalar(val)
                metric_grid(i,j) = val;
            end
        end
    end

    % Store region-wise means
    binned_metrics_three(1).(name) = mean_vals_from_indices(metric_grid, corners);
    binned_metrics_three(2).(name) = mean_vals_from_indices(metric_grid, edges);
    binned_metrics_three(3).(name) = mean_vals_from_indices(metric_grid, center);

    % Compute variability (absolute deviation from normalized mean)
    if all(isnan(metric_grid(:))) || nnz(~isnan(metric_grid)) < 3
        binned_metrics_variability(1).(name) = NaN;
        binned_metrics_variability(2).(name) = NaN;
        binned_metrics_variability(3).(name) = NaN;
    else
        norm_grid = metric_grid / mean(metric_grid(:), 'omitnan');
        var_grid = abs(norm_grid - 1);

        binned_metrics_variability(1).(name) = mean_vals_from_indices(var_grid, corners);
        binned_metrics_variability(2).(name) = mean_vals_from_indices(var_grid, edges);
        binned_metrics_variability(3).(name) = mean_vals_from_indices(var_grid, center);
    end
end
end

%% ===== Helper Function =====
function avg = mean_vals_from_indices(mat, idx_pairs)
    vals = cellfun(@(idx) mat(idx(1), idx(2)), idx_pairs);
    avg = mean(vals, 'omitnan');
end
