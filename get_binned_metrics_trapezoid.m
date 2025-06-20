function [binned_sac, binned_metrics_two, binned_metrics_variability] = get_binned_metrics_trapezoid(full_map, shape)
% Divides trapezoid map into left/right halves and computes grid metrics

% [1] Split full map into 2 vertical halves (left vs right)
cols = size(full_map, 2);
mid_col = floor(cols / 2);
left_map  = full_map(:, 1:mid_col);
right_map = full_map(:, mid_col+1:end);
mini_maps = {left_map, right_map};

% [2] Get SACs
binned_sac = cellfun(@xPearson, mini_maps, 'UniformOutput', false);

% [3] Metric names
metric_names = {'stGrd_h', 'expGrd_h', 'scale_h', ...
                'stGrd_h_el', 'expGrd_h_el', 'scale_h_el', ...
                'eccent', 'orient', 'xyScale', 'abScale', 'ellipicity'};

% Init output structs
binned_metrics_two = struct();
binned_metrics_variability = struct();

% [4] Extract metrics
metric_grid = nan(1, 2);  % left and right

for b = 1:2
    map = mini_maps{b};
    sac = binned_sac{b};
    try
        [stGrd, expGrd, scale] = multiGridness(sac, shape, false, map);
        [xyScale, eccent, orient, abScale, ellip] = gridEllipse_fit(sac, true);
        regSac = gridEllipse_correct(sac, abScale, orient);
        [stGrd_el, expGrd_el, scale_el] = multiGridness(regSac, shape, true, map);

        bin_struct = struct( ...
            'stGrd_h',     stGrd, ...
            'expGrd_h',    expGrd, ...
            'scale_h',     scale, ...
            'stGrd_h_el',  stGrd_el, ...
            'expGrd_h_el', expGrd_el, ...
            'scale_h_el',  scale_el, ...
            'eccent',      eccent, ...
            'orient',      rad2deg(orient), ...
            'xyScale',     xyScale, ...
            'abScale',     abScale, ...
            'ellipicity',  ellip ...
        );

        binned_metrics_two(b) = bin_struct;
    catch
        nan_struct = cell2struct(repmat({NaN}, 1, numel(metric_names)), metric_names, 2);
        binned_metrics_two(b) = nan_struct;
    end
end

% [5] Compute variability: abs(norm(x) - 1)
for m = 1:numel(metric_names)
    name = metric_names{m};
    vals = [binned_metrics_two(1).(name), binned_metrics_two(2).(name)];

    if all(isnan(vals)) || nnz(~isnan(vals)) < 2
        binned_metrics_variability(1).(name) = NaN;
        binned_metrics_variability(2).(name) = NaN;
    else
        mean_val = mean(vals, 'omitnan');
        norm_vals = vals / mean_val;
        var_vals = abs(norm_vals - 1);

        binned_metrics_variability(1).(name) = var_vals(1);
        binned_metrics_variability(2).(name) = var_vals(2);
    end
end
end
