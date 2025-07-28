function [binned_sac, binned_metrics_three, binned_metrics_nine, binned_metrics_variability] = get_binned_metrics(full_map)
% GET_BINNED_METRICS Computes grid metrics in 3x3 spatial bins.
% Outputs:
%   - binned_sac:            3x3 cell array of SACs (spatial autocorrelograms)
%   - binned_metrics_nine:   Struct of grid metrics for each bin
%   - binned_metrics_three:  Aggregated metrics for corners, edges, center
%   - binned_metrics_variability: Normalized deviation from bin mean
% Note: for the gridness metrics, gridness is included for every Grid cell
% regardless f wether it has enough peaks, with the other metrics hwoeer,
% grid cells are only included if every region has enough peaks (6 - from
% stensola methods)

% --- 1. Divide map into blocks ---
[rows, cols] = size(full_map);
row_blocks = round([rows/3, rows/3, rows - 2*round(rows/3)]);
col_blocks = round([cols/3, cols/3, cols - 2*round(cols/3)]);
mini_maps = mat2cell(full_map, row_blocks, col_blocks);

% --- 2. Compute SACs for each bin ---
binned_sac = cellfun(@xPearson, mini_maps, 'UniformOutput', false);

% --- 3. Define metrics to compute ---
metric_names = {'stGrd_h', 'expGrd_h', 'scale_h', ...
                'stGrd_h_el', 'expGrd_h_el', 'scale_h_el', ...
                'eccent', 'orient', 'xyScale', 'abScale', 'ellipicity', 'orient_peak', ...
                'stGrd_s', 'expGrd_s', 'scale_s', ...
                'stGrd_s_el', 'expGrd_s_el', 'scale_s_el', 'sixPeaks'};

binned_metrics_nine = struct();

% --- 4. Compute metrics for each bin ---
for i = 1:3
    for j = 1:3
        map = mini_maps{i,j};
        sac = binned_sac{i,j};

        smoothed_map = smoothdata(map);
        peaks = FastPeakFind(smoothed_map);
        has_enough_peaks = numel(peaks)/2 > 6;

        for m = 1:numel(metric_names)
            binned_metrics_nine(i,j).(metric_names{m}) = NaN;
        end

        % Always compute basic gridness metrics (outside try block)
        [stGrd, expGrd, scale, orient_peak] = multiGridness(sac, 'hexagon', false, map);
        [stGrd_sq, expGrd_sq, scale_sq, ~] = multiGridness(sac, 'square', false, map);

        % Pre-set default values for metrics that depend on ellipse fitting
        [stGrd_el, expGrd_el, scale_el] = deal(NaN);
        [stGrd_sq_el, expGrd_sq_el, scale_sq_el] = deal(NaN);
        [eccent, orient, xyScale, abScale, ellip] = deal(NaN);

        % Attempt ellipse-dependent metrics
        try
            [xyScale, eccent, orient, abScale, ellip] = gridEllipse_fit(sac, false);
            regSac = gridEllipse_correct(sac, abScale, orient);
            [stGrd_el, expGrd_el, scale_el] = multiGridness(regSac, 'hexagon', true, map);

            regSac_sq = gridEllipse_correct(sac, abScale, orient);
            [stGrd_sq_el, expGrd_sq_el, scale_sq_el] = multiGridness(regSac_sq, 'square', true, map);
        catch
            % If ellipse fails, fallback values remain NaN
        end

        % Save all values regardless of try success
        values = {
            stGrd, expGrd, scale, ...
            stGrd_el, expGrd_el, scale_el, ...
            eccent, rad2deg(orient), xyScale, abScale, ellip, orient_peak, ...
            stGrd_sq, expGrd_sq, scale_sq, ...
            stGrd_sq_el, expGrd_sq_el, scale_sq_el, has_enough_peaks
        };

        for m = 1:numel(metric_names)
            binned_metrics_nine(i,j).(metric_names{m}) = values{m};
        end

    end
end

% --- 5. Identify non-gridness metrics and validate 6-peak coverage ---
gridness_metrics = metric_names(contains(metric_names, 'Grd'));
non_grid_metrics = setdiff(metric_names, gridness_metrics);

% Check if any bin has insufficient peaks
six_peaks_mask = false(3,3);
for i = 1:3
    for j = 1:3
        val = binned_metrics_nine(i,j).sixPeaks;
        six_peaks_mask(i,j) = islogical(val) && val;
    end
end

if ~all(six_peaks_mask(:))
    for i = 1:3
        for j = 1:3
            for m = 1:numel(non_grid_metrics)
                binned_metrics_nine(i,j).(non_grid_metrics{m}) = NaN;
            end
        end
    end
end

% --- 6. Define region groupings ---
corners = {[1,1], [1,3], [3,1], [3,3]};
edges   = {[1,2], [2,1], [2,3], [3,2]};
center  = {[2,2]};
regions = {corners, edges, center};

binned_metrics_three = struct();
binned_metrics_variability = struct();

for m = 1:length(metric_names)
    name = metric_names{m};
    metric_grid = nan(3,3);
    for i = 1:3
        for j = 1:3
            val = binned_metrics_nine(i,j).(name);
            if isnumeric(val) && isscalar(val)
                metric_grid(i,j) = val;
            end
        end
    end

    for r = 1:3
        binned_metrics_three(r).(name) = mean_vals_from_indices(metric_grid, regions{r});
    end

    if all(isnan(metric_grid(:))) || nnz(~isnan(metric_grid)) < 3
        for r = 1:3
            binned_metrics_variability(r).(name) = NaN;
        end
    else
        norm_grid = metric_grid / mean(metric_grid(:), 'omitnan');
        var_grid = abs(norm_grid - 1);
        for r = 1:3
            binned_metrics_variability(r).(name) = mean_vals_from_indices(var_grid, regions{r});
        end
    end
end
end

function avg = mean_vals_from_indices(mat, idx_pairs)
    vals = cellfun(@(idx) mat(idx(1), idx(2)), idx_pairs);
    avg = mean(vals, 'omitnan');
end
