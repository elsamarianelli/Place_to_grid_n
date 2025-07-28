function [binned_sac, binned_metrics_two] = get_binned_metrics_halves(full_map)
% GET_BINNED_METRICS Computes grid metrics in halve of env area spatial bins.
% Outputs:
%   - binned_sac:            cell array of SACs (spatial autocorrelograms)
%   - binned_metrics_two:   Struct of grid metrics for each bin

% --- 1. Divide map into blocks ---
rotated_map = rot90(full_map, 1);

% Remove median values
map_median = median(rotated_map(~isnan(rotated_map)), 'all');
rotated_map(rotated_map == map_median) = NaN;

% Count non-NaN values and compute half
non_nan_count = sum(~isnan(rotated_map), 'all');
half_count = non_nan_count / 2;

% Find vertical split (column) for equal non-NaNs
[rows, cols] = size(rotated_map);
cumulative_col_count = zeros(cols, 1);
for j = 1:cols
    cumulative_col_count(j) = sum(~isnan(rotated_map(:, 1:j)), 'all');
end
[~, split_col] = min(abs(cumulative_col_count - half_count));

% Split rotated_map into two vertical blocks (left and right)
col_blocks = [split_col, cols - split_col];
mini_maps = mat2cell(rotated_map, rows, col_blocks);

% % Optional: plot map with vertical split
% figure; imagesc(rotated_map);
% colormap('jet'); colorbar; axis image;
% title('Rotated Map with Vertical Split'); hold on;
% plot([split_col, split_col], [1, rows], 'r-', 'LineWidth', 2); hold off;

% --- 2. Compute SACs for each bin ---
binned_sac = cellfun(@xPearson, mini_maps, 'UniformOutput', false);

% --- 3. Define metrics to compute ---
metric_names = {'stGrd_h', 'expGrd_h', 'scale_h', ...
                'stGrd_h_el', 'expGrd_h_el', 'scale_h_el', ...
                'eccent', 'orient', 'xyScale', 'abScale', 'ellipicity', 'orient_peak', ...
                'stGrd_s', 'expGrd_s', 'scale_s', ...
                'stGrd_s_el', 'expGrd_s_el', 'scale_s_el', 'sixPeaks'};

binned_metrics_two = struct();

% --- 4. Compute metrics for each bin ---
for i = 1:2
    map = mini_maps{i};
    sac = binned_sac{i};

    smoothed_map = smoothdata(map);
    peaks = FastPeakFind(smoothed_map);
    has_enough_peaks = numel(peaks)/2 > 6;

    for m = 1:numel(metric_names)
        binned_metrics_two(i).(metric_names{m}) = NaN;
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
        binned_metrics_two(i).(metric_names{m}) = values{m};
    end
end

% --- 5. Identify non-gridness metrics and validate 6-peak coverage ---
gridness_metrics = metric_names(contains(metric_names, 'Grd'));
non_grid_metrics = setdiff(metric_names, gridness_metrics);

% Check if any bin has insufficient peaks
six_peaks_mask = false(1, 2);
for i = 1:2
    val = binned_metrics_two(i).sixPeaks;
    six_peaks_mask(i) = islogical(val) && val;
end

if ~all(six_peaks_mask(:))
    for i = 1:2
        for m = 1:numel(non_grid_metrics)
            binned_metrics_two(i).(non_grid_metrics{m}) = NaN;
        end
    end
end

