function [binned_sac, binned_metrics_three, binned_metrics_nine] = get_binned_metrics(full_map, shape)
% Function to compute gridness metrics in 3x3 environment bins
% Applies ellipse fitting and correction to SACs before metric extraction

% [1] Divide map into 9 blocks
rows_per_block = round(size(full_map,1) / 3);  
cols_per_block = round(size(full_map,2) / 3);  
row_blocks = [rows_per_block, rows_per_block, size(full_map,1) - 2 * rows_per_block];
col_blocks = [cols_per_block, cols_per_block, size(full_map,2) - 2 * cols_per_block];
mini_maps = mat2cell(full_map, row_blocks, col_blocks);

% [2] Get SACs for each mini-map
binned_sac = cellfun(@xPearson, mini_maps, 'UniformOutput', false);

% [3] Prepare storage for metrics
metric_names = {'stGrd_h', 'expGrd_h', 'scale_h','stGrd_h_el', 'expGrd_h_el', 'scale_h_el', 'eccent', 'orient', 'xyScale', 'abScale', 'ellipicity'};
binned_metrics_nine = struct();
binned_metrics_three = struct();


% [4] Process each bin: ellipse fit, correct SAC, and compute metrics


for i = 1:3
    for j = 1:3
        sac = binned_sac{i,j};
        map = mini_maps{i,j};

        try
            % normal version
            [stGrd, expGrd, scale] = multiGridness(sac, shape, false, map);
            % elipse version - essentially fits ellipse, rescales 
            % figure;
            [xyScale, eccent, orient, abScale, ellip] = gridEllipse_fit(sac, false);
            % saveas(gcf, fullfile('grids_figures/ellipse_fit_example', [num2str(i),'i_', num2str(j), 'j_ellipseSac', '.fig']));
            % saveas(map, fullfile('grids_figures/ellipse_fit_example', [num2str(i),'i_', num2str(j), 'j_map', '.fig']));
            % % 
            % debugged this funcction so should get circle now!
            regSac = gridEllipse_correct(sac, abScale, orient); 
            [stGrd_el, expGrd_el, scale_el] = multiGridness(regSac, shape, false, map);
            % title([num2str(i), num2str(j)])
            
       
            % save
            binned_metrics_nine(i,j).stGrd_h = stGrd;
            binned_metrics_nine(i,j).expGrd_h = expGrd;
            binned_metrics_nine(i,j).scale_h = scale;
            binned_metrics_nine(i,j).stGrd_h_el = stGrd_el;
            binned_metrics_nine(i,j).expGrd_h_el = expGrd_el;
            binned_metrics_nine(i,j).scale_h_el = scale_el;
            binned_metrics_nine(i,j).eccent = eccent;
            binned_metrics_nine(i,j).orient = orient;
            binned_metrics_nine(i,j).xyScale = xyScale;
            binned_metrics_nine(i,j).abScale = abScale;
            binned_metrics_nine(i,j).ellipicity = ellip;
        catch
            % In case of any failure, assign NaNs
            for m = 1:length(metric_names)
                binned_metrics_nine(i,j).(metric_names{m}) = NaN;
            end
        end
    end
end

% [5] Aggregate into 3-bin structure (corners, edges, center)
corners = {[1,1], [1,3], [3,1], [3,3]};
edges   = {[2,1], [2,3], [1,2], [3,2]};
center  = {[2,2]};

for m = 1:length(metric_names)
    name = metric_names{m};
    val_mat = cellfun(@(idx) binned_metrics_nine(idx(1), idx(2)).(name), ...
                      [corners, edges, center], 'UniformOutput', false);

    % 1 = corners
    corner_vals_raw = val_mat(1:4);
    corner_vals = cell2mat(corner_vals_raw(cellfun(@(x) isnumeric(x) && isscalar(x) && ~isnan(x), corner_vals_raw)));
    binned_metrics_three(1).(name) = mean_or_nan(corner_vals);

    % 2 = edges
    edge_vals_raw = val_mat(5:8);
    edge_vals = cell2mat(edge_vals_raw(cellfun(@(x) isnumeric(x) && isscalar(x) && ~isnan(x), edge_vals_raw)));
    binned_metrics_three(2).(name) = mean_or_nan(edge_vals);

    % 3 = center
    center_val_raw = val_mat{9};
    if isnumeric(center_val_raw) && isscalar(center_val_raw) && ~isnan(center_val_raw)
        center_val = center_val_raw;
    else
        center_val = NaN;
    end    
    binned_metrics_three(3).(name) = mean_or_nan(center_val);
end
end

