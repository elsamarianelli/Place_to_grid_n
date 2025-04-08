%% Function to get metrics for binned environment 
 
%  Takes:    full_map - GC fr map
%            shape - 'hexagon' or 'square' depending on gridness metric 
%            bin_structe - 'nine' or 'three' depending on wether you want 3
%            by 3 environment bins returned, or for 'three' option, return
%            1 = corners 2 = sides 3 = centre
%  Returns:  binned_metrics - 3x3 cell array with gridness metrics in each cell structure
%            binned_sac - 3x3 sac array  

function [binned_sac, binned_metrics] = get_binned_metrics(full_map, shape, bin_struct)

% [1] Divide map into 9 blocks
rows_per_block = round(size(full_map,1) / 3);  
cols_per_block = round(size(full_map,2) / 3);  

row_blocks = [rows_per_block, rows_per_block, size(full_map,1) - 2 * rows_per_block];
col_blocks = [cols_per_block, cols_per_block, size(full_map,2) - 2 * cols_per_block];

mini_maps = mat2cell(full_map, row_blocks, col_blocks);

% [2] Take SACs of each mini-map
binned_sac = cellfun(@xPearson, mini_maps, 'UniformOutput', false);

% [3] Define a wrapper function to collect multiple outputs from 'multiGridness'
applyMultiGridness = @(sac, map) multiGridness(sac, shape, map, "plot");

% [4] Apply the function using cellfun and manually collect the outputs
stGrd_h = cell(3,3);
expGrd_h = cell(3,3);
scale_h = cell(3,3);

for i = 1:3
    for j = 1:3
        [stGrd_h{i,j}, expGrd_h{i,j}, scale_h{i,j}] = applyMultiGridness(binned_sac{i,j}, mini_maps{i,j});
    end
end

% [5] Organize the results into a structured format

metric_names = {'stGrd_h', 'expGrd_h', 'scale_h'};

% Initialize structure
binned_metrics = struct();

if strcmp(bin_struct, 'nine')
    for i = 1:3
        for j = 1:3
            for m = 1:length(metric_names)
                metric_name = metric_names{m};
                binned_metrics(i,j).(metric_name) = eval([metric_name '{i,j}']);
            end
        end
    end

elseif strcmp(bin_struct, 'three')
    % Define index groups
    corners = {[1,1], [1,3], [3,1], [3,3]};
    edges = {[2,1], [2,3], [1,2], [3,2]};

    for m = 1:length(metric_names)
        metric_name = metric_names{m};
        
        % Extract the relevant metric array
        metric_array = eval(metric_name);

        % Compute mean for corners
        corners_data = cellfun(@(idx) metric_array{idx(1), idx(2)}, corners, 'UniformOutput', false);
        binned_metrics(1).(metric_name) = mean(corners_data);

        % Compute mean for edges
        edges_data = cellfun(@(idx) metric_array{idx(1), idx(2)}, edges);
        binned_metrics(2).(metric_name) = mean(edges_data);

        % Assign the center metric 
        binned_metrics(3).(metric_name) = metric_array{2,2};
    end
end
end