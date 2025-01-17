function [grid_scales, orientations, sacs_stored, stats_store] = get_binned_scales(cells, n_cells, n_iter, binned, shape)

% depending on whether you want the full environemnt or the 9 binned
% environemnt (aka 'binned' or 'not_binned')

if strcmp (binned, 'binned')
 
% Analyze grid scales within defined sub-regions of the environment
n_sub_regions = 3;
grid_scales = zeros(n_sub_regions, n_sub_regions, n_cells, n_iter);
orientations = zeros(n_sub_regions, n_sub_regions, n_cells, n_iter);
sacs_stored = zeros(165, 231, n_sub_regions, n_sub_regions, n_cells, n_iter);
stats_store = cell(n_sub_regions, n_sub_regions, n_cells, n_iter);


x_lims = [3, 118; 119, 234; 235, 350]; % X-axis limits for sub-regions
y_lims = [3, 85; 86, 168; 169, 251]; % Y-axis limits for sub-regions
% x_lims = [3, 89; 90, 175; 176, 261; 262, 347];  % X-axis limits for 4 sub-regions
% y_lims = [3, 65; 66, 127; 128, 189; 190, 251];  % Y-axis limits for 4 sub-regions 

% Loop over all cells and sub-regions to calculate grid scales
for iter = 1:n_iter
    cells_iter = cells{iter};
    for i = 1:n_cells
        for j = 1:n_sub_regions
            for k = 1:n_sub_regions
                sac = xPearson(cells_iter{i}.grid(y_lims(j, 1):y_lims(j, 2), x_lims(k, 1):x_lims(k, 2)));
                sacs_stored(:,:, j, k, i, iter) = sac;
                try
                    stats = sacProps(sac, shape); % Extract properties from autocorrelogram
                    stats_store(j, k, i, iter) = {stats};
                    grid_scales(j, k, i, iter) = stats.scale;
                    orientations(j, k, i, iter) = stats.peakOrient(1);
                catch
                    grid_scales(j, k, i, iter) = NaN; % Handle cases where scale cannot be calculated
                    orientations(j, k, i, iter) = NaN;
                end
            end
        end
    end
    disp(iter)
end


elseif strcmp(binned, 'not_binned')

    % Analyze grid scales within full environment  
    grid_scales = zeros(n_cells, n_iter);
    orientations = zeros(n_cells, n_iter);
    sacs_stored = zeros(505, 703, n_cells, n_iter);
    stats_store = cell(n_cells, n_iter);
    
    % Loop over all cells and sub-regions to calculate grid scales
    for iter = 1:n_iter
        cells_iter = cells{iter};
        for i = 1:n_cells
            sac = xPearson(cells_iter{i}.grid);
            sacs_stored(:,:,i,iter) = sac;
            try
                stats = sacProps(sac, shape); % Extract properties from autocorrelogram
                stats_store( i, iter) = {stats};
                grid_scales(i, iter) = stats.scale;
                orientations( i, iter) = stats.peakOrient(1);
            catch
                grid_scales(i, iter) = NaN; % Handle cases where scale cannot be calculated
                orientations( i, iter) = NaN;
            end
            disp(i)
        end
        disp(iter)
    end
end
end
% figure
% subplot(2, 1, 1)
% imagesc(grid)
% title("uniform")
% subplot(2, 1, 2)
% imagesc(grids_T)
% title("sanders")