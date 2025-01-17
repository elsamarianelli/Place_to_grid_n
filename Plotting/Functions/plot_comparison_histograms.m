function [row_T, col_T, row_U, col_U] = plot_comparison_histograms(U_grid_scales, T_grid_scales, n_bins, measure, axis_handle)
    % Function to plot side-by-side histograms comparing U_grid_scales and T_grid_scales
    % across a specified number of bins.
    
    % Number of cells and iterations
    [~, ~, n_cells, n_iter] = size(T_grid_scales);
    
    % Preallocate arrays for variability and mean scale for both U and T grid scales
    variability_U = zeros(n_cells, n_iter);
    Mean_scale_U = zeros(n_cells, n_iter);

    variability_T = zeros(n_cells, n_iter);
    Mean_scale_T = zeros(n_cells, n_iter);
    


    % Compute variability and mean scale for U_grid_scales and T_grid_scales
    for iter = 1:n_iter
        for i = 1:n_cells
            % For U grid scales
            Cell_scales_U = U_grid_scales(:, :, i, iter);
            Cell_scales_U = reshape(Cell_scales_U, 1, []);
            Cell_scales_SD_U = std(Cell_scales_U / (sum(Cell_scales_U)));
            variability_U(i, iter) = Cell_scales_SD_U;
            Mean_scale_U(i, iter) = mean(Cell_scales_U);
            
            % For T grid scales
            Cell_scales_T = T_grid_scales(:, :, i, iter);
            Cell_scales_T = reshape(Cell_scales_T, 1, []);
            Cell_scales_SD_T = std(Cell_scales_T / (sum(Cell_scales_T)));
            variability_T(i, iter) = Cell_scales_SD_T;
            Mean_scale_T(i, iter) = mean(Cell_scales_T);
        end
    end

    % find cell indexes with max and minimum variability in
    % scale/orientation measure 
    check_variability = variability_T;
    check_variability(1:20, :) = NaN; check_variability(120:end, :) = NaN;
    [~, linear_idx_max] = max(check_variability(:));
    [row_T, col_T] = ind2sub(size(variability_T), linear_idx_max);
    check_variability = variability_U;
    check_variability(1:20, :) = NaN; check_variability(120:end, :) = NaN;
    [~, linear_idx_min] = min(check_variability(:));
    [row_U, col_U] = ind2sub(size(variability_U), linear_idx_min);

    % Determine the global minimum and maximum across both Mean_scale matrices
    min_val = min([Mean_scale_U(:); Mean_scale_T(:)]);
    max_val = max([Mean_scale_U(:); Mean_scale_T(:)]);

    % Round the min and max values for easier binning
    min_val = floor(min_val);
    max_val = ceil(max_val);
    % Create bin edges (with rounded numbers) and calculate the bin centers
    %bin_edges = linspace(min_val, max_val, n_bins+1);
    if strcmp(measure, 'scale')
        bin_edges = linspace(40, 130, n_bins+1);
    elseif strcmp(measure, 'orientation')
        bin_edges = linspace(0 , 150, n_bins+1);
    end
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    % Preallocate arrays for storing average variability per bin
    avg_variability_U = zeros(1, n_bins);
    avg_variability_T = zeros(1, n_bins);

    % Calculate average variability for U_grid_scales and T_grid_scales
    for bin = 1:n_bins
        % Logical mask for the current bin (U grid scales)
        bin_mask_U = Mean_scale_U >= bin_edges(bin) & Mean_scale_U < bin_edges(bin+1);
        variability_in_bin_U = variability_U(bin_mask_U);
        avg_variability_U(bin) = mean(variability_in_bin_U, 'omitnan');
        
        % Logical mask for the current bin (T grid scales)
        bin_mask_T = Mean_scale_T >= bin_edges(bin) & Mean_scale_T < bin_edges(bin+1);
        variability_in_bin_T = variability_T(bin_mask_T);
        avg_variability_T(bin) = mean(variability_in_bin_T, 'omitnan');
    end

    % Plot the side-by-side histograms
    axis(axis_handle);
    hold on;

    % Width of the bars
    bar_width = 0.3; 

    % Plot the histogram for U_grid_scales (shifted to the left)
    bar(bin_centers - bar_width*5, avg_variability_U, bar_width,'FaceColor', [0 0.7 0.7], 'FaceAlpha', 0.6, 'DisplayName', 'Uniform input');

    % Plot the histogram for T_grid_scales (shifted to the right)
    bar(bin_centers + bar_width*5, avg_variability_T, bar_width,'FaceColor', 'b', 'FaceAlpha', 0.6, 'DisplayName', 'Boundary Compressed input');

    % Add labels, title, and legend
    if strcmp(measure, 'scale')
        xlabel('mean scale across environment');
        ylabel('Average Variability (normalised scales)');
        title('Uniform vs Boundary Compressed: Grid scale variability');
    elseif strcmp(measure, 'orientation')
        xlabel('mean orientation angle (from x axis) across environment');
        ylabel('Average Variability (normalised orientation)');
        title('Uniform vs Boundary Compressed: Grid orientation variability');
    end
    legend('Location', 'northeast');

    % Set the x-ticks to the bin centers
    set(gca, 'XTick', bin_centers);

    % % Make sure the x-axis labels are visible
    % xtickangle(45); % Rotate the x-axis labels if necessary

    hold off;
end
