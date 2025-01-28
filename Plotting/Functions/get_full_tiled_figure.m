function tiled_fig = get_full_tiled_figure(Place_cell_type, place_sf, grid, place_columns, grid_scales, gridnessData, n, list_grid, n_cells, n_iter)
    % Create figure and tiled layout
    figure()
    tiled_fig = tiledlayout(6, n);
    title(tiled_fig, Place_cell_type, 'FontSize', 24);  % Larger font size for the main title
    
    % Set 'TileSpacing' and 'Padding' to 'compact' to reduce spacing
    tiled_fig.TileSpacing = 'tight';  % Compresses spacing between tiles
    tiled_fig.Padding = 'tight';      % Compresses the padding around the layout

    % Plot place cells used for successor features 
    for i = 1:n
        nexttile
        imagesc(place_sf{i*23}.fmap); axis off; colormap(jet);
    end
    nexttile(3)
    title('Place cells for successor features', 'FontSize', 16);  % Standard font size for headings

    % Plot some random grids 
    for i = 1:n
        nexttile
        imagesc(grid{list_grid(i)}.grid); axis off; colormap(jet);
    end
    nexttile(8)
    title('Grid cells from PCA', 'FontSize', 16);  % Standard font size for headings

    % Plot spatial autocorrelograms of those grid cells
    for i = 1:n
        nexttile 
        sac = xPearson(grid{list_grid(i)}.grid);
        imagesc(sac); colormap(jet);
        axis off;
    end
    nexttile(13)
    title('Spatial Autocorrelograms of Grid Cells', 'FontSize', 16);  % Standard font size for headings

    % Plot the place cells taken as columns 
    for i = 1:n
        nexttile
        imagesc(place_columns{i*10}.place); axis off; colormap(jet);
    end
    nexttile(18)
    title('Place Cells from SR Columns', 'FontSize', 16);  % Standard font size for headings

    % Plot grid scale distribution 
    nexttile(21, [2 3])
    mean_g = nanmean(grid_scales, 3); 
    mean_g = (squeeze(mean_g))/71.71;
    std_g = std(mean_g, 0, 3); 
    mean_g = mean(mean_g, 3);
    bar(reshape(mean_g, 1, [])); % Average across cells
    hold on;
    errorbar(reshape(mean_g, 1, []), reshape(std_g, 1, [])./sqrt(n_iter), 'k.');
    
    % Set axis properties
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 14);  % Standardized font size for axis numbers
    set(gcf, 'color', 'w');
    set(gca, 'FontSize', 14);  % Make axis numbers more legible
    xlabel('Environment Bin', 'FontSize', 12);  % Add x-axis label with smaller font size than the title
    ylim([0, 1.2]);  % Set y-axis limit to 1.2
    title('Grid scale', 'FontSize', 16);  % Standard font size for headings
    ylim([0.8 1.2])
    % Gridness plots
    
    % 1) square standard
    nexttile((5*n)-1)
    edges_grid = linspace(-1.5, 3, 21); 
    customColor = [0, 153, 153] / 255;
    bin_centers = (edges_grid(1:end-1) + edges_grid(2:end)) / 2;
    bar(bin_centers, gridnessData.gridness_square_std_mean, 'FaceColor', customColor);
    hold on;
    errorbar(bin_centers,  gridnessData.gridness_square_std_mean,  gridnessData.gridness_square_std_std, 'k', 'linestyle', 'none');
    title('Square standard', 'FontSize', 12);  % Standard font size for headings
    ylabel('Number of cells', 'FontSize', 12);  % Add y-axis label with smaller font size than title
    
    % 2) square expanded 
    nexttile((5*n))
    bar(bin_centers, gridnessData.gridness_square_exp_mean,'FaceColor', customColor);
    hold on;
    errorbar(bin_centers,  gridnessData.gridness_square_exp_mean,  gridnessData.gridness_square_exp_std, 'k', 'linestyle', 'none');
    title('Square expanded', 'FontSize', 12);  % Standard font size for headings

    % 3) hexagonal standard
    nexttile((6*n)-1)
    bar(bin_centers, gridnessData.gridness_hex_std_mean, 'FaceColor', customColor);
    hold on;
    errorbar(bin_centers,  gridnessData.gridness_hex_std_mean,  gridnessData.gridness_hex_std_std, 'k', 'linestyle', 'none');
    title('Hex standard', 'FontSize', 12);  % Standard font size for headings
    ylabel('Number of cells', 'FontSize', 12);  % Add y-axis label with smaller font size than title

    % 4) hexagonal expanded
    nexttile((6*n))
    bar(bin_centers, gridnessData.gridness_hex_exp_mean, 'FaceColor', customColor);
    hold on;
    errorbar(bin_centers,  gridnessData.gridness_hex_exp_mean,  gridnessData.gridness_hex_exp_std, 'k', 'linestyle', 'none');
    title('Hex expanded', 'FontSize', 12);  % Standard font size for headings
    xlabel('Gridness score', 'FontSize', 12);  % Add x-axis label with smaller font size than title

    % Add xlabel for the bottom two gridness plots
    nexttile((5*n)); 
    xlabel('Gridness score', 'FontSize', 12);  % Add x-axis label with smaller font size than title

end
