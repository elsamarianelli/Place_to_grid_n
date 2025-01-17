function proportion = get_good_grids(gridness, cells, threshold)
ids = find(gridness(:,1) > threshold)'; % Find cells with gridness > gridness threshold 
proportion = length(ids)/length(cells);
for i = 1:length(ids)
    figure;
    subplot(1, 2, 1)
    imagesc(cells{ids(i)}.grid);
    colormap jet;
    set(gca, 'LineWidth', 2);
    set(gcf, 'color', 'w');
    set(gca, 'FontSize', 16);
    axis equal tight;
    title('Grid Cell Firing Rate Map');

    % and visualise autocorrelogram for that cell
    subplot(1, 2, 2)
    imagesc(cells{i}.sac);
    axis equal tight;
    colorbar;
    title('Spatial Autocorrelogram of grid cell Cell');
end
