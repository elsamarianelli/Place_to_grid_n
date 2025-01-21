function spatialAutocorrelogram(rate_map, cell_type)
    % rate_map: 2D matrix representing the firing rate map of a neuron
    % cell_type: String indicating whether it's a 'place' or 'grid' cell
    % Ensure the rate map is normalized
    rate_map = rate_map / max(rate_map(:));
    rate_map(isnan(rate_map)) = 0;

    % Compute the 2D autocorrelation of the firing rate map
    auto_corr = xcorr2(rate_map);

    % Normalize the autocorrelogram
    auto_corr = auto_corr / max(auto_corr(:));

    % Plot the original rate map
    figure;
    subplot(1,2,1);
    imagesc(rate_map);
    axis equal tight;
    colorbar;
    title([cell_type ' Cell Firing Rate Map']);
    
    % Plot the autocorrelogram
    subplot(1,2,2);
    imagesc(auto_corr);
    axis equal tight;
    colorbar;
    title(['Spatial Autocorrelogram of ' cell_type ' Cell']);
end

