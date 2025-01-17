% Function to process and visualize subregions
function plotSubregions(grid, x_lims, y_lims, index)
    % Set up a tiled layout for the combined map
    t = tiledlayout(3, 3, 'TileSpacing', 'none', 'Padding', 'none');
    
    % Loop through each subregion
    for j = 1:size(y_lims, 1)
        for k = 1:size(x_lims, 1)
            try
                % Get the subregion
                subregion = grid(y_lims(j, 1):y_lims(j, 2), x_lims(k, 1):x_lims(k, 2));
                
                % Compute SAC (Spatial Autocorrelogram) for the subregion
                sac = xPearson(subregion);
                
                % Compute metrics using autoCorrProps
                in.sac = sac;
                metrics = autoCorrProps(in); % This assumes autoCorrProps is well-defined
                
                % Plot the result
                ax = nexttile;
                imagesc(sac); % Visualize the SAC
                ax.XTick = [];
                ax.YTick = [];
                ax.XColor = 'none';
                ax.YColor = 'none';
                
            catch ME
                % Handle subregion-specific errors
                fprintf('Error in Subregion (%d, %d): %s\n', j, k, ME.message);
                continue;
            end
        end
    end
    
    % Set a title for the tiled layout
    title(t, sprintf('Subregions for Index %d', index));
end
