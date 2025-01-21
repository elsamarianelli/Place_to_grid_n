function grid_map = get_plot_grid(cells, eigvectors, eigvalues, env)

% Loop through each cell to compute its grid cell representation
for i = 1:length(cells)
    % Initialize the grid map as a zero matrix with the same size as the firing map
    grid_map = zeros(size(cells{1}.fmap));
    
    % Get the eigenvector corresponding to the current cell, use only the real part
    % this is not a completely symmetrical matrix, some Eig Vectors will be complex as a result, 
    % and so real(V) of all eig vectors are take to generate grid cells
    ev = real(V(:,i));
    
    % Sum place cell inputs weighted by the eigenvector components
    for j = 1:length(cells)
        grid_map = grid_map + ev(j) * cells{j}.fmap;
    end
    
    % Normalize the grid map to a maximum value of 1
    grid_map = grid_map ./ max(grid_map(:));
    
    % Apply Gaussian smoothing to the grid map, adjusting for environment mask
    grid_map = imgaussfilt(grid_map, smth_sig) ./ imgaussfilt(1 * (env.L == 2), smth_sig);
    
    % Set areas outside the environment (where L <= 1) to NaN
    grid_map(env.L <= 1) = NaN;
    
    % Threshold the grid map to ensure no negative values
    grid_map = max(grid_map - 0, 0);
    
    % Store the computed grid map in the cell structure
    cells{i}.grid = grid_map;
end

figure
hold on
for i = 1:n*m
    subplot(n,m,i)
    h = imagesc(cells{i}.grid); colormap jet;
    set(h,'AlphaData',env.L > 1)
    axis off
    title(sprintf('grid %i',i))
    pbaspect([env.dim_x, env.dim_y, 1])
end
set(gcf,'Position',[100 100 m*100 n*100])