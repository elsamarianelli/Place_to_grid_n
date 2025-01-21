function cells = get_grid(cells, M, env)
% get_grid - Function to project eigenvectors from eigendecomposition back onto
% the environment to generate grid-like firing maps from place cells.

% Outputs:
%   cells - The input cell array updated with a new field 'grid' for each cell

% Smoothing parameter 
smth_sig = 3;
n = 10; 
m = 10;
% eigendecomposition
[V, ~] = eig(M);
V = real(V);
for i = 1:length(cells)
    % Initialize the grid map as a zero matrix with the same size as the firing map
    grid_map = zeros(size(cells{1}.fmap));
    
    ev = V(:, i);
    % Sum place cell inputs weighted by the eigenvector components
    for j = 1:length(cells)
        grid_map = grid_map + ev(j) * cells{j}.fmap;
    end
    
    % Normalize grid map by dividing by max value
    %  Gaussian smoothing
    grid_map = imgaussfilt(grid_map, smth_sig) ./ imgaussfilt(1 * (env.L == 2), smth_sig);
   
    % Normalize 
    grid_map = grid_map ./ max(grid_map(:)); 

    % Set areas outside the environment (where L <= 1) to NaN
    grid_map(env.L <= 1) = NaN;
    
    % % Ensure no negative values by thresholding at zero
    % grid_map = max(grid_map - 0, 0);
    
    % Store the computed grid map in the cell structure
    cells{i}.grid = grid_map;
end

% Create a new figure for visualizing the grid cell representations
figure;
hold on;

% Loop through each cell to plot its grid map
for i = 1:n*m
    subplot(n, m, i);  
    h = imagesc(cells{i}.grid); 
    colormap jet; 
    set(h, 'AlphaData', env.L > 1);
    axis off; 
    title(sprintf('grid %i', i)); % Title for each subplot
    pbaspect([env.dim_x, env.dim_y, 1]); % Set the plot box aspect ratio
end

% Adjust figure size and position on screen
set(gcf, 'Position', [100, 100, m*100, n*100]);

end
