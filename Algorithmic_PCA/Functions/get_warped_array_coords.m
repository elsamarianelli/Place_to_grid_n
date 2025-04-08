function [coordinates, env] = get_warped_array_coords(env, dim_x, dim_y, n_cells, type)
    
    % get distances to environment walls 
    x_dists = env.L == 0; % Identify wall locations
    x_dists([2, dim_x], :) = 0; % Set boundary rows to zero
    x_dists = bwdist(x_dists); % Compute distance to the nearest wall in x direction
    
    y_dists = env.L == 0; % Identify wall locations
    y_dists(:, [2, dim_y]) = 0; % Set boundary columns to zero
    y_dists = bwdist(y_dists); % Compute distance to the nearest wall in y direction
    
    xy_dist = env.dwmap; % Distance to walls (combined distance map)
    
    % Set distances outside the defined environment to NaN
    x_dists(env.L ~= 2) = NaN;
    y_dists(env.L ~= 2) = NaN;
    xy_dist(env.L ~= 2) = NaN;
    
    % Store the calculated distances in the environment structure
    env.x_dists = x_dists;
    env.y_dists = y_dists;

    % Generate evenly spaced points in [0, dim_x] to control density
    t_x = (linspace(0+2, dim_x-2, ceil(sqrt(n_cells))));
    t_y = (linspace(0+2, dim_y-2, ceil(sqrt(n_cells))));
    
    if strcmp(type, 'warped') % need to update this because old way wasnt working
        %%  Apply the boundary warping need to add
        values_x = t_x;
        values_y = t_y;
    elseif strcmp(type, 'uniform')
        values_x = t_x;
        values_y = t_y;
    end
    
    % Use these values as the locaitons in x and y get warped array coords
    [X, Y] = meshgrid(values_x, values_y);
    % Reshape into a list of coordinates
    coordinates = [X(:), Y(:)];
    coordinates = floor(coordinates);
end