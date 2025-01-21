function [xy_field, env, bin_prob] = getPlaceFieldCentres_opt(env, n_cells, dim_x, dim_y, boundary_effect, distribution_type)
    % Generate place field centers with options for random or array-like
    % distribution
    %
    % Inputs:
    % - env: Environment structure with boundary data
    % - n_cells: Number of place cells
    % - dim_x, dim_y: Dimensions of the environment
    % - boundary_effect: Controls weighting by distance to boundaries
    % - distribution_type: 'random' or 'array' for the place field distribution
    %
    % Outputs:
    % - xy_field: Coordinates of place field centers
    % - env: Updated environment structure with distance data
    % - bin_prob: Sampling probabilities for each position


    if strcmp(distribution_type, 'random')
            % Calculate distances to walls
        x_dists = bwdist(env.L == 0); % Distance to walls in x direction
        y_dists = bwdist(env.L == 0); % Distance to walls in y direction
        xy_dist = env.dwmap; % Combined distance map
    
        % Set distances outside the environment to NaN
        x_dists(env.L ~= 2) = NaN;
        y_dists(env.L ~= 2) = NaN;
        xy_dist(env.L ~= 2) = NaN;
    
        % Store calculated distances in the environment structure
        env.x_dists = x_dists;
        env.y_dists = y_dists;
    
        % Calculate sampling weights for field centers (based on inverse distance)
        poss_id = find(~isnan(x_dists)); % Identify valid locations
        bin_prob = (xy_dist).^(-1 / boundary_effect); % Inverse distance weighting
        id_prob = bin_prob(poss_id); % Get probabilities for sampling
        
        % Random distribution (original implementation)
        id_field = randsample(poss_id, n_cells, true, id_prob);
        [y_field, x_field] = ind2sub(size(x_dists), id_field);
        xy_field = [x_field, y_field];

    elseif strcmp(distribution_type, 'array')
        xy_field = get_warped_array_coords(dim_x, dim_y, n_cells, boundary_effect);
    end


    % Sort place fields for visualization (optional)
    num_coords = size(xy_field, 1);
    sorted_coords = zeros(num_coords, 2);
    current_index = 1;
    sorted_coords(1, :) = xy_field(current_index, :);
    visited = false(num_coords, 1);
    visited(current_index) = true;

    for i = 2:num_coords
        % Calculate distances from the current point to all unvisited points
        dists = sqrt(sum((xy_field(~visited, :) - xy_field(current_index, :)).^2, 2));
        
        % Get the indices of unvisited points
        unvisited_indices = find(~visited);
        
        % Find the index of the closest unvisited point
        [~, min_idx] = min(dists);
        next_index = unvisited_indices(min_idx);
        
        % Update the sorted coordinates and mark the point as visited
        sorted_coords(i, :) = xy_field(next_index, :);
        visited(next_index) = true;
        
        % Move to the next point
        current_index = next_index;
    end

    xy_field = sorted_coords; % Final sorted place field centers

end

    
