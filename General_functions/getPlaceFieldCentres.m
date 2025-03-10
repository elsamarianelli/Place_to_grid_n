function [xy_field, env, bin_prob] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, boundary_effect)
    % Identify wall locations
    wall_mask = (env.L == 0);

    % Compute Euclidean distance to the nearest wall (single shortest distance)
    xy_dist = bwdist(wall_mask);

    % Ensure distances outside the defined environment are NaN
    xy_dist(env.L ~= 2) = NaN;

    % Store calculated distances in the environment structure
    env.xy_dist = xy_dist; % Store full distance map

    % Compute sampling weights using inverse shortest distance (no separate x/y weighting)
    poss_id = find(~isnan(xy_dist)); % Identify valid locations
    bin_prob = (xy_dist).^(-1/boundary_effect); % Inverse distance weighting (closer = higher probability)
    
    % Normalize the probability weights
    bin_prob = bin_prob / max(bin_prob(:)); 
    id_prob = bin_prob(poss_id); % Get probabilities for sampling

    % Sample place field center locations based on corrected probabilities
    id_field = randsample(poss_id, n_cells, true, id_prob);
    [x_field, y_field] = ind2sub(size(xy_dist), id_field);
    xy_field = [x_field, y_field];

    % Order place fields for visualization
    num_coords = size(xy_field, 1);
    sorted_coords = zeros(num_coords, 2);
    current_index = 1;
    sorted_coords(1, :) = xy_field(current_index, :);
    visited = false(num_coords, 1);
    visited(current_index) = true;

    % Resort the coordinates so they are ordered to be close to each other
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

    xy_field = sorted_coords;
end
% figure; plot(xy_field(:,1), xy_field(:,2),'.')

