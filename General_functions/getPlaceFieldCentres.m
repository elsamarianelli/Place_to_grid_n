function [xy_field, env, bin_prob]= getPlaceFieldCentres(env, n_cells, dim_x, dim_y, boundary_effect)
    % Calculate the minimum distance to walls in both x and y directions

    x_dists = env.L == 0; % Identify wall locations
    x_dists([2, dim_y], :) = 0; % Set boundary rows to zero
    x_dists = bwdist(x_dists); % Compute distance to the nearest wall in x direction
    
    y_dists = env.L == 0; % Identify wall locations
    y_dists(:, [2, dim_x]) = 0; % Set boundary columns to zero
    y_dists = bwdist(y_dists); % Compute distance to the nearest wall in y direction
   
    xy_dist = env.dwmap; % Distance to walls (combined distance map)
    
    % Set distances outside the defined environment to NaN
    x_dists(env.L ~= 2) = NaN;
    y_dists(env.L ~= 2) = NaN;
    xy_dist(env.L ~= 2) = NaN;
    
    % Store the calculated distances in the environment structure
    env.x_dists = x_dists;
    env.y_dists = y_dists;
    
    % Calculate sampling weights for field centers (based on inverse distance)
    poss_id = find(~isnan(x_dists)); % Identify valid locations
    bin_prob = (xy_dist).^(-1/boundary_effect); % Inverse distance weighting
    id_prob = bin_prob(poss_id); % Get probabilities for sampling
    
    % Sample field center locations based on the calculated probabilities
    id_field = randsample(poss_id, n_cells, true, id_prob);
    [y_field, x_field] = ind2sub(size(x_dists), id_field); % Convert to subscripts
    xy_field = [x_field, y_field]; % Store field centers
   
    % Resort the coordinates so they are ordered to be close to each other -
    % for visualising covar matrix
    dists = sqrt(sum((xy_field - repmat(xy_field(1,:),size(xy_field,1),1)).^2,2));
    [~,ord] = sort(dists); clear dists
    xy_field = xy_field(ord,:); clear ord
    
end

