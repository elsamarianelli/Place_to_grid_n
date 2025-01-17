function xy_sorted = get_sorted_coords(xy_field)
% put cells next to each other
num_coords = size(xy_field, 1);
sorted_coords = zeros(num_coords, 2);
current_index = 1;
sorted_coords(1, :) = xy_field(current_index, :);
visited = false(num_coords, 1);
visited(current_index) = true;

% Resort the coordinates so they are ordered to be close to each other -
% for visualising covar matrix
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

xy_sorted = sorted_coords;
end