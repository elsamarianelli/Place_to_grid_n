function cells = generateSanderPCs(env, n_cells, xy_field)
%GENERATESANDERPCS Generates place cell activity maps (rate maps) for each cell.
%   Inputs:
%       env: Environment structure containing maps and dimensions.
%       n_cells: Number of place cells to generate.
%       xy_field: Matrix of [x, y] coordinates for the center of each place field.
%   Outputs:
%       cells: A cell array where each cell contains a structure with the place field map (fmap).

% Initialize an empty cell array to hold place cell data
cells = cell(1, n_cells);

% Find all valid (non-wall) positions in the environment
bin_id = find(env.L == 2);

% Loop over each cell to generate its place field map
for n = 1:n_cells
    % Center of the place field for the current cell
    mean_x = xy_field(n, 1); 
    mean_y = xy_field(n, 2);

    % Standard deviations in x and y directions (~ width of field / 3)
    sig_x = fieldWidth(env.x_dists(mean_y, mean_x)) / 3; 
    sig_y = fieldWidth(env.y_dists(mean_y, mean_x)) / 3;

    % Initialize the rate map (place field map) as a zero matrix
    place_map = zeros(size(env.L));
    
    % Loop over all valid positions in the environment to compute the firing rate
    for i = 1:length(bin_id)
        [y, x] = ind2sub(size(place_map), bin_id(i)); 
        % Calculate firing rate
        place_map(y, x) = exp(-(x-mean_x)^2 / (2*sig_x^2)) * exp(-(y-mean_y)^2 / (2*sig_y^2)) / (2*pi*sig_x*sig_y);
    end
    
    % Normalize the place map by its maximum value
    place_map = place_map / max(place_map(:));
    
    % Store the normalized place map in the cells array
    cells{n}.fmap = place_map;
end

end

% % Helper function to calculate the firing rate at a given position (x, y)
% % based on a Gaussian distribution centered at (mean_x, mean_y) with standard
% % deviations sig_x and sig_y.
% function fr = firingRate(x, y, mean_x, mean_y, sig_x, sig_y)
%     % Calculate the firing rate using a 2D Gaussian function
%     fr = exp(-(x-mean_x)^2 / (2*sig_x^2)) * exp(-(y-mean_y)^2 / (2*sig_y^2)) / (2*pi*sig_x*sig_y);
% end

