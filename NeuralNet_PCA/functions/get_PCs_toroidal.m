function ToroidalPlaceCellMaps = get_PCs_toroidal(env, n_cells, xy_field,size_fact, boundary_compressed)
% function returns place cell firing maps which map onto toroidal space
% instead of normal 2d environemnt, as used in dordek in their NN, works with
% periodic trajectory where mouse runs through boudnaries and ends up on
% other side. Note, this doesn't really make sense to do, but just to check
% the conditions for the PCA NN being used by dordek. 

% INPUT
% - n_cells: number of place cells 
% - dim_x / dim_y: environment sizes
% - bounds_compressed: 1 or 0 depending on wether you want to have uniform
%                      or sanders type place cells
% OUTPUT 
% - ToroidalPlaceCellMaps: matrix with dimensions??--

%%
% Initialize an empty cell array to hold place cell data
cells = cell(1, n_cells);

% recalculate environment distances to boundaries without wall conditions
% get distances to environment walls 
x_dists = env.L == 1; % Identify wall locations as edge of env
env.x_dists = bwdist(x_dists); % Compute distance to the nearest wall in x direction
y_dists = env.L == 1; 
env.y_dists = bwdist(y_dists); 

% Find all valid positions in environment, (wall included for periodic
% boudnay condition)
bin_id = find(ones(env.dim_x, env.dim_y));

% Fixed standard deviations in x and y directions based on the fixed width
av_bound_dist = nanmean(env.dwmap, 'all');
fw = fieldWidth(av_bound_dist) / size_fact; % field width for uniform
sig_x = fw;
sig_y = fw;

% Loop over each cell to generate its place field map
for n = 1:n_cells
    % Center of the place field for the current cell
    mean_x = xy_field(n, 1); 
    mean_y = xy_field(n, 2);

    % take field widths relative to distance from boundary if sanders PC
    % type
    if boundary_compressed 
        sig_x = fieldWidth(env.x_dists(mean_y, mean_x)) / size_fact; 
        sig_y = fieldWidth(env.y_dists(mean_y, mean_x)) / size_fact;
    end

    % Initialize the rate map (place field map) as a zero matrix
    place_map = zeros(size(env.L));
    
    % Loop over all valid positions in the environment to compute the firing rate
    for i = 1:length(bin_id)
        [y, x] = ind2sub(size(place_map), bin_id(i)); 
        % calculate distance to PF centre with periodic boundary conditions
        diff_x=min(mod(x-mean_x,env.dim_x),mod(mean_x-x, env.dim_x));
        diff_y=min(mod(y-mean_y,env.dim_y),mod(mean_y-y, env.dim_y));
        fr = exp(-(diff_x)^2 / (2*sig_x^2)) * exp(-(diff_y)^2 / (2*sig_y^2)) / (2*pi*sig_x*sig_y);
        place_map(y, x) = fr;
    end
    
    % Normalize the place map by its maximum value
    place_map = place_map / max(place_map(:));
    % Store the normalized place map in the cells array
    cells{n}.fmap = place_map;

end

ToroidalPlaceCellMaps = cells; 

end

