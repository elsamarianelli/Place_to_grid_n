function cells = compressed_PCs_find_alpha(env, n_cells, xy_field)
% Goal: find a way to generate place cells which are compressed at
% boundaries with widths based on distance to boundaries, but where the
% relationship between boundary distance and width is modulated by an alpha
% paramater which is optimised to reduce diring rate variability across
% enviroment(averaged over all cells) with a constrained mean firing rate
% Based on idea that no mater compression to boundaries the averaged firing
% rate should be kept constant 

%   Inputs:
%       env: Environment structure containing maps and dimensions.
%       n_cells: Number of place cells to generate.
%       xy_field: Matrix of [x, y] coordinates for the center of each place field.
%       avgFiringRateFn: Function to compute the average firing rate across all rate maps.
%   Outputs:
%       cells: A cell array where each cell contains a structure with the place field map (fmap).
    
    % Define the range for the mean firing rate for penalising cost
    % funciton
    lower_threshold = 0.05;
    upper_threshold = 0.3; 
    
    % Omptimisation funciton which takes variable alpha values to adjust
    % relationship between width and place field centre to boundary
    optFun = @(alpha) constrainedVariability(env, n_cells, xy_field, alpha, lower_threshold, upper_threshold);
    
    % find alpha that minimises variability in firing rate ac
    optimalAlpha = fminsearch(optFun, 1);  % Initial guess for alpha = 1

    % generate final set of place cells
    cells = cell(1, n_cells);

    % Find all valid (non-wall) positions in the environment
    bin_id = find(env.L == 2);

    for n = 1:n_cells
        % Center of the place field for the current cell
        mean_x = xy_field(n, 1); 
        mean_y = xy_field(n, 2);

        % % Standard deviations in x and y directions scaled by alpha
        sig_x = fieldWidth(env.x_dists(mean_y, mean_x)) / (optimalAlpha); 
        sig_y = fieldWidth(env.y_dists(mean_y, mean_x)) / (optimalAlpha);
        % alternative linear relationship between SD and field width
        % sig_x = (env.x_dists(mean_y, mean_x)) .* (optimalAlpha); 
        % sig_y = (env.y_dists(mean_y, mean_x)) .* (optimalAlpha);

        % Initialize the rate map (place field map) as a zero matrix
        place_map = zeros(size(env.L));

        % Loop over all valid positions in the environment to compute the firing rate
        for i = 1:length(bin_id)
            [y, x] = ind2sub(size(place_map), bin_id(i)); 
            % Calculate firing rate using a 2D Gaussian function
            place_map(y, x) = exp(-(x-mean_x)^2 / (2*sig_x^2)) * exp(-(y-mean_y)^2 / (2*sig_y^2)) / (2*pi*sig_x*sig_y);
        end

        % Normalize the place map by its maximum value
        place_map = place_map / max(place_map(:));

        % Store the normalized place map in the cells array
        cells{n}.fmap = place_map;
    end

    fmap_avg = get_fmap_avg(cells);
    figure; imagesc(fmap_avg);
    figure; plotPlace(cells, 10, 10, env)

end
