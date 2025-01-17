function [env] = GenerateEnv(polys, dim_x, dim_y, shape)
%GENERATEENV Generates an environment structure based on input polygons
%   Inputs:
%       polys: Cell array containing polygons that define boundaries in the environment
%       dim_x: Width of the environment (in pixels)
%       dim_y: Height of the environment (in pixels)
%       shape:  Can specify circle or trapezoid depending on environment
%       boundary type
%   Outputs:
%       env: Structure containing the environment map and various properties

if strcmp(shape, 'circle')
    % Initialize an empty map of size (dim_y+1) x (dim_x+1)
    map = zeros(dim_y+1, dim_x+1);
    
    % Define the center and radius of the circle
    center_x = (dim_x + 1) / 2;  % Center the circle horizontally
    center_y = (dim_y + 1) / 2;  % Center the circle vertically
    radius = min(dim_x, dim_y) / 2 - 2;  % Radius based on the map size, with some padding
    
    % Insert the circular boundary as a line into the map
    map = insertShape(map, 'Circle', [center_x, center_y, radius], 'Color', 'r');
    
    % Store the circle information in the environment structure 'env'
    env.circle = struct('center', [center_x, center_y], 'radius', radius);
    
    % Convert the map to a logical matrix where true represents the presence of a boundary
    map = logical(map(:,:,1));
elseif strcmp(shape, 'trapezoid')
    % Initialize an empty map of size (dim_y+1) x (dim_x+1)
    map = zeros(dim_y+1, dim_x+1);
    
    % Loop through each polygon in the input
    for i = 1:length(polys)
        % Insert the polygon as a line into the map
        % The polygon is defined as a sequence of (x,y) points, and the line is drawn in red
        map = insertShape(map, 'Line', polys{i}, 'Color', 'r');
    
        % Store the polygon points in the environment structure 'env'
        % The polygon points are reshaped into an N-by-2 matrix, where N is the number of points
        env.polys{i} = reshape(polys{i}, 2, length(polys{i}) / 2)';
    end
    
    % Convert the map to a logical matrix where true represents the presence of a boundary
    map = logical(map(:,:,1));
end
% Label connected components in the map
% 'L' is a matrix of the same size as 'map', where each connected component (region) is given a unique label
L = bwlabel(~map, 4); % 4-connected neighborhood

% Create a distance-to-wall map using the binary map
dwmap = map;
dwmap = bwdist(dwmap); % Compute the distance to the nearest boundary (wall)

% Set distances within the walls (label 1) to NaN
dwmap(L == 1) = NaN;

% Store the generated map, distance map, and dimensions in the output structure 'env'
env.map = map;       % The binary map of the environment
env.dwmap = dwmap;   % The distance-to-wall map
env.dim_x = dim_x + 1; % Width of the environment
env.dim_y = dim_y + 1; % Height of the environment
env.L = L;           % Labeled regions in the environment

% Plot the generated environment
figure;
imagesc(~map);       % Display the environment, flipping black/white for visibility
colormap gray;       % Use grayscale colormap
axis xy on;          % Turn on the axis and ensure correct orientation
title('Image of environment', 'FontWeight', 'normal');

end
