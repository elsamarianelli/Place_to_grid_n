function [In] = GenerateEnv(In)
%GENERATEIn.env Generates an In.environment structure based on input polygons
%   Inputs:
%       polys: Cell array containing polygons that define boundaries in the In.environment
%       In.dim_x: Width of the In.environment (in pixels)
%       In.dim_y: Height of the In.environment (in pixels)
%       In.shape:  Can specify circle or trapezoid depending on In.environment
%       boundary type
%   Outputs:
%       In.env: Structure containing the In.environment map and various properties

if strcmp(In.shape, 'circle')
    % Initialize an empty map of size (In.dim_y+1) x (In.dim_x+1)
    map = zeros(In.dim_y+1, In.dim_x+1);
    
    % Define the center and radius of the circle
    center_x = (In.dim_x + 1) / 2;  % Center the circle horizontally
    center_y = (In.dim_y + 1) / 2;  % Center the circle vertically
    radius = min(In.dim_x, In.dim_y) / 2 - 2;  % Radius based on the map size, with some padding
    
    % Insert the circular boundary as a line into the map
    map = insertShape(map, 'Circle', [center_x, center_y, radius], 'Color', 'r');
    
    % Store the circle information in the In.environment structure 'In.env'
    In.env.circle = struct('center', [center_x, center_y], 'radius', radius);
    
    % Convert the map to a logical matrix where true represents the presence of a boundary
    map = logical(map(:,:,1));
elseif strcmp(In.shape, 'trapezoid')
    % Initialize an empty map of size (In.dim_y+1) x (In.dim_x+1)
    map = zeros(In.dim_y+1, In.dim_x+1);
    
    % Loop through each polygon in the input
    for i = 1:length(In.polys)
        % Insert the polygon as a line into the map
        % The polygon is defined as a sequence of (x,y) points, and the line is drawn in red
        map = insertShape(map, 'Line', In.polys{i}, 'Color', 'r');
    
        % Store the polygon points in the In.environment structure 'In.env'
        % The polygon points are reshaped into an N-by-2 matrix, where N is the number of points
        In.polys{i} = reshape(In.polys{i}, 2, length(In.polys{i}) / 2)';
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

% Store the generated map, distance map, and dimensions in the output structure 'In.env'
In.env.map = map;       % The binary map of the In.environment
In.env.dwmap = dwmap;   % The distance-to-wall map
In.env.In.dim_x = In.dim_x + 1; % Width of the In.environment
In.env.In.dim_y = In.dim_y + 1; % Height of the In.environment
In.env.L = L;           % Labeled regions in the In.environment

% Plot the generated In.environment
% figure;
% imagesc(~map);       % Display the In.environment, flipping black/white for visibility
% colormap gray;       % Use grayscale colormap
% axis xy on;          % Turn on the axis and ensure correct orientation
% title('Image of In.environment', 'FontWeight', 'normal');

end
