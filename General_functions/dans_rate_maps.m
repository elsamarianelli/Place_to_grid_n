%% Computing wall distance and firing rate maps using the meshgrid function
%  Daniel Bush, UCL (2025) drdanielbush@gmail.com

function [PlaceCellsUni, PlaceCellsTanni] = dans_rate_maps(env, 

%  Set the environment size
envSizeX        = env.dim_x;  % Environment width (columns)
envSizeY        = env.dim_y;  % Environment height (rows)

%  Generate a meshgrid
[x, y]         = meshgrid(1:envSizeX, 1:envSizeY);

%  Compute distance to each wall
w_dist(:,:,1)   = x;                % Distance to left wall
w_dist(:,:,2)   = envSizeY - y + 1; % Distance to bottom wall
w_dist(:,:,3)   = envSizeX - x + 1; % Distance to right wall
w_dist(:,:,4)   = y;                % Distance to top wall
w_dist          = min(w_dist,[],3); % Distance to closest wall

%  Set a place field centre and SD
pf_centre       = [20,30];          % x and y coordinates
pf_sd           = 10;               % SD of field size
peak_rate       = 10;
pf_dist         = sqrt((x - pf_centre(1)).^2 + (y - pf_centre(2)).^2);
rate_map        = peak_rate * exp(-pf_dist.^2 / (2 * pf_sd.^2));

figure, imagesc(rate_map), axis square

%  Version with separate SD along the x and y axes
x_dist(:,:,1)   = x;                % Distance to left wall
x_dist(:,:,2)   = envSizeX - x + 1; % Distance to right wall
y_dist(:,:,1)   = envSizeY - y + 1; % Distance to bottom wall
y_dist(:,:,2)   = y;                % Distance to top wall
x_dist          = min(x_dist,[],3); % Distance to closest horizontal wall
y_dist          = min(y_dist,[],3); % Distance to closest vertical wall

pf_centre       = [20,30];          % x and y coordinates
pf_x_sd         = 10;               % SD of field size in x
pf_y_sd         = 5;                % SD of field size in y
peak_rate       = 10;
rate_map        = peak_rate * exp(-((x - pf_centre(1)).^2 / (2 * pf_x_sd.^2) + ...
                                    (y - pf_centre(2)).^2 / (2 * pf_y_sd.^2)));

figure, imagesc(rate_map), axis square
