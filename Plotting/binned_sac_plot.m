%% Code to generate some plots looking at grid metrics in different environment bins 

% [1] Load results
% [2] Looking at Hexagonality thresholds
% [3] Combined Sac maps and generated metric data
% [4] Generate data for orientation and ellipsis tilt plots
% [4b] Plot with mean and SD of orienation and ellipses tilt
% [5] Boxplots
% [6] Figure showing B/W heat maps of the environment bins with meaned metrics

addpath('C:\Users\Elsa Marianelli\Documents\GitHub\boundary_warped_place2grid\analysis-matlab-master\GridAnalysis')
addpath('C:\Users\Elsa Marianelli\Documents\GitHub\boundary_warped_place2grid\BasicFunctions')
addpath('C:\Users\Elsa Marianelli\Documents\GitHub\boundary_warped_place2grid\gridness_results_rect')
addpath('C:\Users\Elsa Marianelli\Documents\GitHub\boundary_warped_place2grid\analysis-matlab-master\Miscellaneous')
addpath('C:\Users\Elsa Marianelli\Documents\GitHub\boundary_warped_place2grid\results_grid_maps_rect')
addpath('C:\Users\Elsa Marianelli\OneDrive - University College London\Documents\MATLAB\PCA on Sanders PC\gridness_results_rect\')

% on mac
addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/analysis-matlab-master')
addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid')
addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/BasicFunctions/')

%% [1] Load results

output_dir = 'gridness_results_rect'; 
n_iterations = 1;  % number of iterations to load 
results = load_results(output_dir,n_iterations);

%% [2] find cells with high hexagonality thresholds
% square_gridness = results.gridness_square_U{1}; square_gridness = abs(square_gridness(:,1));
% hex_gridness = results.gridness_hex_U{1}; hex_gridness = abs(hex_gridness(:,1));
% figure; plot(square_gridness,hex_gridness, 'o'); hold on;
% plot([-0 2], [0 2]);
% xlabel('square gridness score')
% ylabel('hexagonal gridness score')
% 
% figure; 
% plot(1:200, abs(hex_gridness), 'o')
% [~, overThreshold_T] = get_binned_average_gridness(n_iterations, results, '_T');
[~, overThreshold_U] = get_binned_average_gridness(n_iterations, results, '_U');
% idx_hex_U = overThreshold_U.gridness_hex_U{2};
% idx_hex_T = overThreshold_T.gridness_hex_T{1};

%% [3] plotting binned combined sac maps for those cells and getting grid metrics for later plots

n_iter = 1;   % iterations of GC generation to run through 
idx = 50:150;  % GC PCs that you want to use (note that the lower and higher ones of the 200 have problems and metrics are harder to assess)

% initialise matrices to store grid metrics
orientation_stored = zeros(3, 3, length(idx), n_iter); % ORIENTATION: degrees anti clockwise from horizontal line where first SAC peak is
gridness_stored = zeros(3, 3, length(idx), n_iter);  % GRIDNESS: hexagonal gridness score as assessed in UCL-hippocmapus code 
                                                            %   (need to use eliptical/square measures in newer code)
scale_stored = zeros(3, 3, length(idx), n_iter); % SCALE: distance to closest peak
ellipticity_stored = zeros(3, 3, length(idx), n_iter); 
phi_stored = zeros(3, 3, length(idx), n_iter);

% initialise matrices to store grid metrics

orientation_stored_U = orientation_stored;
gridness_stored_U =gridness_stored;                                                       
scale_stored_U = scale_stored ;
ellipticity_stored_U= ellipticity_stored;
phi_stored_U = phi_stored;

orientation_stored = orientation_stored_U;
gridness_stored =gridness_stored_U;                                                       
scale_stored = scale_stored_U ;
ellipticity_stored= ellipticity_stored_U;
phi_stored = phi_stored_U;

% set limits for dividing environment into 9 subregions
x_lims = [3, 118; 119, 234; 235, 350]; % X-axis limits for sub-regions
y_lims = [3, 85; 86, 168; 169, 251]; % Y-axis limits for sub-regions

% Get metrics for each iteration
for iteration = 1:5
    disp(iteration);

    for p = 1:length(idx)

        % Pick an index of a grid cell map which has a high hexagonal gridness score
        i = idx(p);
        disp(i)
        cells = results.cells_U{iteration}; %% CHANGE
        
        % % uncomment to plot
        % figure;
        % t = tiledlayout(3, 3, 'TileSpacing', 'none', 'Padding', 'none');

        % Loop over j (vertical) and k (horizontal)
        for j = 1:3
            for k = 1:3
                try
                    % Get the sac for the current subregion
                   
                    % % uncomment to plot
                    % ax = nexttile;

                    sac = xPearson(cells{i}.grid(y_lims(j, 1):y_lims(j, 2), x_lims(k, 1):x_lims(k, 2)));
                    
                    % Pass sac to autoCorrProps
                    in.sac = sac;
                    in.PLOT_ON = true; hold on;
                    in.PLOT_Ellipse_ON = true;
                    metrics = autoCorrProps(in); % should be one with EM comments, which had added ellipse function

                    % Storing metrics for subregions
                    orientation_stored(j, k, p, iteration) = metrics.orientation;
                    gridness_stored(j, k, p, iteration) = metrics.gridness;
                    scale_stored(j, k, p, iteration) = metrics.scale;
                    ellipticity_stored(j, k, p, iteration) = metrics.ellipticity;
                    phi_stored(j, k, p, iteration) = metrics.ellipse_metrics.phi;

                catch ME
                    % Handle the error: store NaN in the output arrays
                    % fprintf('Error at Region (%d, %d, p=%d, iter=%d): %s\n', j, k, p, iteration, ME.message);
                    orientation_stored(j, k, p, iteration) = NaN;
                    gridness_stored(j, k, p, iteration) = NaN;
                    scale_stored(j, k, p, iteration) = NaN;
                    ellipticity_stored(j, k, p, iteration) = NaN;
                    phi_stored(j, k, p, iteration) = NaN;
                end

                % % Adjusting plot appearance
                % ax.XTick = [];   
                % ax.YTick = []; 
                % ax.XColor = 'none'; 
                % ax.YColor = 'none';

            end
        end
    end
end


%% 
%% [4] looking at meaned gridness scores

% take mean across iterations and across subgroups of environemnt (corner, sides and central)

%start point becaue lower PCs arent as good 
SP =1; % not sure what to set this as for now,
        % lower PCs have very poor gridness to begin with expecially whe you chunck them into smaller environemnts

mean_orientation = nanmean(orientation_stored(:,:,SP:end,1), 3);
mean_gridness = nanmean((gridness_stored(:,:,SP:end,1)), 3);
mean_scale = nanmean(scale_stored(:,:,SP:end,1), 3);
mean_ellip = nanmean(ellipticity_stored(:,:,SP:end,1), 3);
sd_orientation = nanstd(orientation_stored(:,:,SP:end,1),1, 3);
mean_phi = nanmean(phi_stored(:,:,SP:end,1), 3);
sd_phi = nanstd(phi_stored(:,:,SP:end,1),1, 3);

%% [4b] plot as in figure I , bu showing orienation and ellipses tilt seperately
% plots for each of the 9 subregions, the mean orientation/elipsis tilt,
% and their sds), the orientation is the angle from horizontal to the
% nearest anticlockwise peak, and the ellipsis tilt is taken as the angle
% to the semimajor axis? (nb should be converted from radians to degrees)

% plot orientation mean with SD for each subregion
figure;
ax = gca;
plot_mean_sd_orientations_binned(ax, mean_orientation, sd_orientation, 'Mean Orientation with SD');

% plot elipsis tilt mean with SD for each subregion
figure;
ax = gca;
plot_mean_sd_orientations_binned(ax, rad2deg(mean_phi), rad2deg(sd_phi), 'Mean Elipsis tilt with SD');

%% [5] boxplots of metrics
% as seen in stensola moser paper, taking corner, side, and central boxes of environment as 3 seperate categories
% currently having an issue with the extreme outliers coming from poor
% metrics, think especially with the binned environments autoCorr funciton
% struggles to get clean metrics

figure;

% Orientation metric
subplot(2, 2, 1);
[corner_values, side_values, central_values] = get_environment_sections(orientation_stored(:,:,:,1));
box_plot_metric(corner_values, side_values, central_values, gca);
title('Orientation');

% Gridness metric
subplot(2, 2, 2);
[corner_values, side_values, central_values] = get_environment_sections(gridness_stored(:,:,:,1));
box_plot_metric(corner_values, side_values, central_values, gca);
title('Gridness');

% Scale metric
subplot(2, 2,3);
[corner_values, side_values, central_values] = get_environment_sections(scale_stored(:,:,:,1));
box_plot_metric(corner_values, side_values, central_values, gca);
title('Scale');

% Scale metric
subplot(2, 2, 4);
[corner_values, side_values, central_values] = get_environment_sections(ellipticity_stored(:,:,:,1));
box_plot_metric(corner_values, side_values, central_values, gca);
title('ellipticity');

%% [6] Figure showing B/W heat maps of the environment bins with meaned metrics

% Create figure
figure;

% Plot Mean Orientation
subplot(2, 2, 1);
imagesc(mean_orientation);
colormap(gray); % Black and white colormap
colorbar; % Add colorbar if needed
title('Mean Orientation');
xlabel('X-axis');
ylabel('Y-axis');
axis equal tight; % Ensures aspect ratio and axis limits are tight

% Plot Mean Gridness
subplot(2, 2, 2);imagesc(mean_gridness);
colormap(gray); % Black and white colormap
colorbar; % Add colorbar if needed
title('Mean Gridness');
xlabel('X-axis');
ylabel('Y-axis');
axis equal tight; % Ensures aspect ratio and axis limits are tight

% Plot Mean Scale
subplot(2, 2, 3);
imagesc(mean_scale);
colormap(gray); % Black and white colormap
colorbar; % Add colorbar if needed
title('Mean Scale');
xlabel('X-axis');
ylabel('Y-axis');
axis equal tight; % Ensures aspect ratio and axis limits are tight

% Plot Mean Elllipticity
subplot(2, 2, 4);
imagesc(mean_ellip);
colormap(gray); % Black and white colormap
colorbar; % Add colorbar if needed
title('Mean ellipticity');
xlabel('X-axis');
ylabel('Y-axis');
axis equal tight; % Ensures aspect ratio and axis limits are tight

% Adjust layout to avoid overlap
set(gcf, 'Position', [100, 100, 1200, 400]); % Resize figure window




%%
% % same thing but for whole environemnt (not_binned)
% 
% idx = overThreshold_U.gridness_square_U{1};
% 
% 
% for p = 20:length(idx)
%     try
%         i = idx(p);
%         cells = results.cells_U{1};
% 
%         figure;
% 
%         sac = xPearson(cells{i}.grid);
% 
%         % Pass sac to autoCorrProps
%         in.sac = sac;
%         metrics = autoCorrProps(in); 
% 
%         ax.XTick = []; 
%         ax.YTick = []; 
%         ax.XColor = 'none'; 
%         ax.YColor = 'none';
% 
%         title(t, num2str(p)); % Set to 'off' if no title is needed
% 
%     catch ME
%         % Handle the error and move on to the next iteration
%         fprintf('Error')
%     end
% 
% end
