%% Plotting saved results (from SR_eigdecomp.m)
%  Code to plots
% 1st row = Place cells which form SR successor features generated
% uniformly and using Sanders style boundary compression 
% 2nd row = Grid cells generated from PCA on trained SR (examples of range
% of module sizes)
% 3rd row = Spatial autocorrelogram for those grid cells
% 4th row Place cells from columns
% 5th row left = Grid scale averaged across iterations for each environmnet
% bin
% 5th row right = Gridness according to hexagonal/square measures and
% standard/expanded measures 

% path to ucl matlab grid funcitons here 
addpath('C:\Users\Elsa Marianelli\Desktop\analysis-matlab-master\analysis-matlab-master\GridAnalysis')
addpath('C:\Users\Elsa Marianelli\OneDrive - University College London\Documents\MATLAB\PCA on Sanders PC\gridness_results_rect\')

% Specify the directory where the results are saved
output_dir = 'gridness_results_rect'; 
n_iterations = 10;  % Number of iterations
% load results 
results = load_results(output_dir,n_iterations);

% set parameters for plotting
n_cells = 200;
it = 1;
n = 5;
n_iter = 10;
list_grid = [24, 32, 76, 98, 132];

% Tanni style results for plotting for 1 iteration
T_place_sf = results.PlaceCellsTanni{it} ;
T_grid = results.cells_T{it} ;
T_place_columns = results.cells_T{it} ;
[T_grid_scales, T_orientations] = get_binned_scales(results.cells_T, n_cells, n_iter, 'binned', 'hexagon');
[gridnessData_T, overThreshold_T] = get_binned_average_gridness(n_iter, results, '_T');


% Uniform style results for plotting
U_place_sf = results.PlaceCellsUni{it} ;
U_grid = results.cells_U{it} ;
U_place_columns = results.cells_U{it} ;
[U_grid_scales, U_orientations, U_sacs_stored, U_stats_stored] = get_binned_scales(results.cells_U, n_cells, n_iter, 'binned', 'hexagon');
[gridnessData_U, overThreshold_U] = get_binned_average_gridness(n_iter, results, '_U');


% [4] Autocorrelogram metric plots 
figure;
cell_idx = 30;
in.sac = U_sacs_stored(:,:,cell_idx,1);
subplot(1, 2, 1)
imagesc(U_sacs_stored(:,:,cell_idx,1))
in.FIND_CENTROID = true;
in.GET_PERIM_FIELDS = false;
in.PLOT_ON = true;
in.FULL_GRIDNESS_RANGE = true; 
in.GET_MEAN_R_AND_SCALE = true;
in.FIELD_EXTENT_METHOD = true;
in.GET_PERIM_GRIDNESS_MASK  = true;
hold
subplot(1, 2, 2)
metrics =  autoCorrProps(in);


% checking sac for high orientatios grids 
idx = find(U_orientations>120);
hi = results.cells_U{1};
for i = 1:10:100
    figure;
    imagesc(xPearson(hi{i}.grid)); 
end

%% plotting? 
% Dimensions of sac map
idx = overThreshold_U.gridness_hex_U{1}; 
for p =  1:length(idx)

[m, n, ~, ~, ~, ~] = size(U_sacs_stored);

% Initialize combined map for j = 1:3, k = 1:3 (3x3 grid of sac maps)
combined_map = zeros(3*m, 3*n);

% Specify i and iter values
i = idx(p); % Desired i
iter = 1; % Iteration to plot

% Loop over j (horizontal) and k (vertical)
for j = 1:3
    for k = 1:3
        % Extract the sac map for given i, iter, j, and k
        sac_map = U_sacs_stored(:, :, j, k, i, iter);
        %include stats info on map 
        stats = U_stats_stored(j,k, i, iter);% get stats for subregion 
        peakCoord = stats{1}.peakCoord;
        % sac_map(peakCoord(:,1), peakCoord(:,2)) = 1;
        % Determine subregion indices in combined map
        row_start = (j-1)*m + 1; % Vertical subregion start
        row_end = j*m;           % Vertical subregion end
        col_start = (k-1)*n + 1; % Horizontal subregion start
        col_end = k*n;           % Horizontal subregion end
        
        % Place sac map in the corresponding subregion
        % in.sac = U_sacs_stored(:,:,i,1);
        % metrics =  autoCorrProps(in);

        combined_map(row_start:row_end, col_start:col_end) = sac_map;
    end
end

% Plot the combined map
figure;
imagesc(combined_map);
axis image; % Keep aspect ratio
colorbar; % Add colorbar for reference
xlabel('Horizontal Position');
ylabel('Vertical Position');
title(num2str(idx(p)))

end

%% 
% % plotting tiled layout 
% tanni_figure = get_full_tiled_figure('Boundary warped Place Cells', T_place_sf, T_grid, ...
%     T_place_columns, T_grid_scales, gridnessData_T,  n, list_grid, n_cells, n_iter);
% 
% uni_figure = get_full_tiled_figure('Uniform Place Cells', U_place_sf, U_grid, ...
%     U_place_columns, U_grid_scales, gridnessData_U, n, list_grid, n_cells, n_iter);

% plot comparison of variability in normalised scales and orientations in each bin

%filter out grid cells below either square or hexagonal threshold
idx_U = [overThreshold_U.gridness_square_U{1}; overThreshold_U.gridness_hex_U{1}];
idx_T = [overThreshold_T.gridness_square_T{1}; overThreshold_T.gridness_hex_T{1}];
idx_U = unique(sort(idx_U, 'ascend'));
idx_T = unique(sort(idx_T, 'ascend'));
U_grid_scales_thresholded = U_grid_scales(:,:, idx_U);
T_grid_scales_thresholded = T_grid_scales(:,:, idx_T);
U_orientations_thresholded = U_orientations(:,:, idx_U);
T_orientations_thresholded = T_orientations(:,:, idx_T);

% plot
figure; 

tiledlayout(4, 4)
axis1 = nexttile(1, [2 2]);
[row_T, col_T, row_U, col_U]= plot_comparison_histograms(U_grid_scales_thresholded, T_grid_scales_thresholded, 9, 'scale', axis1);
title('Grid Scale', 'FontSize', 16)
legend('off')

nexttile(9)
imagesc(xPearson(results.cells_T{col_T}{row_T}.grid));
title('high variation')
xticks(''); yticks('');

nexttile(10)
imagesc(xPearson(results.cells_T{col_U}{row_U}.grid));
title('low variation')
xticks(''); yticks('');

axis2 = nexttile([2 2]);
[row_T, col_T, row_U, col_U] = plot_comparison_histograms(U_orientations_thresholded, T_orientations_thresholded, 15, 'orientation', axis2);
title('Grid Orientation', 'FontSize', 16)

nexttile(11)
imagesc(xPearson(results.cells_T{col_T}{row_T}.grid)); 
title('high variation')
xticks(''); yticks('');

nexttile(12)
imagesc(xPearson(results.cells_T{col_U}{row_U}.grid));
title('low variation')
xticks(''); yticks('');

% Plot grid scale distribution 
nexttile(13)
scales = T_grid_scales;
mean_g = nanmean(scales, 3); 
mean_g = (squeeze(mean_g))/71.71;
std_g = std(mean_g, 0, 3); 
mean_g = mean(mean_g, 3);
hold on;
bar(reshape(mean_g, 1, [])); % Average across cells
hold on;
errorbar(reshape(mean_g, 1, []), reshape(std_g, 1, [])./sqrt(n_iter), 'k.');
ylim([0.9 1.1])
box off

% Set axis properties
set(gca, 'FontSize', 10);  % Standardized font size for axis numbers
set(gcf, 'color', [0 0.7 0.7]);
set(gca, 'FontSize', 10);  % Make axis numbers more legible
xlabel('Environment Bin', 'FontSize', 12);  % Add x-axis label with smaller font size than the title
title('BC', 'FontSize', 16);  % Standard font size for headings

nexttile(14)
scales =U_grid_scales;
mean_g = nanmean(scales, 3); 
mean_g = (squeeze(mean_g))/71.71;
std_g = std(mean_g, 0, 3); 
mean_g = mean(mean_g, 3);
hold on;
bar(reshape(mean_g, 1, []), 'FaceColor', [0 0.7 0.7]); % Average across cells
hold on;
errorbar(reshape(mean_g, 1, []), reshape(std_g, 1, [])./sqrt(n_iter), 'k.');

% Set axis properties
set(gca, 'FontSize', 10);  % Standardized font size for axis numbers
set(gcf, 'color', 'w');
set(gca, 'FontSize', 10);  % Make axis numbers more legible
xlabel('Environment Bin', 'FontSize', 12);  % Add x-axis label with smaller font size than the title
title('Uni', 'FontSize', 16);  % Standard font size for headings
ylim([0.9 1.1])
box off

% Plot grid scale distribution 
nexttile(15)
orientations = T_orientations;
mean_g = nanmean(orientations, 3); 
mean_g = (squeeze(mean_g));
std_g = std(mean_g, 0, 3); 
mean_g = mean(mean_g, 3);

bar(reshape(mean_g, 1, [])); % Average across cells
hold on;
errorbar(reshape(mean_g, 1, []), reshape(std_g, 1, [])./sqrt(n_iter), 'k.');

% Set axis properties
set(gca, 'FontSize', 10);  % Standardized font size for axis numbers
set(gcf, 'color', 'w');
set(gca, 'FontSize', 10);  % Make axis numbers more legible
xlabel('Environment Bin', 'FontSize', 12);  % Add x-axis label with smaller font size than the title
title('BC', 'FontSize', 16);  % Standard font size for headings
ylim([35 45])
box off

nexttile(16)
orientations =U_orientations;mean_g = nanmean(orientations, 3); 
mean_g = (squeeze(mean_g));
std_g = std(mean_g, 0, 3); 
mean_g = mean(mean_g, 3);
bar(reshape(mean_g, 1, []), 'FaceColor', [0 0.7 0.7]); % Average across cells
hold on;
errorbar(reshape(mean_g, 1, []), reshape(std_g, 1, [])./sqrt(n_iter), 'k.');

% Set axis properties
set(gca, 'FontSize', 10);  % Standardized font size for axis numbers
set(gcf, 'color', 'w');
set(gca, 'FontSize', 10);  % Make axis numbers more legible
xlabel('Environment Bin', 'FontSize', 12);  % Add x-axis label with smaller font size than the title
title('Uni', 'FontSize', 16);  % Standard font size for headings
ylim([35 45])
box off