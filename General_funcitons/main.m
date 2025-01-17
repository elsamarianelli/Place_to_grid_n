%% eiegndecomp on Tanni style place cell covariance matrices
% (place cell density and gaussian width varying with distance to enclosure
% walls)

% [1] Generate environment  

% Define environment dimensions and the number of cells
dim_x = 252;
dim_y = 252;
% Define the environment boundaries with a slight offset
n_polys = 1;
polys = cell(n_polys, 1);
polys{1} = [0 0, (dim_x-2) 0, (dim_x-2) (dim_y-2), 0 (dim_y-2), 0 0] + 2;
env = GenerateEnv(polys, dim_x, dim_y, 'circle');

% [2] populate with place cells 
n_cells = 200;

% a) classic Place cell --> using code which ignores boundary inverse
% effect on sampling probablility

x_coords = randi([2, dim_x], n_cells, 1);
y_coords = randi([2, dim_y], n_cells, 1);
plot(x_coords, y_coords, '.')
xy_field = [x_coords, y_coords];
xy_field = get_sorted_coords(xy_field);
[xy_field, env, bin_prob] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, 2);

% take average distance from boundary as field width 
av_bound_dist = nanmean(env.dwmap, 'all');
fw = fieldWidth(av_bound_dist) / 3; % taken from sanders field width in x and y direction calc
PlaceCellsUni = generateUniformPCs(env, n_cells, xy_field, fw);
idx = floor((rand(25, 1)).*200);
plotPlace(PlaceCellsUni(idx'),5,5,env)

fmap_avg = get_fmap_avg(PlaceCellsUni);
[cells, fmap_avg] = my_compressed_PCs(env, n_cells, xy_field, bin_prob, 1, 1);

% Now fmap_avg contains the element-wise average across all 200 fmaps

% b) Tanni style Place cells, more frequent and with smaller firing fields
% near boundaries 

% sampling probability as 1/sqr root(x) distance from nearest wall
% boundary_effects = [0.5, 0.8, 1, 2, 100]; 
% % Create a figure for the plots
% figure;
% hold on
% % Loop through each boundary effect and generate the plots
% for i = 1:length(boundary_effects)
%     % Get the boundary effect
%     boundary_effect = boundary_effects(i);
% 
%     % Call the function to get place field centers and related data
%     [xy_field, env, ~] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, boundary_effect);
% 
%     % Create a subplot for each boundary effect
%     subplot(5, 1, i);  % Adjust the layout according to the number of boundary effects
%     plot(xy_field(:, 1), xy_field(:, 2), '.'); hold on % Plot the place field centers
%     title(['Boundary Effect: ', num2str(boundary_effect)]);
%     xlim([0 dim_x]);
%     ylim([0 dim_y]);
%     % xlabel('X');
%     % ylabel('Y');
%     hold on
% end

[xy_field, env, bin_prob] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, 2); 
figure; plot(xy_field(:,1), xy_field(:,2),'.')
% note xy_field given reorganised to have place cells in list closest to
% % each other for now
x = xy_field(:, 1); 
y = xy_field(:, 2);
gcf 
hold on
plot(x, y, '.', 'MarkerSize', 7)
PlaceCells = generateSanderPCs(env, n_cells, xy_field);
plotPlace(PlaceCells,10, 10,env)
fmap_avg = get_fmap_avg(PlaceCells);

cells = my_compressed_PCs(env, n_cells, xy_field);
plotPlace(cells,10, 10,env)

% [3] Generate trajectory 

n_steps = 360000;
traj = HasselmoForage(env, n_steps);
traj = round(traj);

% [4] Generate PC activity x time matrix 
PlaceCells = PlaceCellsUni;

x_indices = traj(:, 1);
y_indices = traj(:, 2);
fmaps = cellfun(@(PC) PC.fmap, PlaceCells, 'UniformOutput', false);

% Populate PCxTime with activity data
for i = 1:n_steps
    x = x_indices(i);
    y = y_indices(i);
    for j = 1:n_cells
        PC_fmap = fmaps{j};
        PCxTime(j, i) = PC_fmap(x, y);
    end
end

figure
imagesc(PCxTime)
hold on 
xlabel('steps in foraging trajectory')
ylabel('place cell activity')

% [5] use NeuronxTime matrix to generate Place cell activity correlation
% matrix and covariance matrix

% i) Correlation matrix
figure
subplot(1, 2, 1)
corr_M = corrcoef(PCxTime');
imagesc(corr_M)
hold on
xlabel('Place cells activity')
ylabel('Place cells activity')
title('Correlation')

% ii) Covariance matrix
subplot(1, 2, 2)
cov_M = cov(PCxTime');
imagesc(cov_M)
hold on
xlabel('Place cells activity')
ylabel('Place cells activity')
title('Covariance')

% [4/5 b] alternativ route: use trajectory to train SR matrix 
% Initialize the model parameters for training
M = ones(n_cells); % Initialize connection matrix
R = zeros(n_cells, 1); % Initialize reward vector
time_lag = 1; % Time lag of 100ms (a theta cycle)

% Train the model using the generated trajectory
[M, R] = trainModel(cells, M, R, traj, time_lag);
imagesc(M)
cov_M = M;
gridcells= getGrid(cells, M, env, 'off');

% [6] eigen decomposition on correlation matrix and sort output

[V, D] = eig(cov_M);
eigenvalues = diag(D);
% normalise eigenvalues according to largest 
eigenvalues = eigenvalues./max(eigenvalues);
% sort eigenvalues in D diagonal and reorder V matrix 
[sorted_eigenvalues, indices] = sort(eigenvalues, 'descend');
sorted_eigenvectors = V(:, indices);
sorted_eigenvectors = real(sorted_eigenvectors);

% [7] Plot sorted eigenvalues and variance explained by eigenvectors 

% Plot eigenvalues against their corresponding eigenvectors (4a in dordek paper)
figure
subplot(1, 2, 1)
threshold = 100;
plot(indices(end-threshold:end), eigenvalues(end-threshold:end), '.', 'MarkerSize', 10)
xlabel('# of Eigenvector')
ylabel('Eigenvalues')
hold on

% Calculate and plot the cumulative variance explained (4b in dordek paper)
subplot(1, 2, 2)
variance_explained = sorted_eigenvalues / sum(sorted_eigenvalues);
cumulative_variance = cumsum(variance_explained);
plot(cumulative_variance(1:threshold), '-');
xlabel('# of Eigenvector');
ylabel('Variance Explained');

% [8] Projecting high variance eigenvectors back into 2D space
cells_cov = get_grid(PlaceCells, cov_M, env);
% cells_corr = get_grid(PlaceCells, corr_M, env);

% [9] look at gridness for hexagonal grid cells and square grid cells

[cells_cov_square, gridness_square, scales_square] = analGrids(cells_cov, 'square');
[cells_cov_hex, gridness_hex, scales_hex] = analGrids(cells_cov, 'hexagon');

% [10] look at proportion of grid cells with ovre threshold gridnes for
% hexagonal and square gridness

proportion_square_tanni = get_good_grids(gridness_square, cells_cov_square, 0.7);
proportion_hex_tanni = get_good_grids(gridness_hex, cells_cov_hex, 0.7);
gridness_square(1,:)
figure;
x = ["Std square" "Std Hexagon" "tanni square" "tanni hexagon"];
y = [proportion_square proportion_hex, proportion_square_tanni proportion_hex_tanni];
bar(x, y)
% Visualize grid cells with high square scores
gridness = gridness_square;
gridness_threshold = 0.5;
cells = cells_cov_square;
ids = find(gridness(:,1) > gridness_threshold)'; % Find cells with gridness > gridness threshold 
for i = 1:length(ids)
    figure;
    subplot(1, 2, 1)
    imagesc(cells{ids(i)}.grid);
    colormap jet;
    set(gca, 'LineWidth', 2);
    set(gcf, 'color', 'w');
    set(gca, 'FontSize', 16);
    axis equal tight;
    title('Grid Cell Firing Rate Map');

    % and visualise autocorrelogram for that cell
    subplot(1, 2, 2)
    imagesc(cells{i}.sac);
    axis equal tight;
    colorbar;
    title('Spatial Autocorrelogram of grid cell Cell');
end

% Visualize grid cells with high hexagonal scores
gridness = gridness_hex;
gridness_threshold = 0.7;
ids = find(gridness(:,1) > gridness_threshold)'; % Find cells with gridness > gridness threshold 
cells = cells_cov_hex;
for i = 1:length(ids)
    figure;
    subplot(1, 2, 1)
    imagesc(cells{ids(i)}.grid);
    colormap jet;
    set(gca, 'LineWidth', 2);
    set(gcf, 'color', 'w');
    set(gca, 'FontSize', 16);
    axis equal tight;
    title('Grid Cell Firing Rate Map');

    % and visualise autocorrelogram for that cell
    subplot(1, 2, 2)
    imagesc(cells{i}.sac);
    axis equal tight;
    colorbar;
    title('Spatial Autocorrelogram of grid Cell');
end

% [11] looking at sclae scores in different portions enironemnt 
% Analyze grid scales within defined sub-regions of the environment
n_sub_regions = 3;
grid_scales = zeros(n_sub_regions, n_sub_regions, n_cells);
% X-axis limits for the sub-regions
x_lims = [1, 84; 85, 168; 169, 252];
% Y-axis limits for the sub-regions
y_lims = [1, 84; 85, 168; 169, 252];

% Loop over all cells and sub-regions to calculate grid scales
for i = 1:n_cells
    for j = 1:n_sub_regions
        for k = 1:n_sub_regions
            sac = xPearson(cells{i}.grid(y_lims(j, 1):y_lims(j, 2), x_lims(k, 1):x_lims(k, 2)));
            try
                stats = sacProps(sac); % Extract properties from autocorrelogram
                grid_scales(j, k, i) = stats.scale;
            catch
                grid_scales(j, k, i) = NaN; % Handle cases where scale cannot be calculated
            end
        end
    end
end

% Plot the average grid scales with error bars
figure;
bar(reshape(nanmean(grid_scales/71.71, 3), 1, [])); % Average across cells
hold on;
errorbar(reshape(nanmean(grid_scales/71.71, 3), 1, []), reshape(nanstd(grid_scales/71.71, [], 3), 1, [])./sqrt(n_cells), 'k.');
set(gca, 'LineWidth', 2);
set(gcf, 'color', 'w');
set(gca, 'FontSize', 16);
