%% retry plotting from scratch 
%% [1] Load Data & Initialize Parameters
addpath('C:\Users\Elsa Marianelli\Documents\GitHub')
addpath('C:\Users\Elsa Marianelli\Documents\GitHub\Results_Grids')

% Base path for results
base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Results_Grids\results_grid_maps_rect_2\';

% Settings for Grid Analysis
in = struct('GET_PERIM_FIELDS', false, ...
            'PLOT_ON', true, ...
            'FULL_GRIDNESS_RANGE', true, ...
            'GET_MEAN_R_AND_SCALE', true, ...
            'FIELD_EXTENT_METHOD', true, ...
            'GET_PERIM_GRIDNESS_MASK', true, ...
            'FIND_CENTROID', true);

n_iterations = 1;

for i = 1:n_iterations
    % Load iteration data
    file = sprintf('iteration_%d.0', i);
    folder_path = fullfile(base_path, file);
    files = dir(fullfile(folder_path, '*.mat'));
    
    if isempty(files)
        warning('No .mat files found in %s', folder_path);
        continue;
    end
    
    % Load first file
    data = load(fullfile(folder_path, files(1).name));
    n_cells = numel(data.cells_U);
    [x_dim, y_dim] = size(data.cells_U{1}.fmap);
    
    % Compute bin boundaries
    x_split = round(linspace(3, x_dim-2, 4));  % Round to nearest integer
    y_split = round(linspace(3, y_dim-2, 4));  
   
    % Initialize data structures
    gridness_U_h = zeros(n_cells, 9, 2);
    gridness_T_h = zeros(n_cells, 9, 2);
    gridness_U_s = zeros(n_cells, 9, 2);
    gridness_T_s = zeros(n_cells, 9, 2);
    
    metrics_U = cell(n_cells, 9);
    metrics_T = cell(n_cells, 9);

    % Loop over cells and process gridness
    for cell_idx = 1:200
        fprintf('Processing PC Grid %d (Iteration %d)\n', cell_idx, i);

        for x = 1:3
            for y = 1:3
                box = (x-1)*3 + y;
                x_idx = x_split(x):x_split(x+1);
                y_idx = y_split(y):y_split(y+1);
                
                % Compute spatial autocorrelation
                sac_U = xPearson(data.cells_U{cell_idx}.grid(x_idx, y_idx));
                sac_T = xPearson(data.cells_T{cell_idx}.grid(x_idx, y_idx));

                % Compute gridness scores
                gridness_U_h(cell_idx, box, :) = computeGridness(sac_U, 'hexagon');
                gridness_U_s(cell_idx, box, :) = computeGridness(sac_U, 'square');
                gridness_T_h(cell_idx, box, :) = computeGridness(sac_T, 'hexagon');
                gridness_T_s(cell_idx, box, :) = computeGridness(sac_T, 'square');
                
                % Compute additional metrics if gridness is above threshold
                metrics_U{cell_idx, box} = computeMetrics(sac_U, in, gridness_U_s(cell_idx, box, 2), threshold);
                metrics_T{cell_idx, box} = computeMetrics(sac_T, in, gridness_T_s(cell_idx, box, 2), threshold);
            end
        end
    end
end

%% [2] Plot Gridness Measures Across Cells
figure;
titleFontSize = 16;  
axisFontSize = 14;   
lineWidth = 2;       
y_limit = [-0.5, 1.7];  
labels = {'Uni Hexagonal', 'BW Hexagonal', 'Uni Square', 'BW Square'};

grid_data = {gridness_U_h, gridness_T_h, gridness_U_s, gridness_T_s};

for idx = 1:4
    subplot(2,2,idx);
    plot(1:200, mean(grid_data{idx}(:,:,1), 2), 'b', 'LineWidth', lineWidth); hold on;
    plot(1:200, mean(grid_data{idx}(:,:,2), 2), 'r', 'LineWidth', lineWidth);
    title(labels{idx}, 'FontSize', titleFontSize, 'FontWeight', 'bold');
    ylim(y_limit);
    set(gca, 'FontSize', axisFontSize, 'Box', 'off');
end

% Create a single legend for all subplots
legend({'std', 'exp'}, 'FontSize', axisFontSize, 'Location', 'north', 'Orientation', 'Horizontal');
sgtitle('Gridness Scores Across PC Grids', 'FontSize', 18, 'FontWeight', 'bold');



% old code same thing...
% %% [1] load data 
% addpath('C:\Users\Elsa Marianelli\Documents\GitHub')
% addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Results_Grids'
% % base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\results_grid_maps_rect_1\';
% n_iterations = 1;
% base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Algorithmic_PCA\Functions\results_grid_maps_rect_2\';
% base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Results_Grids\results_grid_maps_rect_2\';
% % metricss settings for Grid analysis
% in.GET_PERIM_FIELDS = false;
% in.PLOT_ON = true;
% in.FULL_GRIDNESS_RANGE = true; 
% in.GET_MEAN_R_AND_SCALE = true;
% in.FIELD_EXTENT_METHOD = true;
% in.GET_PERIM_GRIDNESS_MASK  = true;
% in.FIND_CENTROID = true;
% 
% for i = 1:n_iterations 
%     % data directory
%     file = ['iteration_', num2str(i), '.0'];
%     folder_path = [base_path file];
%     files = dir(fullfile(folder_path, '*.mat'));
% 
%     % Load data and get dims
%     file_path = fullfile(folder_path, files(1).name);
%     data = load(file_path);
%     n_cells = size(data.cells_U, 2);
%     x_dim = size(data.cells_U{1}.fmap, 1);
%     y_dim = size(data.cells_U{1}.fmap, 2);
%     x_split = [3; 3+(x_dim-4)/3; 3+(2*((x_dim-4)/3)); x_dim-2];
%     y_split = [3; 3+(y_dim-4)/3; 3+(2*((y_dim-4)/3)); y_dim-2];  
% 
%     % reformat data - and get sac metrics for each bin
%     grids_U = zeros(n_cells, 9, ((x_dim-4)/3)+1, ((y_dim-4)/3)+1);
%     grids_T =  zeros(n_cells, 9, ((x_dim-4)/3)+1, ((y_dim-4)/3)+1);
%     sacs_U = zeros(n_cells, 9, 2*((x_dim-4)/3)-1, 2*((y_dim-4)/3)-1);
%     sacs_T = zeros(n_cells, 9, 2*((x_dim-4)/3)-1, 2*((y_dim-4)/3)-1);
%     metrics_U = cell(n_cells, 9);
%     metrics_T = cell(n_cells, 9);
%     gridness_U_h = zeros(n_cells, 9, 2);
%     gridness_T_h = zeros(n_cells, 9, 2);
%     gridness_U_s = zeros(n_cells, 9, 2);
%     gridness_T_s = zeros(n_cells, 9, 2);
% 
%     for cell_idx = 1:200
% 
%         fprintf('  Assessing Principle Component Grid %d Iteration %d \n', cell_idx, i);
% 
%         for x = 1:3
%             for y = 1:3
% 
%                 % get environment bin details
%                 x_idx = (x_split(x):x_split(x+1));
%                 y_idx = (y_split(y):y_split(y+1));
%                 box = ((x*3)-3)+(y);
% 
%                 % bin maps
%                 grids_U(cell_idx, box, 1:length(x_idx),1:length(y_idx))...
%                         = data.cells_U{cell_idx}.grid(x_idx, y_idx);
%                 grids_T(cell_idx, box, 1:length(x_idx),1:length(y_idx))...
%                         = data.cells_T{cell_idx}.grid(x_idx, y_idx);
% 
%                 % bin sacs
%                 sac_U = xPearson(data.cells_U{cell_idx}.grid(x_idx, y_idx));
%                 sac_T = xPearson(data.cells_T{cell_idx}.grid(x_idx, y_idx));
%                 sacs_U(cell_idx, box, 1:size(sac_U, 1), 1:size(sac_U, 2)) = sac_U;
%                 sacs_T(cell_idx, box, 1:size(sac_U, 1), 1:size(sac_U, 2)) = sac_T;
% 
%                 % UNIFORM metrics
%                 try
%                    % Compute square gridness (standard and expanded)
%                     [stGrd, expGrd, ~] = multiGridness(sac_U, 'square');
%                     gridness_U_s(cell_idx, box, 1) = stGrd;
%                     gridness_U_s(cell_idx, box, 2) = expGrd;
% 
%                     % Compute hexagonal gridness (standard and expanded)
%                     [stGrd, expGrd, ~] = multiGridness(sac_U, 'hexagon');
%                     gridness_U_h(cell_idx, box, 1) = stGrd;
%                     gridness_U_h(cell_idx, box, 2) = expGrd;
% 
%                     % Check threshold for expGrd
%                     if abs(expGrd) > 0.5
%                         % Compute metrics
%                         in.sac = sac_U;
%                         metrics = autoCorrProps(in);
%                         metrics_U{cell_idx, box} = metrics; 
%                     else
%                         metrics_U{cell_idx, box} = NaN; 
%                     end
% 
%                 catch ME  % Catch any error
%                     warning('Error encountered at PC %d, box %d, Uniform: %s', cell_idx, box, ME.message);
%                     metrics_U{cell_idx, box} = NaN; 
%                 end
% 
%                 % TANNI metrics
%                 try 
% 
%                     % Compute square gridness (standard and expanded)
%                     [stGrd, expGrd, ~] = multiGridness(sac_T, 'square');
%                     gridness_T_s(cell_idx, box, 1) = stGrd;
%                     gridness_T_s(cell_idx, box, 2) = expGrd;
% 
%                     % Compute hexagonal gridness (standard and expanded)
%                     [stGrd, expGrd, ~] = multiGridness(sac_T, 'hexagon');
%                     gridness_T_h(cell_idx, box, 1) = stGrd;
%                     gridness_T_h(cell_idx, box, 2) = expGrd;
% 
%                     % Check threshold for expGrd
%                     if abs(expGrd) > 0.5
%                         % Compute metrics
%                         in.sac = sac_T;
%                         metrics = autoCorrProps(in);
%                         metrics_T{cell_idx, box} = metrics;
%                     else
%                         metrics_T{cell_idx, box} = NaN;
%                     end
% 
%                 catch ME  % Catch any error
%                     warning('Error encountered at PC %d, box %d, Tanni: %s', cell_idx, box, ME.message);
%                     metrics_T{cell_idx, box} = NaN;
%                 end
% 
%             end
%         end
%     end
% 
% end
% 
% % checking mean sac maps (like in stensola grid shearing paper within module)
% meaned = squeeze(mean(sacs_U(cell_idx,:,:,:), 1));
% disp(size(meaned))
% figure;
% for p = 1:9
%     subplot(3, 3, p)
%     imagesc(squeeze(meaned(p, :, :)));
% end
% % no common orientations so no pattern occurs when averaging 
% 
% % find cells and boxes with expanded gridness bove threshold of .8
% gridness_U = gridness_U_s;
% threshold = .8;
% mask = gridness_U(:,:,2) > threshold;
% [cell_idx, box_idx] = find(mask);
% above_threshold = gridness_U(mask);
% cell_idx = unique(cell_idx);
% 
% % visualise bins across sac maps 
% for i = 1:length(cell_idx)
%     idx = cell_idx(i);
%     sac = squeeze(sacs_U(idx, :, :,:));
%     figure;
%     for p = 1:9
%         subplot(3, 3, p)
%         imagesc(squeeze(sac(p, :, :)));
%     end
% end
% 
% % convert metric format and return metrics of interest
% metrics = metrics_U;
% [scale, orientation, gridness, ellipticity] = convert_metrics(metrics);
% 
% scale_mean = cellfun(@mean, scale); scale_std = cellfun(@std, scale);
% orientation_mean = cellfun(@mean, orientation); orientation_std = cellfun(@std, orientation);
% gridness_mean = cellfun(@mean, gridness); gridness_std = cellfun(@std, gridness);
% % ellipticity_mean = cellfun(@mean, ellipticity); ellipticity_std =
% % cellfun(@std, ellipticity); 
% 
% %% plot gridness measures across cells
% figure;
% 
% titleFontSize = 16;  
% axisFontSize = 14;   
% lineWidth = 2;       
% y_limit = [-.5, 1.7];  
% 
% % Subplot 1: Gridness score across PC Grids (Uni Hexagonal)
% ax1 = subplot(2,2,1);
% plot(1:200, mean(gridness_U_h(:,:,1), 2), 'b', 'LineWidth', lineWidth); hold on;
% plot(1:200, mean(gridness_U_h(:,:,2), 2), 'r', 'LineWidth', lineWidth);
% title('Uni Hexagonal', 'FontSize', titleFontSize, 'FontWeight', 'bold');
% ylim(y_limit);
% set(gca, 'FontSize', axisFontSize, 'Box', 'off');
% 
% % Subplot 2: Gridness score across PC Grids (Tanni Hexagonal)
% ax2 = subplot(2,2,2);
% plot(1:200, mean(gridness_T_h(:,:,1), 2), 'b', 'LineWidth', lineWidth); hold on;
% plot(1:200, mean(gridness_T_h(:,:,2), 2), 'r', 'LineWidth', lineWidth);
% title('BW Hexagonal', 'FontSize', titleFontSize, 'FontWeight', 'bold');
% ylim(y_limit);
% set(gca, 'FontSize', axisFontSize, 'Box', 'off');
% 
% % Subplot 3: Gridness score across PC Grids (Uni Square)
% ax3 = subplot(2,2,3);
% plot(1:200, mean(gridness_U_s(:,:,1), 2), 'b', 'LineWidth', lineWidth); hold on;
% plot(1:200, mean(gridness_U_s(:,:,2), 2), 'r', 'LineWidth', lineWidth);
% title('Uni Square', 'FontSize', titleFontSize, 'FontWeight', 'bold');
% ylim(y_limit);
% set(gca, 'FontSize', axisFontSize, 'Box', 'off');
% 
% % Subplot 4: Gridness score across PC Grids (Tanni Square)
% ax4 = subplot(2,2,4);
% plot(1:200, mean(gridness_T_s(:,:,1), 2), 'b', 'LineWidth', lineWidth); hold on;
% plot(1:200, mean(gridness_T_s(:,:,2), 2), 'r', 'LineWidth', lineWidth);
% title('BW Square', 'FontSize', titleFontSize, 'FontWeight', 'bold');
% ylim(y_limit);
% set(gca, 'FontSize', axisFontSize, 'Box', 'off');
% 
% % Create a single legend for all subplots
% legend({'std', 'exp'}, 'FontSize', axisFontSize, 'Location', 'north', 'Orientation', 'Horizontal');
% 
% % Set global figure properties
% sgtitle('Gridness Scores Across PC Grids', 'FontSize', 18, 'FontWeight', 'bold');
% 
% % Remove top and right subplot borders
% set(ax1, 'Box', 'off'); set(ax2, 'Box', 'off');
% set(ax3, 'Box', 'off'); set(ax4, 'Box', 'off');
