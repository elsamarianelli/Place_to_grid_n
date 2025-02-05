%% retry plotting from scratch 

%% [1] load data 
addpath('C:\Users\Elsa Marianelli\Documents\GitHub')
% base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\results_grid_maps_rect_1\';
n_iterations = 1;
base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Algorithmic_PCA\Functions\results_grid_maps_rect_2\';

% metricss settings for Grid analysis
in.GET_PERIM_FIELDS = false;
in.PLOT_ON = true;
in.FULL_GRIDNESS_RANGE = true; 
in.GET_MEAN_R_AND_SCALE = true;
in.FIELD_EXTENT_METHOD = true;
in.GET_PERIM_GRIDNESS_MASK  = true;
in.FIND_CENTROID = true;
    
for i = 1:n_iterations 
    % data directory
    file = ['iteration_', num2str(i), '.0'];
    folder_path = [base_path file];
    files = dir(fullfile(folder_path, '*.mat'));

    % Load data and get dims
    file_path = fullfile(folder_path, files(1).name);
    data = load(file_path);
    n_cells = size(data.cells_U, 2);
    x_dim = size(data.cells_U{1}.fmap, 1);
    y_dim = size(data.cells_U{1}.fmap, 2);
    x_split = [3; 3+(x_dim-4)/3; 3+(2*((x_dim-4)/3)); x_dim-2];
    y_split = [3; 3+(y_dim-4)/3; 3+(2*((y_dim-4)/3)); y_dim-2];  
    
    % reformat data - and get sac metrics for each bin
    grids_U = zeros(n_cells, 9, ((x_dim-4)/3)+1, ((y_dim-4)/3)+1);
    grids_T =  zeros(n_cells, 9, ((x_dim-4)/3)+1, ((y_dim-4)/3)+1);
    sacs_U = zeros(n_cells, 9, 2*((x_dim-4)/3)-1, 2*((y_dim-4)/3)-1);
    sacs_T = zeros(n_cells, 9, 2*((x_dim-4)/3)-1, 2*((y_dim-4)/3)-1);
    metrics_U = cell(n_cells, 9);
    metrics_T = cell(n_cells, 9);
    gridness_U_h = zeros(n_cells, 9, 2);
    gridness_T_h = zeros(n_cells, 9, 2);
    gridness_U_s = zeros(n_cells, 9, 2);
    gridness_T_s = zeros(n_cells, 9, 2);

    for cell_idx = 1:200

        fprintf('  Assessing Principle Component Grid %d Iteration %d \n', cell_idx, i);

        for x = 1:3
            for y = 1:3

                % get environment bin details
                x_idx = (x_split(x):x_split(x+1));
                y_idx = (y_split(y):y_split(y+1));
                box = ((x*3)-3)+(y);
                
                % bin maps
                grids_U(cell_idx, box, 1:length(x_idx),1:length(y_idx))...
                        = data.cells_U{cell_idx}.grid(x_idx, y_idx);
                grids_T(cell_idx, box, 1:length(x_idx),1:length(y_idx))...
                        = data.cells_T{cell_idx}.grid(x_idx, y_idx);
                
                % bin sacs
                sac_U = xPearson(data.cells_U{cell_idx}.grid(x_idx, y_idx));
                sac_T = xPearson(data.cells_T{cell_idx}.grid(x_idx, y_idx));
                sacs_U(cell_idx, box, 1:size(sac_U, 1), 1:size(sac_U, 2)) = sac_U;
                sacs_T(cell_idx, box, 1:size(sac_U, 1), 1:size(sac_U, 2)) = sac_T;
                
                % UNIFORM metrics
                try
                   % Compute square gridness (standard and expanded)
                    [stGrd, expGrd, ~] = multiGridness(sac_U, 'square');
                    gridness_U_s(cell_idx, box, 1) = stGrd;
                    gridness_U_s(cell_idx, box, 2) = expGrd;

                    % Compute hexagonal gridness (standard and expanded)
                    [stGrd, expGrd, ~] = multiGridness(sac_U, 'hexagon');
                    gridness_U_h(cell_idx, box, 1) = stGrd;
                    gridness_U_h(cell_idx, box, 2) = expGrd;
                   
                    % Check threshold for expGrd
                    if abs(expGrd) > 0.5
                        % Compute metrics
                        in.sac = sac_U;
                        metrics = autoCorrProps(in);
                        metrics_U{cell_idx, box} = metrics; 
                    else
                        metrics_U{cell_idx, box} = NaN; 
                    end

                catch ME  % Catch any error
                    warning('Error encountered at PC %d, box %d, Uniform: %s', cell_idx, box, ME.message);
                    metrics_U{cell_idx, box} = NaN; 
                end

                % TANNI metrics
                try 
                    
                    % Compute square gridness (standard and expanded)
                    [stGrd, expGrd, ~] = multiGridness(sac_T, 'square');
                    gridness_T_s(cell_idx, box, 1) = stGrd;
                    gridness_T_s(cell_idx, box, 2) = expGrd;

                    % Compute hexagonal gridness (standard and expanded)
                    [stGrd, expGrd, ~] = multiGridness(sac_T, 'hexagon');
                    gridness_T_h(cell_idx, box, 1) = stGrd;
                    gridness_T_h(cell_idx, box, 2) = expGrd;

                    % Check threshold for expGrd
                    if abs(expGrd) > 0.5
                        % Compute metrics
                        in.sac = sac_T;
                        metrics = autoCorrProps(in);
                        metrics_T{cell_idx, box} = metrics;
                    else
                        metrics_T{cell_idx, box} = NaN;
                    end

                catch ME  % Catch any error
                    warning('Error encountered at PC %d, box %d, Tanni: %s', cell_idx, box, ME.message);
                    metrics_T{cell_idx, box} = NaN;
                end

            end
        end
    end
    
end

% checking mean sac maps (like in stensola grid shearing paper within module)
meaned = squeeze(mean(sacs_U(cell_idx,:,:,:), 1));
disp(size(meaned))
figure;
for p = 1:9
    subplot(3, 3, p)
    imagesc(squeeze(meaned(p, :, :)));
end
% no common orientations so no pattern occurs when averaging 

% find cells and boxes with expanded gridness bove threshold of .8
threshold = 1.2;
mask = gridness_U(:,:,2) > threshold;
[cell_idx, box_idx] = find(mask);
above_threshold = gridness_U(mask);
cell_idx = unique(cell_idx);

% visualise bins across sac maps 
for i = 1:length(cell_idx)
    idx = cell_idx(i);
    sac = squeeze(sacs_U(idx, :, :,:));
    figure;
    for p = 1:9
        subplot(3, 3, p)
        imagesc(squeeze(sac(p, :, :)));
    end
end

% convert metric format and return metrics of interest
metrics = metrics_U;
[scale, orientation, gridness, ellipticity] = convert_metrics(metrics);

scale_mean = cellfun(@mean, scale); scale_std = cellfun(@std, scale);
orientation_mean = cellfun(@mean, orientation); orientation_std = cellfun(@std, orientation);
gridness_mean = cellfun(@mean, gridness); gridness_std = cellfun(@std, gridness);
% ellipticity_mean = cellfun(@mean, ellipticity); ellipticity_std =
% cellfun(@std, ellipticity); 

%% gridness measures
figure; plot(1:200, mean(gridness_U_h(:,:,1), 2)); hold on;
plot(1:200, mean(gridness_U_h(:,:,2), 2)); hold on; 
legend('std', 'exp'); title('Gridness score across PC Grids (Uni Hexagonal)')

figure; plot(1:200, mean(gridness_T_h(:,:,1), 2)); hold on;
plot(1:200, mean(gridness_T_h(:,:,2), 2)); hold on; 
legend('std', 'exp'); title('Gridness score across PC Grids (Tanni Hexagonal)')

figure; plot(1:200, mean(gridness_U_s(:,:,1), 2)); hold on;
plot(1:200, mean(gridness_U_s(:,:,2), 2)); hold on; 
legend('std', 'exp'); title('Gridness score across PC Grids (Uni Square)')

figure; plot(1:200, mean(gridness_T_s(:,:,1), 2)); hold on;
plot(1:200, mean(gridness_T_s(:,:,2), 2)); hold on; 
legend('std', 'exp'); title('Gridness score across PC Grids (Tanni Square)')
