%% retry plotting from scratch 

%% [1] load data 
addpath('C:\Users\Elsa Marianelli\Documents\GitHub')
% base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\results_grid_maps_rect_1\';
n_iterations = 1;
base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Algorithmic_PCA\Functions\results_grid_maps_rect_1\';

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
    
    % reformat data - and get spatial autocorrelograms an dmetrics for each
    % bin
    grids_U = zeros(n_cells, 9, ((x_dim-4)/3)+1, ((y_dim-4)/3)+1);
    grids_T =  zeros(n_cells, 9, ((x_dim-4)/3)+1, ((y_dim-4)/3)+1);
    sacs_U = zeros(n_cells, 9, 2*((x_dim-4)/3)-1, 2*((y_dim-4)/3)-1);
    sacs_T = zeros(n_cells, 9, 2*((x_dim-4)/3)-1, 2*((y_dim-4)/3)-1);
    metrics_U = cell(n_cells, 9);
    metrics_T = cell(n_cells, 9);
    gridness_U = zeros(n_cells, 9, 2);
    gridness_T = zeros(n_cells, 9, 2);

    for cell_idx = 1:n_cells
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
                
                % check gridness (standard and expanded)
                [ stGrd, expGrd, ~ ] = multiGridness( sac_U, 'hexagon' );
                gridness_U(cell_idx, box, 1) = stGrd;
                gridness_U(cell_idx, box, 2) = expGrd;
                [ stGrd, expGrd, ~ ] = multiGridness( sac_T, 'hexagon' );
                gridness_T(cell_idx, box, 1) = stGrd;
                gridness_T(cell_idx, box, 2) = expGrd;

                % metrics 
                in.sac = sac_U;
                metrics_U{cell_idx, box}= autoCorrProps(in);
                
                in.sac = sac_T;
                metrics_T(cell_idx, box) = autoCorrProps(in);
                
                

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

    % % visualise bins across sac maps 
    % for i = 60:10:200
    %     sac = squeeze(sacs_U(i, :, :,:));
    %     figure;
    %     for p = 1:9
    %         subplot(3, 3, p)
    %         imagesc(squeeze(sac(p, :, :)));
    %     end
    % end

end


    % % old data
    % % files(:).name
    % % 1 = T_results (gridness for different measures)
    % % 2 = U_results ("")
    % % 3 = grid_place_cells_SR (cells_T and cells_U)
    % % 4 = place_cells (
    % 
    % imagesc(data_U.PlaceCellsUni{1, 100}.fmap)
    % imagesc(data_U.cells_U{1, 100}.fmap)
    % imagesc(data_U.hexagon_cells_U{1, 100}.fmap)
