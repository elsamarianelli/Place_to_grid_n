%% Grid Cell Generation and Parameter Sweep Analysis
% Author: Elsa Marianelli, UCL (2025) – zcbtetm@ucl.ac.uk
% Adapted from: Will de Cothi (2018) – Successor Representation code
%
% Description:
% This script runs a single parameter setting to generate grid cells from
% eigendecomposition of one of the following:
%   - Covariance matrix of place cell activity over time
%   - Learned Successor Representation (SR) matrix
%
% Using Place Cells which can be either:
%   - Uniformly distributed place cells
%   - "Tanni-style" place cells: smaller, more densely distributed near boundaries
%
% Outputs (saved per iteration):
%   - Place cell configurations (PlaceCellsUni / PlaceCellsTanni)
%   - Grid cell metrics:
%       - Hexagonal and square gridness
%       - Spatial scale
%       - Spatial autocorrelograms (SACs)
%       - 2D grid cell rate maps

%% [0] SETUP & PARAMETERS

% Add paths
addpath('/Users/elsamarianelli/Documents/GitHub')
addpath('/Users/elsamarianelli/Documents/GitHub/bound_warped_grids_new/Algorithmic_PCA/Functions/')
addpath('/Users/elsamarianelli/Documents/GitHub/bound_warped_grids_new/General_functions/')
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Algorithmic_PCA\Functions'
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n'
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\General_functions'

% Parameters to set....
use_SR  = false;             % true = SR matrix, false = Covariance matrix
use_traj = 'hasselmo';       % uniform = every bin sampled evenly (for covar)
                             % hasselmo = standard one with wall avoidance
                             % and speed angle changes 
                             % thigmotaxis = '#
trap_add = 0;                % set environment warping - 0 = normal rectangle, use 80 for trapezoid?
In.shape = 'trapezoid';     % environemtn shape - 'trapezoid' (rectangle or trapexoid) OR can be 'circle'

% Place Cell controls - both true = tanni, both false = uniform, can vary
% independantly % to run - true true, false false, true flase
pc_density = false;           % true = density varies with distance to boundary
pc_size    = false;          % true = size also varies relatively

% Additional parameters and environemnt details...should remain constant...
In.pf_width_cntrl = 2;          % Field width divisor (2 = narrower PCs)
n_iterations = 5;
In.n_cells = 250;               % number of place cells
In.n_steps = 150000;             % trajectory length
In.dim_x = 351;                 % environment dimensions
In.dim_y = 252;
In.n_polys = 1;
In.NumberOfPC = In.n_cells; % number of princ comps - should be the same as the number of place cells
In.bound_ctrl = 2;
            
% Folder naming tags - vary according to settings
base_dir = 'grids_data';
method_tag = 'SR'; if ~use_SR; method_tag = 'covar'; end
traj_tag = use_traj; 
density_tag = 'densityVaried' ; if ~ pc_density; density_tag = 'densityConstant'; end
size_tag = 'sizeVaried' ; if ~ pc_size; size_tag = 'sizeConstant'; end

% Saving 
output_dir = fullfile(base_dir, ['data_' method_tag, ...
    '_', traj_tag, '_', density_tag, '_', size_tag]);
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% [1] Generate Environment (same for all iterations)

fprintf('Generating Environment...\n');
In.polys{1} = [0 trap_add, 349 0, 349 250, 0 250-trap_add, 0 0] + 2;  
In = GenerateEnv(In);   % returns In strucutre now with environemnt info (In.env)

%% MAIN LOOP 
for iter = 1:n_iterations

    fprintf('--- Iteration %d/%d ---\n', iter, n_iterations);

    % Make folder
    subfolder = fullfile(output_dir, sprintf('iteration_%.1f', iter));
    if ~exist(subfolder, 'dir'); mkdir(subfolder); end

    %% [2] Load Trajectory

    fprintf('Load Trajectory...\n')

    if strcmp(traj_tag, 'uniform')
        % in order to completely remove time effects from the covariance
        % matrix since that is the goal "traj" input should be uniform
        % sampling of th environment...
        traj_cov = unique(combinations((1:In.dim_x)', (1:In.dim_y)'));  
        % double up so the number of steps is roughtly the same as whats used in SR
        traj = table2array([traj_cov; traj_cov]);
    elseif strcmp(traj_tag, 'hasselmo')    
        traj = load_premade_traj(iter);
        traj = traj(1:In.n_steps, :);  % trim to desired length
    elseif strcmp(traj_tag, 'thigmotaxis')%
            %To DO - generate different type of 
    end

    %% [3] Generate Place Cells
    fprintf('Generating Place Cells...\n');
    [Cells, In] = generate_place_cells(In, pc_density, pc_size);
    save(fullfile(subfolder, 'orig_place_cells.mat'), 'Cells');

    %% [4] Generate SR matrix or Covariance matrix and do PCA

    fprintf('Creating Grid Cells via %s\n', method_tag);

    if use_SR % SR Matrix
        M = ones(In.n_cells); R = ones(In.n_cells);
        [M, ~] = trainModel(Cells, M, R, traj, 1);
        cells = getPlace(Cells, M, In.env);
        cells = getGrid(cells, M, In.env, 'off');
    else % Covariance Matrix 
        [NeuronxEnvMat, NeuronxTimeMat] = reformat_firing_maps(Cells, traj);
        NeuronxTimeMat = NeuronxTimeMat - mean(NeuronxTimeMat, 2);
        eigvec = pca(NeuronxTimeMat', 'Algorithm', 'eig', 'Centered', false);
    end

    %% [5] Compute Grid Metrics

    fprintf('Computing Gridness Metrics...\n');
    GC_metrics = cell(In.n_cells, 1);
    h = waitbar(0, 'Processing PCs...');

    for PC = 1:In.n_cells
        waitbar(PC / In.n_cells, h);

        if use_SR
            map = cells{PC}.grid; % grid already calculated before
        else
        % in case of Covariance matrix, combine place fields according
        % to Principle component vectors and average
            map = comb_fields(NeuronxEnvMat, eigvec(:, PC)); 
        end

        map = map(3:end-2, 3:end-2);  % remove boundaried sections so you dont get line sin the sac
        sac = xPearson(map); % spatial autocorrelogram

        if all(isnan(sac(:)))
            GC_metrics{PC} = NaN; % somethimes happens for the first map, stops loop breaking
        else
            [metrics.stGrd_s, metrics.expGrd_s, metrics.scale_s] = multiGridness(sac, 'square', map, "off");
            [metrics.stGrd_h, metrics.expGrd_h, metrics.scale_h] = multiGridness(sac, 'hexagon', map, "off");
        
            metrics.map = map; metrics.sac = sac;
            GC_metrics{PC} = metrics; % save metrics
        end
    end

    close(h);

    %% [6] Save Results

    fprintf('Saving...\n');
    save(fullfile(subfolder, 'metrics_and_maps.mat'), 'In', 'GC_metrics', '-v7.3');

end

fprintf('\nAll done.\n');
