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
%   - You can modulate the size and density seperately
%
% Outputs (saved per iteration):
%   - Place cell configurations (PlaceCellsUni / PlaceCellsTanni)
%   - Grid cell metrics:
%       - Hexagonal and square gridness
%       - Spatial scale
%       - Spatial autocorrelograms (SACs)
%       - 2D grid cell rate maps
%% [0] SETUP & PARAMETERS

% Parameters to set....
use_SR  = true;              % true = SR matrix, false = Covariance matrix
use_traj = 'generate';       % uniform = every bin sampled evenly (for covar)
                             % hasselmo = standard one with wall avoidance
                             % and speed angle changes 
% lengths of the shorter and longer parallel walls 0.2m and 
% 0.9m respectively with angled walls equal to 1.9m; 0.5 m height - krupic
% dimensions - dimensions here match this 
trap_add = 0;%(351-78)/2;   % set environment warping - 0 = normal rectangle, use 136 for trapezoid?
In.shape = 'trapezoid';        % environemtn shape - 'trapezoid' (rectangle or trapexoid) OR can be 'circle'
PCA_type = 'Standard';       % FISTA or sharp_Asymptotics or Non negative

% Place Cell controls - both true = tanni, both false = uniform, can vary
% independantly % to run - true true, false false, true flase
pc_density = true  ;       % true = density varies with distance to boundary
pc_size    = false   ;      % true = size also varies relatively

mean_firing_match = true;   % true = the size of generated place fields is set such that 
                            % the mean firing rate matches that of the varied setting

% Additional parameters and environemnt details...should remain constant...
In.pf_width_cntrl = 2;      % Field width divisor (2 = narrower PCs)
n_iterations = 5;
In.n_cells = 500;           %%CHANGE BACK T 500?     % number of place cells - set higher when NN = true
In.n_steps = 360000;        % trajectory length
In.dim_x = 351;%728;        % environment dimensions
In.dim_y = 351;             % y = 351 for square, 195 for trapezoid;
In.n_polys = 1;
In.NumberOfPC = 250;        % number of grids to generate
In.bound_ctrl = 2;
            
% Folder naming tags - vary according to settings
base_dir = 'grids_data_trap_env_SR_500_placecells'; 
method_tag = 'SR'; if ~use_SR; method_tag = 'covar'; end
traj_tag = use_traj;  
density_tag = 'densityVaried' ; if ~ pc_density; density_tag = 'densityConstant'; end
size_tag = 'sizeVaried' ; if ~ pc_size; size_tag = 'sizeConstant'; end
uniform_ctrl_tag = 'WidthControled' ; if ~ mean_firing_match; uniform_ctrl_tag = 'WidthNormal'; end

% saving 
user_docs = fullfile(getenv('USERPROFILE'), 'Documents');  % On Windows
dir = fullfile(user_docs, base_dir);
output_dir = fullfile(dir, ['data_' method_tag, '_', traj_tag, '_', density_tag, '_', size_tag, '_',uniform_ctrl_tag ]);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% pregenerate iteration files to populate in parfor loop
for iter = 1:n_iterations
    subfolder = fullfile(output_dir, sprintf('iteration_%.1f', iter));
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end
end

%% [1] Generate Environment (same for all iterations)

fprintf('Generating Environment...\n');
In.polys{1} = [0 trap_add, (In.dim_x-2) 0, (In.dim_x-2) (In.dim_y-2), 0 ((In.dim_y-2)-trap_add), 0 0] + 2;  
In = GenerateEnv(In);   % returns In strucutre now with environemnt info (In.env)
In.env.dim_x = In.dim_x; In.env.dim_y = In.dim_y;

%% MAIN LOOP  
% parpool();  % Start parallel pool
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) disp(msg));  % Print each message from workers

parfor iter = 1:n_iterations
    send(dq, sprintf('--- Iteration %d/%d started ---', iter, n_iterations));
    
    % Initialise so parfor doesn't get confused
    In_local = In;
    traj = []; cells = []; NeuronxEnvMat = []; NeuronxTimeMat = []; eigvec = [];

    % [2] Load Trajectory
    send(dq, sprintf('Iteration %d: Load/generate trajectory...', iter));
    if strcmp(traj_tag, 'uniform')
        traj_cov = unique(combinations((1:In_local.dim_x)', (1:In_local.dim_y)'));  
        % replicate so that the the traj is roughly the same length as 2
        % hour of hasselm trajectory (353808 steps vs 360000)
        traj = table2array([traj_cov; traj_cov;traj_cov; traj_cov ]);
    elseif strcmp(traj_tag, 'hasselmo')    
        traj = load_premade_traj(iter);
        traj = traj(1:In_local.n_steps, :);  %% need to regenerate as    env dims changed
    elseif strcmp(traj_tag, 'generate')
        traj = HasselmoForage(In_local.env, In_local.n_steps);
        traj = floor(traj);
    end
    
    % [3] Generate Place Cells
    send(dq, sprintf('Iteration %d: Generate place cells...', iter));

    if pc_density && pc_size
        [Cells, In_local] = generate_place_cells(In_local, pc_density, pc_size);
    else
        if mean_firing_match
        [TanniCells, In_local] = safe_generate_place_cells(In_local, true, true);
        [best_width, UniCells, ~] = tune_pf_width_to_match_activity( ...
            In_local, TanniCells, 30, 0.01, pc_density, pc_size);
        In_local.pf_width_cntrl = best_width;
        end
        [Cells, In_local] = safe_generate_place_cells(In_local, pc_density, pc_size);
    end

    output_dir_local = fullfile(getenv('USERPROFILE'), 'Documents', base_dir, ...
        ['data_' method_tag, '_', traj_tag, '_', density_tag, '_', size_tag, '_',uniform_ctrl_tag]);
    subfolder = fullfile(output_dir_local, sprintf('iteration_%.1f', iter));
    parsave2(fullfile(subfolder, 'orig_place_cells.mat'), Cells, In_local);
    
    % [4] Generate Matrix and PCA
    send(dq, sprintf('Iteration %d: Compute eigenvectors (%s)...', iter, method_tag));

    if strcmp(PCA_type, 'Standard')
        if use_SR
             M = ones(In_local.n_cells); R = ones(In_local.n_cells);
            [M, ~] = trainModel(Cells, M, R, traj, 1);
            cells = getPlace(Cells, M, In_local.env);
            cells = getGrid(cells, M, In_local.env); %, 'off');
        else % Covar
            [NeuronxEnvMat, NeuronxTimeMat] = reformat_firing_maps(Cells, traj);
            NeuronxTimeMat = NeuronxTimeMat - mean(NeuronxTimeMat, 2);
            eigvec = pca(NeuronxTimeMat', 'Algorithm', 'eig', 'Centered', false, ...
                'NumComponents', In_local.NumberOfPC);
        end
C = cov(NeuronxTimeMat');
imagesc(C)

    elseif strcmp(PCA_type, 'sharp_assymptotics') % Motanari and Richards 2014 algorithm
         %quite slow
        [NeuronxEnvMat, NeuronxTimeMat] = reformat_firing_maps(Cells, traj);
        zero_mean = 'spatial';
        eigvec = runNNPCA(NeuronxTimeMat, In_local.NumberOfPC, zero_mean);

    elseif strcmp(PCA_type, 'FISTA_DOG') %% not really working 

        % 1- Dordek DOG version for place cell size paramater
        % sweep, doesnt actually do PCA on a covar matrix 
        sigma_vector = [5, 10];  % tunable inner and outer sigmas for gaussian generation
        dims = [In_local.dim_y   In_local.dim_x]
        iterations = 2000;      
        nCells = 10;           % n grid cells to generate
        boundaries ='replicate';
        [NeuronxEnvMat, NeuronxTimeMat] = reformat_firing_maps(Cells, traj);
        [eigvec, GC_maps, Energy_array] = runFISTA_PCs(NeuronxEnvMat, dims, sigma_vector, nCells, boundaries, iterations)    
        
    elseif strcmp(PCA_type, 'FISTA_merged') 

        % 2- wang and wang FISTA version wih penalty that does run on covariance matrix
        % Assume you have a data matrix X (neurons x time)
        C = cov(NeuronxTimeMat');        
        K = 30;                % number of grid cells
        lambda = 0.05;         % contribution of DOG filter penalty term 
                               % can adgust gaussian sizes within the
                               % funciton
        [V, energy] = nn_pca_fista(C, K, max_iter, NeuronxEnvMat, lambda);
        
    end
    
    % [5] Compute Grid Metrics

    fprintf('Computing Gridness Metrics...\n');
    GC_metrics = cell(In_local.NumberOfPC, 1);
    % figure; 
    for PC = 1:In_local.NumberOfPC
        disp(PC);

        if use_SR
            map = cells{PC}.grid;
            % imagesc(map)
        else
            map = comb_fields(NeuronxEnvMat, eigvec(:, PC));
        end

        map = map(3:end-2, 3:end-2);
        sac = xPearson(map);

        metrics = struct();
        metrics.map = map;
        metrics.sac = sac;

        GC_metrics{PC} = metrics; 
    end

    % [6] Save Results
    send(dq, sprintf('Iteration %d: Saving results...', iter));
    parsave(fullfile(subfolder, 'metrics_and_maps.mat'), GC_metrics);
end
   
disp('All done.');
