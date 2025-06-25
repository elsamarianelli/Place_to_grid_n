function [V, energy_log] = nn_pca_fista(C, K, max_iter, NeuronxEnvMat)
% NN_PCA_FISTA
% Performs Non-Negative PCA using FISTA on a covariance matrix of NeuronxTime
% Implements the optimization used by Wang & Wang (2021) and Dordek (2016).
%
% INPUTS:
%   C         : [n x n] symmetric covariance matrix
%   K         : number of non-negative components to extract
%   max_iter  : number of FISTA iterations per component
%
% OUTPUTS:
%   V         : [n x K] matrix of non-negative eigenvectors
%   energy_log: [max_iter x K] energy of each component over time


    n = size(C,1);
    V = zeros(n, K);
    energy_log = nan(max_iter, K);

    L = 2 * max(eig(C));  % Lipschitz constant of -vᵀCv gradient

    % FISTA options
    opts.lambda = 0;         % no extra penalty
    opts.max_iter = max_iter;
    opts.tol = 1e-4;        % tolerance for breaing search 
    opts.verbose = false;
    tlo = tiledlayout(ceil(K/10), 10.*2, ...
    'TileSpacing', 'none', ...
    'Padding', 'compact');

    for k = 1:K
        fprintf(['Generating Principle Component...' num2str(k) '\n'])
    
        % --- Define optimization problem ---
        grad = @(v) -2 * C * v;  % gradient of -vᵀCv
        proj = @(v, ~) max(0, v) / norm(max(0, v));  % non-neg + normalize
        cost = @(v) -v' * C * v;

        % --- Initialize ---
        x0 = abs(randn(n, 1));
        x0 = x0 / norm(x0);

        % --- Run FISTA ---
        Xhist = abs(randn(n, 5));  % 5 prior inits
        [vk, ~, energy] = fista_general(grad, proj, Xhist, L, opts, cost);
        V(:, k) = vk;
        energy_log(1:length(energy), k) = energy;

        % --- Deflation: remove contribution from current component ---% rank-1 deflation (Wang & Wang style)
        Cv = C * vk;
        C = C - Cv * vk';  

        % Visualise
        ax = nexttile(tlo, k*2-1);
        imagesc(ax, comb_fields(NeuronxEnvMat, V(:, k)));
        axis(ax, 'image', 'off');
        drawnow;  
        ax = nexttile(tlo, k*2);
        imagesc(ax, xPearson(comb_fields(NeuronxEnvMat, V(:, k))));
        axis(ax, 'image', 'off');
        drawnow;  
        
    end
end
