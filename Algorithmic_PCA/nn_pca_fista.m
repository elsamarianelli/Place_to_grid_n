function [V, energy_log] = nn_pca_fista(C, K, max_iter, NeuronxEnvMat, lambda)
% NN_PCA_FISTA
% Performs Non-Negative PCA using FISTA on a covariance matrix of NeuronxTime
% Implements the optimization used by Wang & Wang (2021) and Dordek (2016).
%
% INPUTS:
%   C         : [n x n] symmetric covariance matrix
%   K         : number of non-negative components to extract
%   max_iter  : number of FISTA iterations per component
%   lambda    : contribution of penalty term
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
    opts.verbose = true;  
    opts.adaptive_L = true;
    opts.random_init = true; % random initionalisation or average of last 5 maps as initialition 
    % (was an attempt to make faster but not sure it makes sense)
    tlo = tiledlayout(ceil(K/10), 10.*2, ...
    'TileSpacing', 'none', ...
    'Padding', 'compact');

    for k = 1:K
        fprintf(['Generating Principle Component...' num2str(k) '\n'])
    
        % --- Define optimization problem ---
        grad = @(v) -2 * C * v;  % gradient of -vᵀCv
        proj = @(v, ~) max(0, v) / norm(max(0, v));  % non-neg + normalize
        % cost = @(v) -v' * C * v;
        sz = 10;             % kernel size (e.g., 15x15)
        sigma_exc = 2;       % std dev of the excitatory Gaussian
        sigma_inh = 4;       % std dev of the inhibitory Gaussian
        dims = [ 253   352]; % EM - hardcoded for now fix
        strength = 0.8;      % inhibition strength relative to excitation
        periodicity_kernel = mexican_hat_kernel(sz, sigma_exc, sigma_inh, strength);
        cost = @(v) -v' * C * v + ...
            lambda * norm(conv2(reshape(comb_fields(NeuronxEnvMat, v), dims), periodicity_kernel, 'same'), 'fro')^2;

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
