function [X, iter, min_cost] = fista_general(grad, proj, Xhist, L, opts, calc_F)   
% function [X, iter, min_cost] = fista_general(grad,proj, Xinit, L, opts, calc_F)   
% * A Fast Iterative Shrinkage-Thresholding Algorithm for 
% Linear Inverse Problems.
% * Solve the problem: X = arg min_X F(X) = f(X) + lambda*g(X) where:
%   - X: variable, can be a matrix.
%   - f(X): a smooth convex function with continuously differentiable 
%       with Lipschitz continuous gradient `L(f)` (Lipschitz constant of 
%       the gradient of `f`).
%  INPUT:
%       grad   : a function calculating gradient of f(X) given X.
%       proj   : a function calculating pL(x) -- projection
%       Xinit  : a matrix -- initial guess.
%       L      : a scalar the Lipschitz constant of the gradient of f(X).
%       opts   : a struct
%           opts.lambda  : a regularization parameter, can be either a scalar or
%                           a weighted matrix.
%           opts.max_iter: maximum iterations of the algorithm. 
%                           Default 300.
%           opts.tol     : a tolerance, the algorithm will stop if difference 
%                           between two successive X is smaller than this value. 
%                           Default 1e-8.
%           opts.verbose : showing F(X) after each iteration or not. 
%                           Default false. 
%       calc_F: optional, a function calculating value of F at X 
%               via feval(calc_F, X). 
%  OUTPUT:
%      X        : solution
%      iter     : number of run iterations
%      min_cost : the achieved cost
% Modifications:
% 06/17/2016: set default value for opts.pos = false
% -------------------------------------
% Author: Tiep Vu, thv102, 4/6/2016
% (http://www.personal.psu.edu/thv102/)
% -------------------------------------
    Xinit_hist = Xhist;  % 5 prior inits

    % ---- Default opts ----
    if ~isfield(opts, 'max_iter'), opts.max_iter = 500; end
    if ~isfield(opts, 'tol'), opts.tol = 1e-6; end
    if ~isfield(opts, 'verbose'), opts.verbose = false; end
    if ~isfield(opts, 'lambda'), opts.lambda = 0; end
    if ~isfield(opts, 'pos'), opts.pos = false; end

    Linv = 1 / L;
    lambdaLiv = opts.lambda * Linv;

    % ---- Smart Initialization ----
    if size(Xinit_hist, 2) > 1
        Xinit = mean(Xinit_hist, 2);  % use mean of prior inits
    else
        Xinit = Xinit_hist;
    end
    Xinit = Xinit / norm(Xinit);     % normalize

    % ---- FISTA Initialization ----
    x_old = Xinit;
    y_old = Xinit;
    t_old = 1;
    iter = 0;
    cost_old = 1e10;
    opts_proj = opts;
    opts_proj.lambda = lambdaLiv;

    % ---- Main FISTA Loop ----
    while iter < opts.max_iter
        iter = iter + 1;

        % Gradient and proximal step
        grad_step = y_old - Linv * grad(y_old);
        x_new = proj(grad_step, opts_proj);

        % FISTA momentum update
        t_new = 0.5 * (1 + sqrt(1 + 4 * t_old^2));
        y_new = x_new + ((t_old - 1)/t_new) * (x_new - x_old);

        % ---- Adaptive Restarting (Beck & Teboulle) ----
        if (x_new - x_old)' * (y_old - x_new) > 0
            t_new = 1;
            y_new = x_new;
        end

        % ---- Convergence check ----
        err = norm(x_new - x_old, 1) / numel(x_new);
        if err < opts.tol
            break;
        end

        % ---- Early stopping on cost plateau ----
        if nargin >= 6
            cost_new = calc_F(x_new);
            if abs(cost_old - cost_new) < 1e-7
                break;
            end
            cost_old = cost_new;
        end

        % ---- Update for next iteration ----
        x_old = x_new;
        y_old = y_new;
        t_old = t_new;

        % ---- Verbose logging ----
        if opts.verbose
            if nargin >= 6
                fprintf('iter = %3d | cost = %.6f\n', iter, cost_old);
            else
                fprintf('iter = %3d\n', iter);
            end
        end
    end

    X = x_new;

    % Final cost
    if nargout == 3 && nargin >= 6
        min_cost = calc_F(X);
    elseif nargout == 3
        min_cost = NaN;
    end
end
