function [best_width, UniCells, scaling_factor_history] = tune_pf_width_to_match_activity(In, TanniCells, max_iters, tolerance, pc_density, pc_size)
% TUNE_PF_WIDTH_TO_MATCH_ACTIVITY Adjusts pf_width_cntrl to equalize mean population activity
%
%   [best_width, UniCells, scaling_factor_history] = tune_pf_width_to_match_activity(In, TanniCells, max_iters, tolerance)
%
%   This function iteratively adjusts the 'In.pf_width_cntrl' parameter used in
%   uniform place field generation to match the mean population firing rate
%   of the given boundary-modulated (Tanni-style) place cells.
%
%   Inputs:
%       In         - Environment and parameter structure
%       TanniCells - Precomputed boundary-modulated place cells (from generate_place_cells)
%       max_iters  - Maximum number of adjustment iterations (e.g. 20)
%       tolerance  - Tolerance threshold for convergence (e.g. 0.01)
%
%   Outputs:
%       best_width             - The optimal value of In.pf_width_cntrl
%       UniCells               - Final uniform place cells with matched activity
%       scaling_factor_history - Vector of scaling factors per iteration (for plotting/debugging)

    % Initialization
    pf_width = In.pf_width_cntrl;  % Starting value
    step = 0.5;
    scaling_factor_history = [];
    best_width = NaN;

    for iter = 1:max_iters
        In.pf_width_cntrl = pf_width;

        % Generate uniform place cells
        [UniCells, In] = generate_place_cells(In, pc_density, pc_size);

        % Calculate scaling factor
        scaling_factor = balancing_population_activity(TanniCells, UniCells);
        scaling_factor_history(end+1) = scaling_factor;

        fprintf('Iter %d: pf_width_cntrl = %.4f → scaling factor = %.4f\n', iter, pf_width, scaling_factor);

        % Check convergence
        if abs(scaling_factor - 1) < tolerance
            best_width = pf_width;
            fprintf('✓ Converged at pf_width_cntrl = %.4f\n', pf_width);
            return;
        end

        % Adjust width based on whether uniform is too strong or too weak
        if scaling_factor < 1
            pf_width = pf_width + step;
        else
            pf_width = pf_width - step;
        end

        % Shrink step for more precise convergence
        step = step * 0.8;
    end

    warning('! Maximum iterations reached without full convergence. Closest width: %.4f (scaling factor: %.4f)', pf_width, scaling_factor);
    best_width = pf_width;
end
