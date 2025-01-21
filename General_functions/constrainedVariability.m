
function [cost] = constrainedVariability(env, n_cells, xy_field, alpha, lower_threshold, upper_threshold)
    % variability and mean firing rate of average fmaps for all place cells
    [variability, mean_fr] = variabilityOfAveragedFiringRate(env, n_cells, xy_field, alpha);

    % Apply a penalty if mean firing rate is outside range
    if mean_fr < lower_threshold || mean_fr > upper_threshold
        penalty = 1e6 * ((max(0, mean_fr - upper_threshold))^2 + (max(0, lower_threshold - mean_fr))^2);
    else
        penalty = 0;  
    end
    %cost which you will adjust alpha to try and minimise
    cost = variability + penalty;
end
