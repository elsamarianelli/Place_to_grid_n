function [corner_values, side_values, central_values] = get_environment_sections(metric_stored)
% function to seperate metrics matrix for each subregion into corner, side
% and central matrix 
    data = metric_stored;
    
    % Preallocate matrices for results
    corner_values = zeros(4, size(data, 3)); % For corner values
    side_values = zeros(4, size(data, 3));   % For side values
    
    % Extract corner values
    corner_indices = [1, 1; 1, 3; 3, 1; 3, 3]; % Indices for corners (row, col)
    for i = 1:4
        corner_values(i, :) = squeeze(data(corner_indices(i, 1), corner_indices(i, 2), :));
    end

    % Extract central value
    central_values = squeeze(data(2, 2, :))'; % Center of 3x3 matrix

    % Extract side values
    side_indices = [1, 2; 2, 1; 2, 3; 3, 2]; % Indices for sides (row, col)
    for i = 1:4
        side_values(i, :) = squeeze(data(side_indices(i, 1), side_indices(i, 2), :));
    end

    corner_values = mean(corner_values, 1);
    side_values = mean(side_values, 1);
end
