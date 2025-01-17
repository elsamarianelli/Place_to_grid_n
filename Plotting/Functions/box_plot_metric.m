function box_plot_metric(corner_values, side_values, central_values, ax)
    if nargin < 4
        figure; % Create a new figure if no axis handle is provided
        ax = gca; % Use current axis
    end

    % Combine data into a single vector
    all_values = [corner_values(:); side_values(:); central_values(:)];

    % Create grouping variable
    num_corner = numel(corner_values);
    num_side = numel(side_values);
    num_centre = numel(central_values);
    group = [repmat({'Corner'}, num_corner, 1); 
             repmat({'Side'}, num_side, 1); 
             repmat({'Centre'}, num_centre, 1)];

    % Plot into the provided axis
    axes(ax); % Set the target axis
    boxplot(all_values, group);
    ylabel('Values');
end
