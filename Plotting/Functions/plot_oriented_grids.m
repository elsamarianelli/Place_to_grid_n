function plot_oriented_grids(angleMatrix)
    % angleMatrix is a 3x3 matrix where each element represents an angle in degrees
    % Check if the matrix is 3x3
    if ~isequal(size(angleMatrix), [3, 3])
        error('Input matrix must be 3x3.');
    end
    
    % Create a figure
    figure;
    hold on;
    
    % Define the size of each bin in the grid (3x3 grid)
    bin_size = 10;
    
    % Define the number of lines to plot in each bin
    num_lines_per_bin = 10;
    
    % Define the range of x-coordinates for the lines
    x_range = linspace(-bin_size, bin_size, 100);
    
    % Loop through each bin
    for row = 1:3
        for col = 1:3
            % Get the corresponding angle from the matrix
            angle = angleMatrix(row, col);
            
            % Convert angle to radians for rotation
            theta = deg2rad(angle);
            
            % Calculate the rotation matrix
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            
            % Determine the center of the current bin
            x_center = (col - 1) * bin_size - bin_size;
            y_center = (3 - row) * bin_size - bin_size;
            
            % Generate and rotate multiple parallel lines
            for line_num = 1:num_lines_per_bin
                % Generate a set of parallel lines by offsetting them in the y-direction
                y_offset = linspace(-bin_size / 2, bin_size / 2, num_lines_per_bin);
                y_line = y_offset(line_num) * ones(size(x_range));
                
                % Rotate the lines
                rotated_line = R * [x_range; y_line];
                X_rot = rotated_line(1, :);
                Y_rot = rotated_line(2, :);
                
                % Shift the rotated lines to the center of the current bin
                X_shifted = X_rot + x_center;
                Y_shifted = Y_rot + y_center;
                
                % Plot the rotated and shifted lines
                plot(X_shifted, Y_shifted, 'k', 'LineWidth', 1.5);
            end
        end
    end
    
    % Set axis properties
    axis equal;
    axis([-1.5*bin_size 1.5*bin_size -1.5*bin_size 1.5*bin_size]);
    xlabel('X');
    ylabel('Y');
    title('Continuous Oriented Lines in a 3x3 Grid');
    hold off;
end


