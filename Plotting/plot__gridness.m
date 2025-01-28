%% Gridness Plot 
% Set the threshold value
threshold = 0.8;

% Directory containing the results
output_dir = 'Results\gridness_results_rect';

% Load all .mat files from the directory
n_files = 7;

% Initialize arrays to store the percentages for each gridness type
percent_square_U = zeros(n_files, 1);
percent_square_T = zeros(n_files, 1);
percent_hex_U = zeros(n_files, 1);
percent_hex_T = zeros(n_files, 1);

% Base path
base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Results\gridness_results_rect\iteration_';

% Loop through iteration numbers
for iteration = 1:n_files 
    % Create the folder path dynamically
    folder_path = [base_path, num2str(iteration), '.0'];
    
    % Check if the folder exists
    if isfolder(folder_path)
        % List all .mat files in the folder
        files = dir(fullfile(folder_path, '*.mat'));
        
        % Load each .mat file
        for i = 1:length(files)
            file_path = fullfile(folder_path, files(i).name);
            fprintf('Loading file: %s\n', file_path); 
            data = load(file_path);
            
           % Calculate the percentage of values in column 1 exceeding the threshold
            percent_square_U(i) = 100 * sum(data.gridness_square_U(:,1) > threshold) / length(data.gridness_square_U(:,1));
            percent_square_T(i) = 100 * sum(data.gridness_square_T(:,1) > threshold) / length(data.gridness_square_T(:,1));
            percent_hex_U(i) = 100 * sum(data.gridness_hex_U(:,1) > threshold) / length(data.gridness_hex_U(:,1));
            percent_hex_T(i) = 100 * sum(data.gridness_hex_T(:,1) > threshold) / length(data.gridness_hex_T(:,1));
        end
    else
        fprintf('Folder not found: %s\n', folder_path);
    end
end

% Plot the results
figure;
hold on;
plot(1:n_files, percent_square_U, '-bo', 'LineWidth', 1.5, 'DisplayName', 'Square_U');
plot(1:n_files, percent_square_T, '-ro', 'LineWidth', 1.5, 'DisplayName', 'Square_T');
plot(1:n_files, percent_hex_U, '-go', 'LineWidth', 1.5, 'DisplayName', 'Hex_U');
plot(1:n_files, percent_hex_T, '-mo', 'LineWidth', 1.5, 'DisplayName', 'Hex_T');
hold off;

% Add labels and legend
xlabel('Boundary Effect Index');
ylabel('Percentage of Values > Threshold');
legend('show');
title(sprintf('Percentage of Gridness Values Exceeding Threshold = %.2f', threshold));
