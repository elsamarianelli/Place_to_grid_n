%% Gridness Plot 
% Set the threshold value (you can modify this as needed)
threshold = input('Enter the threshold value: ');

% Directory containing the results
output_dir = 'gridness_result';

% Load all .mat files from the directory
files = dir(fullfile(output_dir, '*.mat'));
n_files = length(files);

% Initialize arrays to store the percentages for each gridness type
percent_square_U = zeros(n_files, 1);
percent_square_T = zeros(n_files, 1);
percent_hex_U = zeros(n_files, 1);
percent_hex_T = zeros(n_files, 1);

% Process each file
for i = 1:n_files
    % Load the data
    data = load(fullfile(output_dir, files(i).name));
    
    % Calculate the percentage of values in column 1 exceeding the threshold
    percent_square_U(i) = 100 * sum(data.gridness_square_U(:,1) > threshold) / length(data.gridness_square_U(:,1));
    percent_square_T(i) = 100 * sum(data.gridness_square_T(:,1) > threshold) / length(data.gridness_square_T(:,1));
    percent_hex_U(i) = 100 * sum(data.gridness_hex_U(:,1) > threshold) / length(data.gridness_hex_U(:,1));
    percent_hex_T(i) = 100 * sum(data.gridness_hex_T(:,1) > threshold) / length(data.gridness_hex_T(:,1));
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
