% Load the saved data (replace 'gridness_data.mat' and 'scale_data.mat' with your file names)
% Assuming the data is saved in arrays: gridness_scores and scale_scores
% Each array has dimensions: (n_bins, 3), where 3 represents the three iterations

load('gridness_data.mat');  % Load the binned gridness data
load('scale_data.mat');     % Load the binned scale data

% Assume gridness_scores and scale_scores are matrices of size [n_bins x 3]

% Number of bins (assuming both datasets have the same number of bins)
n_bins = size(gridness_scores, 1);

% Calculate the mean across the three iterations for gridness and scale
mean_gridness = mean(gridness_scores, 2);
mean_scale = mean(scale_scores, 2);

% Create a new figure for the gridness scores
figure;
hold on;
% Plot each individual iteration (gridness)
for i = 1:3
    plot(1:n_bins, gridness_scores(:, i), '--o', 'DisplayName', ['Iteration ' num2str(i)]);
end
% Plot the mean gridness scores
plot(1:n_bins, mean_gridness, '-ko', 'LineWidth', 2, 'DisplayName', 'Mean Gridness');
hold off;

% Add labels and title for gridness
xlabel('Bin');
ylabel('Square Gridness');
title('Binned Square Gridness (Mean and Iterations)');
legend('show');
grid on;

% Create a new figure for the scale scores
figure;
hold on;
% Plot each individual iteration (scale)
for i = 1:3
    plot(1:n_bins, scale_scores(:, i), '--o', 'DisplayName', ['Iteration ' num2str(i)]);
end
% Plot the mean scale scores
plot(1:n_bins, mean_scale, '-ko', 'LineWidth', 2, 'DisplayName', 'Mean Scale');
hold off;

% Add labels and title for scale
xlabel('Bin');
ylabel('Square Scale');
title('Binned Square Scale (Mean and Iterations)');
legend('show');
grid on;
