%% Gridness Plot 
% Set the threshold value
threshold = 0.8;

% Directory containing the results
% output_dir = 'Results\gridness_results_rect';

% Load all .mat files from the directory
n_files = 20;

% Initialize arrays to store the percentages for each gridness type
% (uniform and tanni and expanded and standard gridness measures)
percent_square_U_e = zeros(n_files, 1); % expanded
percent_square_T_e = zeros(n_files, 1);
percent_hex_U_e = zeros(n_files, 1);
percent_hex_T_e = zeros(n_files, 1);

percent_square_U_s = zeros(n_files, 1); % standard
percent_square_T_s = zeros(n_files, 1);
percent_hex_U_s = zeros(n_files, 1);
percent_hex_T_s = zeros(n_files, 1);

% Base path
% base_path = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Results\gridness_results_rect\iteration_';
base_path = 'C:\Users\Elsa Marianelli\OneDrive - University College London\Documents\MATLAB\PCA on Sanders PC\gridness_results_rect\';

for i = 1:n_files 
    % Create the folder path dynamically
    file = ['iteration_', num2str(i), '.0'];
    folder_path = [base_path file];
    
    % Check if the folder exists
    if isfolder(folder_path)
        % List all .mat files in the folder
        files = dir(fullfile(folder_path, '*.mat'));
        
        % Load data
        % Uniform data
        file_path = fullfile(folder_path, files(2).name);
        fprintf('Loading file: %s\n', file_path); 
        data_U = load(file_path);

        % Tanni data
        file_path = fullfile(folder_path, files(1).name);
        fprintf('Loading file: %s\n', file_path); 
        data_T = load(file_path);
        
        % Calculate the percentage of values in column 1 (standard gridness) exceeding the threshold
        type = 2; % 2 for expanding 1 for standard gridness measure
        percent_square_U_e(i) = 100 * sum(data_U.gridness_square_U(:,type) > threshold) / length(data_U.gridness_square_U(:,type));
        percent_square_T_e(i) = 100 * sum(data_T.gridness_square_T(:,type) > threshold) / length(data_T.gridness_square_T(:,type));
        percent_hex_U_e(i) = 100 * sum(data_U.gridness_hex_U(:,type) > threshold) / length(data_U.gridness_hex_U(:,type));
        percent_hex_T_e(i) = 100 * sum(data_T.gridness_hex_T(:,type) > threshold) / length(data_T.gridness_hex_T(:,type));
        type = 1; % 2 for expanding 1 for standard gridness measure
        percent_square_U_s(i) = 100 * sum(data_U.gridness_square_U(:,type) > threshold) / length(data_U.gridness_square_U(:,type));
        percent_square_T_s(i) = 100 * sum(data_T.gridness_square_T(:,type) > threshold) / length(data_T.gridness_square_T(:,type));
        percent_hex_U_s(i) = 100 * sum(data_U.gridness_hex_U(:,type) > threshold) / length(data_U.gridness_hex_U(:,type));
        percent_hex_T_s(i) = 100 * sum(data_T.gridness_hex_T(:,type) > threshold) / length(data_T.gridness_hex_T(:,type));
    else
        fprintf('Folder not found: %s\n', folder_path);
    end
end

% Plot
list = {percent_square_U_e,percent_square_T_e, percent_hex_U_e, percent_hex_T_e};
list_2 = {percent_square_U_s,percent_square_T_s, percent_hex_U_s, percent_hex_T_s};

figure;
for i=1:length(list)
    mean_val = mean(list{i});
    bar(i*2, mean_val, 'k'); 
    hold on
    mean_val = mean(list_2{i});
    bar((i*2)-1, mean_val,'w'); 
    hold on
end

xticks(1:2:8)
xticklabels({'square uniform','square tanni', 'hex uniform', 'hex tanni'});
hold on;
ylabel('Percentage of Values > Threshold');
title(sprintf('Percentage of Expanded Gridness Values Over %.2f', threshold));