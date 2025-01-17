% Set parameters 
dim_x = 252;
dim_y = 252;
n_cells = 200;
boundary_effects = [0.5, 0.8, 1, 2, 100]; 
output_dir = 'gridness_results_trap';                    
n_iterations = 3; % Number of iterations for each boundary effect

% Load the data for all boundary effects and iterations at the beginning
all_data = struct();
for i = 1:length(boundary_effects)
    b_effect = boundary_effects(i);
    all_data(i).boundary_effect = b_effect;
    all_data(i).square_gridness_all = [];
    all_data(i).square_gridness_exp_all = [];
    all_data(i).square_scale_all = [];

    for iter = 1:n_iterations
        data_file = fullfile(output_dir, sprintf('boundary_effect_%.1f/results_iter_%d.mat', b_effect, iter));
        data = load(data_file);

        % Collect gridness and scale scores
        all_data(i).square_gridness_all = [all_data(i).square_gridness_all, data.gridness_square_T(:, 1)];
        all_data(i).square_gridness_exp_all = [all_data(i).square_gridness_exp_all, data.gridness_square_T(:, 2)];
        all_data(i).square_scale_all = [all_data(i).square_scale_all, data.scale_square_T];
    end
end

% Load uniform data 
uniform_data = struct();
uniform_data.square_gridness_all = [];
uniform_data.square_gridness_exp_all = [];
uniform_data.square_scale_all = [];

for iter = 1:n_iterations
    data_file = fullfile(output_dir, sprintf('uniform/results_iter_%d.mat', iter));
    data = load(data_file);

    % Collect gridness and scale scores for the uniform data
    uniform_data.square_gridness_all = [uniform_data.square_gridness_all, data.gridness_square_U(:, 1)];
    uniform_data.square_gridness_exp_all = [uniform_data.square_gridness_exp_all, data.gridness_square_U(:, 2)];
    uniform_data.square_scale_all = [uniform_data.square_scale_all, data.scale_square_U];
end

% Create a new figure for the merged plots
figure;
hold on;

% Plot settings for gridness and scale
edges_grid = linspace(-1.5, 1.5, 21); 
edges_scale = linspace(0, 300, 21);
colors = lines(length(boundary_effects));  % Generate distinct colors for each boundary effect

% Iterate through each boundary effect and plot the first three columns (Gridness, Expanded Gridness, Scale)
for i = 1:length(boundary_effects)
    b_effect = boundary_effects(i);
    
    % If boundary effect is 100, use uniform data
    if b_effect == 100
        % Gridness for uniform data
        subplot(6, 4, ((i*4)-3));  % For gridness
        histogram(mean(uniform_data.square_gridness_all, 2), edges_grid, 'LineWidth', 1.5);
        ylim([0 100])
        hold on;

        % Expanded Gridness for uniform data
        subplot(6, 4, ((i*4)-2));  % For expanded gridness
        histogram(mean(uniform_data.square_gridness_exp_all, 2), edges_grid, 'LineWidth', 1.5);
        ylim([0 100])
        hold on;

        % Scale for uniform data
        subplot(6, 4, ((i*4)-1));  % For scale
        histogram(mean(uniform_data.square_scale_all, 2), edges_scale, 'LineWidth', 1.5);
        ylim([0 70])
        hold on;
        
        % Right column for place field centers (Uniform)
        % Call the function to get place field centers for the uniform condition
        [xy_field, env, ~] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, b_effect); % Uniform version
        subplot(6, 4, i*4);  % Adjust the layout for place field centers
        plot(xy_field(:, 1), xy_field(:, 2), '.'); hold on;
        title('Uniform Data');
        xlim([0 dim_x]);
        ylim([0 dim_y]);
        
    else
        % Gridness
        subplot(6, 4, ((i*4)-3));  % For gridness
        histogram(mean(all_data(i).square_gridness_all, 2), edges_grid, 'LineWidth', 1.5);
        ylim([0 100])
        hold on;

        % Expanded Gridness
        subplot(6, 4, ((i*4)-2));  % For expanded gridness
        histogram(mean(all_data(i).square_gridness_exp_all, 2), edges_grid, 'LineWidth', 1.5);
        ylim([0 100])
        hold on;

        % Scale
        subplot(6, 4, ((i*4)-1));  % For scale
        histogram(mean(all_data(i).square_scale_all, 2), edges_scale, 'LineWidth', 1.5);
        ylim([0 70])
        hold on;

        % Right column for place field centers
        [xy_field, env, ~] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, b_effect);
        subplot(6, 4, i*4);  % Adjust the layout for place field centers
        plot(xy_field(:, 1), xy_field(:, 2), '.'); hold on;
        title(['Boundary Effect: ', num2str(b_effect)]);
        xlim([0 dim_x]);
        ylim([0 dim_y]);
    end
end


% ** Final row: k-density for Gridness, Expanded Gridness, and Scale **
for i = 1:length(boundary_effects)
    b_effect = boundary_effects(i);

    % Gridness - k-density
    subplot(6, 4, ((6*4)-3));  % Final row, column 1 for Gridness
    [f_grid, xi_grid] = ksdensity(mean(all_data(i).square_gridness_all, 2));
    plot(xi_grid, f_grid, 'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;

    % Expanded Gridness - k-density
    subplot(6, 4, ((6*4)-2));  % Final row, column 2 for Expanded Gridness
    [f_grid_exp, xi_grid_exp] = ksdensity(mean(all_data(i).square_gridness_exp_all, 2));
    plot(xi_grid_exp, f_grid_exp, 'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;

    % Scale - k-density
    subplot(6, 4, ((6*4)-1));  % Final row, column 3 for Scale
    [f_scale, xi_scale] = ksdensity(mean(all_data(i).square_scale_all, 2));
    plot(xi_scale, f_scale, 'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;
end

% Set the titles for the first row of the figure
subplot(6, 4, 1);  
title('Gridness');
subplot(6, 4, 2);  
title('Expanded Gridness');
subplot(6, 4, 3);  
title('Scale');
subplot(6, 4, 4);  
title('Place Field Centers');

% Set the title for the final row (k-density)
subplot(6, 4, ((6*4)-3));
xlabel('Gridness Score');
ylabel('Density');
% title('K-Density: Gridness');

subplot(6, 4, ((6*4)-2));
xlabel('Expanded Gridness Score');
ylabel('Density');
% title('K-Density: Expanded Gridness');

subplot(6, 4, ((6*4)-1));
xlabel('Scale Score');
ylabel('Density');
% title('K-Density: Scale');

% Add labels for each row (boundary effect value) on the left side
for i = 1:length(boundary_effects)
    effect = boundary_effects(i);
    if effect == 100
        ylabel('Uniform');
    else
        subplot(6, 4, ((i*4)-3));
        ylabel(num2str(effect));
    end
end

% Adjust layout
hold off;


% % Set parameters 
% dim_x = 252;
% dim_y = 252;
% n_cells = 200;
% n_steps = 200000;
% boundary_effects = [0.5, 0.8, 1, 2, 100]; 
% output_dir = 'gridness_results';                    
% n_iterations = 3; % Number of iterations for each boundary effect
% 
% % Load the data for all boundary effects and iterations at the beginning
% all_data = struct();
% for i = 1:length(boundary_effects)
%     b_effect = boundary_effects(i);
%     all_data(i).boundary_effect = b_effect;
%     all_data(i).square_gridness_all = [];
%     all_data(i).square_gridness_exp_all = [];
%     all_data(i).square_scale_all = [];
% 
%     for iter = 1:n_iterations
%         data_file = fullfile(output_dir, sprintf('boundary_effect_%.1f/results_iter_%d.mat', b_effect, iter));
%         data = load(data_file);
% 
%         % Collect gridness and scale scores
%         all_data(i).square_gridness_all = [all_data(i).square_gridness_all, data.gridness_square_T(:, 1)];
%         all_data(i).square_gridness_exp_all = [all_data(i).square_gridness_exp_all, data.gridness_square_T(:, 2)];
%         all_data(i).square_scale_all = [all_data(i).square_scale_all, data.scale_square_T];
%     end
% end
% 
% % Load uniform data 
% uniform_data = struct();
% uniform_data.square_gridness_all = [];
% uniform_data.square_gridness_exp_all = [];
% uniform_data.square_scale_all = [];
% 
% for iter = 1:n_iterations
%     data_file = fullfile(output_dir, sprintf('uniform/results_iter_%d.mat', iter));
%     data = load(data_file);
% 
%     % Collect gridness and scale scores for the uniform data
%     uniform_data.square_gridness_all = [uniform_data.square_gridness_all, data.gridness_square_U(:, 1)];
%     uniform_data.square_gridness_exp_all = [uniform_data.square_gridness_exp_all, data.gridness_square_U(:, 2)];
%     uniform_data.square_scale_all = [uniform_data.square_scale_all, data.scale_square_U];
% end
% 
% % Create a new figure for the merged plots
% figure;
% hold on;
% 
% % Plot settings for gridness and scale
% edges_grid = linspace(-1.5, 1.5, 21); 
% edges_scale = linspace(0, 300, 21);
% colors = lines(length(boundary_effects));  % Generate distinct colors for each boundary effect
% 
% % Iterate through each boundary effect and plot the first three columns (Gridness, Expanded Gridness, Scale)
% for i = 1:length(boundary_effects)
%     b_effect = boundary_effects(i);
% 
%     % If boundary effect is 100, use uniform data
%     if b_effect == 100
%         % Gridness for uniform data
%         subplot(6, 4, ((i*4)-3));  % For gridness
%         histogram(mean(uniform_data.square_gridness_all, 2), edges_grid, 'LineWidth', 1.5);
%         ylim([0 100])
%         hold on;
% 
%         % Expanded Gridness for uniform data
%         subplot(6, 4, ((i*4)-2));  % For expanded gridness
%         histogram(mean(uniform_data.square_gridness_exp_all, 2), edges_grid, 'LineWidth', 1.5);
%         ylim([0 100])
%         hold on;
% 
%         % Scale for uniform data
%         subplot(6, 4, ((i*4)-1));  % For scale
%         histogram(mean(uniform_data.square_scale_all, 2), edges_scale, 'LineWidth', 1.5);
%         ylim([0 70])
%         hold on;
% 
%         % Right-most column for place field centers (Uniform)
%         % Call the function to get place field centers for the uniform condition
%         [xy_field, env, ~] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, b_effect); % Uniform version
%         subplot(6, 4, i*4);  % Adjust the layout for place field centers
%         plot(xy_field(:, 1), xy_field(:, 2), '.'); hold on;
%         title('Uniform Data');
%         xlim([0 dim_x]);
%         ylim([0 dim_y]);
% 
%     else
%         % Gridness
%         subplot(6, 4, ((i*4)-3));  % For gridness
%         histogram(mean(all_data(i).square_gridness_all, 2), edges_grid, 'LineWidth', 1.5);
%         ylim([0 100])
%         hold on;
% 
%         % Expanded Gridness
%         subplot(6, 4, ((i*4)-2));  % For expanded gridness
%         histogram(mean(all_data(i).square_gridness_exp_all, 2), edges_grid, 'LineWidth', 1.5);
%         ylim([0 100])
%         hold on;
% 
%         % Scale
%         subplot(6, 4, ((i*4)-1));  % For scale
%         histogram(mean(all_data(i).square_scale_all, 2), edges_scale, 'LineWidth', 1.5);
%         ylim([0 70])
%         hold on;
% 
%         % Right-most column for place field centers
%         [xy_field, env, ~] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, b_effect);
%         subplot(6, 4, i*4);  % Adjust the layout for place field centers
%         plot(xy_field(:, 1), xy_field(:, 2), '.'); hold on;
%         title(['Boundary Effect: ', num2str(b_effect)]);
%         xlim([0 dim_x]);
%         ylim([0 dim_y]);
%     end
% end
% 
% 
% % ** Final row: k-density for Gridness, Expanded Gridness, and Scale **
% for i = 1:length(boundary_effects)
%     b_effect = boundary_effects(i);
% 
%     % Gridness - k-density
%     subplot(6, 4, ((6*4)-3));  % Final row, column 1 for Gridness
%     [f_grid, xi_grid] = ksdensity(mean(all_data(i).square_gridness_all, 2));
%     plot(xi_grid, f_grid, 'Color', colors(i,:), 'LineWidth', 1.5);
%     hold on;
% 
%     % Expanded Gridness - k-density
%     subplot(6, 4, ((6*4)-2));  % Final row, column 2 for Expanded Gridness
%     [f_grid_exp, xi_grid_exp] = ksdensity(mean(all_data(i).square_gridness_exp_all, 2));
%     plot(xi_grid_exp, f_grid_exp, 'Color', colors(i,:), 'LineWidth', 1.5);
%     hold on;
% 
%     % Scale - k-density
%     subplot(6, 4, ((6*4)-1));  % Final row, column 3 for Scale
%     [f_scale, xi_scale] = ksdensity(mean(all_data(i).square_scale_all, 2));
%     plot(xi_scale, f_scale, 'Color', colors(i,:), 'LineWidth', 1.5);
%     hold on;
% end
% 
% % Set the titles for the first row of the figure
% subplot(6, 4, 1);  
% title('Gridness');
% subplot(6, 4, 2);  
% title('Expanded Gridness');
% subplot(6, 4, 3);  
% title('Scale');
% subplot(6, 4, 4);  
% title('Place Field Centers');
% 
% % Set the title for the final row (k-density)
% subplot(6, 4, ((6*4)-3));
% xlabel('Gridness Score');
% ylabel('Density');
% % title('K-Density: Gridness');
% 
% subplot(6, 4, ((6*4)-2));
% xlabel('Expanded Gridness Score');
% ylabel('Density');
% % title('K-Density: Expanded Gridness');
% 
% subplot(6, 4, ((6*4)-1));
% xlabel('Scale Score');
% ylabel('Density');
% % title('K-Density: Scale');
% 
% % % Add labels for each row (boundary effect value) on the left side
% % for i = 1:length(boundary_effects)
% %     effect = boundary_effects(i);
% %     if effect == 100
% %         ylabel('Uniform');
% %     else
% %         subplot(6, 4, ((i*4)-3));
% %         ylabel(num2str(effect));
% %     end
% % end
% 
% % Adjust layout
% hold off;
