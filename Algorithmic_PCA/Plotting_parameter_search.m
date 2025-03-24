%% Plotting Parameter search output from covariance_eigendecomp.m 
%  Shows how different numbers of trajectory steps affect different gridness metrics.

% -------- 1) SETTINGS --------
baseName     = 'Tanni_Covar_ED_steps_';
stepsList    = 6000:6000:30000;
nIterations  = 5;

% List of metrics to plot
metricNames = {'stGrd_s',  'stGrd_h'};
colors = lines(numel(metricNames));  % use distinguishable colors

% -------- 2) INIT STORAGE --------
meanVals = zeros(length(stepsList), numel(metricNames));
stdVals  = zeros(length(stepsList), numel(metricNames));

% -------- 3) MAIN LOOP --------
for m = 1:length(metricNames)
    metricName = metricNames{m};
    allMetrics = cell(length(stepsList), 1);

    for s = 1:length(stepsList)
        stepVal = stepsList(s);
        stepMetrics = zeros(nIterations, 1);

        for iter = 1:nIterations
            fname = fullfile([baseName num2str(stepVal)], sprintf('iteration_%.1f', iter), 'metrics_and_maps.mat');
            if isfile(fname)
                data = load(fname);
                % Use arrayfun if metrics is a struct array
                gridness_vals = arrayfun(@(x) x.(metricName), data.metrics);
                stepMetrics(iter) = mean(gridness_vals); % or max(gridness_vals)
            else
                warning('Missing file: %s', fname);
                stepMetrics(iter) = NaN;
            end
        end

        allMetrics{s} = stepMetrics;
    end

    % Compute mean and std for this metric across iterations
    meanVals(:,m) = cellfun(@(x) mean(x, 'omitnan'), allMetrics);
    stdVals(:,m)  = cellfun(@(x) std(x, 'omitnan'), allMetrics);
end

% -------- 4) PLOT --------
figure; hold on;
for m = 1:length(metricNames)
    errorbar(stepsList, meanVals(:,m), stdVals(:,m), '-o', ...
        'LineWidth', 1.5, 'DisplayName', strrep(metricNames{m}, '_', '\_'), ...
        'Color', colors(m,:));
end
xlabel('Number of Steps');
ylabel('Gridness (mean ± SD)');
title('Gridness Metrics Across Trajectory Lengths');
legend('Location', 'best');
grid on;
