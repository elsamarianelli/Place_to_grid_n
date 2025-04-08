%% Plotting Parameter search output from covariance_eigendecomp.m 
% Shows how different numbers of trajectory steps affect different gridness metrics.

% (1) SETTINGS 
baseName     = 'param_sweeps\Tanni_Covar_ED_runningSpeed_';    % Change depending on varied parameter
settings     = 5:5:35;                                    % Parameter values used in search
nIterations  = 5;                                        % Number of iterations for each parameter
metricNames  = {'expGrd_h', 'expGrd_s', 'stGrd_h', 'stGrd_s'}; % List of metrics to plot
threshold    = 0.8;                                      % Threshold for "gridness"
nPCs         = 200;                                       % Number of principal components
parameter    = 'Running Speeds';                          % for plotting

% (2) SETUP
colors = lines(numel(metricNames));       

% Preallocate scale arrays
meanProps = nan(numel(metricNames), numel(settings));
stdProps  = nan(numel(metricNames), numel(settings));
allScale_h = nan(length(settings), nIterations);  
allScale_s = nan(length(settings), nIterations);  

% (3) LOOP OVER METRICS
for m = 1:numel(metricNames)

    metricName = metricNames{m}; disp(metricName)
    allProps = nan(numel(settings), nIterations);

    for s = 1:numel(settings)

        stepVal = settings(s); disp(stepVal)

        for iter = 1:nIterations
            fname = fullfile([baseName num2str(stepVal)], ...
                             sprintf('iteration_%.1f', iter), ...
                             'metrics_and_maps.mat');
            if isfile(fname)
                data = load(fname);

                % Extract metric values across PCs
                metric_vals = cellfun(@(x) x.(metricName), data.GC_metrics, 'UniformOutput', false);
                metric_vals = cell2mat(metric_vals);

                % Proportion of PCs exceeding threshold
                allProps(s, iter) = sum(metric_vals > threshold) / size(metric_vals, 1); %% ---- change back stepVal to nPCs!!!!!!!!!

                % ALSO extract scale_h and scale_s
                scale_h_vals = cellfun(@(x) x.scale_h, data.GC_metrics, 'UniformOutput', false);
                scale_s_vals = cellfun(@(x) x.scale_s, data.GC_metrics, 'UniformOutput', false);
                scale_h_vals = cell2mat(scale_h_vals);
                scale_s_vals = cell2mat(scale_s_vals);
                
                % Save mean field scale (per iteration)
                allScale_h(s, iter) = mean(scale_h_vals, 'omitnan');
                allScale_s(s, iter) = mean(scale_s_vals, 'omitnan');

            else
                warning('Missing file: %s', fname);
            end
        end
    end

    % Compute mean and std over iterations
    meanProps(m, :) = mean(allProps, 2, 'omitnan');
    stdProps(m, :)  = std(allProps, 0, 2, 'omitnan');
    meanScale_h = mean(allScale_h, 2, 'omitnan');
    stdScale_h  = std(allScale_h, 0, 2, 'omitnan');
    meanScale_s = mean(allScale_s, 2, 'omitnan');
    stdScale_s  = std(allScale_s, 0, 2, 'omitnan');
    
end

% (4) PLOT: Gridness and Field Scale Comparison
figure('Color', 'w', 'Position', [100, 100, 1200, 550]);

% === SUBPLOT 1: Gridness Metrics ===
subplot(1,2,1); hold on;
metricNames = {'expanded hexagonal', 'expanded square', 'standard hexagonal', 'standard square'};
colors = lines(numel(metricNames));  
markers = {'o', 's', 'd', '^'};      

for m = 1:numel(metricNames)
    errorbar(settings, meanProps(m,:), stdProps(m,:), ...
        '-o', ...
        'Color', colors(m,:), ...
        'Marker', markers{m}, ...
        'MarkerFaceColor', colors(m,:), ...
        'LineWidth', 2, ...
        'MarkerSize', 8, ...
        'DisplayName', metricNames{m});
end

xlabel(parameter, 'FontSize', 16, 'FontWeight', 'bold');
ylabel(sprintf('Proportion of PCs > %.2f Gridness', threshold), 'FontSize', 16, 'FontWeight', 'bold');
title('Gridness Metrics', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off', ...
    'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');
grid off;

% === SUBPLOT 2: Scale Metrics ===
subplot(1,2,2); hold on;
scaleMetrics = {'scale_h', 'scale_s'};
meanScale = nan(numel(scaleMetrics), numel(settings));
stdScale  = nan(numel(scaleMetrics), numel(settings));

for m = 1:numel(scaleMetrics)
    metricName = scaleMetrics{m};
    allVals = nan(numel(settings), nIterations);

    for s = 1:numel(settings)

        stepVal = settings(s);
        disp(stepVal)
        for iter = 1:nIterations
            fname = fullfile([baseName num2str(stepVal)], ...
                             sprintf('iteration_%.1f', iter), ...
                             'metrics_and_maps.mat');
            if isfile(fname)
                data = load(fname);
                metric_vals = cellfun(@(x) x.(metricName), data.GC_metrics, 'UniformOutput', false);
                metric_vals = cell2mat(metric_vals);

                % Mean across PCs
                allVals(s, iter) = mean(metric_vals, 'omitnan');
            else
                warning('Missing file: %s', fname);
            end
        end
    end

    meanScale(m,:) = mean(allVals, 2, 'omitnan');
    stdScale(m,:)  = std(allVals, 0, 2, 'omitnan');
end

% Plotting scale data
scaleNames = {'Hexagonal Field Scale', 'Square Field Scale'};
for m = 1:2
    errorbar(settings, meanScale(m,:), stdScale(m,:), ...
        '-o', ...
        'Color', colors(m,:), ...
        'Marker', markers{m}, ...
        'MarkerFaceColor', colors(m,:), ...
        'LineWidth', 2, ...
        'MarkerSize', 8, ...
        'DisplayName', scaleNames{m});
end

xlabel(parameter, 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Mean Field Scale (pixels)', 'FontSize', 16, 'FontWeight', 'bold');
title('Field Scale Comparison', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off', ...
    'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');
grid off;

% %% for comparison between 2 different setting use this to change x axis ticks...
% tick_positions = sort(settings);  % ensures increasing order
% tick_labels = {'Random', 'Arrayed'};  % adjust label order to match
% tick_labels = {'5', '10', '15', '20', '25'};  % adjust label order to match
% 
% subplot(1,2,1);
% xticks(tick_positions);
% xticklabels(tick_labels);
% 
% subplot(1,2,2);
% xticks(tick_positions);
% xticklabels(tick_labels);

%% visualise gridness scores over threshold for each PC
baseName = 'param_sweeps\Tanni_Covar_ED_boundaryEffect_2.5';  % Change depending on varied parameter
nIterations = 5;
nPCs = 200;
threshold = .8;
gridness_visualisation_allPC(baseName, nIterations, nPCs, threshold)
%% List your 4 .fig files
% === Open each figure and store left/right axes handles ===
f1 = openfig('bound_effect.fig', 'invisible');
a1 = findall(f1, 'Type', 'axes'); a1 = flipud(a1);

f2 = openfig('number_place_cells.fig', 'invisible');
a2 = findall(f2, 'Type', 'axes'); a2 = flipud(a2);

f3 = openfig('traj_length.fig', 'invisible');
a3 = findall(f3, 'Type', 'axes'); a3 = flipud(a3);

f4 = openfig('random_vs_arrayed.fig', 'invisible');
a4 = findall(f4, 'Type', 'axes'); a4 = flipud(a4);

% === Create the combined figure ===
figure;
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% === Row 1 (with legends) ===
nexttile;
copyobj(a1(1).Children, gca);
legend(findobj(gca, '-regexp', 'Type', 'line|errorbar'), 'Location', 'best');
title('');
xlabel(a1(1).XLabel.String);
ylabel(a1(1).YLabel.String);
set(gca, 'FontSize', 12);

nexttile;
copyobj(a1(2).Children, gca);
legend(findobj(gca, '-regexp', 'Type', 'line|errorbar'), 'Location', 'best');
title('');
xlabel(a1(2).XLabel.String);
ylabel(a1(2).YLabel.String);
set(gca, 'FontSize', 12);

% === Row 2 ===
nexttile; copyobj(a2(1).Children, gca); title(''); xlabel(a2(1).XLabel.String); ylabel(a2(1).YLabel.String); set(gca, 'FontSize', 12);
nexttile; copyobj(a2(2).Children, gca); title(''); xlabel(a2(2).XLabel.String); ylabel(a2(2).YLabel.String); set(gca, 'FontSize', 12);

% === Row 3 ===
nexttile; copyobj(a3(1).Children, gca); title(''); xlabel(a3(1).XLabel.String); ylabel(a3(1).YLabel.String); set(gca, 'FontSize', 12);
nexttile; copyobj(a3(2).Children, gca); title(''); xlabel(a3(2).XLabel.String); ylabel(a3(2).YLabel.String); set(gca, 'FontSize', 12);

% === Row 4 ===
nexttile; copyobj(a4(1).Children, gca); title(''); xlabel(a4(1).XLabel.String); ylabel(a4(1).YLabel.String); set(gca, 'FontSize', 12);
nexttile; copyobj(a4(2).Children, gca); title(''); xlabel(a4(2).XLabel.String); ylabel(a4(2).YLabel.String); set(gca, 'FontSize', 12);

% === Close source figures ===
close(f1); close(f2); close(f3); close(f4);

% === Your custom column titles and y-axis overrides ===
colTitles = {'Gridness', 'Scale'};
customYLabels = {'% PCs > Threshold', 'Field Scale'};

% Column Titles (slightly up and to the left)
annotation(gcf, 'textbox', [0.22, 0.965, 0.2, 0.03], ...
    'String', colTitles{1}, 'EdgeColor', 'none', ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

annotation(gcf, 'textbox', [0.70, 0.965, 0.2, 0.03], ...
    'String', colTitles{2}, 'EdgeColor', 'none', ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Override Y labels for consistency
for i = 1:4
    nexttile((i - 1) * 2 + 1); ylabel(customYLabels{1}, 'FontSize', 12);
    nexttile((i - 1) * 2 + 2); ylabel(customYLabels{2}, 'FontSize', 12);
end

% === Custom x-ticks for Random vs Arrayed ===
tick_positions = sort([0, 4.5]);
tick_labels = {'Random', 'Arrayed'};

nexttile(7);
xticks(tick_positions);
xticklabels(tick_labels);

nexttile(8);
xticks(tick_positions);
xticklabels(tick_labels);

% === Custom x-axis label for Trajectory Length ===
nexttile(5); xlabel('Number of Steps');
nexttile(6); xlabel('Number of Steps');



% (6) LOOKING AT BINNED STATS 

-----