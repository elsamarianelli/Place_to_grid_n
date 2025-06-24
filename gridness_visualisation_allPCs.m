function gridness_visualisation_allPCs(metrics_all, threshold)
% Visualizes hex vs square gridness across all PCs and iterations
% Inputs:
%   metrics_all - {nIterations x nPCs} struct array with grid metrics
%   threshold    - scalar threshold for gridness

[nIterations, nPCs] = size(metrics_all);
totalPCs = nIterations * nPCs;

% Preallocate
allData_standard = nan(totalPCs, 3);  % [PC_index, stGrd_h, stGrd_s]
allData_expanded = nan(totalPCs, 3);  % [PC_index, expGrd_h, expGrd_s]

% === Extract Metrics ===
for iter = 1:nIterations
    for pc = 1:nPCs
        idx = (iter-1)*nPCs + pc;
        metrics = metrics_all{iter, pc};

        stH = get_field(metrics, 'stGrd_h');
        stS = get_field(metrics, 'stGrd_s');
        exH = get_field(metrics, 'expGrd_h');
        exS = get_field(metrics, 'expGrd_s');

        allData_standard(idx, :) = [idx, stH, stS];
        allData_expanded(idx, :) = [idx, exH, exS];
    end
end

%% === PLOT ===
figure('Color', 'w', 'Position', [100 100 1200 500]);

for i = 1:2
    subplot(1,2,i); hold on;

    if i == 1
        data = allData_standard;
        titleText = 'Standard Gridness';
    else
        data = allData_expanded;
        titleText = 'Expanded Gridness';
    end

    hex = data(:,2);
    square = data(:,3);

    % Define logicals
    both = hex > threshold & square > threshold;
    neither = hex <= threshold & square <= threshold;
    hexOnly = hex > threshold & square <= threshold;
    squareOnly = hex <= threshold & square > threshold;

    scatter(hex(neither), square(neither), 40, [0.6 0.6 0.6], 'filled');
    scatter(hex(both), square(both), 40, [0.2 0.6 1], 'filled');
    scatter(hex(hexOnly), square(hexOnly), 40, [0.2 0.8 0.4], 'filled');
    scatter(hex(squareOnly), square(squareOnly), 40, [1 0.6 0.2], 'filled');

    xline(threshold, '--k');
    yline(threshold, '--k');

    xlabel('Hexagonal Gridness', 'FontSize', 14);
    ylabel('Square Gridness', 'FontSize', 14);
    title(titleText, 'FontSize', 16);
    axis square;
    box off;
    legend({'Neither > t', 'Both > t', 'Hex only', 'Square only'}, 'Location', 'bestoutside');
end
end


