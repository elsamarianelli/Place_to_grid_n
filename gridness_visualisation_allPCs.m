%% cells over threhsold check 
function gridness_visualisation_allPC(basePath, nIerations, nPCs, threshold)
% Initialize preallocated arrays with NaNs
allData_standard = nan(nPCs * nIterations, 3);  % [PC_index, stGrd_h, stGrd_s]
allData_expanded = nan(nPCs * nIterations, 3);  % [PC_index, expGrd_h, expGrd_s]

% going through all iterations of  some covariance eigendecomp generated
% data and pulling out gridness metrics for each PC
for iter = 1:nIterations
    fname = fullfile(basePath, sprintf('iteration_%.1f', iter), 'metrics_and_maps.mat');
    if ~isfile(fname)
        warning('Missing: %s', fname);
        continue;
    end

    data = load(fname);
    GC = data.GC_metrics;

    for pc = 1:nPCs
        idx = (iter-1)*nPCs + pc;

        metrics = GC{pc};

        % allocate (check if empty)
        stH = NaN; if isfield(metrics, 'stGrd_h') && ~isempty(metrics.stGrd_h), stH = metrics.stGrd_h; end
        stS = NaN; if isfield(metrics, 'stGrd_s') && ~isempty(metrics.stGrd_s), stS = metrics.stGrd_s; end
        exH = NaN; if isfield(metrics, 'expGrd_h') && ~isempty(metrics.expGrd_h), exH = metrics.expGrd_h; end
        exS = NaN; if isfield(metrics, 'expGrd_s') && ~isempty(metrics.expGrd_s), exS = metrics.expGrd_s; end

        % sore values
        allData_standard(idx, :) = [idx, stH, stS];
        allData_expanded(idx, :) = [idx, exH, exS];
    end
end


%% PLOTTING

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

    % Define regions
    both = hex > threshold & square > threshold;
    neither = hex <= threshold & square <= threshold;
    hexOnly = hex > threshold & square <= threshold;
    squareOnly = hex <= threshold & square > threshold;

    % Plot regions with different colors
    scatter(hex(neither), square(neither), 40, [0.6 0.6 0.6], 'filled');
    scatter(hex(both), square(both), 40, [0.2 0.6 1], 'filled');
    scatter(hex(hexOnly), square(hexOnly), 40, [0.2 0.8 0.4], 'filled');
    scatter(hex(squareOnly), square(squareOnly), 40, [1 0.6 0.2], 'filled');

    % Threshold lines
    xline(threshold, '--k');
    yline(threshold, '--k');

    xlabel('Hexagonal Gridness', 'FontSize', 14);
    ylabel('Square Gridness', 'FontSize', 14);
    title(titleText, 'FontSize', 16);
    axis square;
    box off;

    legend({'Neither > t', 'Both > t', 'Hex only', 'Square only'}, ...
        'Location', 'bestoutside');
end

%% Plot Rate Maps for Cells That Are Grid-Like in Both Expanded Measures

% Get indices of qualifying PCs
allData = allData_expanded;
threshold = 1.8;
both = allData(:,2) > threshold & allData(:,3) > threshold;
gridLikePCs = allData(bothExpanded, 1);  % This gives global indices (1 to 1000)

fprintf('Found %d PCs with expanded hex and square gridness > %.2f\n', numel(gridLikePCs), threshold);

% Settings for plotting
nCols = 5;
nToShow = min(numel(gridLikePCs), 20);  % Limit number shown
nRows = ceil(nToShow / nCols);
figure('Color', 'w', 'Name', 'Grid-like PCs', 'Position', [100 100 1000 200 + 180*nRows]);

for i = 1:nToShow
    idx = gridLikePCs(i);
    
    pc = mod(idx-1, nPCs) + 1;  % Convert back to per-iteration PC index

    % Load the correct metrics file
    fname = fullfile(basePath, sprintf('iteration_%.1f', iter), 'metrics_and_maps.mat');
    if ~isfile(fname)
        warning('Missing file: %s', fname);
        continue;
    end
    data = load(fname);
    GC = data.GC_metrics;

    % Plot the rate map
    subplot(nRows, nCols, i);
    imagesc(GC{pc}.sac);
    axis image off;
    title(sprintf('PC %d ', pc), 'FontSize', 12);
end

ax = findall(gcf, 'Type', 'axes');
for i = 1:length(ax)
    axis(ax(i), 'tight');
end

end