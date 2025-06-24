function plot_gridness_by_region(metrics_all, threshold, subfolder)
% Visualize hex vs square gridness across Corners, Edges, Center
% metrics_all must contain .three.{expGrd_h, expGrd_s}

[nIterations, nPCs] = size(metrics_all);
region_labels = {'Corners', 'Edges', 'Center'};
nRegions = numel(region_labels);

% Preallocate
allData = nan(nPCs * nIterations, 5, nRegions);  % [idx, hex, square, iter, pc]

for iter = 1:nIterations
    for pc = 1:nPCs
        idx = (iter - 1) * nPCs + pc;
        metrics = metrics_all{iter, pc};
        if isempty(metrics) || ~isfield(metrics, 'three'), continue; end

        tbl_hex = struct2table(metrics.three);
        if ~all(ismember({'expGrd_h_el', 'expGrd_s_el'}, tbl_hex.Properties.VariableNames))
            continue;
        end

        for r = 1:nRegions
            hex = tbl_hex.expGrd_h(r);
            sq  = tbl_hex.expGrd_s(r);
            allData(idx,:,r) = [idx, hex, sq, iter, pc];
        end
    end
end

% === Plotting ===
colors = struct( ...
    'neither', [0.6 0.6 0.6], ...
    'both',    [0.1 0.6 0.6], ...
    'hexOnly', [0.4 0.7 1], ...
    'sqOnly',  [0.5 0.9 0.5]);

fig = figure('Color', 'w', 'Position', [100 100 800 250]);
t = tiledlayout(1, nRegions, 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:nRegions
    nexttile(r); hold on;
    data = squeeze(allData(:,:,r));
    hex = data(:,2); square = data(:,3);

    % Categories
    both     = hex > threshold & square > threshold;
    neither  = hex <= threshold & square <= threshold;
    hexOnly  = hex > threshold & square <= threshold;
    sqOnly   = hex <= threshold & square > threshold;

    scatter(hex(neither), square(neither), 50, colors.neither, 'filled');
    scatter(hex(both),    square(both),    50, colors.both,    'filled');
    scatter(hex(hexOnly), square(hexOnly), 50, colors.hexOnly, 'filled');
    scatter(hex(sqOnly),  square(sqOnly),  50, colors.sqOnly,  'filled');

    xline(threshold, '--k', 'LineWidth', 1.5);
    yline(threshold, '--k', 'LineWidth', 1.5);

    xlabel('Hex Gridness', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Square Gridness', 'FontSize', 16, 'FontWeight', 'bold');
    title(region_labels{r}, 'FontSize', 18, 'FontWeight', 'bold');
    axis square;
    set(gca, 'FontSize', 16, 'LineWidth', 2, 'Box', 'off');
end

% Add legend manually outside of tiles
legend_entries = {'Neither > t', 'Both > t', 'Hex only', 'Square only'};
lgd = legend(legend_entries, ...
    'Position', [0.92 0.3 0.05 0.4], ...
    'FontSize', 14, 'Box', 'off');
set(gcf, 'Renderer', 'painters');  % For vector output

% Save
if nargin > 2 && isfolder(subfolder)
    saveas(fig, fullfile(subfolder, 'Gridness_By_Region.svg'));
end
end


