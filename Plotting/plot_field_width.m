function plot_field_width(Place_cells_all, In, subfolder, pc_size)
% PLOT_FIELD_WIDTH
% Plots standard deviations (σx, σy) and nearest-neighbor distances of place cells
% grouped by environment regions (corners, edges, center), using field width logic only.
%
% INPUTS:
%   Place_cells_all : [1 x nPCs] array of cells, each with .fmap
%   In              : struct containing env.dwmap, env.L, pf_width_cntrl, etc.
%   subfolder       : folder path to save figure
%   pc_size         : true for boundary-scaled widths, false for fixed width

% === [1] Get PC Centers ===
fmaps = cellfun(@(pc) pc.fmap, Place_cells_all, 'UniformOutput', false);
linear_idx = cellfun(@(f) find(f == max(f(:)), 1), fmaps);
[rows, cols] = cellfun(@(f, idx) ind2sub(size(f), idx), fmaps, num2cell(linear_idx));
xy_field = [cols', rows'];  % x, y

% === [2] Compute σx and σy ===
nPCs = size(xy_field, 1);
sig_x = nan(nPCs, 1);
sig_y = nan(nPCs, 1);
% recalculate x and y dists to use
x_dists = In.env.L == 0; % Identify wall locations
x_dists([2, In.env.dim_y], :) = 0; % Set boundary rows to zero
x_dists = bwdist(x_dists); % Compute distance to the nearest wall in x direction

y_dists = In.env.L == 0; % Identify wall locations
y_dists(:, [2, In.env.dim_x]) = 0; % Set boundary columns to zero
y_dists = bwdist(y_dists); % Compute distance to the nearest wall in y direction
In.env.x_dists = x_dists; 
In.env.y_dists = y_dists; 

if pc_size
    for i = 1:nPCs
        cx = xy_field(i, 1); 
        cy = xy_field(i, 2);
        if isfinite(cx) && isfinite(cy) && ...
           cx > 0 && cy > 0 && ...
           cy <= size(In.env.x_dists,1) && cx <= size(In.env.x_dists,2)

            dx = In.env.x_dists(cy, cx);
            dy = In.env.y_dists(cy, cx);
            sig_x(i) = fieldWidth(dx) / 4;
            sig_y(i) = fieldWidth(dy) / 4;
        end
    end
else
    av_bound_dist = nanmean(In.env.dwmap(:));
    fixed_width = fieldWidth(av_bound_dist) / In.pf_width_cntrl;
    sig_x(:) = fixed_width;
    sig_y(:) = fixed_width;
end

% === [3] Assign bins to 3x3 regions ===
region_labels = {'Corners', 'Edges Ver.', 'Edges Hor.', 'Center'};
region_ids = { [1 3 7 9], [4 6], [2 8], 5 };
colors = {
    [0.75 0.85 0.95],  % light blue
    [0.45 0.65 0.85],  % medium light blue
    [0.30 0.50 0.75],  % medium blue
    [0.15 0.30 0.55]   % dark blue
};
[nRows, nCols] = size(In.env.map);
row_thirds = floor(nRows / 3);
col_thirds = floor(nCols / 3);
bin_map = zeros(nRows, nCols);
bin_id = reshape(1:9, [3 3])';
for i = 1:3
    for j = 1:3
        r_idx = (1 + (i-1)*row_thirds) : min(i*row_thirds, nRows);
        c_idx = (1 + (j-1)*col_thirds) : min(j*col_thirds, nCols);
        bin_map(r_idx, c_idx) = bin_id(i,j);
    end
end

% === [4] Group σx, σy, and coords by region ===
field_widths_x = cell(1, 4); field_widths_y = cell(1, 4); coords = cell(1, 4);
for i = 1:nPCs
    row = xy_field(i,2); col = xy_field(i,1);
    if row <= 0 || col <= 0 || row > nRows || col > nCols, continue; end
    bin = bin_map(row, col);
    for r = 1:4
        if ismember(bin, region_ids{r})
            field_widths_x{r}(end+1) = sig_x(i);
            field_widths_y{r}(end+1) = sig_y(i);
            coords{r}(end+1, :) = [col, row];
            break
        end
    end
end

% === [5] Nearest Neighbor Distances ===
closest_dists = cell(1, 4);
for r = 1:4
    if size(coords{r},1) >= 2
        D = squareform(pdist(coords{r}));
        D(D == 0) = NaN;
        minD = min(D, [], 2);
        closest_dists{r} = minD;
    else
        closest_dists{r} = NaN;
    end
end

% === [6] PLOT ===
fig = figure('Color', 'w', 'Position', [200 200 800 350]);  % Smaller figure
tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');  % Tighter layout

titles = {'X Field Width (σₓ)', 'Y Field Width (σᵧ)', 'Nearest PC Distance'};
ylabels = {'bins', '', ''};
data = {field_widths_x, field_widths_y, closest_dists};

for p = 1:3
    nexttile; hold on;
    for r = 1:4
        vals = data{p}{r};
        if isempty(vals), continue; end
        vals = vals(~isnan(vals));
        m = mean(vals); s = std(vals);
        bar(r, m, 'FaceColor', colors{r}, 'EdgeColor', 'k', 'LineWidth', 2);  % Thicker bars
        errorbar(r, m, s, 'k', 'LineWidth', 2);  % Thicker error bars
    end
    title(titles{p}, 'FontSize', 20, 'FontWeight', 'bold');
    ylabel(ylabels{p}, 'FontSize', 18);
    xticks(1:4); xticklabels(region_labels);
    set(gca, 'FontSize', 16, 'Box', 'off', 'LineWidth', 2);  % Thicker axes
    ylim padded;
end

sgtitle('Place Field Statistics by Region', 'FontSize', 22, 'FontWeight', 'bold');
saveas(fig, fullfile(subfolder, 'FieldStats_ByRegion_StdOnly.svg'));

end  