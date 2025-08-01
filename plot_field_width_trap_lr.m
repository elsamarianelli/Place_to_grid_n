function plot_field_width_trap_lr(Place_cells_all, In, subfolder, pc_size, include_centerline)
% PLOT_FIELD_WIDTH_TRAP_LR
% Compare place-field widths and nearest-neighbour distances between
% Left vs Right halves of a trapezoidal environment using a row-wise center line.
%
% INPUTS:
%   Place_cells_all  : [1 x nPCs] array of cells, each with .fmap
%   In               : struct with fields:
%                      - env.map (matrix with NaN or 0 outside environment, non-NaN inside)
%                      - env.L, env.dim_x, env.dim_y, env.dwmap, pf_width_cntrl, etc.
%   subfolder        : output folder path for saving figure(s)
%   pc_size          : true → boundary-scaled widths; false → fixed width
%   include_centerline (optional, default=false) :
%                      if true, assign centerline pixels to Left (≤ center)
%
% OUTPUTS:
%   Saves SVG figure: 'FieldStats_LR_Trap_StdOnly.svg'
%
% ASSUMPTIONS:
%   - fieldWidth() is available on the path (your existing helper).
%   - In.env.map marks valid in-bounds cells (preferably non-NaN inside).
%   - Place-cell "center" is the argmax of each .fmap (as in your original).

if nargin < 5 || isempty(include_centerline), include_centerline = false; end

% === [1] Get PC Centers ===
fmaps = cellfun(@(pc) pc.fmap, Place_cells_all, 'UniformOutput', false);
linear_idx = cellfun(@(f) find(f == max(f(:)), 1), fmaps);
[rows, cols] = cellfun(@(f, idx) ind2sub(size(f), idx), fmaps, num2cell(linear_idx));
xy_field = [cols', rows'];  % x, y
nPCs = size(xy_field, 1);

% === [2] Compute σx and σy (same logic as your square function) ===
sig_x = nan(nPCs, 1);
sig_y = nan(nPCs, 1);

% Build distance maps to walls
x_dists = In.env.L == 0;            % Identify wall locations
x_dists([2, In.env.dim_y], :) = 0;  % Set boundary rows to zero
x_dists = bwdist(x_dists);

y_dists = In.env.L == 0;
y_dists(:, [2, In.env.dim_x]) = 0;
y_dists = bwdist(y_dists);

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

% === [3] Build environment mask and row-wise center line ===
map = In.env.map;
[nRows, nCols] = size(map);

% Prefer non-NaN as "inside". If map has no NaNs, fall back to >0.
env_mask = ~isnan(map);
if ~any(env_mask(:))
    env_mask = map ~= 0;
end

left_mask  = false(nRows, nCols);
right_mask = false(nRows, nCols);

for r = 1:nRows
    cols_in = find(env_mask(r, :));
    if isempty(cols_in), continue; end
    c_min = cols_in(1);
    c_max = cols_in(end);
    c_center = floor( (c_min + c_max) / 2 );

    if include_centerline
        % Include centerline with Left
        left_mask(r, c_min : c_center) = env_mask(r, c_min : c_center);
        if c_center + 1 <= c_max
            right_mask(r, c_center+1 : c_max) = env_mask(r, c_center+1 : c_max);
        end
    else
        % Exclude centerline column from both halves
        if c_center - 1 >= c_min
            left_mask(r, c_min : c_center-1) = env_mask(r, c_min : c_center-1);
        end
        if c_center + 1 <= c_max
            right_mask(r, c_center+1 : c_max) = env_mask(r, c_center+1 : c_max);
        end
    end
end

% === [4] Assign each PC to Left or Right ===
region_labels = {'Left','Right'};
colors = {
    [0.75 0.85 0.95],  % light blue
    [0.30 0.50 0.75],  % medium blue
};

field_widths_x = {[], []};
field_widths_y = {[], []};
coords         = {[], []};

for i = 1:nPCs
    row = xy_field(i,2); col = xy_field(i,1);
    if row <= 0 || col <= 0 || row > nRows || col > nCols, continue; end

    if left_mask(row, col)
        ridx = 1;
    elseif right_mask(row, col)
        ridx = 2;
    else
        % not inside either half (centerline or outside env) → skip
        continue;
    end

    field_widths_x{ridx}(end+1) = sig_x(i);
    field_widths_y{ridx}(end+1) = sig_y(i);
    coords{ridx}(end+1, :)      = [col, row];
end

% === [5] Nearest-neighbour distances within each half ===
closest_dists = cell(1, 2);
for r = 1:2
    if size(coords{r},1) >= 2
        D = squareform(pdist(coords{r}));
        D(D == 0) = NaN;
        closest_dists{r} = min(D, [], 2);
    else
        closest_dists{r} = NaN;
    end
end

% === [6] Plot ===
fig = figure('Color', 'w', 'Position', [200 200 800 250]);

tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

titles  = {'X Field Width (σₓ)', 'Y Field Width (σᵧ)', 'Nearest PC Distance'};
ylabels = {'bins', '', ''};
data    = {field_widths_x, field_widths_y, closest_dists};

for p = 1:3
    nexttile; hold on;
    for r = 1:2
        vals = data{p}{r};
        if isempty(vals), continue; end
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        m = mean(vals);
        s = std(vals);
        bar(r, m, 'FaceColor', colors{r}, 'EdgeColor', 'k', 'LineWidth', 2);
        errorbar(r, m, s, 'k', 'LineWidth', 2);
    end
    title(titles{p}, 'FontSize', 20, 'FontWeight', 'bold');
    ylabel(ylabels{p}, 'FontSize', 18);
    xticks(1:2); xticklabels(region_labels);
    set(gca, 'FontSize', 16, 'Box', 'off', 'LineWidth', 2);
    ylim padded;
end

% sgtitle(sprintf('Place Field Statistics: Left vs Right (Trapezoid)%s', ...
        % ternary(include_centerline,' — centerline included','')), ...
        % 'FontSize', 22, 'FontWeight', 'bold');

if ~exist(subfolder, 'dir'), mkdir(subfolder); end
saveas(fig, fullfile(subfolder, 'FieldStats_LR_Trap_StdOnly.svg'));

end  % function

% --- tiny helper for nicer title text ---
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
