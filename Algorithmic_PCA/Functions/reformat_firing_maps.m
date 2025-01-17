% Reformat Place cell for correct use
function [format_1, format_2] = reformat_firing_maps(fm, traj)

y_dim = size(fm{1}.fmap, 1);
x_dim = size(fm{1}.fmap, 2);
n_cells = size(fm, 2);

% format 1 - PC fr each x and y position, for all cells (dim x by dim y by
% ncells)
format_1 = zeros(x_dim, y_dim, n_cells);
for i = 1:n_cells
    format_1(:, :, i) = fm{i}.fmap'; % reflip x and y
end

% format 2 - each PC firing for each step in trajectoy (ncells by nsteps)
format_2 = zeros(n_cells, length(traj));
for ii = 1:length(traj)
    x = traj(ii, 1);
    y = traj(ii, 2);
    fr = format_1(x, y, :);
    format_2(:, ii) = fr;
end

end
