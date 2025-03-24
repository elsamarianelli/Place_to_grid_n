%% Reformat Place Cell Firing Maps for Analysis
function [format_1, format_2] = reformat_firing_maps(fm, traj)

% This function reformats firing maps of place cells into two different structures:
% 1. A 3D matrix representing spatial firing rates across x and y positions for all cells.
% 2. A 2D matrix representing the firing rates of each cell at each step along a given trajectory.
%
% Inputs:
%   fm   - A cell array where each element contains a structure with a firing map (fmap).
%          fmap is a 2D matrix representing the spatial firing rate of a place cell.
%   traj - A matrix of size [n_steps x 2], where each row represents a position (x, y) in the trajectory.
%
% Outputs:
%   format_1 - A 3D matrix (x_dim x y_dim x n_cells), where each slice contains the firing
%              rate map of a different place cell.
%   format_2 - A 2D matrix (n_cells x n_steps), where each column corresponds to the firing 
%              rates of all place cells at a specific step in the trajectory.

% Get dimensions of the firing map from the first cell
y_dim = size(fm{1}.fmap, 1);  % Number of y positions
x_dim = size(fm{1}.fmap, 2);  % Number of x positions
n_cells = size(fm, 2);        % Number of place cells

% Initialize format_1: a 3D matrix storing firing maps for all place cells
format_1 = zeros(x_dim, y_dim, n_cells);
for i = 1:n_cells
    format_1(:, :, i) = fm{i}.fmap'; % Transpose fmap to swap x and y dimensions
end

% Initialize format_2: a 2D matrix storing firing rates along the trajectory
format_2 = zeros(n_cells, length(traj));
for ii = 1:length(traj)
    x = traj(ii, 1);  % X position from trajectory
    y = traj(ii, 2);  % Y position from trajectory
    fr = format_1(x, y, :);  % Extract firing rate for all cells at (x, y)
    format_2(:, ii) = fr;     % Store firing rates for this time step
end

end

