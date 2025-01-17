function [cells] = getGrid(cells, M, env, nn)
% getGrid - Computes grid cell representations from place cell data using
% eigenvectors of a SR matrix
%
% Inputs:
%   cells - A cell array where each element is a structure representing a place cell,
%           containing at least a field 'fmap' which is the firing rate map.
%   M     - A matrix representing spatial relationships, used for eigen decomposition.
%   env   - A structure representing the environment, containing at least a field 'L'
%           which is a binary matrix used for masking and smoothing operations.
%
% Outputs:
%   cells - The input cell array updated with a new field 'grid' for each cell, which
%           is the computed grid cell representation.

% Set the standard deviation for Gaussian smoothing
smth_sig = 3;

% eigen decmompostion of matrix
if strcmp(nn, 'off')
    % Perform eigen decomposition of matrix M
    [V, ~] = eig(M);
else   
    % non negative eigen decomp/PCA methods
    V = NNPCA2014(M, 200);
    % [W, H] = nnmf(M, 10);
    % [V] = nonNegativeEig(M, 200);
end

% Loop through each cell to compute its grid cell representation
for i = 1:length(cells)
    % Initialize the grid map as a zero matrix with the same size as the firing map
    grid_map = zeros(size(cells{1}.fmap));
    
    % Get the eigenvector corresponding to the current cell, use only the real part
    % this is not a completely symmetrical matrix, some Eig Vectors will be complex as a result, 
    % and so real(V) of all eig vectors are take to generate grid cells
    ev = real(V(:,i));
    
    % Sum place cell inputs weighted by the eigenvector components
    for j = 1:length(cells)
        grid_map = grid_map + ev(j) * cells{j}.fmap;
    end
    
    % Normalize the grid map to a maximum value of 1
    grid_map = grid_map ./ max(grid_map(:));
    
    % Apply Gaussian smoothing to the grid map, adjusting for environment mask
    grid_map = imgaussfilt(grid_map, smth_sig) ./ imgaussfilt(1 * (env.L == 2), smth_sig);
    
    % Set areas outside the environment (where L <= 1) to NaN
    grid_map(env.L <= 1) = NaN;
    
    % Threshold the grid map to ensure no negative values
    grid_map = max(grid_map - 0, 0);
    
    % Store the computed grid map in the cell structure
    cells{i}.grid = grid_map;
end
end
