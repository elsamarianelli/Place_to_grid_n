%% vectorised version -COMB_FIELDS Combines multiple spatial maps using a weighted sum.
% This function takes a set of 2D spatial maps and combines them into a 
% single map using a weighted sum, then normalizes the result.
% Inputs:
%   maps  - A 3D matrix of size [x_dim, y_dim, n_maps], where each slice 
%           maps(:,:,ii) represents an individual spatial field.
%   combo - A vector of size [n_maps, 1] containing weights for each map.
%
% Output:
%   map   - A 2D matrix [x_dim, y_dim] representing the weighted combination 
%           of all input maps, normalized by the number of maps.
function map = comb_fields(maps, combo)
    % Validate dimensions
    if size(maps, 3) ~= numel(combo)
        error('Number of maps does not match length of combo vector');
    end

    % Vectorized weighted sum
    combo = reshape(combo, 1, 1, []);  % Ensure broadcasting shape
    map = sum(maps .* combo, 3);      % Weighted sum across third dimension

    % Normalize central region (excluding 2-cell boundary)
    inner = map(3:end-2, 3:end-2);
    inner = (inner - min(inner(:))) / max(inner(:) - min(inner(:)) + eps);  % +eps to avoid /0
    map(3:end-2, 3:end-2) = inner;
end




%% old version

% function [map] = comb_fields(maps, combo)
% 
% % Initialize the output map as a zero matrix with the same x and y dimensions
% map = zeros(size(maps, [1 2]));
% 
% % Loop through each map in the third dimension and apply the weight from combo
% for ii = 1:size(maps,3)
%     map = map + combo(ii) * maps(:,:,ii); % Weighted sum of maps
% end
% 
% inner = map(3:end-2, 3:end-2);  
% map(3:end-2, 3:end-2) = (inner - min(inner(:))) / (max(inner(:)) - min(inner(:)));
% 
% end
