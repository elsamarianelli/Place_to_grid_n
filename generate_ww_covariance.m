function [C, NeuronxEnvMat] = generate_ww_covariance(n_cells, dims, field_radius)

% Parameters
grid_size = prod(dims);            % total environment bins
NeuronxEnvMat = zeros(n_cells, grid_size);  % [n_cells × env_bins]

% Generate random place field centers
[x, y] = meshgrid(1:dims(2), 1:dims(1));
xy = [x(:), y(:)];

for i = 1:n_cells
    % Random center
    cx = randi(dims(2));
    cy = randi(dims(1));

    % Gaussian field
    d2 = (xy(:,1) - cx).^2 + (xy(:,2) - cy).^2;
    place_field = exp(-d2 / (2 * field_radius^2));
    
    % Normalize field and store
    place_field = place_field / max(place_field);
    NeuronxEnvMat(i, :) = place_field';
end

% Simulate "activity over time" via pseudo random trajectory
n_timepoints = 10000;
indices = randi(grid_size, [1, n_timepoints]);
activity = NeuronxEnvMat(:, indices);  % [n_cells × time]

% Zero-mean each neuron (as Wang & Wang do)
activity = activity - mean(activity, 2);

% Covariance of neuron responses
C = cov(activity');  % [n_cells × n_cells]

end
