function [M, R] = trainModel(cells, M, R, traj, time_lag)
% trainModel - Trains a model using temporal difference learning to update 
% the successor representation matrix (M) and the reward vector (R).
%
% Inputs:
%   cells     - A structure or matrix representing spatial cells or locations.
%   M         - The successor representation matrix to be updated.
%   R         - The reward vector associated with states.
%   traj      - A matrix containing trajectory data with each row representing a state (x, y).
%   time_lag  - The time step difference used to predict future states.
%
% Outputs:
%   M         - Updated successor representation matrix.
%   R         - Updated reward vector.

% Initialize a progress bar for monitoring the training process
h = waitbar(0, 'Training model');

% Set up a figure to display the successor matrix during training
figure
hold on
imagesc(M) % Display the initial successor matrix
title('Successor Matrix')

% Calculate the total number of update steps
total_steps = length(traj) - time_lag;

% Initialize learning rate parameters
start_alpha = 1e-4; % Starting learning rate
end_alpha = 1e-4;   % Ending learning rate

% Linearly interpolate learning rates over the total steps
alpha = (1:total_steps) * (end_alpha - start_alpha) / total_steps + start_alpha;

% Initialize reward accumulator (not used in this function but might be used in other parts)
r = 0;

% Main loop to update the successor representation matrix (M) and reward vector (R)
for t = 1:total_steps
    % Update the progress bar
    waitbar(t / total_steps)
    
    % Store the current state of the successor matrix
    old_M = M;
    
    % Extract the current state (x, y) from the trajectory
    x = traj(t, 1);
    y = traj(t, 2);
    
    % Extract the next state (n_x, n_y) from the trajectory based on the time lag
    n_x = traj(t + time_lag, 1);
    n_y = traj(t + time_lag, 2);
    
    % Get the activity rates or features for the current state
    phi = GetRates(ceil(x), ceil(y), cells);
    
    % Update the successor matrix (M) and reward vector (R) using a temporal difference method
    [M, R] = SR_Update(r, phi, GetRates(ceil(n_x), ceil(n_y), cells), M, R, alpha(t));
    
    % Periodically update the displayed successor matrix to visualize training progress
    if mod(t, 800) == 0
        imagesc(M)
    end
end

% Add a color bar to the figure for better visualization
colorbar

% Close the progress bar after training is complete
close(h)

end

