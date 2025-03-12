function [Position, Direction] = HasselmoForage(env, n_steps)
%HASSELMOFORAGE Generates a foraging trajectory using Hasselmo's model.
%
%   This function simulates an agent foraging in an environment using a 
%   velocity-based movement model inspired by Hasselmo et al. The agent's 
%   movement is influenced by random velocity fluctuations and wall avoidance.
%
%   INPUTS:
%       env        - Structure containing the environment, including:
%                    - env.L: A matrix representing the environment layout
%                    - env.dim_x, env.dim_y: Dimensions of the environment
%       n_steps  - Number of time steps in the foraging trajectory
%
%   OUTPUTS:
%       Position   - (n_steps x 2) matrix of x, y coordinates over time
%       Direction  - (n_steps x 1) vector of agent's heading directions
%
%   MOVEMENT MODEL:
%   - Velocity follows a Rayleigh distribution (biased towards lower speeds).
%   - The agent turns randomly with a Gaussian distribution.
%   - When near a wall, the agent slows down and turns away from it.
%   - The movement updates at a rate of 50 Hz (dt = 1/50 sec per step).
%

%% **Initialize Parameters for Velocity and Turning**
mu = 0;            % Mean of random turn angle (Gaussian)
sigma = 0.2;       % Standard deviation of random turn angle
b = 16;            % Scale parameter for Rayleigh-distributed velocity
v = 20;            % Initial velocity
Dir = rand * 2 * pi; % Initial random heading direction
dt = 1 / 50;       % Time step (50 Hz update rate)

% Generate random velocity and turning samples
RandomTurn = normrnd(mu, sigma, [1, n_steps]);   % Gaussian turn angles
RandomVel = min(raylrnd(b, [1, n_steps]), 49);   % Rayleigh-distributed speeds

%% **Initialize Position and Direction Storage**
Position = zeros(n_steps, 2);  % Stores (x, y) positions
Direction = zeros(n_steps, 1); % Stores direction angles

% Select a random valid starting position
start = datasample(find(env.L == 2), 1); % Choose a position where env.L == 2 (valid locations)
[j, i] = ind2sub([env.dim_y, env.dim_x], start); % Convert index to (y, x)
Position(1, :) = [i, j]; % Store initial position
Direction(1) = Dir; % Store initial direction

%% **Foraging Loop**
for step = 2:n_steps
    % Get distance and angle to the nearest wall
    [dWall, aWall] = minDistAngle(env, round(Position(step-1, :)), Dir);
    
    % Update speed from the precomputed random velocity sample
    v = RandomVel(step);
    
    % **Wall Avoidance Behavior**
    if (dWall < 7 && abs(aWall) < pi/2) % If near a wall and facing it
        angle = sign(aWall) * (pi/2 - abs(aWall)) + RandomTurn(step); % Turn away
        v = v - 0.5 * max([0, (v - 5)]); % Reduce speed when near a wall
    else
        angle = RandomTurn(step); % Otherwise, turn randomly
    end
    
    % **Move Agent**
    Position(step, :) = Position(step-1, :) + v * dt * [cos(-Dir), sin(-Dir)];
    
    % **Boundary Check: Keep Position Inside Valid Region**
    if (env.L(round(Position(step, 2)), round(Position(step, 1))) ~= 2)
        Position(step, :) = Position(step-1, :); % Prevent moving into walls
    end
    
    % **Update Direction**
    Dir = Dir + angle;
    Direction(step) = Dir;
end

%% **Plot the Generated Trajectory**
figure;
imagesc(~env.map'); % Display environment map (flipped for visibility)
colormap gray;
hold on;
plot(Position(:, 2), Position(:, 1), 'LineWidth', 2); % Plot trajectory
title(sprintf('Foraging Path (Time = %.2f mins)', n_steps / (60 * 50)));

% Set plot appearance
axis tight;
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 16);
set(gcf, 'Position', [100, 100, 800, 800]);

end
