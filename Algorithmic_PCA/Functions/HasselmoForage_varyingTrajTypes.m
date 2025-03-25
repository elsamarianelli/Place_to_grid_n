function [Position, Direction] = HasselmoForage_varyingTrajTypes(In)
%% Generates a foraging trajectory with adjustable mean speed.
%
%   INPUTS:
%       In.env           - Structure containing the environment (env.L, env.dim_x, env.dim_y)
%       In.n_steps       - Number of time steps
%       In.b             - Rayleigh scale parameter controlling average speed
%       In.wallHugBias   - tendency to go towards environment boundaries (0.01 -
%       0.05 range)     ----set to 0 to remove bias
%       In.wallAvoidDist - how far from al before you tilt away from it
%       (standard = 7)
%       In.wallHugDist   - how far from wall before you start tilting back
%       towards it (standard = 20)

%   OUTPUTS:
%       Position      - (n_steps x 2) matrix of x, y positions over time
%       Direction     - (n_steps x 1) vector of agent's heading directions

%% Initialize movement parameters
mu    = 0;        % Mean of turn angle
sigma = 0.2;      % Std dev of turning
dt    = 1 / 50;   % Time step (50 Hz update rate)
Dir   = rand * 2 * pi;  % Initial random direction

% Pre-generate movement randomness
RandomTurn = normrnd(mu, sigma, [1, In.n_steps]);    % Turning samples
RandomVel  = min(raylrnd(In.b, [1, In.n_steps]), 49);    % Velocity samples (Rayleigh-distributed)

%% Initialize output
Position  = zeros(In.n_steps, 2);
Direction = zeros(In.n_steps, 1);

% Start at a valid random location
start_idx = datasample(find(In.env.L == 2), 1);
[j, i]     = ind2sub([In.env.dim_y, In.env.dim_x], start_idx);
Position(1,:) = [i, j];
Direction(1)  = Dir;

%% Foraging Loop
for step = 2:In.n_steps
    [dWall, aWall] = minDistAngle(In.env, round(Position(step-1,:)), Dir);
    v = RandomVel(step);

    % Wall avoidance
    if dWall < In.wallAvoidDist && abs(aWall) < pi/2
        angle = sign(aWall) * (pi/2 - abs(aWall)) + RandomTurn(step);
        v = v - 0.5 * max([0, (v - 5)]);

    % Boundary hugging 
    elseif dWall > In.wallHugDist
        % if far from wall agent tends to turn back towards walls
        angle = -In.wallHugBias * aWall + RandomTurn(step);  % 0.1 is a gentle bias factor

    else
        angle = RandomTurn(step);

    end

    % Update position
    Position(step,:) = Position(step-1,:) + v * dt * [cos(-Dir), sin(-Dir)];

    % Boundary check
    if In.env.L(round(Position(step,2)), round(Position(step,1))) ~= 2
        Position(step,:) = Position(step-1,:);
    end

    % Update direction
    Dir = Dir + angle;
    Direction(step) = Dir;
end

%% Plot result
figure;
imagesc(~In.env.map'); colormap gray; axis xy; hold on;
plot(Position(:,2), Position(:,1), 'r-', 'LineWidth', 2);
title(sprintf('Foraging Path (Time = %.2f mins, b = %.1f)', In.n_steps / (60 * 50), In.b));
axis tight;
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 800]);

end
