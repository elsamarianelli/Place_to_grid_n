% Simplified MATLAB code to simulate grid cell firing

% Parameters
envSize = 100;          % Size of the environment (100x100)
numCells = 20;          % Number of grid cells
gridSpacing = 30;       % Grid spacing (approx. in cm)
gridOrientation = pi/6; % Orientation of the grid (30 degrees)
numSteps = 40000;         % Number of steps in random trajectory

% generating hasselmo trajectory
polys = cell(1, 1);
polys{1} = [0 0, envSize-2, 0, envSize-2 envSize-2, 0 envSize-2, 0 0] + 2; 
env = GenerateEnv(polys, envSize, envSize, 'trapezoid');
trajectory = generate_trajectory(env, numSteps);

% Generate grid cell properties
theta = gridOrientation + [0, 2*pi/3, -2*pi/3]; % Axes of grid
gridOffsets = envSize * rand(numCells, 2); % Random offsets for each cell

% Generate grid cell firing fields
[X, Y] = meshgrid(1:envSize, 1:envSize); % Environment grid
gridFiringFields = zeros(envSize, envSize, numCells); % Store firing fields
for cellIdx = 1:numCells
    % Compute grid firing for the cell
    phaseOffsets = [cos(theta); sin(theta)]' * (gridOffsets(cellIdx, :) ./ gridSpacing)';
    gridFiringFields(:, :, cellIdx) = ...
        cos(2*pi/gridSpacing * (cos(theta(1)) * X + sin(theta(1)) * Y + phaseOffsets(1))) + ...
        cos(2*pi/gridSpacing * (cos(theta(2)) * X + sin(theta(2)) * Y + phaseOffsets(2))) + ...
        cos(2*pi/gridSpacing * (cos(theta(3)) * X + sin(theta(3)) * Y + phaseOffsets(3)));
    gridFiringFields(:, :, cellIdx) = max(0, gridFiringFields(:, :, cellIdx)); % Rectify firing
end

% Simulate grid cell activations along the trajectory
activations = zeros(numSteps, numCells);
for stepIdx = 1:numSteps
    x = round(trajectory(stepIdx, 1));
    y = round(trajectory(stepIdx, 2));
    x = max(min(x, envSize), 1); % Bound x within the environment
    y = max(min(y, envSize), 1); % Bound y within the environment
    for cellIdx = 1:numCells
        activations(stepIdx, cellIdx) = gridFiringFields(y, x, cellIdx);
    end
end

% Plot results
figure;
subplot(2, 2, 1);
plot(trajectory(:, 1), trajectory(:, 2), 'r', 'LineWidth', 1.5);
title('Random Trajectory');
xlim([0 envSize]); ylim([0 envSize]);
xlabel('X'); ylabel('Y');

subplot(2, 2, 2);
imagesc(gridFiringFields(:, :, 1));
title('Example Grid Cell Firing Field');
colorbar; axis equal tight;

subplot(2, 2, [3, 4]);
imagesc(activations');
title('Grid Cell Activations Along Trajectory');
xlabel('Step'); ylabel('Grid Cell');
colorbar;

% Plot all grid cell firing fields
figure;
numRows = ceil(sqrt(numCells)); % Number of rows and columns in the grid
for cellIdx = 1:numCells
    subplot(numRows, numRows, cellIdx);
    imagesc(gridFiringFields(:, :, cellIdx));
    title(['Cell ' num2str(cellIdx)]);
    axis equal tight;
    colorbar;
    xlabel('X'); ylabel('Y');
end
sgtitle('Grid Cell Firing Fields'); % Title for the entire figure
