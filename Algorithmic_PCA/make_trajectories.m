%% Generate and save 20 trajectories for later use

% Define environment dimensions
dim_x = 351;
dim_y = 252;
n_polys = 1;
polys = cell(n_polys, 1);
polys{1} = [0 0, 349 0, 349 250, 0 250, 0 0] + 2; % Rectangular boundary

n_steps = 360000; % Number of steps in the trajectory
n_iterations = 20; % Number of trajectories to generate

output_dir = 'premade_traj_and_env'; % Output folder for storing data

% Ensure the output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for n = 1:n_iterations
    
    % [1] Generate Environment
    fprintf('Generating Environment for trajectory %d/%d\n', n, n_iterations);
    env = GenerateEnv(polys, dim_x, dim_y, 'trapezoid'); % Modify if needed
    
    % [2] Generate Random Trajectory
    fprintf('Generating Trajectory %d/%d\n', n, n_iterations);
    traj = generate_trajectory(env, n_steps); 
    
    % [3] Store extra info
    details = struct();
    details.thigmotaxis = false;  % Change to `true` if wall-following behavior is enabled
    
    % [4] Save Environment & Trajectory
    filename = fullfile(output_dir, sprintf('trajectory_%02d.mat', n));
    save(filename, 'env', 'traj', 'details');
    
    fprintf('Saved trajectory %d/%d to %s\n', n, n_iterations, filename);
end

fprintf('All trajectories successfully generated and saved!\n');

