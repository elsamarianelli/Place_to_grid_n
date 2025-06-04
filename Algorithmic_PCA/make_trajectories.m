%% Generate and save trajectories for later use

% Define environment dimensions
In.dim_x = 351.*3;
In.dim_y = 252.*3;
n_polys = 1;
In.polys = cell(n_polys, 1);
trap_add = 0;
In.polys{1} = [0 trap_add, (In.dim_x-2) 0, (In.dim_x-2) (In.dim_y-2), 0 ((In.dim_y-2)-trap_add), 0 0] + 2;  
In.shape = 'trapezoid'; 

In.n_steps = 360000.*9; % Number of steps in the trajectory
n_iterations = 5; % Number of trajectories to generate

output_dir = 'premade_traj_and_env_9xEnv'; % Output folder for storing data

% Ensure the output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for n = 1:n_iterations
    
    % [1] Generate Environment
    fprintf('Generating Environment for trajectory %d/%d\n', n, n_iterations);
    In = GenerateEnv(In); % Modify if needed
    In.env.dim_x = In.dim_x; In.env.dim_y = In.dim_y; 

    % [2] Generate Random Trajectory
    fprintf('Generating Trajectory %d/%d\n', n, n_iterations);
    traj = generate_trajectory(In.env, In.n_steps); 
    
    % [3] Store extra info
    details = struct();
    details.thigmotaxis = false;  % Change to `true` if wall-following behavior is enabled
    
    % [4] Save Environment & Trajectory
    filename = fullfile(output_dir, sprintf('trajectory_%02d.mat', n));
    save(filename, 'env', 'traj', 'details');
    
    fprintf('Saved trajectory %d/%d to %s\n', n, n_iterations, filename);
end

fprintf('All trajectories successfully generated and saved!\n');

