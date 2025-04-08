%% Generate_premade_trajectories.m
%  Generate and save multiple trajectories with varying speed and wall-hug bias
%  Set parameters as used in de cothi code

% Generate environment
In.dim_x = 351;
In.dim_y = 252;
In.n_polys = 1;
polys = cell(In.n_polys, 1);
In.polys{1} = [0 0, 349 0, 349 250, 0 250, 0 0] + 2; % rectangular, change for warped trapezoid
env = GenerateEnv(In.polys, In.dim_x, In.dim_y, 'trapezoid');
save(fullfile(base_dir, 'environment.mat'), 'env');  % Save once

% Generate premade trajectories with wall bias or speed variation
base_dir = fullfile(pwd, 'premade_traj_and_env');
if ~exist(base_dir, 'dir'); mkdir(base_dir); end

n_steps = 150000;
n_iterations = 5;

% [1] Vary wall hug bias 

fixed_b = 16;
hug_bias_values = 0:0.05:0.2;

for hug_bias = hug_bias_values
    disp(hug_bias)
    for iter = 1:n_iterations
        disp(iter)
        In.env = env;
        In.n_steps = n_steps;
        In.b = fixed_b;
        In.wallHugBias = hug_bias;
        In.wallAvoidDist = 7;
        In.wallHugDist = 20;

        [Position, Direction] = HasselmoForage_varyingTrajTypes(In);

        hug_str = strrep(sprintf('%.2f', hug_bias), '.', 'p');
        fname = sprintf('trajectory_b%02d_hug%s_iter%02d.mat', fixed_b, hug_str, iter);
        fpath = fullfile(base_dir, fname);

        traj.traj = Position;
        traj.dir = Direction;
        traj.b = fixed_b;
        traj.wallHugBias = hug_bias;
        traj.iter = iter;
        save(fpath, 'traj');

        fprintf('Saved (hug bias sweep): %s\n', fname);
    end
end

% [2] Vary speed 
speed_values = 30:5:45;
fixed_hug_bias = 0;

for b = speed_values
    for iter = 1:n_iterations
        In.env = env;
        In.n_steps = n_steps;
        In.b = b;
        In.wallHugBias = fixed_hug_bias;
        In.wallAvoidDist = 7;
        In.wallHugDist = 20;

        [Position, Direction] = HasselmoForage_varyingTrajTypes(In);

        hug_str = strrep(sprintf('%.2f', fixed_hug_bias), '.', 'p');
        fname = sprintf('trajectory_b%02d_hug%s_iter%02d.mat', b, hug_str, iter);
        fpath = fullfile(base_dir, fname);

        traj.traj = Position;
        traj.dir = Direction;
        traj.b = b;
        traj.wallHugBias = fixed_hug_bias;
        traj.iter = iter;
        save(fpath, 'traj');

        fprintf('Saved (speed sweep): %s\n', fname);
    end
end


%% Plotting
% Plot example trajectories: 2 rows (hug bias and speed), 5 columns each
figure('Color', 'w', 'Position', [100 100 1200 500]);
n_cols = 5;
base_dir = 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\premade_traj_and_env';
% --- First row: vary wall hug bias ---
hug_bias_values = 0:0.05:0.2;
fixed_b = 16;
for i = 1:n_cols
    hug_bias = hug_bias_values(i);
    hug_str = strrep(sprintf('%.2f', hug_bias), '.', 'p');
    fname = sprintf('trajectory_b%02d_hug%s_iter01.mat', fixed_b, hug_str);  % show iter 1
    fpath = fullfile(base_dir, fname);
    if isfile(fpath)
        data = load(fpath);
        subplot(2, n_cols, i); hold on;
        plot(data.traj.traj(:,1), data.traj.traj(:,2), 'r', 'LineWidth', 1.5);
        title(sprintf('Wall Bias = %.2f', hug_bias), 'FontSize', 12);
    end
end

% --- Second row: vary speed ---
speed_values = 5:5:25;
fixed_hug_bias = 0;
for i = 1:n_cols
    b = speed_values(i);
    hug_str = strrep(sprintf('%.2f', fixed_hug_bias), '.', 'p');
    fname = sprintf('trajectory_b%02d_hug%s_iter01.mat', b, hug_str);
    fpath = fullfile(base_dir, fname);
    if isfile(fpath)
        data = load(fpath);
        subplot(2, n_cols, i + n_cols); hold on;
        plot(data.traj.traj(:,1), data.traj.traj(:,2), '.', 'MarkerSize',3);
        title(sprintf('Speed b = %d', b), 'FontSize', 12);
    end
end

sgtitle('Example Foraging Trajectories: Wall Bias (Top) vs Speed (Bottom)', 'FontSize', 16, 'FontWeight', 'bold');