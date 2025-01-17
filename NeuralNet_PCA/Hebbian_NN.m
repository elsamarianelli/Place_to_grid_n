%% Hebbian learning version
%  code to replicate non - negative weighted hebbian learning 
%  netowork from Dordek et al. 2015

%% Paths
% addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/Misc')
% addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/BasicFunctions')

%% [1] simulation Parameters 
In.x = 300;         % environment dimensions
In.y = 300;
n_polys = 1;        % enironemnt shape
polys = cell(n_polys, 1);
polys{1} = [0 0, (In.x - 2), 0, (In.x -2) (In.y - 2), 0 (In.y-2), 0 0] + 2; 
In.n_steps = 600000;% number of steps in trajectory
In.n_PCs = 15^2;    % number of place cells, need more if doing random distribution
In.n_GCs = 1;       % number of grid cells

% generating environment
env = GenerateEnv(polys, In.x, In.y, 'trapezoid');

% if using hasselmo trajectory, otherwise not needed...
% In.traj = generate_trajectory(env, In.n_steps);

% get place cell centres and fmaps which wrap around boundaries
distribution_type ='random';
[~, ~ , ~, PlaceCenters_u, PlaceCenters_t] = generate_place_cells(env, In.n_PCs , In.x, In.y, 2, distribution_type);
In.PC_centers = PlaceCenters_u'; 
xy_field = PlaceCenters_u;
boundary_compressed = 0;
size_fact = 4; % controls how much to divide fw by, larger size_fact = smaller fields
PlaceCells_uf = get_PCs_toroidal(env, In.n_PCs, xy_field, size_fact, boundary_compressed);

% reformat PC fmaps for later use in Network (into dimx x dimy x n_cells
% instead of stored in structures with each fmap)
PlaceCells = zeros(In.x+1, In.y+1, In.n_PCs);
for i = 1:In.n_PCs
    PlaceCells(:, :, i) = PlaceCells_uf{i}.fmap;
end
In.PCs = PlaceCells;

% visualiing some place cells
for i = 1:5:In.n_PCs
figure; 
imagesc(PlaceCells(:,:,i))
end

%% [2] network parameters 
In.epsilon = 1e6;     % real epsilon in function is 1/(delta+epsilon), bigger = slower
In.delta = 1;
In.maxWeight = 0.04; 
In.slow_learning = 1; % set wether you want to slow learning over time
                      % (delta is scaled by current step, so real epsilon
                      % reduces over time)

In.slope = 100;       % activation slope parameter to amplify difference in output
In.NonNeg = 1;        % Non negativity constraint
In.Visualise = 1;     % plotting
In.Lateral = 0;       % lateral connecitons in grid cell layer
In.rho = 0.4;         % Importance of lateral connections to GC firing 
In.Hasselmo = 0;      % trajectory to be used 
                      % (either hasselmo 1, or one which runs through boundaries onto other side 0)

% for grids to be acheived mean 0 parameter is also required                      
In.mean_zero = 'off'; % 'diff' or 'adapt'
In.adapt_rate = 0.5./3; % growth rate for inactivation of adaptation??
In.Sat = 30; % used when calculating gaussian fr in their version, although 
                % here the fr maps are calculate before hand then passed in

%% [2] 
[W, J] = Network(In); % W = lateral connections between GCs, J = PCs to GCs

%% visualise periodic boundaries trajectory 
n_steps = 5000;
x = 1; y = 1; direction = 0;
traj = zeros(2, n_steps);
figure;
for i = 1:n_steps
    [x, y, direction] = get_next_time(x, y, In, direction);
    traj(1,i) = x; 
    traj(2,i) = y; 
    plot(traj(1, :), traj(2, :), '.')
    drawnow;
end
