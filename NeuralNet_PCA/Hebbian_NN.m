%% Network version (Ojas)
%  code to replicate non - negative weighted hebbian learning 
%  netowork from Dordek et al. 2015

% note: problem with the algorithmic NN PCA was the number and size of
% grids was way too low, see if this helps here
%% Paths
% addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/Misc')
% addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/BasicFunctions')

%% [1] simulation Parameters 
In.x = 300;         % environment dimensions
In.y = 300;
n_polys = 1;        % enironemnt shape
polys = cell(n_polys, 1);
polys{1} = [0 0, (In.x - 2), 0, (In.x -2) (In.y - 2), 0 (In.y-2), 0 0] + 2; 
In.n_steps = 1000000;% number of steps in trajectory
In.n_PCs = 50^2;     % number of place cells, need more if doing random distribution
In.n_GCs = 3;        % number of grid cells

% generating environment
env = GenerateEnv(polys, In.x, In.y, 'trapezoid');

% if using hasselmo trajectory, otherwise not needed...
% In.traj = generate_trajectory(env, In.n_steps);

% OR simpler random walk trajectory (ken's) % also goes through bounds
speed = 5;
steps = (cumsum(randn(In.n_steps,2))*speed)+1;
x = ceil(mod(steps(:,1),In.x)); 
y = ceil(mod(steps(:,2),In.y)); 
In.traj_s = [x,y];

% Populating Place Cells
distribution_type = 'array'; 
[PlaceCellsUni, PlaceCellsTanni, env, xy_field_u, xy_field_t] = ...
    generate_place_cells(env, In.n_PCs, In.x, In.y, 2, distribution_type);

% options for periodic boundaries with PC fr wrapping around
size_control = 7; % bigger for smaller PC firing fields
ToroidalPlaceCellMaps = get_PCs_toroidal(env, In.n_PCs, xy_field_u, size_control, 0); 

% put into formats needed for network and plotting 
[format_1, format_2] = reformat_firing_maps(PlaceCellsUni, In.traj_s);
PlaceCells = format_1;

% visualising some place cells
% for i = 1:5:In.n_PCs
% figure; 
% imagesc(PlaceCells(:,:,i))
% end

%% [2] network parameters 
In.epsilon = 5e7;     % real epsilon in function is 1/(delta+epsilon), bigger = slower
In.delta = 1;
In.maxWeight = 0.04; 
In.slow_learning = 1; % set wether you want to slow learning over time
                      % (delta is scaled by current step, so real epsilon
                      % reduces over time)

In.slope = 100;        % activation slope parameter to amplify difference in output
In.NonNeg = 1;        % Non negativity constraint
In.Visualise = 1;     % plotting
In.Lateral = 1;       % lateral connecitons in grid cell layer
In.rho = 0.4;         % Importance of lateral connections to GC firing 
In.trajectory = 'rand_walk';    % trajectory to be used 
                      % (either hasselmo, p_bounds, rand_walk)

% for grids to be acheived mean 0 parameter is also required                      
In.mean_zero = 'diff';  % 'diff' or 'adapt'
In.adapt_rate = 0.5./3; % growth rate for inactivation of adaptation??
In.Sat = 30;            % used when calculating gaussian fr in their version, although 
                        % here the fr maps are calculate before hand then passed in

%% [2] 
[W, J] = Network(In); % W = lateral connections between GCs
                      % J = PCs to GCs

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

%% visualise other trajectories
figure;
for i = 1:n_steps
    plot(In.traj_s(i, 1), In.traj_s(i, 2), '.', 'Color', 'b')
    drawnow
    hold on
end


%%  testing these PCs using dordeks PCA Network - 
% %  is it a problem with the network or the parameters?
% Struct.simdur = 5000000; 
% Struct.Resolution =  75; 
% Struct.NmEC =  1; 
% Struct.N =  50 ; 
% Struct.N1 =  50^2; 
% Struct.epsilon =  1000000; 
% Struct.delta =  1; 
% Struct.maxWeight =  0.0400; 
% Struct.minx =  0 ;
% Struct.maxx =  300 ;
% Struct.miny =  0; 
% Struct.maxy =  300; 
% Struct.vel =  0.2867 ;
% Struct.angular =  0.2000 ;
% Struct.placeCellType =  'Gaussian'; 
% Struct.psiSat =  30; 
% Struct.meanZaro =  'diff' ;
% Struct.Arc =  'single';
% Struct.output =  'linear'; 
% Struct.NonNegativity =  1;
% Struct.saveInput =  0; 
% 
% % % generating environment
% % env = GenerateEnv(polys, In.x, In.y, 'trapezoid');
% % 
% % % Populating Place Cells
% % distribution_type = 'array'; 
% % [PlaceCellsUni, PlaceCellsTanni, env, xy_field_u, xy_field_t] = ...
% %     generate_place_cells(env, In.n_PCs, In.x, In.y, 2, distribution_type);
% % 
% % % options for periodic boundaries with PC fr wrapping around
% % size_control = 7; % bigger for smaller PC firing fields
% % ToroidalPlaceCellMaps = get_PCs_toroidal(env, In.n_PCs, xy_field_u, size_control, 0); 
% % 
% % % put into formats needed for network and plotting 
% % [format_1, format_2] = reformat_firing_maps(PlaceCellsUni, In.traj_s);
% % 
% % Struct.TemporalInput =  format_2; 
% Struct.placeCenters; 
% [Struct] = PCANetworkFunc(Struct, 1);
