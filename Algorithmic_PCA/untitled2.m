%% NN algorithmic PCA whole 

%% [1] Generating fr for pcs across whole trajectory (reformated)
 
% Define environment dimensions and the number of cells
dim_x = 351; dim_y = 252;
n_polys = 1; polys = cell(n_polys, 1);
polys{1} = [0 0, 349, 0, 349 250, 0 250, 0 0] + 2; 
n_cells = 300;
n_steps = 360000;

% Generating Environment
env = GenerateEnv(polys, dim_x, dim_y, 'trapezoid');
    
% Generating Trajectory
traj = generate_trajectory(env, n_steps);

% Populating Place Cells
distribution_type = 'random';
[PlaceCellsUni, PlaceCellsTanni, env, xy_field_u, xy_field_t] = generate_place_cells(env, n_cells, dim_x, dim_y, 2, distribution_type);

% reformat for later use (want format_2 in this case)
[format_1, format_2] = reformat_firing_maps(PlaceCellsUni, traj);

%% [2] Generate grids
addpath '/Users/elsamarianelli/Documents/GitHub/place2grid/matlab' % to dans

method = 'nnpca';

NumberOfPC = 25; % number of principle components
elec = format_2; % firing for each step
pfields = format_1; % firing rate maps for whole environment

switch method
    case 'pca'
        u = pca(elec' - mean(elec)');
    case 'nnpca'
        r_saved = elec;
        fprintf('Calculating Non-Negative PCA...\n')
        L = length(elec);
        for z=1:NumberOfPC %per PC
            % subtract the mean from every feature
            meanInputMat = mean(r_saved,2);
            InputMatFixed = r_saved(:,1:L) - repmat(meanInputMat,1,L);

            % send the fixed data to the Monatrani algorithm, see 
            % (Montanari A, Richard E. Non-negative principal component 
            % analysis: Message passing algorithms and sharp asymptotics. 
            % arXiv preprint arXiv:14064775. 2014.

            u(:,z) = NNPCA2014( InputMatFixed);
            % check if there's Nan's in the eigenvectors!
            if sum( isnan(u(:,z))) >1
                fprintf('NaNs!! breaking...\n');
                break;
            end

            % caclulate the corresponding eigenvalue
            alpha1(:,z) = u(:,z)'*(InputMatFixed*InputMatFixed' )*u(:,z);

            % subtract the projection of the data on the 1st eigenvector 
            % from the data, and reuse the "new data" for the next PC calc.
            r_saved = InputMatFixed - u(:,z) * (u(:,z)'*InputMatFixed);

        end

    otherwise
        error('not a valid method, choose pca or nnpca')
end

%% try and project the first eigenvector (column of v) back to real world
% first take a column of v and multiply it by the place fields from the OG
% data to make the grid cell representation

figure(2); % GC firing
for ii = 1:25
    subplot(5,5,ii)
    map = p2g_combine_placefields(pfields, u(:,ii));
    imagesc(map)
    colormap jet
    drawnow
end

figure(3); % cross corr
for iii = 1:25
    subplot(5,5,iii)
    map = p2g_combine_placefields(pfields, u(:,iii));
    map = xPearson(map);
    imagesc(map)
    colormap jet
    drawnow
end

