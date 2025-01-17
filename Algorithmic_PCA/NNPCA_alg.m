%% NN algorithmic PCA whole 

%% [1] Generating fr for pcs across whole trajectory (reformated)
 
% Define environment dimensions and the number of cells
dim_x = 300; dim_y = 300;
n_polys = 1; polys = cell(n_polys, 1);
polys{1} = [0 0, dim_x-2, 0, dim_x-2 dim_y-2, 0 dim_y-2, 0 0] + 2; 
n_cells = 1000; % seems like a lot of cells needed to get good grids
n_steps = 360000;

% Generating Environment
env = GenerateEnv(polys, dim_x, dim_y, 'trapezoid');
    
% Generating Trajectory
traj = generate_trajectory(env, n_steps);

% Populating Place Cells
distribution_type = 'array';
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

grids_fmap = cell(NumberOfPC, 1);
grids_sac = cell(NumberOfPC, 1);

figure(3); % cross corr
for iii = 1:NumberOfPC
    subplot(5,5,iii)
    map = p2g_combine_placefields(pfields, u(:,iii));
    grids_fmap{iii}.map = map;
    cross_corr = xPearson(map);
    grids_sac{iii}.sac = cross_corr;
    imagesc(cross_corr); 
    colormap jet
    drawnow
end

%% [3] analyse resulting Grids (not binned)
addpath '/Users/elsamarianelli/GitHub/boundary_warped_place2grid/analysis-matlab-master/GridAnalysis/'
addpath '/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/Misc'
addpath '/Users/elsamarianelli/GitHub/boundary_warped_place2grid/analysis-matlab-master/Miscellaneous/'

figure(4);
grid_metrics = cell(n_cells, 1);
% Loop over j (vertical) and k (horizontal)
for GC_PC = 1:length(grids_fmap)
    try
        subplot(5,5,GC_PC)
        sac = xPearson(grids_fmap{GC_PC}.map);
        % Pass sac to autoCorrProps
        in.sac = sac;
        in.PLOT_ON = true; hold on;
        in.PLOT_Ellipse_ON = true;
        metrics = autoCorrProps(in); % should be one with EM comments, which had added ellipse function
        grid_metrics{GC_PC} = metrics;
    catch
        fprintf(strcat('error ', num2str(GC_PC)));
    end
end

%% [4] analyse resulting Grids binned
% set limits for dividing environment into 9 subregions
x_lims = [0, dim_x/3; dim_x/3, 2*(dim_x/3); 2*(dim_x/3), dim_x-0]; % X-axis limits for sub-regions
y_lims = [0, dim_y/3; dim_y/3, 2*(dim_y/3); 2*(dim_y/3), dim_y-0]; % Y-axis limits for sub-regions
grid_metrics_binned = cell(NumberOfPC, 3, 3);
% Loop over j (vertical) and k (horizontal)
for GC_PC = 1:length(grids_fmap)
    for j = 1:3
        for k = 1:3
            try
                sac = xPearson(grids_fmap{GC_PC}.grid(y_lims(j, 1):y_lims(j, 2), x_lims(k, 1):x_lims(k, 2)));
                % Pass sac to autoCorrProps
                in.sac = sac;
                in.PLOT_ON = true; hold on;
                in.PLOT_Ellipse_ON = true;
                metrics = autoCorrProps(in); % should be one with EM comments, which had added ellipse function
                grid_metrics_binned{GC_PC, j, k} = metrics;
            catch
                fprintf('error...\n')
            end
        end
    end
end

