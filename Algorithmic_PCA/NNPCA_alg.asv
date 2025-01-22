%% Algorithmic zero-meaned NNPCA implementation

%  code to implemetnt PCA which takes firing rates of PCs along whole
%  trajectory (zero meaned) as direct input (instead of eigdecomp on
%  SR matrix learnt using TDLR)

%  Using "Montanari A, Richard E. Non-negative principal component analysis  
%  Message Passing Algorithms and Sharp Asymptotics" approach. (FISTA would
%  be faster though).

%% Parameter considerations:

% Dordek uses "place cells in a 
% rectangular grid, such that a place cell is centered at each pixel of the 
% image (that is – number of place cells equals the number of image
% pixels)" - Here, need around 1000 Place cells or more, 
% with an arrayed distribution across the environment, instead of randomly 
% distributed place centres, can do random distribution but doesnt work as
% well.

% Additionaly, as in the dordek paper, periodic boudnary condition firing 
% fields are a requirement (aka the space they map onto toroidal space)
% (although trajectory here is not generated with periodic boundary 
% condition as in neural net version).

%% [1] Generating fr for pcs across whole trajectory (reformated)

% clear existing variables and add paths
clear all; close all; 
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\General_functions'
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Algorithmic_PCA\Functions'
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\boundary_warped_place2grid\Neural_net'
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\Algorithmic_PCA'

% Define environment dimensions and the number of cells
dim_x = 300; dim_y = 250;
n_polys = 1; polys = cell(n_polys, 1);
polys{1} = [0 0, dim_x-2, 0, dim_x-2 dim_y-2, 0 dim_y-2, 0 0] + 2; 
n_cells = 2500;
n_steps = 500000;
speed = 10;

% Generating Environment
env = GenerateEnv(polys, dim_x, dim_y, 'trapezoid');
    
% Generating Trajectory
% [a] Hasselmo trajectory - turns away from boundaries, velocity is altered
% traj = generate_trajectory(env, n_steps);

% [b] OR simpler random walk trajectory (ken's)
steps = (cumsum(randn(n_steps,2))*speed)+1;
x = ceil(mod(steps(:,1),dim_x)); 
y = ceil(mod(steps(:,2),dim_y)); 
traj = [x,y];

% Populating Place Cells
distribution_type = 'array'; 
[PlaceCellsUni, PlaceCellsTanni, env, xy_field_u, xy_field_t] = ...
    generate_place_cells(env, n_cells, dim_x, dim_y, 2, distribution_type);

% options for periodic boundaries with PC fr wrapping around
size_control = 7; % bigger for smaller PC firing fields
ToroidalPlaceCellMaps = get_PCs_toroidal(env, n_cells, xy_field_u, size_control, 0); 

% put into formats needed for PCA functions and plotting 
[format_1, format_2] = reformat_firing_maps(ToroidalPlaceCellMaps, traj);

% visualise
figure(1); subplot(1,2,1); plot(xy_field_u(:,1), xy_field_u(:,2),'.'); hold on;
           subplot(1,2,2); imagesc(format_1(:,:,100))

%% [2] Generate grids - with pca and nn pca (zero mean input)
%  NNPCA subtracts projection of each PC ouput from data before calculating next PC

NumberOfPC = 30;       % number of principle components to be generated (heirarchical structure)
zero_mean = 'spatial'; % option for how to zero-mean input data

% [a] standard PCA (zero mean input data)
PC_Standard = pca(format_2' - mean(format_2)');

% [b] NN PCA (hierarchical) 
r_saved = format_2;
L = length(format_2);

for z=1:NumberOfPC 
    
    % zero-mean input data
    if strcmp(zero_mean, 'spatial')
        meanInputMat = mean(r_saved,2);
        InputMatFixed = r_saved(:,1:L) - repmat(meanInputMat,1,L);
    elseif strcmp(zero_mean, 'temporal')
        meanInputMat = mean(r_saved,1);
        InputMatFixed = r_saved(:,1:L) - repmat(meanInputMat,size(format_1,1),1);
    end

    % Using Montanari A, Richard E. Non-negative principal component analysis 
    PC_NN(:,z) = NNPCA2014(InputMatFixed);

    % break if nans
    if sum( isnan(PC_NN(:,z))) > 1
        fprintf('NaNs!! breaking...\n');
        break;
    end

    % calculate the corresponding eigenvalue
    alpha1(:,z) = PC_NN(:,z)'*(InputMatFixed*InputMatFixed' )*PC_NN(:,z);

    % subtract the projection of the data on the 1st eigenvector 
    % from the data, and reuse the "new data" for the next PC calc.
    r_saved = InputMatFixed - PC_NN(:,z) * (PC_NN(:,z)'*InputMatFixed);

end
 
%% project the first eigenvector back to real world
% first take a column of v and multiply it by the place fields from the OG
% data to make the grid cell representation
 
PC = PC_NN; % switch from PC_Standard

figure; % GC firing
for ii = 1:NumberOfPC
    subplot(5,NumberOfPC/5,ii)
    map = comb_fields(format_1, PC(:,ii));
    imagesc(map)
    colormap jet
    drawnow
end

grids_fmap = cell(NumberOfPC, 1);
grids_sac = cell(NumberOfPC, 1);

figure; % sac
for iii = 1:NumberOfPC
    subplot(5,NumberOfPC/5,iii)
    map = comb_fields(format_1, PC(:,iii));
    grids_fmap{iii}.map = map;
    cross_corr = xPearson(map);
    grids_sac{iii}.sac = cross_corr;
    imagesc(cross_corr); 
    colormap jet
    drawnow
end

%% [3] analyse resulting Grids (not binned)
% addpath '/Users/elsamarianelli/GitHub/boundary_warped_place2grid/analysis-matlab-master/GridAnalysis/'
% addpath '/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/Misc'
% addpath '/Users/elsamarianelli/GitHub/boundary_warped_place2grid/analysis-matlab-master/Miscellaneous/'
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\analysis-matlab-master\GridAnalysis'
addpath 'C:\Users\Elsa Marianelli\Documents\GitHub\Place_to_grid_n\analysis-matlab-master\Miscellaneous'
figure;
grid_metrics = cell(n_cells, 1);
% Loop over j (vertical) and k (horizontal)
for GC_PC = 1:length(grids_fmap)
    try
        subplot(5,NumberOfPC/5,GC_PC)
        sac = xPearson(grids_fmap{GC_PC}.map);
        % Pass sac to autoCorrProps
        in.sac = sac;
        in.PLOT_ON = true; hold on;
        in.PLOT_Ellipse_ON = true;
        metrics = autoCorrProps(in); % should be one with EM comments, which has added ellipse function
        grid_metrics{GC_PC} = metrics;
    catch
        fprintf(strcat('error ', num2str(GC_PC)));
    end
end

%% [4] analyse resulting Grids binned --> grid scale too large to split bins for at least first 30 PCs
% set limits for dividing environment into 9 subregions
x_lims = [1, dim_x/3; dim_x/3, 2*(dim_x/3); 2*(dim_x/3), dim_x-1]; % X-axis limits for sub-regions
y_lims = [1, dim_y/3; dim_y/3, 2*(dim_y/3); 2*(dim_y/3), dim_y-1]; % Y-axis limits for sub-regions
grid_metrics_binned = cell(NumberOfPC, 3, 3);
figure;
j_l = 0:3:6;
% Loop over j (vertical) and k (horizontal)
for GC_PC = 28:length(grids_fmap)
    figure;
    tiledlayout(3, 3, 'TileSpacing','compact', 'Padding','tight')
    for j = 1:3
        for k = 1:3
            try
                nexttile
                % imagesc(grids_fmap{GC_PC}.map(y_lims(j, 1):y_lims(j, 2), x_lims(k, 1):x_lims(k, 2)))
                % disp((j_l(j) + (k)))
                sac = xPearson(grids_fmap{GC_PC}.map(y_lims(j, 1):y_lims(j, 2), x_lims(k, 1):x_lims(k, 2)));
                % Pass sac to autoCorrProps
                in.sac = sac;
                in.PLOT_ON = true; hold on;
                in.PLOT_Ellipse_ON = true;
                metrics = autoCorrProps(in);
                grid_metrics_binned{GC_PC, j, k} = metrics;
            catch
                fprintf('error...\n')
            end           
        end
    end
end

%% trying original 
% addpath 'Dordek-et-al.-Matlab-code\'Yedidyah''s code'\functions\PCAEvecCalc'
% PC_hold = struct;
% r_saved = format_2;
% NumberOfPC = 15;
% 
% [PC_hold] = PCAEvecCalc(PC_hold,r_saved,NumberOfPC);
% 
% map = comb_fields(format_1, PC_hold.PCANNEvec(:,1));
% figure; % GC firing
% for ii = 1:NumberOfPC
%     subplot(5,NumberOfPC/5,ii)
%     map = comb_fields(format_1, PC_hold.PCANNEvec(:,ii));
%     imagesc(map)
%     colormap jet
%     drawnow
% end