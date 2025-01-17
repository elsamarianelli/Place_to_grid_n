function [W, J] = Network(In)
% code to emualte nn PCA using a hebbian network with non-negative weights
% Simplified from Dordek version, using Tanni warped PCs and Hasselmo trajectory instead

% addpath('/Users/elsamarianelli/Documents/GitHub/boundary_warped_place2grid/Neural_net')

% Takes In structure with...
% - x and y dimension sizes for environmnet 
% - n_steps in trajectory 
% - n_PCs and n_GCs in each layer 
% - PCs : PC firing maps (dim x by dim y by n_PCs)
% - epsilon and delta : learning rate parameters 
    % learning rate = 1/(delta+epsilon)
% - maxWeight : weight limit for feedforward connections (J)
% - slowlearning : delta scaled by step number in trajectory if 1
% - slope : scales output (h) for bigger difference
% - NonNeg : sets all PC to GC weights below 0 to 0 for each iteration if 1
% - Visaulise : plotting if 1
% - Lateral : lateral connecitons between GCs if 1, with ...
% - rho : paramater which adjusts how much the lateral inputs influence
%   activity over PC inputs
% - Hasselmo : hasselmo trajectory which avoids boundaries (1) or
%   trajectory which goes through walls and reappears on opposite side
%   (space is toroidal) (0)
% - mean_zero : 'off' if normal PC firing at location used to update weights
%               'adapt' (Kropff and Treves 2008) method, subtracting weighted
% sum of past firing rates
%               'diff' which uses derivatives as rate

% Returns...
% - J = PCs to GCs weights
% - W = GC lateral connections weights

% [0] Initialise weights and trajectory starting point 

% PC to GC layer weights ( normalised ) 
J = rand(In.n_GCs, In.n_PCs); 
normFactor = repmat(sum(J,2),1,In.n_PCs);
J = In.maxWeight*J./normFactor;

% Lateral connection weights in GC layer 
W = zeros(In.n_GCs);

% initial input from place to grid cells ( for lateral architecture )
h = zeros(1,In.n_GCs);

% initialise for zero-mean condition
h_avg = zeros(1,In.n_GCs);

%initilize derivatives
r_prev = 0;
r_prev_prev = 0;
r_gaussian =0;

% store trajectory for plot later 
traj_log = zeros(2,In.n_steps);

% initial traj parameters
x = randi(In.x, 1);
y = randi(In.y, 1);
direction = 0;

% simulation Loop, update weights in network over traectory 
numcell = 1; % GC for later plotting
fig1 = figure;

% simulate Ojas rule over the course of a trajectory 
for step = 1:In.n_steps    

    disp(step)

    % update trajectory log
    traj_log(1, step) = x;
    traj_log(2, step) = y;

    % [2] get array of PC activity in that location
    PC_activity = squeeze(In.PCs(x, y, :))';
    
    % zero mean constraint 
    r_prev_prev = r_prev;
    r_prev = r_gaussian;
    r_gaussian = PC_activity.*In.Sat;

    if strcmp(In.mean_zero,'diff') % mean zero using derivatives
        delta_r = r_gaussian - r_prev;
        delta_r2 = r_prev - r_prev_prev;
        delta2_r = (delta_r - delta_r2);
        r = delta2_r;
    elseif strcmp(In.mean_zero,'adapt') % subtract mean for each iteration
        h_adapt = h-h_avg;
        h_avg = h_avg + In.adapt_rate*(h-h_avg);
        r = r_gaussian;
    elseif strcmp(In.mean_zero, 'off')
        r = PC_activity;
    end

    % [3] update learning rule, slows as learning progresses to stabilise
    % towards end
    if In.slow_learning 
        epsilon = (1/(step*In.delta+In.epsilon));
    else 
        epsilon = 1/(In.delta+In.epsilon);
    end

    % [4] calculate PC layer output based on current weights
    % output is linear with high activation funciton slope     
    if In.Lateral
        h = (1-In.rho)*(J*r')' + In.rho*(W*(h+eps)')';
    else
        h = (J*r')' ;
    end
    
    if strcmp(In.mean_zero,'adapt') %% not too sure about their adapt 
        h = h_adapt;
    end

    % increase slope activation function slope
    h = In.slope*h;
    
    % [5a] update PC to GC weights (Ojas rule)
    normalisation_term = eye(In.n_GCs,In.n_GCs).*(h'*h)*J;
    deltaJ =  h'*r - normalisation_term;
    J = J + epsilon*deltaJ;
    
    % [b] non negativity constraint
    if In.NonNeg 
        J(J<0) = 0;
    end

    % [6] update GC layer lateral weights
    if In.Lateral
    % update GC lateral weights if In.Lateral = 1
        W = W - 0.001*epsilon*(h'*h); % slower learning rate
        % removing recurrent inputs
        W = (ones(size(W)) - eye(size(W))).*W; 
    end 
   
    % [7] PLotting
     if step == 1
        % Initialize plots
        subplot1 = subplot(4, 2, 1:2);
        plotW = plot(subplot1, W(numcell, :), 'b');
        title(subplot1, 'Lateral weights');
        xlabel('Grid cells');
    
        subplot2 = subplot(4, 2, 3:4);
        plotJ = plot(subplot2, J(numcell, :), 'b');
        ylim([0 .03]);
        title(subplot2, 'Inter-layer weights');
        xlabel('Place cells');
    
        subplot3 = subplot(4, 2, 5);
        map = reweighted_fr_maps(In, numcell, J, W);
        img1 = imagesc(map, 'Parent', subplot3);
        axis(subplot3, 'equal');  % Ensure square scaling
        title(subplot3, 'Reweighted Firing Rate Map');
    
        subplot4 = subplot(4, 2, 6);
        cross_map = xPearson(map);
        img2 = imagesc(cross_map, 'Parent', subplot4);
        axis(subplot4, 'equal');  % Ensure square scaling
        title(subplot4, 'Autocorrelation Map');

        % subplot5 = subplot(4, 2, 7);
        % traj_plot = plot(traj_log(1,:), traj_log(2,:),'.');
        % ylim([0 In.y]); xlim([0 In.x]);
        % axis(subplot5, 'equal');  % Ensure square scaling
        % title(subplot5, 'Trajectory');

    end
    
    % Update data periodically
    if mod(step, In.n_steps / 100) == 0
        % Update plot data
        set(plotW, 'YData', W(numcell, :));
        set(plotJ, 'YData', J(numcell, :));
    
        % Update images
        map = reweighted_fr_maps(In, numcell, J, W);
        set(img1, 'CData', map);
        set(img2, 'CData', xPearson(map));
    
        % Fix axis limits and scaling
        axis(subplot3, [1 size(map, 2) 1 size(map, 1)]);
        axis(subplot4, [1 size(cross_map, 2) 1 size(cross_map, 1)]);
    
        % Update visualization
        drawnow limitrate;
    end
    % [8] Take step 
    if In.Hasselmo % HASSELMO THAT AVOIDS BOUNDARIES DOESNT WORK HERE BECASUE NEED PERIODIC BOUNDARY CONDITIONS
        x = In.traj(step+1, 1);
        y = In.traj(step+1, 1);
    else
        [x, y, direction] = get_next_time(x,y,In,direction);
        x = x+1;
        y = y+1;  % cant be 0s
    end

end

end




