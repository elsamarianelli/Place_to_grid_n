function [maps] = reweighted_fr_maps(In, GC_idx, J, W)
% Reweight and sum spatial place cell activity to visualize grid cell input

    % get reweighted input to each grid cell first

    GC_maps = zeros(In.x+1, In.y+1, In.n_GCs); % initialise maps
    
    for GC = 1:In.n_GCs
        % get weights
        weights = J(GC_idx, :);
            
        % Initialize the summed map
        PC_summed = zeros(In.x + 1, In.y + 1);
    
        % Loop through each place cell
        for PC_idx = 1:In.n_PCs
    
            % get firing rate map for this place cell
            PC_map = In.PCs(:, :, PC_idx);
    
            % Reweight the map 
            weighted_map = PC_map * weights(PC_idx);
    
            % Add the weighted map to the summed map
            PC_summed = PC_summed + weighted_map;
        end

        GC_maps(:,:,GC) = PC_summed;

    end
    
    % reweight GC maps according to lateral weights 
    % get weights
    lat_weights = W(GC_idx, :);
        
    % Initialize the summed map
    GC_summed = zeros(In.x + 1, In.y + 1);

    % Loop through each grid cell
    for GC = 1:In.n_GCs

        % get firing rate map for this place cell
        GC_map = GC_maps(:, :, GC);

        % Reweight the map 
        weighted_map = GC_map * lat_weights(GC);

        % Add the weighted map to the summed map
        GC_summed = GC_summed + weighted_map;
    end
    
    % take the final map as the input from the lateral connections with
    % rho controlling how much it contributes, summed with the PC input
    % for that cell
    if In.Lateral 
        maps = (GC_summed.*In.rho) + ((1-In.rho)*(GC_maps(:,:,GC_idx)));
    else 
        maps = GC_maps(:,:, GC_idx);
    end
end
