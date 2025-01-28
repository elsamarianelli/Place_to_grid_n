function [gridnessData, overThreshold] = get_binned_average_gridness(n_iter, results, type)

    n_bins = 21;
    edges_grid = linspace(-1.5, 3, n_bins); 
    
    % Fields to process (you will vary between square and hex)
    gridness_types = {'gridness_square', 'gridness_hex'};
    
    % Initialize a structure to store the results
    gridnessData = struct();
    overThreshold = struct();

    % Loop through each gridness type (square and hex)
    for g = 1:length(gridness_types)
        gridness_type_base = gridness_types{g};
        field_name = [gridness_type_base, type];
        
        % Initialize arrays to store binned data for each gridness type
        binnedData_exp = zeros(n_iter, n_bins-1);
        binnedData_std = zeros(n_iter, n_bins-1);
        
        % Loop through each iteration and extract the gridness measure
        for i = 1:n_iter
            % Extract the 200x2 matrix from the ith cell for both exp and std (column 1 and 2)
            gridness_iter = results.(field_name){i};
            
            % Extract column 1 for standard and column 2 for expanded 
            gridness_measure_std = gridness_iter(:,1); 
            gridness_measure_exp = gridness_iter(:,2);  
            
            %store high gridness index 
            index_over_thresh = find(gridness_measure_exp > 1);
            overThreshold.(field_name){i} = index_over_thresh;
            
            % binning for standard data
            [counts_std, ~] = histcounts(gridness_measure_std, edges_grid);
            binnedData_std(i, :) = counts_std;
            
            % binning for expanded data
            [counts_exp, ~] = histcounts(gridness_measure_exp, edges_grid);
            binnedData_exp(i, :) = counts_exp;
        end
        
        % Compute mean and std for the standard and experimental data
        gridnessData.([gridness_type_base, '_std_mean']) = mean(binnedData_std, 1);
        gridnessData.([gridness_type_base, '_std_std']) = std(binnedData_std, 1);
        gridnessData.([gridness_type_base, '_exp_mean']) = mean(binnedData_exp, 1);
        gridnessData.([gridness_type_base, '_exp_std']) = std(binnedData_exp, 1);

    end
end
