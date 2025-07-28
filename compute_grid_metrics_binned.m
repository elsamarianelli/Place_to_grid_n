function metrics_all = compute_grid_metrics_binned(Grid_cells, binning_mode)
% COMPUTE_GRID_METRICS_BINNED - Extracts spatial metrics from grid cell maps using binning.

% Inputs:
%   Grid_cells   - Cell array of grid cell data across iterations.
%   binning_mode - How to bin environment ('three', 'trapezoid_lr')
%
% Output:
%   metrics_all  - Cell array [nIterations x nPCs] containing metrics per bin and PC

nIterations = numel(Grid_cells);
nPCs = numel(Grid_cells{1});
metrics_all = cell(nIterations, nPCs);

parfor iter = 1:nIterations
    disp(iter)
    for pc = 1:nPCs
        
            map = Grid_cells{iter}
            map = map{pc}.map
            % sac = xPearson(map);
            % [~, expGrd, ~, ~] = multiGridness(sac, 'hexagon', false, map);
            % 
            % if expGrd > .4 % arbitrary low threshold to speed up...

            try

                % Select binning function
                switch binning_mode
                    case 'three'
                        [~, binned_three, binned_nine, binned_var] = get_binned_metrics(map);
                        % if anything of the bins fail to have an elipse
                        % fitted exclude the whole pc - change to plotting
                        % because need to have ALL gridness 

                        % if any(arrayfun(@(x) any(isnan(x.eccent)), binned_three))
                        %     metrics_all{iter, pc} = [];
                        % else
                            metrics_all{iter, pc}.three       = binned_three;
                            metrics_all{iter, pc}.nine        = binned_nine;
                            % metrics_all{iter, pc}.variability = binned_var;
                        % end
                        
                    case 'trapezoid_lr' %% need to change
                        % map = rot90(map)
                        [~, binned_two] = get_binned_metrics_halves(map);
                        metrics_all{iter, pc}.trapezoid_lr   = binned_two;
                end
                fprintf("Metrics saved for iter %d, PC %d\n", iter, pc);
            catch
                metrics_all{iter, pc} = [];
            end
     
    end
end
end
