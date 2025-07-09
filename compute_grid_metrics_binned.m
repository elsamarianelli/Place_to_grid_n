function metrics_all = compute_grid_metrics_binned(Grid_cells, binning_mode)
% COMPUTE_GRID_METRICS_BINNED - Extracts spatial metrics from grid cell maps using binning.

% Inputs:
%   Grid_cells   - Cell array of grid cell data across iterations.
%   shape        - Shape used for SAC/gridness ('hexagon', etc.)
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
        map = Grid_cells{iter}{pc}.map;
        smoothed_map = smoothdata(map);
        peaks = FastPeakFind(smoothed_map);  % format: [x1 y1 x2 y2 ...] → divide by 2 for count

        % Set minimum peaks based on binning scheme
        switch binning_mode
            case 'trapezoid_lr'
                num_req_peaks = 10;
            case 'three'
                num_req_peaks = 30;
            otherwise
                error('Unsupported binning_mode: %s', binning_mode);
        end

        if numel(peaks)/2 > num_req_peaks
            try
                % Select binning function
                switch binning_mode
                    case 'three'
                        [~, binned_three, binned_nine, binned_var] = get_binned_metrics(map);
                        if any(arrayfun(@(x) any(isnan(x.eccent)), binned_three))
                            metrics_all{iter, pc} = [];
                        else
                            metrics_all{iter, pc}.three       = binned_three;
                            metrics_all{iter, pc}.nine        = binned_nine;
                            metrics_all{iter, pc}.variability = binned_var;
                        end

                    case 'trapezoid_lr'
                        map = rot90(map)
                        [~, binned_two, binned_var2] = get_binned_metrics_trapezoid(map);
                        if any(arrayfun(@(x) any(isnan(x.scale_h)), binned_two))
                            metrics_all{iter, pc} = [];
                        else
                            metrics_all{iter, pc}.trapezoid_lr           = binned_two;
                            metrics_all{iter, pc}.variability_trapezoid  = binned_var2;
                        end
                end
                fprintf("Metrics saved for iter %d, PC %d\n", iter, pc);
            catch
                metrics_all{iter, pc} = [];
            end
        end
    end
end
end
