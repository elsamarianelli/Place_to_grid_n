%% "Plot a representative distribution of place field locations" 

baseName     = 'param_sweeps\Tanni_Covar_ED_boundaryEffect_';   
settings     = .5:1:4.5;
for s = 1:numel(settings)

    stepVal = settings(s); disp(stepVal)

    for iter = 1:nIterations
        fname = fullfile([baseName num2str(stepVal)], ...
            sprintf('iteration_%.1f', iter), ...
            'orig_place_cells.mat');
        data = load(fname);
        
    end
end

%% note for monday - use dwmap in in, and peak of pc , to get each distance to closest bound across 5 iterations across boundary controls
% and plot distribution.