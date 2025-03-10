function stats = get_figure_and_stats_for_tiled_layout(trial, t, c, in)

    nexttile
    % [1] Rate Map / Polar Plot 
    % spikeData(3) cut==5 as example with firing
    % cell 5, tetrode 3
    sp = trial.spikeData(t);
    in.spikePosInd = sp.pos_sample(sp.cut==c);
    in.posXY = trial.posData.xy;
    in.MAP_TYPE = 'rate'; % 'dir' for polar plot
    in.PLOT_ON = true;
    ratemap = getRatemap(in);
    hold on;
    axis off
    in.spikeTS = sp.timestamp(sp.cut==c);
    in.psr = trial.posData.sample_rate;
    in.x = ratemap.map;
    in.nodwell = ratemap.nodwell;
    in.PLOT_ON = false;
    in.tol = 10^-10;
    sac = autoCorr2D(in);

    % [4] SAC metrics
    nexttile
    in.sac = sac.autocorrelogram;
    in.FIND_CENTROID = true;
    in.ORIENTATION = true;
    in.GET_MEAN_R_AND_SCALE = true;
    in.FULL_GRIDNESS_RANGE = true; 
    in.PLOT_ON = true;
    metrics = autoCorrProps(in);
    axis off
    hold off;

    % save stats 
    stats = metrics;
    
end