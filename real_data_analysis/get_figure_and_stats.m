function stats = get_figure_and_stats(trial, t, c, in)

    figure; 
    fig_title = ['tetrode ', num2str(t), ', cell ', num2str(c)];
    title(fig_title)
    
    % [1] Rate Map / Polar Plot 
    % spikeData(3) cut==5 as example with firing
    % cell 5, tetrode 3
    subplot(2, 2, 1)
    sp = trial.spikeData(t);
    in.spikePosInd = sp.pos_sample(sp.cut==c);
    in.posXY = trial.posData.xy;
    in.MAP_TYPE = 'rate'; % 'dir' for polar plot
    in.PLOT_ON = true;
    ratemap = getRatemap(in);
    hold on;
    
    % [2] Shuffled Gridness 
    % blue - cumulative density plot of shuffled gridness 
    % red  - unshuffled gridness
    subplot(2, 2, 2)
    in.spikeTS = sp.timestamp(sp.cut==c);
    in.SHUFFLES = 1000;
    in.psr = trial.posData.sample_rate;
    in.PLOT_ON = true;
    ret = shuffledGridness(in);
    hold on;
    
    % [3] autocorrelation of 2d maps 
    % using FFT (fast fourier transform)
    subplot(2, 2, 3)
    in.x = ratemap.map;
    in.nodwell = ratemap.nodwell;
    in.PLOT_ON = true;
    in.tol = 10^-10;
    sac = autoCorr2D(in);
    hold on;
    
    % [4] SAC metrics
    subplot(2, 2, 4)
    in.sac = sac.autocorrelogram;
    in.FIND_CENTROID = true;
    in.ORIENTATION = true;
    in.GET_MEAN_R_AND_SCALE = true;
    in.FULL_GRIDNESS_RANGE = true; 
    in.PLOT_ON = true;
    metrics = autoCorrProps(in);
    hold off;

    % save stats 
    stats = metrics;
    
end