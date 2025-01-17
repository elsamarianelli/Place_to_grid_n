% loading results generated from clean_retry.m

% Specify the directory where the results are saved
output_dir = 'gridness_results_rect'; 
n_iterations = 1;  % Number of iterations to load
% load results 
results = load_results(output_dir,n_iterations);
%%

prac = load_trial();


map = results.cells_U{1}{20}.grid;
figure; imagesc(map)
close gcf 

in.x = map;                  
in.nodwell = ones(size(map));
autoCorr = autoCorr2D(in);
figure; imagesc(autoCorr.autocorrelogram)