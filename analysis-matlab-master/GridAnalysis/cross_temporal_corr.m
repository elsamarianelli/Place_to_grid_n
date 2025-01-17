function [ac, bin_left_edge]=cross_temporal_corr(times_a, times_b, bin_size, max_delta, duration)
%
% max_delta = 2; %seconds
% bin_size = 0.05; %seconds.
% duration = 1200; %seconds
%
%---- example dummy data ---
% n_a = 1123;
% n_b = 800;
% times_a = sort(rand(n_a,1))*duration;
% times_b = sort(rand(n_b,1))*duration;
% ------------------
%
% DM, Apr 2016
%

% we are concerned with a temporal window +-max_delta seconds around each spike from cell a.
% we want to find the first and last spike in b that fall in each of these windows.
% We make use of the 2nd output of histc which tells you: "if I try and put
% a point into one of the bins, starting from the left, how many bin edges
% will I cross before I reach the correct bin".  If the point falls outside 
% all the bins the answer given is a meaingless "0", so we avoid this by
% adding in extra wide bins at the start and end. 
times_b_capped = [-max_delta; times_b; duration+max_delta]; 
[~, idx_in_b_window_start] = histc(times_a-max_delta, times_b_capped);
[~, idx_in_b_window_end] = histc(times_a+max_delta, times_b_capped);

% correct for capping..
% Note that (somewhat unsusually) the result of all this is that idx_start is
% exclusive and end_idx is inclusive.
idx_in_b_window_start=idx_in_b_window_start-1;
idx_in_b_window_end=idx_in_b_window_end-1; 

n_b_in_window = idx_in_b_window_end - idx_in_b_window_start;

% now we perpare the indices for our all-important subtraction
idx_a = repeatInd(n_b_in_window);
idx_b = idx_in_b_window_start(idx_a) + countTo(n_b_in_window);

% and we can now do the subtraction 
delta_t = times_b(idx_b) - times_a(idx_a);

% ok, all that's left is to make a histogram from the delta_t's
if mod(max_delta, bin_size) ~= 0
    error('max_delta must be multiple of bin_size '); 
end;
bin_left_edge = -max_delta:bin_size:max_delta;
bin_idx = floor(delta_t/bin_size) + max_delta/bin_size + 1; 
ac = accumarray(bin_idx, 1, [numel(bin_left_edge),1], @sum);
ac = ac(1:end-1);%takes care of delta equals exactly max_delta
bin_left_edge = bin_left_edge(1:end-1)'; % drop extra bin edge, to agree with ac
end
