function ret = repeatInd(n)
% By example:
%    
%        #     1  2  3  4  5  6  7  8
%        n = [ 0, 3, 0, 0, 2, 0, 2, 1]
%        res = [2, 2, 2, 5, 5, 7, 7, 8]
%        
%    That is the input specifies how many times to repeat the given index.
%
% DM, Feb 2015 (originally written in Python, also by DM)
%
    n = n(:);
    n_mask = logical(n);
    n_inds = find(n);
    n_inds(2:end) = n_inds(2:end) - n_inds(1:end-1); % take diff and leave 1st value in place
    n_cumsum = ones('single');
    n_cumsum(2:length(n)+1) = cumsum(n)+1; 
    ret = zeros(n_cumsum(end)-1,1,'single');
    ret(n_cumsum(n_mask)) = n_inds; % note that n_mask is 1 element shorter than n_cumsum
    ret = cumsum(ret);
end
    