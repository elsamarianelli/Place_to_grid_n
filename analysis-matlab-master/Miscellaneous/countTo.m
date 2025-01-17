function ret = countTo(n)
% By example:
%    
%        #     1  2  3  4  5  6  7  8
%        n = [ 0, 3, 0, 0, 2, 0, 2, 1]
%        res = [1, 2, 3, 1, 2, 1, 2, 1]
%        
%    That is we count from 1 to n[i] for each i
%
% DM, Feb 2015
%
    n = n(:);
    n_mask = logical(n);
    n_cumsum = cumsum(n); 
    ret = ones(n_cumsum(end)+1,1,'single');
    ret(n_cumsum(n_mask)+1) = -n(n_mask)+1; 
    ret = cumsum(ret(1:n_cumsum(end)));
end
    