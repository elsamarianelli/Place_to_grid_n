function [ phi ] = GetRates( x, y, cells )
%GETRATES Summary of this function goes here
%   Detailed explanation goes here
n_cells = length(cells);
phi = zeros(1,n_cells);

for i = 1:n_cells
    phi(i) = cells{i}.fmap(y,x);
end

phi(isnan(phi)) = 0;

end

