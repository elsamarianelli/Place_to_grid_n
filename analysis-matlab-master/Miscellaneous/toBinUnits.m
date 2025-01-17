function [ binUnits ] = toBinUnits( nd_data, nd_data_ranges, binsizes)
% See https://github.com/UCL/mTint/tree/master/LoadData for info

[nDataPoints,nDims] = size(nd_data);

% Check for correct input array sizes
if size(nd_data_ranges,1) > 2 || size(nd_data_ranges,1) < 1 || size(nd_data_ranges,2) ~= nDims
    error('Ranges array must have 1 or 2 rows & the same number of columns as the nd_data array.');
end

if size(binsizes,1) ~= 1 || size(binsizes,2) ~= nDims
    error('Binsizes array must have 1 row & the same number of columns as the nd_data array.');
end

% If nd_data_ranges has two rows, the first row is the minima and should be 
% subtracted from the data.
if size(nd_data_ranges,1) == 2
    nd_data = bsxfun(@minus,nd_data,nd_data_ranges(1,:));
    maxBinUnits = diff(nd_data_ranges,1)./binsizes; 
else
    maxBinUnits = nd_data_ranges./binsizes;
end

% Do the division to convert to bin units
binUnits = bsxfun(@rdivide,nd_data,binsizes);

% Cope with both kinds of bad points....
% (a) points are too low
bad_low_points = binUnits < eps; %note we use eps here so than ceil(binUnits) >= 1
binUnits(bad_low_points) = eps; 
% (b) points are too high
bad_high_points = bsxfun(@gt, binUnits, maxBinUnits); % gt="greater than"
[ignore,colInds] = find(bad_high_points);
binUnits(bad_high_points) = maxBinUnits(colInds);

if (sum(bad_low_points) + sum(bad_high_points)) / (nDataPoints*nDims) > 0.1
    warning('MATLAB:mtint:toBinUnits','More than a tenth of data points to be binned by toBinUnits are outside of the expected range.');
end


end
