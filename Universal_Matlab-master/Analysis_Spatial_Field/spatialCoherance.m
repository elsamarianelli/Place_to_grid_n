function [ sc ] = spatialCoherance( rm )
%SPATIALCOHERANCE Calculates spatial coherance of a ratemap
% SC is the correlation of a all occupied bins with the mean of their
% neighbours.
%
% TAKES rm - ratemap, can be smoothed or unsmoothed, though SHOULD BE
% UNSMOOTHED. Unvisted bins should be nans
%
% NB. Tested against results from TW spatail coherance and they both give exactly the same
% values.

%Create filter to calc neighbours
neighbourFilt=ones(3);
neighbourFilt(2,2)=0;

%Deal with unoccupied bins
visitedBins = double(~isnan(rm));
rm(isnan(rm))=0;

%Work out mean firing of neighbouring bins
neighbourFiring=conv2(rm, neighbourFilt, 'same');
nVistedNeigh = conv2(visitedBins, neighbourFilt, 'same');
warning('off', 'MATLAB:divideByZero');
meanNeighFiring = neighbourFiring./nVistedNeigh;
warning('on', 'MATLAB:divideByZero');

%Correlate each occupied bin with the mean of it's occupied neighbours
indVistedBins = find(visitedBins);
%Remove any nans
validMeanNeighFiring=find(~isnan(meanNeighFiring));
valid=intersect(validMeanNeighFiring, indVistedBins);
sc = corr(rm(valid), meanNeighFiring(valid), 'type', 'Pearson');

%Strictly should be Fisher's z-transform of correlation (as defined by
%Muller & Kubie 1989) though tint doesn't do this
% sc=atanh(sc);