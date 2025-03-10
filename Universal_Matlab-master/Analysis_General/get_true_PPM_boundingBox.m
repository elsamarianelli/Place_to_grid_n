function [ppm, iRange, jRange]=get_true_PPM_boundingBox(binPos, binSz, envSz)
%GET_TRUE_PPM_BOUNDINGBOX Find true PPM of square pos map with known bin sz
% Often PPM is incorrectly specicified in the data.pos.header if the 
% environment is square and of a known size then the PPM can be infered
% from the animals path. This code takes the binned pos map for a square
% envioronment - does not have to be perfectly square to the camera - and
% given the bin size in pix as well as the envionrment wall length in m
% retuns an estimate of the true PPM as well as indicies into the i and j
% range of binPos specifying which bins are occuped - these ranges can be
% used for a close crop.
%
% NOTE. This code is only reliable if the environment is well visited - it
% tries to eliminate intermitant sections of the path reaching outside the
% environment but if the environment is poorly sampled will typically
% underestimate the area visited and hence will over estimate ppm
%
% TAKES
% binPos        binned positional data generated using binData.m [2d mat]
% binSz         bin size in camera pixels as provided to binData.m [scalar]
% envSz         side length of square environment in meters [scalar]
%
% RETURNS
% ppm           estimate of true ppm based on square environment [scalar]
% iRange        i axis (i.e. y axis) range of binPos that is occupied (e.g. 2:45)
% jRange        similar to iRange but for j axis (i.e. x axis)

% First find and label visited bins in binPos - labelling is done to find
% regions of 4 edge continuity
visitedBins         =bwlabel(binPos>0,4);

% Second find bounding box that can be used to crop binPos for just visited
% regions
% Note bounding box returns (for a 2d mat) [upper left x coord, uppperleft 
% y coord, x width, y width] there are all using ij coordinates with origin 
% in upper left. But not for a positional matrics in which the bins are 
% centred on the axis integers then the returned bounding box will be for 
% the edge of the bins so a pos map occupying starting in bin (5,6) in ij 
% would be reported as [5.5, 4.5, ...] for the same reason the width will 
% also appear to be one bin bigger as it is around the edges of the bins.
stats               =regionprops(visitedBins, 'boundingbox', 'filledarea');
[~, correctStat]    =max([stats.FilledArea]); %Take largest area returned
filledArea          =stats(correctStat).FilledArea;
boundingBox         =stats(correctStat).BoundingBox;
clear stats
   
%So to select the occupied part of the ratemap need to take from topLeftBinIJ to
%toLeftBinIJ plus widthIJ
ijCorner            =[ceil(boundingBox(2)), ceil(boundingBox(1))];
ijWidth             =[boundingBox(4)-1, boundingBox(3)-1];
iRange              =ijCorner(1):ijCorner(1)+ijWidth(1)-1;
jRange              =ijCorner(2):ijCorner(2)+ijWidth(2)-1;


% Third calculate PPM - do this by taking the sqrt of the binned area to
% get the number of bins per side length, then convert via pix to PPM. This
% should be robust to situations where the environment is tilted relative
% to the camera.
binPerSide          =sqrt(filledArea); %Worked out from filled area
ppm                 =binPerSide * binSz / envSz; 
 


end

