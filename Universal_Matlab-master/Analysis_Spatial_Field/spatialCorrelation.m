function [ spatialCorr, nValidBins ] = spatialCorrelation( rateMap1, rateMap2, stackMaps )
%SPATIAL_CORRELATION Conducts Pearson correlations on paired ratemaps (maps are interporlated
%to the same size linearly). Can also deal with stacked (population ratemaps) i.e.
%provided ratemaps should be 3d with third dim being different cells. Size of third dim
%must be the same between the two stacks
%
%NB ratemaps should have nans for unvisted bins
%
% Takes
% rateMap1 & rateMap2 - rates by bin, don't need to be exactly the same size
% [optional] stackMaps [true or false/ 1 or 0] indicating if ratemaps with multiple cells
%       are conctanated [true] and treated as a single map or if each should be compared
%       individually [false] and a mean spat corr taken. If not supplied defaults to
%       [true] i.e. stack
%
%
% RETURNS
% spatialCorr - spatial correlation (a Pearson coef). If stacked ratemaps are used (3d)
%           then a single spatialCorr is returned when stackMaps is true and multiple when
%           it is false (i.e. one for each comparison)
% nValidBins - number of bins used to calc corr (nb. mutual zeros and nans are
%           removed)



% --------------------------------------------------------------------------------------------------
% --- INTERPPORLATE FILES TO SAME DIMENSIONS -------------------------------------------------------
% --------------------------------------------------------------------------------------------------

%Check if third argument is present
if nargin==2,    stackMaps=true; end

%Check if these are stacked ratemaps and if so is third dim the same size
if size(rateMap1,3)~= size(rateMap2,3)
    error('spatial_correlation: Stacked ratemaps do not have matching size in 3rd dim.\n');
end
stackSize=size(rateMap1,3); %Also find stack size (i.e. size of dim3 that is no of cells)

%Find max size of each of first two dimensions - then scale up to this size
maxDimSize=max([size(rateMap1,1), size(rateMap1,2); size(rateMap2,1), size(rateMap2,2)],[],1); %mn pair


%Note when interporlating nans will spread - this is fine - now loop over dim3 (to deal
%with fact that I can't use interp3 for some reason on 3d mats
[interpRateMap1, interpRateMap2]=deal(zeros([maxDimSize, stackSize])); %pre allocate
for nn=1:stackSize
    interpRateMap1(:,:,nn)=interp2(rateMap1(:,:,nn), ....
        linspace(1, size(rateMap1,2), maxDimSize(2)), ...
        linspace(1, size(rateMap1,1), maxDimSize(1))', '*linear');
    
    interpRateMap2(:,:,nn)=interp2(rateMap2(:,:,nn), ...
        linspace(1, size(rateMap2,2), maxDimSize(2)), ...
        linspace(1, size(rateMap2,1), maxDimSize(1))', '*linear');
end

clear rateMap1  rateMap2

% --------------------------------------------------------------------------------------------------
% --- APPLY CRITERIA AND RUN CORRELATION -----------------------------------------------------------
% --------------------------------------------------------------------------------------------------

%Create masks for mutally zero firing (ignored). Mask also remove bins that were not visted
% and their partner from the other ratemap. Mask indicates bins to exclude
% Note excluding mutual zero is standard for place cells
maskMutualZero=(interpRateMap1==0)& (interpRateMap2==0);
maskNans= isnan(interpRateMap1) | isnan(interpRateMap2);
mask= maskMutualZero | maskNans;
clear maskMutualZero maskNans

nValidBins=sum(~mask(:));

%Do correlation on bins that qualify
if stackMaps %Treat as single large ratemap
    spatialCorr=corr(interpRateMap1(~mask), interpRateMap2(~mask), 'type', 'Pearson');
    
else %Treat as individual ratemaps - correlate each and take mean
    spatialCorr=zeros(stackSize,1);
    for  nn=1:stackSize
        tmpMask1=interpRateMap1(:,:,nn);
        tmpMask2=interpRateMap2(:,:,nn);
        spatialCorr(nn)=corr(tmpMask1(~mask(:,:,nn)), tmpMask2(~mask(:,:,nn)), 'type', 'Pearson');
    end

    %         stdSpatialCorr=std(tmpCorr);
end

