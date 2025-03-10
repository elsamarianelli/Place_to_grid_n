function [ scale, gridness, axisOrient] = sacProps90deg( sac )
% SACPROPS90DEG Return key metrics for 90deg grid cells from SAC.
% Standard gridness type measures look for 6 peaks (i.e. hexagonal grid)
% this code is adapted to only look for 4 peaks (i.e. compatable with a
% square grid as is found in DNNs etc).
%
% NB this code is adapted from the orginal autoCorrProps.m
%
% ARGS
% sac           Smooth SAC [2d mat]
%
% RETURNS
% scale         [5 x1] First value is mean distance of 4 central peaks from 
%               centre measured in bins of the ratemap. Next 4 values are
%               the distances to the 4 peaks found (or less if less found)
%
% gridness      [Scalar] Measure of 90deg rotational symetry vs 45deg
%
% axisOrient    [2x1] Orientaion of 2 main axes in degs anti-clock from 
%               x-axis

%
% Note: In the autocorrelogram it is sensible to exclude bins that were constructed with relatively
% small overlap between the ratemap1 and ratemap2 (Hafting excludes bins with an overlap of 20 or
% less). Set these bins to 0 before passing to this function


% -------------------------------------------------------------------------------------------------
% --- MAIN FUNCTION -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
sac             =real(sac); %Ignore complex


sac(isnan(sac)) =-1; %Sub nans for -1
sacTemp         =sac;        % Do not allow local max that have ...
sacTemp(sac<=0) =-1; %  .. r-value below zero.

peaksAutoCorr   =imregionalmax(sacTemp); %Find local maxima
clear sacTemp

[lableMask, numPeaks]=bwlabel(peaksAutoCorr, 8);

%In case adjacent points share maxima find the centroid of them - NB returns structure array
%stats(n).Centroid containing for each peak the x,y position of the centroid but y is
%counting down from origin at top left
stats=regionprops(lableMask, 'Centroid');
%NL. [n x 2] pairs of x,y coord for max points
xyCoordMaxBin=reshape([stats.Centroid], 2,[])'; %Still x,y pair

% Convert to a new reference frame which as the origin at the centre of the autocorr
% and with y negative at top and y positive at bottom, x negative on left and postive on
% right
centralPoint    =ceil(size(sac)/2); %m,n pair
xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPoint), [size(xyCoordMaxBin,1), 1]);

%Calculate distance of peaks from centre point and find 5 closest (one will be central peak
%disregard this)
distFromCentre  =sum(xyCoordMaxBinCentral.^2,2).^0.5;
[~, orderOfClose]=sort(distFromCentre);

%Get id of closest peaks - note closest peak 1 will be centre
if length(orderOfClose)>=5; closestPeaks=orderOfClose(1:5); %Might be fewer than 5 peaks
else closestPeaks=orderOfClose(1:end);
end


%x,y pairs in cartesian coords with y counting down from top and origin at top left.
closestPeaksCoord=xyCoordMaxBin(closestPeaks,:);
closestPeaksCoord=round(closestPeaksCoord); %Should be integer
%also x,y pairs with y counting down from top and origin in center
closestPeaksCoordCentral=xyCoordMaxBinCentral(closestPeaks,:);


% ----------------------------------------------------------------------------------------
% --- FIND FIELD AROUND EACH PEAK --------------------------------------------------------
% ----------------------------------------------------------------------------------------
%As we've just defined peak of each field need to find the extent of the field around it. Define
%this as the area enclosed within the half-height of the peak. Do this by looping through eack peak:
%find areas of SAC above the half height, then find the patch that includes the peak of interest.
%Could be a problem is peaks overlap as they might.
%Create two results: a labelMatrix with the pixels of each field labeld with diff numbers (i.e. same
%IDs as used above) and a matrix with the perimeter of each peak marked

%Two 3d mats [sizeAutoCorr1, sizeAutoCorr2, numberPeaks plus centralPeak]
peakMasks       =zeros([size(sac), size(closestPeaksCoord,1)+1]);
perimeterMasks  =zeros([size(sac), size(closestPeaksCoord,1)+1]);

for n           =1:size(closestPeaksCoord,1)
    %NB. Have to flip dimensions as these are x,y pairs and need to be m,n
    [peakMasks(:,:,n), perimeterMasks(:,:,n)]=...
        sf_findPeakExtent(sac, closestPeaks(n), [closestPeaksCoord(n,2), closestPeaksCoord(n,1)]);
end
lableMask       =max(peakMasks, [], 3); %Colapse to 2d
clear peakMasks perimeterMasks


% -------------------------------------------------------------------------------------------------
% --- SCALE ----------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%Calculate scale measured in bins
tmp             =distFromCentre(closestPeaks(2:end));
scale           =mean(tmp); %Mean distance to 4 main peaks
scale           =[scale; tmp(:)]; %Return vector of mean then distance to all 4 peaks

% -------------------------------------------------------------------------------------------------
% --- ORIENTATION & NON CIRCULARITY --------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%Calculate orientation of grid - note orientation is measured counterclockwise from horizontal so
%basically polar coordinates.

if numPeaks==1 %If only central peak found
    axisOrient  =nan;
else
    %Nb closestPeaksCoordCentral inlcudes centre as first point
    th          =cart2pol(closestPeaksCoordCentral(:,1), -closestPeaksCoordCentral(:,2)); %Y is negative
    peaksAboveXaxis=find(th>=0); %Remove negative values - can do this as peaks are 180deg radially symetrical
    peaksAboveXaxis=peaksAboveXaxis(peaksAboveXaxis~=1); %Remove central point is present - should always be first in list
    axisOrient  =rad2deg(sort(th(peaksAboveXaxis)));
end



% --------------------------------------------------------------------------------------------------
% ---- GRIDNESS ------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------


[meshX, meshY]  =meshgrid(1:size(sac,2), 1:size(sac,1)); %meshX increase to right, meshy increases down
%Change coordinates so that origin is at centre and to facilitate pythag
meshXCentre     =(meshX-centralPoint(2)).^2;
meshYCentre     =(meshY-centralPoint(1)).^2;
maxDistFromCentre=zeros(length(closestPeaks),1);

for k=1:length(closestPeaks)
    isPresent   =lableMask==closestPeaks(k);
    %Deal with situations where the mask of one peak has bloted out another entirely
    if ~sum(isPresent(:)), maxDistFromCentre(k)=0; continue, end
    maxDistFromCentre(k)=max((meshXCentre(lableMask==closestPeaks(k)) + meshYCentre(lableMask==closestPeaks(k))).^0.5);
end
maxDistFromCentre=max(maxDistFromCentre);


%MaxDistFromCentre can not be greater than the distance from the centre to the nearest edge of
%autoCorr - if it is larger reset it to this lower value
if maxDistFromCentre>min(floor(size(sac)/2))
    maxDistFromCentre=min(floor(size(sac)/2));
end

 
%Create a circular mask using the meshgrids I've created above
gridnessMask    =((meshXCentre+meshYCentre).^0.5)<=maxDistFromCentre; %Mask for points within bound of outer peak
gridnessMask    =gridnessMask & ~(lableMask==closestPeaks(1));
clear maxDistFromCentre meshX meshY

selA             =sac(gridnessMask);
%Now loop round rotating the autoCorr,
%NB rotation is counterclockwise - doesn't make any difference to gridness measure
count=1;
rotatedCorr     =ones(5,1)*nan;
for rotAmnt     =[90, 45, 135] %Change to '0:6:174' to calculate full gridness range
    selB        =imrotate(sac,rotAmnt, 'bilinear', 'crop');
    selB        =selB(gridnessMask);
    valid       =~(isinf(selA) | isinf(selB) | isnan(selA) | isnan(selB));
    rotatedCorr(count)=corr(selA(valid), selB(valid), 'type', 'Pearson'); 
    count       =count+1;
end
clear selA selB

%Finally gridness90 is the difference between the  correlation at 90deg and
%the max of the correlation at 45 or 135deg
gridness= min(rotatedCorr(1)) - max(rotatedCorr([2,3])); %This one for quick verison



% ----------------------------------------------------------------------------------------
% --- SUB FUNCTIONS ----------------------------------------------------------------------
% ----------------------------------------------------------------------------------------

function [peakMask, perimeterMask]=sf_findPeakExtent(autoCorr, peakId, peakCoord)
%Finds extent of field that belongs to each peak - defined as area in half-height and also
%perimieter. NB. peakCoord must by m,n pair in normal matrix coords

peakMask=zeros(size(autoCorr));
perimeterMask=zeros(size(autoCorr));
%Next line defines threshold used to find peak - currently using half height
aboveHalfHeightMask=bwlabel(autoCorr>(autoCorr(peakCoord(1),peakCoord(2)).* (1/2)),8);
peakIdTemp=aboveHalfHeightMask(peakCoord(1),peakCoord(2));
peakMask(aboveHalfHeightMask==peakIdTemp)=peakId;
perimeterMask(bwperim(aboveHalfHeightMask==peakIdTemp))=peakId;
