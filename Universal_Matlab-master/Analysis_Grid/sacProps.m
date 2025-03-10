function stats =sacProps(sac)
%SACPROPS Returns key metrics from SAC
% Based on original function autoCorrProps - note some of this code is old
% (dates from 2006). To do: speed up and refactor code but maintain
% backward compatability.
%
% Returns - as a structure - a range of properties relating to the spatial
% autocorrelation (SAC) of a grid cells (or otherwise).
%
% TAKES
% sac           Smooth SAC of grid cell (normal to calculate based on a
%               smooth ratemap in which case no need to smooth the SAC)
%
% RETURNS
% Note all the following are returned as fields of the strcuture 'stats'
% scale         Grid scale in bins. Calculated as the median distance of 6
%               central peaks from centre of the SAC
%
% gridness      Basic measure of 6 fold rotational symetry. This is the
%               original version based on Sargolini et al (2006) and used
%               in Barry et al (2007) etc.
%
% peakOrient    Orientation of the three main axes of the grid pattern 
%               [sorted] in degs CCW from x-axis. Hence the first value of
%               peakOrient is the overall orientation.
%
% peakCoord     [6x2] position of cloest peaks as xy pairs with origin top 
%               left
%
% perimPeakMask [size(sac)] Bins that define edge of central and six peaks
%
% perimGridMask [size(sac)] Bins that define edge of area used for gridness 
%               calc.
%
% meanROfPeaks  [1] Mean of peak value of the 6 central peaks found in SAC
%
% NB. Unit for wavelength is bins of autocorr
%
% Note: In the autocorrelogram it is sensible to exclude bins that were 
% constructed with relatively small overlap between the ratemap1 and 
% ratemap2 (Hafting excludes bins with an overlap of 20 or
% less). Set these bins to 0 before passing to this function



% -------------------------------------------------------------------------
% --- MAIN FUNCTION -------------------------------------------------------
% -------------------------------------------------------------------------
sac         =real(sac); %Ignore complex


sac(isnan(sac))=-1; %Sub nans for -1
sacTmp      =sac;        % TW. Do not allow local max that have .. lines addopted from TW code
sacTmp(sac<=0) =-1; %  .. r-value below zero.

%Don't consider imaginary components which some times appear in shuffled data
sacTmp      =real(sacTmp);
peaksAutoCorr=imregionalmax(sacTmp); %Find local maxima
clear sacTmp

[lableMask, ~]=bwlabel(peaksAutoCorr, 8);

%In case adjacent points share maxima find the centroid of them - NB returns structure array
%stats(n).Centroid containing for each peak the x,y position of the centroid but y is
%counting down from origin at top left
regStats      =regionprops(lableMask, 'Centroid');
%NL. [n x 2] pairs of x,y coord for max points
xyCoordMaxBin=reshape([regStats.Centroid], 2,[])'; %Still x,y pair

% Convert to a new reference frame which as the origin at the centre of the autocorr
% and with y negative at top and y positive at bottom, x negative on left and postive on
% right
% NB autocorr will always have sides with odd number of bins
centralPoint=ceil(size(sac)/2); %m,n pair
xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPoint), [size(xyCoordMaxBin,1), 1]);

%Calculate distance of peaks from centre point and find seven closest (one will be central peak
%disregard this)
distFromCentre=sum(xyCoordMaxBinCentral.^2,2).^0.5;
[~, orderOfClose]= sort(distFromCentre);

%Get id of closest peaks - note closest peak 1 will be centre
if length(orderOfClose)>=7; closestPeaks=orderOfClose(1:7); %Might be fewer than 7 peaks
else closestPeaks=orderOfClose(1:end);
end


%x,y pairs in cartesian coords with y counting down from top and origin at top left.
peakCoord=xyCoordMaxBin(closestPeaks,:);
peakCoord=round(peakCoord); %Should be integer
%also x,y pairs with y counting down from top and origin in center
closestPeaksCoordCentral=xyCoordMaxBinCentral(closestPeaks,:);

% ----------------------------------------------------------------------------------------
% --- R VALUE OF PEAKS -------------------------------------------------------------------
% ----------------------------------------------------------------------------------------
meanROfPeaks=ones(6,1)*nan;
for nnCP=1:length(closestPeaks)
    meanROfPeaks(nnCP)=sac(find(lableMask==closestPeaks(nnCP),1));
end
meanROfPeaks=mean(meanROfPeaks(~isnan(meanROfPeaks)));

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
peakMasks=zeros([size(sac), size(peakCoord,1)+1]);
perimeterMasks=zeros([size(sac), size(peakCoord,1)+1]);

for n=1:size(peakCoord,1)
    %NB. Have to flip dimensions as these are x,y pairs and need to be m,n
    [peakMasks(:,:,n), perimeterMasks(:,:,n)]=...
        sf_findPeakExtent(sac, closestPeaks(n), [peakCoord(n,2), peakCoord(n,1)]);
end
peakMask=max(peakMasks, [], 3); %Colapse to 2d
perimPeakMask=max(perimeterMasks,[],3);%Colapse to 2d
clear peakMasks perimeterMasks

lableMask=peakMask; %Needed for code below

% -------------------------------------------------------------------------------------------------
% --- SCALE ----------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%Calculate scale measured in bins
closestPeakDistFromCentre=distFromCentre(closestPeaks(2:end));
scale   =median(closestPeakDistFromCentre);

% -------------------------------------------------------------------------------------------------
% --- ORIENTATION & NON CIRCULARITY --------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%Calculate orientation of grid - note orientation is measured counterclockwise from horizontal so
%basically polar coordinates.
%NB Changed after TW spotted error (y component should be negative as y
%axis is increses down in xyCoordMaxBinCentral) this will have affected orientation estimates but
%nothing else. NB2 have checked against ratemaps and this correction was
%needed
% NB3. Just to reiterate: negative y is necessary because DACQ camera coordinates have the
% origin at top left hence y axis increase down (imagesc needs to be dispalyed with axis
% ij as a result). But cart2pol assumes origin is bottom left so have to make correction
% before using cart2pol.

%Nb closestPeaksCoordCentral inlcudes centre as first point
[th, ~]=cart2pol(closestPeaksCoordCentral(:,1), -closestPeaksCoordCentral(:,2)); %Y is negative
peaksAboveXaxis=find(th>=0); %Remove negative values - can do this as peaks are 180deg radially symetrical
peaksAboveXaxis=peaksAboveXaxis(peaksAboveXaxis~=1); %Remove central point is present - should always be first in list
% closestPeaksCoordCentral=closestPeaksCoordCentral(peaksAboveXaxis, :);
peakOrient=rad2deg(sort(th(peaksAboveXaxis)));


% --------------------------------------------------------------------------------------------------
% ---- GRIDNESS ------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
[meshX, meshY]=meshgrid(1:size(sac,2), 1:size(sac,1)); %meshX increase to right, meshy increases down
%Change coordinates so that origin is at centre and to facilitate pythag
meshXCentre=(meshX-centralPoint(2)).^2;
meshYCentre=(meshY-centralPoint(1)).^2;

maxDistFromCentre=zeros(length(closestPeaks),1);



for k=1:length(closestPeaks)
    isPresent=lableMask==closestPeaks(k);
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
gridnessMask=((meshXCentre+meshYCentre).^0.5)<=maxDistFromCentre; %Mask for points within bound of outer peak
gridnessMask=gridnessMask & ~(lableMask==closestPeaks(1));
perimGridMask=bwperim(gridnessMask); %Extent of area used to calc gridness
clear maxDistFromCentre meshX meshY

selA    =sac(gridnessMask);
%Now loop round rotating the autoCorr, selecting central section and calculating correlation. Am
%going to do this for all 3deg intervals although only 5 position are really required to calculate
%gridness. Change if slow.
%NB rotation is counterclockwise - doesn't make any difference to gridness measure
count=1;
rotatedCorr=ones(5,1)*nan;
for rotAmnt     =[60, 120, 30, 90, 150] %Change to '0:6:174' to calculate full gridness range
    selB    =imrotate(sac,rotAmnt, 'bilinear', 'crop');
    selB    =selB(gridnessMask);
    valid   =~(isinf(selA) | isinf(selB) | isnan(selA) | isnan(selB));
    rotatedCorr(count)=corr(selA(valid), selB(valid), 'type', 'Pearson');
    count   =count+1;
end
clear selA selB

%Finally gridness is the difference between the lowest correlation at 60deg and 120deg and the
%highest correlation at 30, 90 and 150deg
% gridness= min(rotatedCorr([11, 21])) - max(rotatedCorr([6, 16, 26])); %Uncomment this one if using full gridness range
gridness= min(rotatedCorr([1, 2])) - max(rotatedCorr([3, 4, 5])); %This one for quick verison


%Combine results into a structure to return
stats.scale     =scale;
stats.gridness  =gridness;
stats.peakOrient=peakOrient;
stats.peakCoord =peakCoord;
stats.perimPeakMask=perimPeakMask;
stats.perimGridMask=perimGridMask;
stats.meanROfPeaks =meanROfPeaks;

end


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
end
