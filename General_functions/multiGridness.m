function [ stGrd, expGrd, scale ] = multiGridness( sac, shape )
%MULTIGRIDNESS Returns standard & expanding annulus version of gridness. NB
%can be combined with code that corrects for eliptical distortion in SAC.
%
% TAKES
% sac       [2d mat] spatial autocorr - should be constructed from smooth 
%           ratemap
%
% RETURNS
% stdGrd    standard gridness. Used in Hafting(2005)
% expGrd    expanding gridness - test increasing radi & take highest. Used
%           in Brandon (2011) & Killian (2012)
%
%Build based on CB's autoCorrProps so uses some preprocessing from that
%function. Note stdGrd should be identical to the value returned for
%gridess by autoCorrProps but is faster. expGrd checks gridness for
%multiple annuli and takes the highest - has been used by Moser group &
%Brandon (2011).
%
%To ipliment the elipse correction versions of gridness (e.g. similar to
%those used by Brandon (2011) and Killian (2012) first generate a
%regularised SAC by running the original SAC through gridEllipse_fit.m
%then gridEllipse_correct.m (which, respecitivly, fit an elipse to the SAC
%and then regualirsed it) and finally put the regularised SAC through 
%this function.



% --- MAIN ----------------------------------------------------------------

% ---Preprocess the sac ---------------------------------------------------
sac             =real(sac); %Ignore complex, shouldn't be present but can get introduced
sac(isnan(sac)) =-1; %Sub nans for -1



% --- Find peaks ----------------------------------------------------------
%Require that peaks - maxima - must have a value > 0
tmp             =sac;
tmp(sac<=0)     =-1;
peakSAC         =imregionalmax(tmp); %Find local maxima
clear tmp

[lableMask]     =bwlabel(peakSAC, 8);

%In case adjacent points share maxima find the centroid of them - NB returns structure array
%stats(n).Centroid containing for each peak the x,y position of the centroid but y is
%counting down from origin at top left
stats           =regionprops(lableMask, 'Centroid');
%NL. [n x 2] pairs of x,y coord for max points
xyCoordMaxBin   =reshape([stats.Centroid], 2,[])'; %Still x,y pair

% Convert to a new reference frame which as the origin at the centre of the autocorr
% and with y negative at top and y positive at bottom, x negative on left and postive on
% right
% NB autocorr will always have sides with odd number of bins
centralPoint    =ceil(size(sac)/2); %m,n pair
xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPoint), [size(xyCoordMaxBin,1), 1]);
x = xyCoordMaxBin(:, 1);
y = xyCoordMaxBin(:, 2);
% figure;
% imagesc(sac)
% hold on; plot(x, y, '.')
%Calculate distance of peaks from centre point and find seven closest (one will be central peak
%disregard this)
distFromCentre  =sum(xyCoordMaxBinCentral.^2,2).^0.5;
[~, orderOfClose]=sort(distFromCentre);

if strcmp(shape, 'hexagon')
    %Get id of closest peaks - note closest peak 1 will be centre
    if length(orderOfClose)>=7; closestPeaks=orderOfClose(1:7); %Might be fewer than 7 peaks
    else closestPeaks=orderOfClose(1:end);
    end
elseif strcmp(shape, 'square') % looking at non-hexagonal gridness
    if length(orderOfClose)>=5; closestPeaks=orderOfClose(1:5); %adjust to 5 peaks (4 surrounding one in centre)
    else closestPeaks=orderOfClose(1:end);
    end
end

%x,y pairs in cartesian coords with origin at top left.
closestPeaksCoord=xyCoordMaxBin(closestPeaks,:);
closestPeaksCoord=round(closestPeaksCoord); %Should be integer
%also x,y pairs with y counting down from top and origin in center
closestPeaksCoordCentral=xyCoordMaxBinCentral(closestPeaks,:);
% 
% figure;
% imagesc(sac); hold on;
% x = closestPeaksCoord(:, 1);
% y = closestPeaksCoord(:, 2);
% plot(x, y, '.', 'MarkerSize', 10)

% --- FIND FIELD AROUND EACH PEAK ----------------------------------------
%As we've just defined peak of each field need to find the extent of the field around it. Define
%this as the area enclosed within the half-height of the peak. Do this by looping through eack peak:
%find areas of SAC above the half height, then find the patch that includes the peak of interest.
%Could be a problem is peaks overlap as they might.
%Create two results: a labelMatrix with the pixels of each field labeld with diff numbers (i.e. same
%IDs as used above) and a matrix with the perimeter of each peak marked

%Two 3d mats [sizeAutoCorr1, sizeAutoCorr2, numberPeaks plus centralPeak]
peakMasks       =zeros([size(sac), size(closestPeaksCoord,1)+1]);
perimeterMasks  =zeros([size(sac), size(closestPeaksCoord,1)+1]);

for n=1:size(closestPeaksCoord,1)
    %NB. Have to flip dimensions as these are x,y pairs and need to be m,n
    [peakMasks(:,:,n), perimeterMasks(:,:,n)]=...
        sf_findPeakExtent(sac, closestPeaks(n), [closestPeaksCoord(n,2), closestPeaksCoord(n,1)]);
end
peakMask        =max(peakMasks, [], 3); %Colapse to 2d
clear peakMasks perimeterMasks

lableMask       =peakMask; %Needed for code below

% --- Gridness time -------------------------------------------------------

% -- Do some pre calculations to use in the different gridness measures
%Create mat of distances from centre at each bin
[meshX, meshY] =meshgrid(1:size(sac,2), 1:size(sac,1)); %meshX increase to right, meshy increases down
%Change coordinates so that origin is at centre and to facilitate pythag
centreRadi      =sqrt((meshX-centralPoint(2)).^2 + (meshY-centralPoint(1)).^2);
clear mesh*

maxAllowRadi    =min(floor(size(sac)./2)); %Max radius is to edge of sac


% -- 1) First do stGrd (i.e. standard gridness)
%Basic idea of standard gridness is to examine each of the fields
%identified in the code above. Find the furthest edge of the fields from
%the centre. Then take a circular area with a max radius equal to this
%value and excluding the central peak. Use this area to calcuate gridness.

outerRadi       =zeros(length(closestPeaks),1);
for k=1:length(closestPeaks)
    isPresent       =lableMask==closestPeaks(k);
    %Deal with situations where the mask of one peak has bloted out another entirely
    if ~sum(isPresent(:)), outerRadi(k)=0; continue, end
    outerRadi(k)    =max(centreRadi(lableMask==closestPeaks(k)));
end
outerRadi       =max(outerRadi);
if outerRadi>maxAllowRadi, outerRadi=maxAllowRadi; end

%Create a circular mask using the radi identified excluding central field
gridnessMask    =centreRadi<=outerRadi & ~(lableMask==closestPeaks(1));


selA            =sac(gridnessMask);
isValidA        =(~isnan(selA) & ~isinf(selA)); %values from selA that are a real num
%Loop to rotate autocorr & calc corr.
count           =1;
% switch for square or hexagon gridness check
if strcmp(shape, 'hexagon')
    rotatedCorr     =zeros(5,1);
    angles_oi = [60, 120, 30, 90, 150];
elseif strcmp(shape, 'square')
    rotatedCorr     =zeros(4,1);
    angles_oi = [ 90, 180, 45, 135];
end

for rotAmnt         =angles_oi
    selB            =imrotate(sac,rotAmnt, 'bilinear', 'crop');
    selB            =selB(gridnessMask);
    isValidB        =(~isnan(selB) & ~isinf(selB));
    isValidAB       =isValidA & isValidB;
    rotatedCorr(count)=corr(selA(isValidAB), selB(isValidAB)); 
    count=count+1;
end
clear selA selB count rotAmnt

%Finally gridness is the difference between the lowest correlation at 60deg and 120deg and the
%highest correlation at 30, 90 and 150deg

if strcmp(shape, 'hexagon')
    stGrd           =min(rotatedCorr([1, 2])) - max(rotatedCorr([3, 4, 5])); %This one for quick verison
elseif strcmp(shape, 'square')
    stGrd           =min(rotatedCorr([1, 2])) - max(rotatedCorr([3, 4])); 
end


% -- 2) Second do expGrd (i.e. expanding radius gridness)    annulus radius
% also allowed to be expanded 
%Taking the definition of this method from Killian (2012). Briefly: define
%inner radius as half of mean distance (in bins) to the 6 central peaks -
%if less than 6 are found just use 4. Then define outer radius starting at
%inner plus 2 bins all the way upto the outer edge of the sac. For each of
%these calculate gridness as normal and take max. [See SI figure 1 of
%Killian (2012) for description.

innerRadi       =round(mean(sqrt(sum(closestPeaksCoordCentral.^2)))/2);
outerRadiRange  =innerRadi+2:maxAllowRadi;
nRadi           =length(outerRadiRange);

%Loop over rotations and radi
if strcmp(shape, 'hexagon')
    num = 5;
elseif strcmp(shape, 'square')
    num = 4;
end

allCorr         =zeros(nRadi,num);
rotAmnt         =angles_oi;
for k=1:num
    rotSac          =imrotate(sac, rotAmnt(k),'bilinear', 'crop'); %Rotate
    
    for n=1:nRadi
        outerRadi   =outerRadiRange(n);
        gridnessMask=centreRadi<=outerRadi & centreRadi>=innerRadi;
        sel         =sac(gridnessMask); %Select annulus from unrotated
        selRot      =rotSac(gridnessMask); %Select from rotate
        isValid     =~isnan(sel) & ~isnan(selRot) & ~isinf(sel) & ~isinf(selRot);
        allCorr(n,k)=corr(sel(isValid), selRot(isValid)); %Do Pearson corr
    end
end
if strcmp(shape, 'hexagon')
    expGrd          =min(allCorr(:,1:2),[],2) - max(allCorr(:,3:5),[],2); %Grd for each radi
elseif strcmp(shape, 'square')
    expGrd          =min(allCorr(:,1:2),[],2) - max(allCorr(:,3:4),[],2); %Grd for each radi
end

expGrd          =max(expGrd,[],1); %Max of all radi
clear innerRadi outerRadiRange nRadi allCorr rotAmnt sel


% -------------------------------------------------------------------------------------------------
% --- SCALE ----------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%Calculate scale measured in bins
% 
closestPeakDistFromCentre = distFromCentre(closestPeaks(2:end));
scale   =median(closestPeakDistFromCentre);

end

% ----------------------------------------------------------------------------------------
% --- SUB FUNCTIONS ----------------------------------------------------------------------
% ----------------------------------------------------------------------------------------

function [peakMask, perimeterMask]=sf_findPeakExtent(autoCorr, peakId, peakCoord)
%Finds extent of field that belongs to each peak - defined as area in half-height and also
%perimieter. NB. peakCoord must by m,n pair in normal matrix coords

peakMask        =zeros(size(autoCorr));
perimeterMask   =zeros(size(autoCorr));
%Next line defines threshold used to find peak - currently using half height
aboveHalfHeightMask=bwlabel(autoCorr>(autoCorr(peakCoord(1),peakCoord(2)).* (1/2)),8);
peakIdTemp      =aboveHalfHeightMask(peakCoord(1),peakCoord(2));
peakMask(aboveHalfHeightMask==peakIdTemp)=peakId;
perimeterMask(bwperim(aboveHalfHeightMask==peakIdTemp))=peakId;

end