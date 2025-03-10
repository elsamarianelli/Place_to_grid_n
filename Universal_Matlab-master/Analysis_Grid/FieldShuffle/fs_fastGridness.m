function [ stdG, expG, gScl] = fs_fastGridness( sac )
%FAST_GET_GRIDNESS Cut down version of autocorr props - just returns
%gridness but does both the standard version (as I've used before) and the
%expanding version as used in Langston et al (2011). Note if I can't find a
%centre field extent then central field radius is set at 5bins (10cm)
%
% This code is called repeatedly in the shuffled grid code so will make
% some variables persistent to avoid having to recalculate.
%
% TAKES
% sac           Spatial autocorrr
%
% RETURNS
% stdG          Standard measure of regularity
% expG          Expanding gridness (i.e. Langston et al 2011)
% gScl          Grid scale in bins - mean of central peaks


% --- VARS ----------------------------------------------------------------
%Vars that control the way expanding gridness is applied - three params
%basically how many bins larger than central peak are taken, how large the
%steps are and how many bins less than shortest env length taken. In
%Langston these are (in bins which are 2.5cm) central area plus 5 bins, in
%1 bin steps until 5 bins short of environment width
expGMin             =5;
expGStep            =1;
expGMax             =5;



% --- PERSISTENT VARIABLES ------------------------------------------------
% Only change if the size of the orginal sac changes. Use persistent since
% this function will be called many times in a shuffled loop on very
% similar data
persistent szSAC centralPt dstFrmCentre hlfMinSz allMsks

if isempty(szSAC) || any(szSAC ~= size(sac)) %Size of sac has changed - recalc
    
    %Set size of sac - used to determine if sac has changed
    szSAC               =size(sac);
    hlfMinSz            =min(floor(szSAC/2));
    
    %Calc central point in arena - used to figure out how far points are
    %from centre
    centralPt           =ceil(szSAC/2); %m,n pair
    
    %Pre calc a grid of points indicating how far each bin is from centre
    %xx increase to right, yy increases down
    [xx,yy]             =meshgrid(1:szSAC(2), 1:szSAC(1));
    %Change coordinates so that origin is at centre and to facilitate pythag
    meshXCentre         =(xx-centralPt(2)).^2;
    meshYCentre         =(yy-centralPt(1)).^2;
    dstFrmCentre        =sqrt(meshXCentre + meshYCentre);
    dstFrmCentre        =round(dstFrmCentre);
    clear xx yy meshXCentre meshYCentre
    
    %Precalc all possible Langston mask starting with radius 1 up to radius
    %equal to half the width the smallest size. Requires a loop and NOTE
    %these do not have the central field subtracted.
    allMsks             = false(szSAC(1), szSAC(2), hlfMinSz);
    for nn              =1:hlfMinSz
        allMsks(:,:,nn)          =dstFrmCentre<=nn;
    end
    
    
end



% --- MAIN FUNCTION -------------------------------------------------------
sac=real(sac); %Ignore complex
sac(isnan(sac))     =-1; %Sub nans for -1
peaksSAC            =imregionalmax(sac);%Find local maxima in x,y coord

% autoCorrTemp        =sac;        % TW. Do not allow local max that have ..
% autoCorrTemp(sac<=0)=-1; %  .. r-value below zero.
%Don't consider imaginary components which some times appear in shuffled data
% sac                 =real(sac);
% peaksSAC       =imregionalmax(autoCorrTemp);
% clear autoCorrTemp

% Following chunk was previously used in case two adjacent peaks were equal
% hight - which is very unlikely - was finding centroid but this is slow
% and not really necessary. Now replace with just a direct find of the
% peaks
% [lableMask]=bwlabel(peaksAutoCorr, 8);
% stats=regionprops(lableMask, 'Centroid');
% %NL. [n x 2] pairs of x,y coord for max points
% xyCoordMaxBin=reshape([stats.Centroid], 2,[])'; %Still x,y pair

[y,x]               =find(peaksSAC); %Location of peaks in mn,n pair
xyCoordMaxBin       =[x,y];
clear x y

% Convert to a new reference frame which as the origin at the centre of the autocorr
% and with y negative at top and y positive at bottom, x negative on left and postive on
% right
% NB autocorr will always have sides with odd number of bins
xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPt), [size(xyCoordMaxBin,1), 1]);

%Get measure of distance from peak  as we just need to rank
[dstFrmCent, closestPeaks]   =sort(sqrt(sum(xyCoordMaxBinCentral.^2,2)));
closestPeaks        =closestPeaks(1:min(7,length(closestPeaks))); %Get id of closests
dstFrmCent          =dstFrmCent(1:min(7,length(closestPeaks))); %Includes centre
gScl                =mean(dstFrmCent(2:end));
if isempty(gScl), gScl=nan; end

%x,y pairs in cartesian coords with y counting down from top and origin at top left.
closestPeaksCoord   =xyCoordMaxBin(closestPeaks,:);
closestPeaksCoord   =round(closestPeaksCoord); %Should be integer



% --- FIND FIELD AROUND EACH PEAK -----------------------------------------
%As we've just defined peak of each field need to find the extent of the field around it. Define
%this as the area enclosed within the half-height of the peak. Do this by looping through eack peak:
%find areas of SAC above the half height, then find the patch that includes the peak of interest.
%Could be a problem is peaks overlap as they might.
%Create two results: a labelMatrix with the pixels of each field labeld with diff numbers (i.e. same
%IDs as used above) and a matrix with the perimeter of each peak marked

%Two 3d mats [sizeAutoCorr1, sizeAutoCorr2, numberPeaks plus centralPeak]
[peakMasks, perimeterMasks]=deal(zeros([size(sac), size(closestPeaksCoord,1)+1]));

for n=1:size(closestPeaksCoord,1)
    %NB. Have to flip dimensions as these are x,y pairs and need to be m,n
    [peakMasks(:,:,n), perimeterMasks(:,:,n)]=...
        sf_findPeakExtent(sac, closestPeaks(n), [closestPeaksCoord(n,2), closestPeaksCoord(n,1)]);
end
peakMask            =max(peakMasks, [], 3); %Colapse to 2d
lableMask           =peakMask; %Needed for code below




% ---- GRIDNESS -----------------------------------------------------------
%First setup for the normal gridness which defines are region bounded by
%the extent of the central six peaks and excludes a ring round the very
%centre peak
fields              =logical(lableMask); %1s for fields
maxDistFromCentre   =max(dstFrmCentre(fields));
maxDistFromCentre   =min(maxDistFromCentre, hlfMinSz); %Don't exceed size of sac
clear fields

%Define mask. Std excludes gentral peak
%and fits round edge of 6th furthest peak.
grdMskStd           =dstFrmCentre<=maxDistFromCentre; %Mask for points within bound of outer peak
grdMskStd           =grdMskStd & ~(lableMask==closestPeaks(1));

%Second setup for the expanding type of gridness as defined by Langston
%2010 SI page 11. Define a central peak radius as either the first local
%minimum in the curve of mean vs radius or first negative corr in same
%plot. Then adance in steps defined in vars above starting at point certain
%number of steps from centre and stopping certain number of points from
%edge (all defined in vars)

%So get mean corr as function of radius (in integer steps)
corrVsRad           =accumarray(dstFrmCentre(dstFrmCentre>0), ...
    sac((dstFrmCentre>0)), [], @mean);
corrVsRad(isnan(corrVsRad))=0;%v.infrequently get a nan - remove
radOfCent           =min(find(imregionalmin(corrVsRad),1), find(corrVsRad<0,1));


%For expanding find the furthest point of central peak from centre
mskExtRad           =radOfCent + expGMin : expGStep : hlfMinSz - expGMax;

%Under some circumstances it's not possible to find a central field or the
%central field is so large that we can't define a mask. maskExtRad will be
%empty. In these situations define the centre at being a radOfCent=5 then
%proceed as before
if isempty(mskExtRad),
    radOfCent=5;
    mskExtRad           =radOfCent + expGMin : expGStep : hlfMinSz - expGMax;
end

%Select the expanding masks to use from the prebuilt ones then subtract
%central field from all
grdMskExp           =allMsks(:,:,mskExtRad);
grdMskExp(repmat(dstFrmCentre<=radOfCent, [1,1, size(grdMskExp,3)]))=0;
clear mskExtRad radOfCent corrVsRad


%Get the rotated SACs and store
%NB rotation is counterclockwise - doesn't make any difference to gridness measure
rotSACs             =zeros([szSAC, 5]);
rotAmnt             =[60, 120, 30, 90, 150];
for  nn             =1:5
    rotSACs(:,:,nn)         =imrotate(sac,rotAmnt(nn), 'bilinear', 'crop');
end   
clear rotAmnt

%Now do the correlations to get gridness - normal first
stdG                =getG(sac, rotSACs, grdMskStd);

%Then do expanding gridness
tmpExpG             =zeros(size(grdMskExp,3),1);
for nn              =1:size(grdMskExp,3)
   tmpExpG(nn)              =getG(sac, rotSACs, grdMskExp(:,:,nn)); 
end
expG                =max(tmpExpG); %Take best



end


% --- SUB FUNCTIONS ----------------------------------------------------------------------
function [peakMask, perimeterMask]=sf_findPeakExtent(sac, peakId, peakCoord)
%Finds extent of field that belongs to each peak - defined as area in half-height and also
%perimieter. NB. peakCoord must by m,n pair in normal matrix coords

[peakMask, perimeterMask]=deal(zeros(size(sac)));

%Next line defines threshold used to find peak - currently using half height
aboveHalfHeightMask     =bwlabel(sac>(sac(peakCoord(1),peakCoord(2)).* (1/2)),8);
peakIdTemp              =aboveHalfHeightMask(peakCoord(1),peakCoord(2));
peakMask(aboveHalfHeightMask==peakIdTemp)=peakId;
perimeterMask(bwperim(aboveHalfHeightMask==peakIdTemp))=peakId;
end


function g          =getG(sac, rotSACs, mask)
%Get gridness for the current mask - does the corr with rotated sacs

allCorrs            =zeros(size(rotSACs,3),1);
for nn              =1:size(rotSACs,3)
   tmp              =rotSACs(:,:,nn);
   allCorrs(nn)     =corr(sac(mask), tmp(mask));
end

g                   =min(allCorrs([1, 2])) - max(allCorrs([3, 4, 5]));
    


end
