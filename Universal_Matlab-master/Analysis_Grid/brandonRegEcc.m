function [ rSac1, rSac2 ] = brandonRegEcc( sac)
%BRANDONREGECC Regularise elliptical SAC using Brandon(2011) method
% To deal with grids that were regular but slightly elliptical - and hence
% which have low gridness - Brandon (2011) introduced a method for
% correcting the sac to be more regular. Basically attempting to put the
% six central peaks back onto a circle.
%
% The method Brandon used and which is implimented here is a bit clunky and
% produces two correct SACs. The basic approach is to find six central
% peaks, then take the closest as set that as the minor axis of the elipse,
% then estimate the major axis of the elipse with respect to the distance
% from the centre of the furtherst point and the angle between it and the
% closest point. This gives one estimate of the major and minor axis. The
% second estimate is obtained by doing the opposite procedure. See Brandon
% (2011) SI for details but note there is at least one error in the
% equations specified - first equation requires a sqrt of a negative value
% to be taken which seems unlikely to be correct - I've just ignored the
% negative sign. Having found the major and minor axis the SAC is
% regularised using my preexisting function regularise_eliptic_grid.m which
% scales the major axis and minor axis to be the same length.
%
% IMPORTANT In general this method doesn't actually work that well -
% certianly not as well as trying to fit an ellipse - using least squares -
% to the 6 central points - but problem with the least squares is that
% sometimes it errors if it encounters multiple minima.
%
% ARGUMENTS
% sac       spatial autocorr made from smooth ratemap
%
% RETURNS
% rSac1     sac regularised by first method
% rSac2     sac regularised by second method
%
% NB. If too few peaks (<4) are found then the untransformed SAC is
% returned (need at least 4 points to fit elipse)
% NB2. If difference between major and minor axis is too great (> 2x) then
% no correction is attempted and original SAC is returned (this is a
% criteria that was used by Brandon)

%
% --- MAIN ----------------------------------------------------------------
% Code is in three chunks - first borrowed from autoCorrProps via
% multiGridness - finds the six central peaks. Second chunk then fits
% ellipses to these. Third does transformations.



% --1. FIND PEAKS ---------------------------------------------------------
% ---Preprocess the sac ---------------------------------------------------
sac=real(sac); %Ignore complex, shouldn't be present but can get introduced
sac(isnan(sac))=-1; %Sub nans for -1


% --- Find peaks ----------------------------------------------------------
sac=real(sac); %Ignore complex
sac(isnan(sac))     =-1; %Sub nans for -1
peaksSAC            =imregionalmax(sac);%Find local maxima in x,y coord

[y,x]               =find(peaksSAC); %Location of peaks in mn,n pair
xyCoordMaxBin       =[x,y];
clear x y


% Convert to a new reference frame which as the origin at the centre of the autocorr
% and with y negative at top and y positive at bottom, x negative on left and postive on
% right
% NB autocorr will always have sides with odd number of bins
szSAC               =size(sac);
centralPt           =ceil(szSAC/2); %m,n pair
xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPt), [size(xyCoordMaxBin,1), 1]);

%Get measure of distance from peak (not euclidian) as we just need to rank
%and exclude the central peaks
[~, closestPeaks]   =sort(sum(xyCoordMaxBinCentral.^2,2));
closestPeaks        =closestPeaks(2:min(7,length(closestPeaks))); %Get id of closests

%x,y pairs in cartesian coords with y counting down from top and origin at top left.
closestPeaksCoord   =xyCoordMaxBin(closestPeaks,:);
closestPeaksCoord   =round(closestPeaksCoord); %Should be integer

if length(closestPeaksCoord)<4
    [rSac1, rSac2]=deal(sac); %Too few peaks to regularise - just return sac
    return
end


% --2. ESTIAMTE AXES OF ELLIPSE -------------------------------------------

%Get angle in rads from x-axis and radius of each point. Notice inversion
%of y axis since sac has origin at top left and cart2pol assumes origin
%bottom left. Value for angle returned is in rads increasing anti-clock
%from x-axis
[pksAng, pksR]      =cart2pol(closestPeaksCoord(:,1), -closestPeaksCoord(:,2));
% [closeR, tmpInd]    =min(pksR); %Closest peak
closeAng            =pksAng(1); %Angle of closest from x-axis
closeR              =pksR(1);
% [distR, tmpInd]     =max(pksR); %Distant peak
distAng             =pksAng(end);
distR               =pksR(end);
clear tmp*

%Now get angle between cloest and furtherest field
theta               =mod(closeAng - distAng, 2*pi);

cosSqTh             =cos(theta)^2;
sinSqTh             =sin(theta)^2;


%Method 1. Minor axis set as vector to closest field - estimate major axis
m1MinR              =closeR; %Lenght minor axis
m1MajTh             =mod(closeAng+(pi/2), 2*pi); %Orient of major axis -at pi/2 to minor
m1MajR              =sqrt( cosSqTh /((sinSqTh/closeR^2)-(distR^(-2))));

%Method 2. Major axis set as vector to furthest field - estimate minor
m2MajR              =distR; %Lenght major axis
m2MajTh             =mod(distAng, 2*pi); %Orient of major axis
m2MinR              =sqrt( sinSqTh /(( closeR^(-2)) - (cosSqTh/distR^2)));


% -- 3. NOW DO TRANSFORMATIONS --------------------------------------------
%Method 1.
if m1MajR>(2*m1MinR) %Major axis is more than twice minor - no correction
    rSac1               =sac;
else
    %Do correction
    [ rSac1 ]       =regularise_eliptic_grid( sac, abs([m1MajR, m1MinR]), rad2deg(m1MajTh));
end

%Method 1.
if m2MajR>(2*m2MinR) %Major axis is more than twice minor - no correction
    rSac2               =sac;
else
    %Do correction
    [ rSac2 ]       =regularise_eliptic_grid( sac, [m2MajR, m2MinR], rad2deg(m2MajTh));
end






end

