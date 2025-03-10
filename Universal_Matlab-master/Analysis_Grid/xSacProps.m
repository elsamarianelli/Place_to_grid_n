function stats =xSacProps(xsac)
%XSACPROPS Key metrics from a cross spatial auto corr
% Grid cells from the same module can be compared using a xSAC. This code
% extracts key metrics from the xSAC, foremost being the relative location
% of the central peak to the main axes.
%
% Code roughly based on orginal autoCorrProps.m
%
% ARGS
% xsac          The cross spatial autocorr (normally calculated using
%               xPearson) between two grid cell ratemaps (should be
%               smoothed and generally would be from the same module).
%
%
% RETURNS THE FOLLOWING AS FIELDS OF stats.xxxx
% phaseMag      The distance in bins of the ratemap of the central peak of
%               the xsac from the origin. The simplest indication of the
%               relative offset between grids. This value is in bins.
%
% phaseAng      Angle of the peak from the centre in radians measured CCW
%               from x-axis. This value will match an image of the SAC if
%               it is ploted with the axis top left.
%
% peakCoord     [1x2] m,n coordinates of the peak closest to the orgin.
%               Units are bins of SAC.
%
% peakCorr      [scalar] correlation value of the peak closest to the
%               origin
%
% axisAngDst    [3x2] angle, distance pairs indicating the position of the
%               three main axis relative to the central peak. Angle is
%               radians measured CCW from x-axis, distance is in bins.


% --- MAIN ---------------------------------------------------------------
% Ignore nans, complex and peaks below 0
xsac            =real(xsac);
xsac(isnan(xsac))=-1;
xsac(xsac<0)    =0;

% Now find peaks and identify the one closest to the centre
peaks           =imregionalmax(xsac); %a 2d logical mat
[peakM, peakN]  =find(peaks); %Peak coord
clear peaks

[sizeM, sizeN]  =size(xsac); %Determine the centre coord
centM           =sizeM/2 +0.5;
centN           =sizeN/2 +0.5;

%Now work out distance and angle of peaks from centre - flipping to measure
%angle (in rads) from CCW from x-axis
[ang, dist]     =cart2pol((peakN - centN), -(peakM - centM));
ang             =wrapTo2Pi(ang);
[dist2Cent, ind]=sort(dist); %Shortest first

stats.phaseMag  =dist2Cent(1); %Absolute distance of central peak from orig
stats.phaseAng  =ang(ind(1)); %Angle to central peak
stats.peakCoord =[peakM(ind(1)), peakN(ind(1))]; %Peak coords in absolute mn
stats.peakCorr  =xsac(stats.peakCoord(1),stats.peakCoord(2)); %Correlation value of the peak


% Now identify the six peaks around the central peak - basically do this by
% repeating the process above but using the central peak instead of the
% centre of the SAC.
[ang, dist]     =cart2pol((peakN - stats.peakCoord(2)), -(peakM - stats.peakCoord(1)));
[dist, ind]     =sort(dist); %Shortest first - shortest shoul be zero the centre
dist            =dist(2:7); %six closest excluding centre
ang             =wrapTo2Pi(ang);
ang             =ang(ind(2:7)); %corresponding angles
[ang, ind]      =sort(ang); %Put them in order - take first 3
stats.axisAngDst=[ang(1:3), dist(ind(1:3))];










