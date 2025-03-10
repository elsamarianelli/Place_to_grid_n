function [meanOffset, peakOffset] = led2v1_dir_offset(led2Dir, led1Dir)
%Compares animal dir vectors from 2led and 1led to determine any inherant bias.
% Animals heading can be determined either from postion of 2LEDs (more accurate but
% potential bias) on head or from the displacement of a single LED (less accurate). But in
% the case of 2LED that position of those relative to the animal's head has to be set
% correctly by the experimenter (pos of large and small led); this might be done
% incorrectly. Difference between 1LED heading and 2LED heading can be used to correct any
% bias (assumption: distribution of 2LED heading is centred on distribution of 1LED
% heading. One way in which this might not be correct is if the animals perferes to turn
% one way and 'looks' into the turn).
%
% NB. To apply correction need to subtract offSet from 2 led heading e.g.
%  tint.pos.dir=mod(tint.pos.dir - meanOffset, 360);
%
% TAKES
% led2Dir - direction vector calcuated from 2LEDs (in deg)
% led1Dir - direction vector calcualated from heading (normally single led) (in deg)
%           NB. values in degs increasing anti-clock from north
%
% RETURNS
% offSet - displacment of 2led distribution from the 1led distribution in deg where 0
%           indicates no displacement and positive indicates 2led distribution is
%           anti-clockwise of the 1 led distribution. Two values returned: first is
%           average offset and second is peak offset.
%
%  [meanOffset, peakOffset] = led2v1_dir_offset(tint.pos.dir, tint.pos.dir_disp)


dirDiff2from1=wrapTo180(led2Dir(:) - led1Dir(:));
meanOffset=rad2deg(circ_mean(deg2rad(dirDiff2from1(~isnan(dirDiff2from1))))); %Circ mean

%To get peak offset bin into 1deg bins and smooth with boxcar then take peak
binCentres=-179.5:179.5;
binnedDif=hist(dirDiff2from1, binCentres);
smthKern = ones(1,7) * (1/7);
binnedDif=imfilter(binnedDif, smthKern, 'circular');
[temp, peakInd]=max(binnedDif);
peakOffset=binCentres(peakInd(1));