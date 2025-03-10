function [ regSac ] = regularise_eliptic_grid( sac, abScale, orient  )
%REGULARISE_ELIPTIC_GRID Reshapes eliptical grid to be regular
% Similar to analyses run by other groups - having found elipticalness of
% grid using grid_xy_elipse_scale.m can then use this function to make grid
% points lie on circle and then redo gridness.
%
% Method is to first rotate the sac so that the major axis is aligned to
% the y-axis. Then use the major minor axis to adjust scale. Finally rotate
% back
%
% ARGS
% sac - spatial autocorr
% abScale - scale of major and minor axis
% orient - orientation in deg of major axis antiC from x-axis

%--- House keeping
abScale=sort(abScale, 'descend'); %Sometimes ab scale is mixed up - major shoudl be first


%--- Main function
%First check if need to do anything
if abScale(1)==abScale(2) %Grid is already regular
     regSac=sac;
    return
end

%Grids aren't regular to start to regularise
%1)First rotate so major axis aligns to x-axis
regSac=imrotate(sac, -orient, 'bilinear');

%2)Decide by how much to resize - major axis is x so work on y axis to
%bring to same scale
sacSize=size(regSac);
xStart=1;
xEnd=sacSize(2);
    tmp=sacSize(1)/2 - (sacSize(1)/2)*(abScale(2)/abScale(1));
    yStart=1+tmp;
    yEnd=sacSize(1)-tmp;

%3) Do the resample the sac to regularise the grid
regSac=interp2(regSac, linspace(xStart,xEnd, sacSize(2)), linspace(yStart, yEnd, sacSize(1))');

%4) Rotate back to original orientation
regSac=imrotate(regSac, orient, 'bilinear');

%5) Finally keep onlyl the central porition of the sac so that it matches
%the original size
sizeDif=round((size(regSac)-size(sac))/2);
regSac=regSac(1+sizeDif(1):end-sizeDif(1), 1+sizeDif(2):end-sizeDif(2));




end

