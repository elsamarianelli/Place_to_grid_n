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