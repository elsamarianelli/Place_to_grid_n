function [ret, ellipse] = autoCorrProps(varargin)
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info

if nargin == 0, ret = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

autoCorr = single(real(in.sac)); %Ignore complex
autoCorr(isnan(autoCorr)) = -1;
autoCorrTemp = autoCorr;
autoCorrTemp(autoCorr<=0) = -1;
sizeAC = size(autoCorr);

% #######################################################################
% [Stage 1] Find peak pixels and identify the 7 closest to the centre

%1a. Find pixels which are local maxima
peaksMask = imregionalmax(autoCorrTemp);
clear autoCorrTemp

%1b. Label each group of contiguous maxima pixels with a unique number
[peaksLabel, numPeaks] = bwlabel(peaksMask, 8);

%1c. Maxima may extend over more than a single pixel, but we want a single point for
%    each peak.  Can either find the centroid (centre of mass) or just ignore all but
%    one of the pixels in each maxima.
if in.FIND_CENTROID
    stats = regionprops(peaksLabel, 'Centroid');
    xyCoordPeaks = reshape([stats.Centroid], 2,[])'; %  [n x 2] for n x-y pairs
else
    [y, x, lb] = find(peaksLabel);
    xyCoordPeaks(lb,:) = [x y];
    clear x y lb
end

%1d. Convert to a new reference frame which has the origin at the centre of the autocorr
%    and with y negative at top and y positive at bottom, x negative on left and postive on
%    right. Note autocorr will always have sides with odd number of bins
centralPoint = ceil(sizeAC/2); %m,n pair
xyCoordPeaksCentral = bsxfun(@minus,xyCoordPeaks,centralPoint([2 1]));

%1e. Calculate distance of peaks from centre point and find seven closest
%    (one will be central peak, but we deal with this later)
peaksDistToCentre = sum(xyCoordPeaksCentral.^2,2).^0.5;
[ignore, orderOfClose] = sort(peaksDistToCentre);

%1f.  Get id and coordinates of closest peaks - note closest peak 1 will be centre
closestPeaks = orderOfClose(1:min(7,end)); %might be less than 7 peaks
%x,y pairs in cartesian coords with y counting down from top and origin at top left.
closestPeaksCoord = xyCoordPeaks(closestPeaks,:);
closestPeaksCoord = round(closestPeaksCoord); %Should be integer
ret.closestPeaksCoordFull = closestPeaksCoord; % EM addition

% #######################################################################
% [Stage 2] Expand peak pixels into the surrounding half-height region

if in.FIELD_EXTENT_METHOD == 1
    
    %2i Prealocate two 3d mats [sizeAutoCorr1, sizeAutoCorr2, numberPeaks]
    peakLabel = zeros([size(autoCorr), numel(closestPeaks)]);
    perimeterLabel = zeros([size(autoCorr), numel(closestPeaks)]);
    
    %2ii For each peak, find the half-height region around it (see subfunction)
    for ii = 1:numel(closestPeaks)
        [peakLabel(:,:,ii), perimeterLabel(:,:,ii)] = ...
            sf_findPeakExtent(autoCorr, closestPeaks(ii), closestPeaksCoord(ii,[2 1]));
    end
    
    %2iii Colapse the 3d mats to 2d
    fieldsLabel = max(peakLabel,[],3);
    fieldsMask = logical(fieldsLabel);
    if in.GET_PERIM_FIELDS
        ret.perimFields = max(perimiterLabel,[],3);
    end
    clear peakLabel perimeterLabel
    
    
elseif in.FIELD_EXTENT_METHOD == 2
    %2a. Find the inverse-drainage-basin for each peak.
    fieldsLabel = watershed_(-autoCorr,peaksMask) + 1; %plus 1, because we want to be able to use label as index
    
    %2b. Work out what threshold to use in each drainage-basin
    nZones = max(fieldsLabel(:));
    closePeakInds = sub2ind(sizeAC,closestPeaksCoord(:,2),closestPeaksCoord(:,1));
    closesPeaksNewId = fieldsLabel(closePeakInds);
    thresholds = Inf(nZones,1);
    thresholds(closesPeaksNewId) = autoCorr(closePeakInds)/2;
    
    %2c. Apply thresholds to get a mask and updated labels
    fieldsMask = autoCorr > thresholds(fieldsLabel);
    fieldsLabel(~fieldsMask) = 0; %note label numbers are the ones from watershed not from finding peaks
    
    %2d. If requested, find the perimiter of the fields
    if in.GET_PERIM_FIELDS
        ret.perimFields = bwperim(fieldsMask); %not actually used in this function
    end
else
    error('FIELD_EXTENT_METHOD must be 1 or 2');
end



if in.GET_MEAN_R_AND_SCALE
    % #######################################################################
    % [Stage 3] Calculate a couple of metrics based on the closest peaks
    
    %3a. Find the (mean) autoCorr value at the closest peak pixles
    sumROfPeaks = accumarray(peaksLabel(peaksMask),autoCorr(peaksMask));
    countPeaksPixels = accumarray(peaksLabel(peaksMask),1);
    ret.meanROfPeaks = sumROfPeaks(closestPeaks)./countPeaksPixels(closestPeaks);
    
    %3b. Calculate scale measured in bins
    closestPeakDistFromCentre = peaksDistToCentre(closestPeaks(2:end));
    ret.scale = median(closestPeakDistFromCentre);
end

if in.GET_ORIENTATION
    % #######################################################################
    % [Stage 4] Calculate orientation
    
    % Note: negative y is necessary because DACQ camera coordinates have the
    % origin at top left hence y axis increase down (imagesc needs to be dispalyed with axis
    % ij as a result). But cart2pol assumes origin is bottom left so have to make correction
    % before using cart2pol.
    
    if numPeaks == 1 %If only central peak found
        ret.orientation = nan;
    else
        %4a. get cartesian coordinates of the closest peaks to centre, excluding centre itself
        closestPeaksCoordCentral = xyCoordPeaksCentral(closestPeaks(2:end),:);
        
        %4b. get angle to centre from peaks
        [th, ignore] = cart2pol(closestPeaksCoordCentral(:,1), closestPeaksCoordCentral(:,2));
        
        %4c. remove negative angles from the list, and return coordaintes
        %    and angles of remaining peaks.
        peaksAboveXaxis = find(th>=0); %Remove negative values - can do this as peaks are 180deg radially symetrical
        ret.closestPeaksCoord = closestPeaksCoordCentral(peaksAboveXaxis, :);
        ret.peaksOrient = rad2deg(sort(th(peaksAboveXaxis)));
        
        %4d. orientation of grid is defined as the angle to the first-peak, anticlockiwse
        ret.orientation = ret.peaksOrient(1);
    end
    
    
end

% #######################################################################
% [Stage 5] Calculate gridness

%5a. build a matrix giving the distance to centre for each pixel in autoCorr
xvals = (1:sizeAC(1))-centralPoint(1);
yvals = (1:sizeAC(2))-centralPoint(2);
distToCentre = bsxfun(@plus,xvals'.^2,yvals.^2).^0.5;
clear xvals yvals

%5b. looking at all the field pixels, find the furthest distance to the centre
maxDistFromCentre = NaN;
if length(closestPeaks) >= 7
    maxDistFromCentre = accumarray(fieldsLabel(fieldsMask),distToCentre(fieldsMask),[],@max);
    maxDistFromCentre = max(maxDistFromCentre);
end
%If neccessary, cap MaxDistFromCentre as being the shortest distance to the edge of the autoCorr
if isnan(maxDistFromCentre) || maxDistFromCentre > min(floor(sizeAC/2))
    maxDistFromCentre = min(floor(sizeAC/2));
end

%5c. Create a solid-circular mask of pixles that are within the MaxDistFromCentre
gridnessMask = (distToCentre <= maxDistFromCentre); %Mask for points within bound of outer peak

%5d. From the solid-circular mask, remove the central field.
centreMask = fieldsLabel == fieldsLabel(centralPoint(1),centralPoint(2));
gridnessMask = gridnessMask & ~centreMask;
clear distToCentre centreMask

%5e. If requested, find the perimiter of the gridness mask
if in.GET_PERIM_GRIDNESS_MASK
    ret.perimGridMask = bwperim(gridnessMask); %not actually used in this function
end

%5f. get the cropped and nan-masked version of the autocorr
W = ceil(maxDistFromCentre);
autoCorrMiddle = autoCorr;
autoCorrMiddle(~gridnessMask) = NaN;
autoCorrMiddle = autoCorrMiddle((-W:W) +centralPoint(1),(-W:W) +centralPoint(2));

%5g. Correlate the autoCorr middle with several rotated versions of itself.
if in.FULL_GRIDNESS_RANGE
    rotList = 0:6:174; % even if we use a longer list we don't actually do anything with the extra correlation values
else
    rotList = [60,120,30,90,150];
end
rotatedCorr = nan(size(rotList));
for ii = 1:numel(rotList)
    autoCorrMiddleRot = imrotate(autoCorrMiddle,rotList(ii), 'bilinear', 'crop');
    rotatedCorr(ii) = nancorr(autoCorrMiddle(:), autoCorrMiddleRot(:), 'Pearson');
end
clear autoCorrMiddle autoCorrMiddleRot

%5h. Calculate gridness, which is defined as:
%   {the lower of the correlations at 60 and 120deg}   minus
%   {the highest of the correlations at 30, 90 and 150deg}
if in.FULL_GRIDNESS_RANGE
    ret.gridness = min(rotatedCorr([11, 21])) - max(rotatedCorr([6, 16, 26]));
else
    ret.gridness = min(rotatedCorr([1, 2])) - max(rotatedCorr([3, 4, 5]));
end

% additional EM, get ellipticity measures 
[ret.ellipticity, ret.ellipse_metrics] = assess_ellipticity(closestPeaksCoord,[],'no');

if in.PLOT_ON
    PlotForDebug(in.sac,gridnessMask,fieldsMask,peaksAboveXaxis,...
        closestPeaksCoord,sizeAC,ret.scale,ret.orientation,centralPoint, ret.gridness);
end;

end

function PlotForDebug(sac,gridnessMask,fieldsMask,peaksAboveXaxis,...
    closestPeaksCoord,sizeAC,scale,orientation,centralPoint, gridness)

%a. get the outer ring of autoCorr in gray, but in an rgb matrix
im = sac;
im = (im+1)/2;
im2 = min(max(im,0),1); %cap any values that were over or under 1
im2(isnan(im)) = 1; %will be white
im2 = repmat(im2,[1 1 3]);

%b. get the region used in gridness calculation in jet colors
cm = jet(60);
im3 = round(im*60.999 + 0.5);
im3 = ind2rgb(im3,cm);
msk = gridnessMask & ~isnan(im);
msk = repmat(msk,[1 1 3]);
im2(msk) = im3(msk);

%c. display the image we've made
image(im2);
axis image
xlabel('bins'); ylabel('bins');
hold on;

%d. plot the field boundaries in a black on white double line
fields = bwboundaries(fieldsMask);
for k = 1:numel(fields)
    plot(fields{k}(:,2), fields{k}(:,1), 'w', 'Linewidth', 3)
    plot(fields{k}(:,2), fields{k}(:,1), 'k', 'Linewidth', 1)
end

%e. plot blue horizontal line to show where angle is measured from
for k=peaksAboveXaxis'
    plot(closestPeaksCoord([1 k+1],1),closestPeaksCoord([1 k+1],2),'Color', [1, 0.5, 0],'Linewidth',2);
end

%f. plot a red curve to show the orientation
plot([closestPeaksCoord(1,1) sizeAC(2)],closestPeaksCoord([1 1],[2 2]),'--b','Linewidth',2);
mag = scale*.75;
th = 0:-0.2:-deg2rad(orientation);
[x,y] = pol2cart(th,mag);
plot(x + centralPoint(2), centralPoint(1)-y,'b','Linewidth',2);

% %g. plot peak centres, starting from the second ring of peaks
% for k=2:8
%     mag = scale*k;
%     th = deg2rad((0:60/k:360)+orientation);
%     [centres_x, centres_y] = pol2cart(th,mag);
%     centres(:,1) = centres_x + centralPoint(2);
%     centres(:,2) = centralPoint(1) - centres_y;
%     %seems that Matlab actually does the clipping for us, but I'm not
%     %sure if it's supposed to, so lets do it anyway
%     clipped = centres(:,2) < 0 | centres(:,2) > sizeAC(1) | centres(:,1) < 0 | centres(:,1) > sizeAC(2);
%     centres(clipped,:) = [];
%     plot(centres(:,1),centres(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','w','MarkerSize',3);
%     clear centres
% end

% %h. plot fitted ellipses with semi major and minor lines
hold on;
[~, ~] = assess_ellipticity(closestPeaksCoord, gca, 'yes');
hold off;


hold off;
axis off;
lastBins = (fliplr(size(im2)));
% text(lastBins(2), lastBins(2),num2str(gridness));
% text(lastBins(2), lastBins(2)/2,num2str(scale));
end


% #######################################################################
%-% SUB FUNCTION, could use Matlab's watershed(im) instead
function L = watershed_(im,im_min)
%same as watershed form image toolbox, but takes im_min as inputs rather than
%calculating it here (would use imregionalmax)
%has to call an imagetoolbox private function, watershed_meyer.
persistent watershed_meyerH doBasic
conn = 8;

if exist('bwconncomp') ~= 2 || (~isempty(doBasic) && doBasic)
    L = watershed(im);
    return
end

try
    if isempty(watershed_meyerH)
        %get a handle to the private function
        old_cd = cd;
        cd([fileparts(which('watershed')) filesep 'private'])
        watershed_meyerH = @watershed_meyer;
        cd(old_cd);
    end
    
    %do the watershedding
    
    
    cc = bwconncomp(im_min, conn);
    L = watershed_meyerH(im,conn,cc);
catch e
    warning('watershed_ function is not improving the speed, might as well just use watershed:\n\n%s',e.message)
    L = watershed(im); %do it the normal way.
    doBasic = true;
end


end


% #######################################################################
%-% SUB FUNCTION
%Finds extent of field that belongs to each peak - defined as area in half-height and also
%perimieter. Note peakCoord must by m,n pair in normal matrix coords
function [peakLabel, perimeterLabel] = sf_findPeakExtent(autoCorr, peakId, peakCoord)

peakLabel = zeros(size(autoCorr));
perimeterLabel = zeros(size(autoCorr));

%Define threshold used to find peak - currently using half height
halfHeight = autoCorr(peakCoord(1),peakCoord(2))*0.5;
aboveHalfHeightLabel = bwlabel(autoCorr > halfHeight,8);

peakIdTemp = aboveHalfHeightLabel(peakCoord(1),peakCoord(2));
peakLabel(aboveHalfHeightLabel == peakIdTemp) = peakId;
perimeterLabel(bwperim(aboveHalfHeightLabel==peakIdTemp)) = peakId;
end




function in = DealWithInputs(varargin)
defaults.sac = [];
defaults.FIND_CENTROID = false;
defaults.GET_PERIM_FIELDS = false;
defaults.GET_PERIM_GRIDNESS_MASK = false;
defaults.GET_ORIENTATION = true;
defaults.GET_MEAN_R_AND_SCALE = true;
defaults.FULL_GRIDNESS_RANGE = false;
defaults.FIELD_EXTENT_METHOD = 2;
defaults.PLOT_ON = true;
defaults.PLOT_Ellipse_ON = true; % EM

VERSION = 1.03;

% Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    in = defaults;
    return;
end

if isstruct(varargin{1})
    if nargin > 1 && VERSION ~= varargin{2}
        error(['%s called with version number %g, but code is version %g.\n' ...
            'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
            'the function without a version number.'],mfilename,varargin{2},VERSION);
    end
    in = ModifyExistingFields(defaults,varargin{1});
else
    in = ModifyExistingFields(defaults,varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
