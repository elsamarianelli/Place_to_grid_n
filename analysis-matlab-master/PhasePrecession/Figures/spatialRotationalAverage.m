function data = spatialRotationalAverage(X,Y,P,I,vars)
%SPATIAL_ROTATIONAL_AVERAGE(X, Y, P,I,  vars)
%
% This function plots to gcf. It creates an image showing the rotational
% average at points in space.  (If that makes any sense!? There's a more
% complete explanation lower down.)
%
% [X] and [Y] give the coordinates at time t. [P] is the phase/angle at time t.
% [I] are indicies corresponding to spikes. (optional)
% [vars] contains the following elements:
% [kernelSize] e.g. 5; it's the first parameter for the gaussian smoothing;
% [kernelSigma] e.g. 1; it's the second paramater for the smoothing.
% [vars.figtitle] e.g. 'exp:1234 tet:2 cell:2' (optional)
% [imageSize] e.g. 25; it's 'basically' the size of the resulting image
% [saturation] e.g. how to set the saturation the hsv colour scale
%
% If you call the function with only one parameter:
% SPATIAL_ROTATIONAL_AVERAGE(vars)
% It produces a legend in the gcf
%
%%    Copyright (C) <2013>  <Daniel Manson> <d.manson@ucl.ac.uk> &
%%                          <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License

%Which mode are we in
if nargin == 1
    legend(X)% X is actually 'vars' in this case
else
    data = doit(X,Y,P,I,vars);
end

% the main body of the function is inside 'doit'.

function data = doit(X,Y,P,I,vars)
% Starting with X, Y and P...
% Convert the phase to an 'a component' and a 'b component'.  Now we take each
% component seperatly: bin the components (sum the values into bins) and
% then smooth with a gaussian kernel.  Then combine the two maps to find the
% average 'angle' (i.e. phase), which becomes the color, and the 'amplitude'
% of the average, which becomes the saturation.  'Average' here means
% the resultant vector.

imageSize = vars.imageSize;
kernelSize = vars.kernelSize;
kernelSigma = vars.kernelSigma;
sat = vars.saturation;

%Decompose phase
Pa = cos(P);
Pb = sin(P);

%Ensure origin (0,0)
X = X - min(X);
X(X<eps) = eps;
Y = Y - min(Y);
Y(Y<eps) = eps;

%Scale xy to size of image
scX = imageSize/max(X);
scY = imageSize/max(Y);

%Get path as an index into image
Xi = ceil(X*scX);
Xi(Xi>imageSize) = imageSize;
Yi = ceil(Y*scY);
Yi(Yi>imageSize) = imageSize;

%Reduce to positions of spike only if needed
if ~isempty(I)
    Xi = Xi(I);
    Yi = Yi(I);
    Pa = Pa(I);
    Pb = Pb(I);
end

%Define circular mask of radius half the size of the image
[x y] = meshgrid(((-imageSize/2)+0.5):((imageSize/2)-0.5),((-imageSize/2)+0.5):((imageSize/2)-0.5));
[t r] = cart2pol(x,y);
circ_mask = r<=imageSize/2;

%Create matrices of indexes into final image
[x y] = meshgrid(1:imageSize, 1:imageSize);

%Bin phase components into among pixels, interpolate missing values and mask
Ba = accumarray([Xi Yi],Pa);  %default accum function is '@sum'
idx = Ba~=0;
Ba = griddata(x(idx),y(idx),Ba(idx),x,y,'nearest');
Ba(~circ_mask) = 0;

%Bin phase components into among pixels, interpolate missing values and mask
Bb = accumarray([Xi Yi],Pb);
idx = Bb~=0;
Bb = griddata(x(idx),y(idx),Bb(idx),x,y,'nearest');
Bb(~circ_mask) = 0;

%Similarly, get occupancy histogram, interpolate missing values and mask
D = accumarray([Xi Yi],1);  %default accum function is '@sum'
idx = D~=0;
D = griddata(x(idx),y(idx),D(idx),x,y,'nearest');
D(~circ_mask) = 0;

%Gaussian smooth phase matrices and mask
% Ba = conv2(Ba,fspecial('gaussian',[kernelSize kernelSize],kernelSigma),'same')';
Ba = imfilter(Ba,fspecial('gaussian',[kernelSize kernelSize],kernelSigma),'symmetric','same','conv')';
Ba(~circ_mask) = 0;
% Bb = conv2(Bb,fspecial('gaussian', [kernelSize kernelSize],kernelSigma),'same')';
Bb = imfilter(Bb,fspecial('gaussian', [kernelSize kernelSize],kernelSigma),'symmetric','same','conv')';
Bb(~circ_mask) = 0;
% D = conv2(D,fspecial('gaussian',[kernelSize kernelSize],kernelSigma),'same')';
D = imfilter(D,fspecial('gaussian',[kernelSize kernelSize],kernelSigma),'symmetric','same','conv')';
D(~circ_mask) = 0;

%Reconstruct estimated phase at each pixel
B = mod(atan2(Bb,Ba),2*pi)/2/pi; %B should be in interval [0,1]

%Get saturation matrix for HSV image in the various cases
switch sat
    case 'legend'
        % For Legend
%         C = 2*r/imageSize;
        C = ones(size(D));
        C(~circ_mask) = 0;
        
    case 'even'
        % Fully saturated
        C = ones(size(D));
        C(~circ_mask) = 0;
        
    case 'var'
        % "Normal"
        C = sqrt(Bb.^2 + Ba.^2)./D;
        %         C = C - min(C(circ_mask));
        C(~circ_mask) = 0;
    case 'variance'
        C = 1-(sqrt(Bb.^2 + Ba.^2)./D);
        %         C = C - min(C(circ_mask));
        %         C(~circ_mask) = 0;
end

% This bit scales and caps the intensity.  Matlab interprets
% anything >= 1 as being 100% intensity.
if isfield(vars,'scale')
    C = C/imageSize;
elseif strcmp(sat,'variance')
    
else
    C = C/max(C(:));
end

fig = gcf;

if strcmp(sat,'variance')
    cmap = [jet(255);[1 1 1]];
    colormap(cmap)
    C(~circ_mask) = 1+(1/255);
    imagesc(C)
    axis square
    box off
    axis off
    colorbar
    data = C;
else
    %         subplot('position', [0 0 1 1]);
    cm = colormap('hsv');
    
    %manually apply the color map
    B = B*(size(cm,1)-2);
    B = round(B) + 1;
    
    %Get RGB image from the phase matrix according to colormap cm (HSV)
    data = ind2rgb(B,cm);
    % Convert to HSV
    data = rgb2hsv(data);
    %Apply appropriate saturation
    data(:,:,2) = C;
    %Value matrix is always 1.
    data(:,:,3) = 1;
    %Convert back to RGB
    data = hsv2rgb(data);
    imshow(data)%,'InitialMagnification','fit')
end
if isfield(vars,'figtitle') && ischar(vars.figtitle)
    set(fig,'NumberTitle','off','Name',vars.figtitle);
end
axis ij


function legend(vars)
t = 0:.01:(2*pi);
x = cos(t)*100 + 101;
y = -sin(t)*100 + 101;
% x = x+ rand(size(x)); %It's not really neccessary to add noise,
% y = y+ rand(size(y)); % but it helps to see what is going on.
% set(gcf,'Toolbar','none','MenuBar','none');
doit(x',y',t,[],vars);