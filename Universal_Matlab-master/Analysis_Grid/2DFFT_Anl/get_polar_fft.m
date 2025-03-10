function [ psPol, th, frq ] = get_polar_fft( ps2D, maxRad, thBin, rBin, sigma  )
%GET_POLAR_FFT Reorganises 2d FFT into polar coords
% Having created a 2D FFT reorganises that data into polar coordinates -
% hence a radiusa and angle (measured anti-clock from x axis). Because 2d
% FFT is symetrical only works on half to save time. Can be limited to only
% look at spetial freqs lower (i.e. wavelengths higher) than some value
% maxRad.
%
% For angles remember that y increases down in ratemap and FFT, but we want
% to measure an angle anti clock from x axis. So need to invert angles. 
%
% NB. Code is fairly fast - slowest area is accumarry which takes 0.2s for 
% a 512N FFT. No point trying to make other vars persistent.
%
% ARGS
% ps2D - the [n x n] 2D FFT produced by get_2d_fft
% maxRad - ignore freq higher (wavelength smaller) than this in FFT points
% thBin - angular bin size to use in deg [typically 1 or 0.5] need to
%       increase number of FFT components to support smaller size
% rBin - radial bin size (i.e. FFT space) best to leave as 1 to match ps2D
% sigma [1] or [1,2] size of sigma for gaus smooth to use on polFFT.
%       Specify one value to use on both or [radialSig, angSig]
%
% RETURNS
% psPol - [nAngBin x nRadBin] 2D FFT reogranised into polar space & smthed
% th - [nAngBins x 1] angs in deg corresponding to dim1 of psPol
% frq - [nRadBins x 1] FFT components for freq (same as in 2d FFT) 

[fftSize]=single(length(ps2D)); %Size of power spec - both dim same size
centCord=floor(fftSize/2)+1;

%1) First get polar coords for all coords in the current 2dFFT
if maxRad %If maxRad doesnt equal 0
    x=-maxRad:+maxRad; %Take x value from centre +- maxRad
    y=-maxRad:0; %Just take half of y from -rad to centre
else
    x=(1:fftSize)-centCord;
    y=-centCord+1:0;
end
ps2D=ps2D(y+centCord, x+centCord); %Select just top half of ps2D and maxRad
%Not y value are -ve - this is correct as y axis inverted might need to add
%0.5 to get bin centres?
[xx,yy]=meshgrid(x,y);
clear x y

[th,r]=cart2pol(xx,-yy); %Get anlge (th) and radius (r)
th=abs(rad2deg(th)); %Convert to deg and change -180 to 180 [no other -ve vals]


%2)Now rebin
rSubs=round(r./rBin)+1;
thSubs=mod(round(th./thBin),180/thBin)+1; %Constrain values to be 0:179 but add 1 for binning
tmpMaxTh=(180/thBin);
tmpMaxR=max(rSubs(:));

%Get mean over bins and fill missing bins with nan
psPol=accumarray(single([thSubs(:),rSubs(:)]), ps2D(:), [tmpMaxTh, tmpMaxR], @(x) mean(x), nan);
psPol=psPol(:,1:maxRad+1);
th=0:thBin:180-thBin;
frq=0:rBin:maxRad;
clear rSubs thSubs tmp*


%3) Fill in missing values
%Ideally would interporlate but that's quite time consuming to do for just
%missing values. Next best is to smooth with gaussian not different sigma
%in x and y.
if sigma %Smooth is sigma not 0
    if length(sigma)==2
        rSig=sigma(1);
        thSig=sigma(2);
    else
        [rSig, thSig]=deal(sigma);
    end
    clear sig
    
    kernSize=ceil([thSig, rSig].*4); %Filter will be 4 times sigma
    kernSize=kernSize+mod(kernSize+1,2);%ensure they are odd
    [xx,yy]=meshgrid((1:kernSize(2))-ceil(kernSize(2)/2),...
        (1:kernSize(1))-ceil(kernSize(1)/2));
    kern=exp(1).^(-((xx.^2)./(2*rSig^2) + (yy.^2)./(2*thSig^2))); %2D gaussian
    clear xx yy
    
    %Now smooth but want to treat angular axis (y axis) as circular so it wraps
    %but not other axis (x-axis). So pad x-axis with nans and discard after
    %filtering - need as many nans as kern is wide.
    psPol=padarray(psPol, [0,kernSize(2)+1], nan, 'pre');
    norm=double(~isnan(psPol));
    psPol(isnan(psPol))=0; %replace nan for filter
    
    psPol=imfilter(psPol, kern, 'circular');
    norm=imfilter(norm, kern, 'circular');
    psPol=psPol./norm; %normalise
    psPol=psPol(:,kernSize(2)+2:end); %lose padded area
end

end


