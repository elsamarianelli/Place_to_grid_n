function [ ps2D ] = get_2d_fft( rm, fftSize, smthSig, maxRad )
%GET_2D_FFT Generate 2d FFT of ratemap
% Ideally takes an unsmoothed ratemap, mean normalises and returns
% frequencies upto some defined cut off. Can apply smoothing if necessary
% and choose max freq to return or 0 for all freq.
%
% ARGS
% rm - normal ratemap, ideally unsmoothed, unvisted bins are nan.
% fftSize - zero pad rm to be squaure with sides this size best if pow2
% smthSig- sigma of gaussian to smth power spec with (0=no smooth).
% maxRad - return centre of power spec <= this many bins radius (0=all)



rm=rm-nanmean(rm(:)); %Mean normalise - remove 0th fft component
rm(isnan(rm))=0; %Replace nans with 0
ps2D=abs(fft2(rm,fftSize,fftSize)); %Do fft with zero padding
ps2D=fftshift(ps2D); %Powerspectra - low freq to centre

if smthSig %Execute if not 0
    %If smoothing required use a guassian kern with sigma specific by
    %smthKern. Actual template is 4 x size of sigma to ensure that all of
    %guassian is in.
    kern=fspecial('gaussian', ceil(smthSig.*4), smthSig);
    %Now smooth - need to take account of number of bins included in
    %smoothing at edge then divide through. Filter2 doesn't do this and
    %will get edge effects otherwise.
    tmp=ones(size(ps2D));
    tmp=filter2(kern, tmp);
    ps2D=filter2(kern, ps2D);
    ps2D=ps2D./tmp; %Smoothed without edge effects
    clear tmp kern
end

if maxRad %Excute if not 0
    %Note for even number of bins centre is at (n/2)+1. So create meshgrid
    %to do Pythag.
    %Meshgrid is slow for large number of loops so make persistent but need
    %to check if it hasn't changed size i.e. if fftSize has not changed
    
    persistent fftSizeTmp maxRadTmp mask
    if isempty(fftSizeTmp) || ... %If persistent don't exist ..
            fftSizeTmp~=fftSize || ... or either change size
            maxRadTmp ~=maxRad%
        centCord=floor(fftSize/2)+1; % ... then create new values
        [xx,yy]=meshgrid(1:fftSize);
        xx=(xx-centCord).^2;
        yy=(yy-centCord).^2;
        mask=sqrt(xx+yy); %and get mask
        mask=mask>maxRad; %Values less than or equal to distance allowed
        
        fftSizeTmp=fftSize;
        maxRadTmp=maxRad;
    end
    
    ps2D(mask)=nan; %Set values outside allowed range to nan
end

end

