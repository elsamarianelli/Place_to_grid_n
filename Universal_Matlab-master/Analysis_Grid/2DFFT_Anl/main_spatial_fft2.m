function [ ps2D, psPol, fCompAngDeg, fCompWaveCm ] = main_spatial_fft2( posData, spkInd )
%Note changed strategy relative to main_spatial_fft. Now do what I think is
%fair though computationally slower. Get 2dFFT, then polar for each
%shuffle. Get max power 95th criteria from these and threshold the polar
%FFT at this level. Then process as before (using secondary thresh_
%MAIN_SPATIAL_FFT Bins, gets 2D FFT, surrogates & sig components
% Master function for analysing FFT components of spatial cells - calls
% other functions to do most of the work. Starts with raw positional data
% and spike indicies, then bins and gets 2D FFT. Also shuffles spikes and
% repeates process to get confidence intervals.
%
% NB Requires that the Git Repo Universal_Matlab is in the path
%
% ARGS
% posData - pos branch of tint global data structure produces by mtint read
%           in functions i.e. tint.pos
% spkInd - [nSpkx1] index that spikes occured at where values are index
%           into posXY. e.g. data.tetrode.pos_sample
% Also requires access to function fftVars which carries a variable struct
%
% RETURNS
% ps2D      [N x N] 2D FFT of ratemap  (low freq at centre)
%
% psPol     [nDeg x maxRad N] ps2D transformed into radial coords (low freq
%           on left) and thresholded at 95% confidence intervals (i.e. by
%           subtracting power obtained from shuffles)
% fCompAngDeg - [nComponents x 1] angle in deg (antiClock from x-axis) of
%       sig FFT components. Note these lie at 90deg to grid axis.
% fCompWaveCm - [nCompoents x 1] wavelength in cm of FFT components will be
%       shorter than grid scale (for reg grid
%       fCompWaveCm=gridScale*cosd(30)


% ------------------------------------------------------------------------
% --- HOUSE KEEPING ------------------------------------------------------
% ------------------------------------------------------------------------
gcp; %Start parallel instances - used to speed shuffling
f           =fftVars; %Load vars structure from this folder

nPosPt      =length(posData.xy);
ppm         =str2double(posData.header{strcmp('pixels_per_metre', posData.header(:,1)),2});
spBin       =ppm/100*f.binSizeCm; %Spatial bin size in pix for a fixed size in cm
clear ppm

% ------------------------------------------------------------------------
% --- MAIN ---------------------------------------------------------------
% ------------------------------------------------------------------------

%1)Start by binning data to get a spatial ratemap (unsmoothed).
binXY       =bin_data('dwell', 'position', spBin, posData);
binSpk      =bin_data('spikes', 'position', spBin, posData, spkInd);
rm          =make_smooth_ratemap(binXY, binSpk, 1);  %ratemap with no smoothing


%2)Now get the real power spectra and convert to polar
ps2D        =get_2d_fft(rm, f.fftSize, f.smth2dSig, f.maxRad );
tmpThSig    =f.thSig/f.thBin; %Smoothing in angular space in bins from degs
tmpRSig     =f.rSig/f.rBin; %Smoothing in freq space in bins from FFT components
[psPol, th] = ...
    get_polar_fft( ps2D, f.maxRad , f.thBin, f.rBin, [tmpRSig, tmpThSig] );



%3)Make surrogates to compare against real power spec.
shufMaxPower=zeros(f.nIt,1); %Preallocate- max from each shuff
% tmpAll=zeros(size(psPol,1), size(psPol,2), f.nIt); %Store all ps  - for debug
minShuf=f.minShuf*50; %Min shift in pos pt from seconds

parfor nn=1:f.nIt
    %First shift spkind by some random amount, make new rm & power spec
    shufSpkInd  =sf_shuf_spk_ind(spkInd, nPosPt, minShuf);
    shufBinSpk  =bin_data('spikes', 'position', spBin, posData, shufSpkInd);
    shufRm      =make_smooth_ratemap(binXY, shufBinSpk, 1);  %ratemap with no smoothing
    %NL. Create 2d FFT & polar version
    shufPs2D    =get_2d_fft( shufRm, f.fftSize, f.smth2dSig, f.maxRad ); %shuf ps
    [ shufPsPol, ~, ~ ] = ...
        get_polar_fft( shufPs2D, f.maxRad , f.thBin, f.rBin, [tmpRSig, tmpThSig] );
    
    shufMaxPower(nn)=max(shufPsPol(:)); %Store max power from each ps
%     tmpAll(:,:,nn)=shufPsPol; %Store all power spec - for debug
end
%Get  threshold 95th
crit95=prctile(shufMaxPower,95);
clear shuf* minShuf

%4) Apply critical values to real ang spec and remove values below zero
psPol       =psPol-crit95;
psPol(psPol<=0)=nan;
angMean=nanmean(psPol,2); %Colapse to angluar space - smoothed in previous step
if ~any(angMean>0) %Check if there are peaks above 95th perc if so continue
    [ps2D, psPol, fCompAngDeg, fCompWaveCm]=deal([]); %Set to all empty and return
    return
end
clear  crit95


%Collapse polar rm to just mean power vs angular and get the peaks in deg
sndPeakThresh=max(angMean).*f.sndPeakCoef;
[angPeakPow,pInd]=... %angPeakPow is not retured - value at peak in ang plot
    findpeaks(angMean, 'minpeakheight', sndPeakThresh, 'minpeakdistance',f.minPeakDist/f.thBin);
fCompAngDeg=th(pInd); %Angle of fourier comp in deg

%Then from polar FFT find for each angle identified the freq componenet
%with the highest power hence get the freq and wavelength. Note wavelength
%will be less than for the grid scale for a totally regular grids
%wavelenght of fourier component will be cosd(30)* grid scale i.e. 0.866
[polPeakVal, peakFComp]=max(psPol(pInd,:),[],2); %polPeakVal not returned - val at peak in 2d
peakFComp=peakFComp-1; %Subtract 1 becuase first is 0th component
fCompWaveCm=f.fftSize./((1/f.binSizeCm).*peakFComp); %Get wavelgh of each FFT comp

%But note grid axis will


%Image for debug
% subplot(1,3,1), imagesc(ps2D), axis square;
% subplot(1,3,2), imagesc(psPol);
% subplot(1,3,3), plot(th,angMean);





end

function spkInd=sf_shuf_spk_ind(spkInd, nPosPt, minShuf)
%Shift all spk ind by some random value not less that +/- some value then
%wrap to be meaningful values
% ARGS
% spkInd - index of spikes in pos points
% nPosPt - total number of pos points - max ind cannot exceed htis
% minShuf - min shuffle allowed in pos pts +/-

randShft=ceil(rand(1).*(nPosPt-minShuf)) + minShuf; % between min shuff and maxind
spkInd=mod(spkInd+randShft, nPosPt); %Add rand shift and wrap
spkInd=max(1,spkInd); %Remve zeros and set to1.
end

